# standard
from pathlib import Path
import subprocess
import shutil
import tempfile
import time
import os

# class
import utility

def tool_call(call, logfile):
    logfile.write(f"Executing: {call}\n")
    rss_label = utility.get_time_rss_label()
    verbose = utility.get_time_verbose_flag()

    # initialize runtime and memory requirements
    runtime = 0
    mem = 0

    call_time = f"/usr/bin/time {verbose} {call}"
    start_time = time.time()
    result = subprocess.run(
        call_time,
        shell=True, # allows handling of '>' redirection
        stdout=subprocess.PIPE, # captures stdout of '/usr/bin/time'
        stderr=subprocess.PIPE, # captures stderr of '/usr/bin/time'
        text=True # decodes stdout/stderr as strings
    )
    end_time = time.time()
    runtime = round(end_time - start_time, 5)
    stderr_output = result.stderr
    logfile.write(stderr_output)
    if result.returncode != 0: # /usr/bin/time passes the exit code of the command on (127 when it is not found)
        raise RuntimeError(f"exit code {result.returncode} from: {call}\n{stderr_output}")
    rss_value = utility.get_rss_from_stderr(stderr_output, rss_label)
    if rss_value > 0: # get_rss_from_stderr returns -1 when the RSS line is not found
        mem = rss_value / utility.get_time_rss_per_mb()

    return runtime, mem

def index_size(path):
    """Bytes of an index: the file, or the files below the directory (giggle and igd write a directory, #61).
    Raises when there is no index: a missing path, or nothing written to it (no index is 0 bytes)."""
    path = Path(path)
    size = sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.is_dir() else path.stat().st_size
    if size == 0:
        raise RuntimeError(f"no index found at {path}")
    return size

def index_call(options, refdirs, label):
    """Tabix creates the index in the same folder as the input file."""
    print(f"Indexing {refdirs['ref'] / f'{label}.bed'} with {options.tool}...")

    runtime = 0
    mem = 0
    idx_size_mb = 0
    def step(call): # a measured step of the index: runtime summed, memory maxed
        nonlocal runtime, mem
        step_rt, step_mem = tool_call(call, options.logfile)
        runtime += step_rt
        mem = max(mem, step_mem)

    if options.tool in ("bedtools_sorted", "bedtk_sorted", "tabix"):
        # the "index" of the sorted variants is the sorted reference, read by the query step; tabix also
        # compresses and indexes it
        step(f"sort -k1,1 -k2,2n -k3,3n {refdirs['ref'] / f'{label}.bed'} > {refdirs['idx'] / f'{label}.bed'}")

    if options.tool == "tabix":
        step(f"bgzip -f {refdirs['idx'] / f'{label}.bed'} > {refdirs['idx'] / f'{label}.bed.gz'}")

        # determine size of the index (in MB) - gzipped and tabixed
        bgzip_size = os.stat(refdirs['idx'] / f'{label}.bed.gz').st_size
        bgzip_size_mb = round(bgzip_size/(1024*1024), 5)
        idx_size_mb += bgzip_size_mb

        # create tabix index
        step(f"tabix -f -C -p bed {refdirs['idx'] / f'{label}.bed'}.gz")

        csi_size = os.stat(refdirs['idx'] / f'{label}.bed.gz.csi').st_size
        csi_size_mb = round(csi_size/(1024*1024), 5)
        idx_size_mb += csi_size_mb

    elif options.tool in ("bedops", "bedmaps"):
        step(f"sort -k1,1 -k2,2n -k3,3n {refdirs['ref'] / f'{label}.bed'} > {refdirs['idx'] / f'{label}.bed'}")

    elif options.tool == "giggle":
        step(f"bash /giggle/scripts/sort_bed {refdirs['ref'] / f'{label}.bed'} {refdirs['idx']} 4")
        step(f"giggle index -i {refdirs['idx'] / f'{label}.bed.gz'} -o {refdirs['idx'] / f'{label}_index'} -f -s")

        indexpath = Path(options.datadir) / "bench" / options.benchname / options.tool
        # for some reason the giggle index is not created in ./giggle/idx/<index> but in ./giggle/<index> - so use this path
        giggle_size = index_size(indexpath / f'{label}_index')
        giggle_size_mb = round(giggle_size/(1024*1024), 5)
        idx_size_mb += giggle_size_mb

    elif options.tool == "igd":
        # copy the reference file to its own index directory
        idxindir = refdirs['idx'] / f'{label}_in'
        idxoutdir = refdirs['idx'] / f'{label}_out'
        idxindir.mkdir(parents=True, exist_ok=True)
        idxoutdir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(refdirs['ref'] / f'{label}.bed', idxindir / f'{label}.bed')

        step(f"igd create {idxindir} {idxoutdir} {label}")

        igd_size = index_size(idxoutdir)
        igd_size_mb = round(igd_size/(1024*1024), 5)
        idx_size_mb += igd_size_mb

    return runtime, mem, idx_size_mb


def sorted_genome(reffiles, label):
    """Genome file in the chromosome order of the sorted data (sort -k1,1), for `bedtools -sorted -g`: bedtools
    then handles chromosomes that only one of the files has. Written once next to the index, not measured."""
    genome = Path(reffiles['idx']).parent / f"{label}.genome"
    if not genome.exists():
        with open(genome, "w") as fh:
            subprocess.run(["sort", "-k1,1", str(reffiles['chromlens'])], stdout=fh, check=True)
    return genome


def query_call(options, label, reffiles, queryfile):
    tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False) # the tool's output, removed by the caller
    if os.stat(queryfile).st_size == 0: # an empty complex bin has no overlaps; granges and gia reject an empty file, tabix -R dumps the reference
        tmpfile.close()
        return 0, 0, tmpfile
    scratch = Path(tempfile.mkdtemp()) # intermediate files of the query, removed at the end

    query_rt = 0
    query_mem = 0
    def step(call): # a measured step of the query: runtime summed, memory maxed
        nonlocal query_rt, query_mem
        step_rt, step_mem = tool_call(call, options.logfile)
        query_rt += step_rt
        query_mem = max(query_mem, step_mem)

    def restore_duplicates(tool_output):
        """Unmeasured. A tool that reports a reference once however many queries hit it (bedtk flt, granges filter,
        bedIntersect -aHitAny) carries no pairing in its output, while the complex score counts one line per (query,
        reference) pair; bedtools prints each reported reference once per query it overlaps. The complex score of these
        tools thus checks the set of references they found, the pairs come from bedtools (#69; bedtk since v0.13)"""
        subprocess.run(f"bedtools intersect -wa -a {tool_output} -b {queryfile} > {tmpfile.name}", shell=True, check=True)

    if options.tool == "tabix":
        step(f"tabix {reffiles['idx']} -R {queryfile} > {tmpfile.name}")

    elif options.tool == "bedtools":
        step(f"bedtools intersect -wa -a {reffiles['ref-unsrt']} -b {queryfile} > {tmpfile.name}")

    elif options.tool == "bedtools_sorted":
        # first sort the query file
        query_sorted = scratch / "query_sorted.bed"
        step(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted}")

        # the sweep algorithm of bedtools (-sorted) on the sorted reference of the index step
        genome = sorted_genome(reffiles, label)
        step(f"bedtools intersect -sorted -g {genome} -wa -a {reffiles['idx']} -b {query_sorted} > {tmpfile.name}")

    elif options.tool == "bedops":
        # first sort the query file
        query_sorted = scratch / "query_sorted.bed"
        step(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted}")

        if "complex" in str(queryfile):
            step(f"bedmap --echo-map --multidelim '\n' {query_sorted} {reffiles['idx']} > {tmpfile.name}")
        else: # basic queries and arbitrary target/query pairs; the reference is the sorted one of the index step
            step(f"bedops --element-of 1 {reffiles['idx']} {query_sorted} > {tmpfile.name}")

    elif options.tool == "bedmaps":
        # first sort the query file
        query_sorted = scratch / "query_sorted.bed"
        step(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted}")
        step(f"bedmap --echo-map --multidelim '\n' {query_sorted} {reffiles['idx']} > {tmpfile.name}")

    elif options.tool == "giggle":
        # giggle treats the indexed and the query intervals as closed on both ends, so a reference that merely touches the
        # query is a hit (#70). A query shrunk to [start+1, end-1] gives exactly the half-open result against the closed
        # reference; a query of 1 bp becomes the point [start, start], which still hits a reference ending at start.
        # Prepared unmeasured; sort_bed needs the .bed suffix
        query_closed = scratch / "query_closed.bed"
        with open(queryfile) as fh, open(query_closed, "w") as out:
            for line in fh:
                chrom, start, end, *rest = line.rstrip("\n").split("\t")
                start, end = (int(start) + 1, int(end) - 1) if int(end) - int(start) >= 2 else (int(start), int(start))
                out.write("\t".join([chrom, str(start), str(end), *rest]) + "\n")
        step(f" bash /giggle/scripts/sort_bed {query_closed} {scratch} 4")

        indexpath = Path(options.datadir) / "bench" / options.benchname / options.tool
        # for some reason the giggle index is not created in ./giggle/idx/<index> but in ./giggle/<index> - so use this path
        step(f"/giggle/bin/giggle search -i {indexpath / f'{label}_index'} -q {scratch / 'query_closed.bed.gz'} -v > {tmpfile.name}")

    elif options.tool == "granges":
        # granges needs the .tsv suffix: copy reference and query (not measured)
        ref_tsv = scratch / "ref.tsv"
        query_tsv = scratch / "query.tsv"
        shutil.copy2(reffiles['ref-srt'], ref_tsv)
        shutil.copy2(queryfile, query_tsv)

        # granges 0.2.2 labels its query trees in natural chromosome order (numbers, X, Y, M) but reads the
        # genome file in file order, so any other order makes it compare the wrong chromosomes (#36); the
        # simulator writes the chromosomes in random order, so rewrite the genome file (not measured)
        genome = scratch / "genome.txt"
        lines = Path(reffiles['chromlens']).read_text().splitlines()
        genome.write_text("".join(line + "\n" for line in sorted(lines, key=lambda line: utility.chrom_sort_key(line.split("\t")[0]))))

        tmpfile2 = scratch / "tool_output.txt"
        step(f"granges filter --genome {genome} --left {ref_tsv} --right {query_tsv} > {tmpfile2}")
        restore_duplicates(tmpfile2) # granges filter keeps each left range once

    elif options.tool == "gia":
        step(f"gia intersect -a {queryfile} -b {reffiles['ref-unsrt']} -t > {tmpfile.name}")

    elif options.tool == "bedtk":
        tmpfile2 = scratch / "tool_output.txt"
        step(f"bedtk flt {queryfile} {reffiles['ref-unsrt']} > {tmpfile2}")
        restore_duplicates(tmpfile2) # bedtk flt reports each reference once

    elif options.tool == "bedtk_sorted":
        query_sorted = scratch / "query_sorted.bed"
        step(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted}")

        tmpfile2 = scratch / "tool_output.txt"
        step(f"bedtk flt {query_sorted} {reffiles['idx']} > {tmpfile2}")
        restore_duplicates(tmpfile2)

    elif options.tool == "awk":
        # determine the path to the awk script
        script_path = Path(__file__).parent / "tools" / "intersect_awk.py"
        step(f"python3 {script_path} -t {queryfile} -q {reffiles['ref-unsrt']} > {tmpfile.name}")

    elif options.tool == "intervaltree":
        # determine the path to the intervaltree script
        script_path = Path(__file__).parent / "tools" / "intersect_intervaltree.py"
        step(f"python3 {script_path} -q {queryfile} -t {reffiles['ref-unsrt']} -r target -o {tmpfile.name}")

    elif options.tool == "igd":
        tmpfile2 = scratch / "tool_output.txt"
        idxpath = Path(options.datadir) / "bench" / options.benchname / options.tool / "idx"
        step(f"igd search {idxpath / f'{label}_out' / f'{label}.igd'} -q {queryfile} -f > {tmpfile2}")

        # process the igd output to match the output of other tools (e.g., BED format)
        chrom = ""
        with open(tmpfile2) as fh, open(tmpfile.name, "w") as out:
            for line in fh:
                if line.startswith("Query"):
                    chrom = line.split(",")[0].split()[1]
                elif line[0].isdigit():
                    # extract start/end positons
                    parts = line.split()
                    start = parts[1].strip()
                    end = parts[2].strip()
                    out.write(f"{chrom}\t{start}\t{end}\n")

    elif options.tool == "ailist":
        tmpfile2 = scratch / "tool_output.txt"
        step(f"ailist {reffiles['ref-unsrt']} {queryfile} > {tmpfile2}")

        # process the ailist output to match the output of other tools (e.g., BED format)
        # extract the lines that contain the overlaps (4th column contains the number of overlaps) - repeat lines
        fh = open(tmpfile2)
        for line in fh:
            count = int(line.split()[3])
            if count != 0:
                for i in range(count):
                    tmpfile.write(line)
        fh.close()

    elif options.tool == "ucsc":
        tmpfile2 = scratch / "tool_output.txt"
        step(f"bedIntersect -aHitAny {reffiles['ref-unsrt']} {queryfile} {tmpfile2}")
        restore_duplicates(tmpfile2) # -aHitAny reports each reference once; without it bedIntersect prints the intersections, not the reference

    tmpfile.close()
    shutil.rmtree(scratch)

    return query_rt, query_mem, tmpfile
