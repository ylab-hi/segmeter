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
    rss_value = utility.get_rss_from_stderr(stderr_output, rss_label)
    if rss_value:
        rss_value_mb = rss_value/(1000000)
        mem = rss_value_mb

    return runtime, mem

def index_call(options, refdirs, label, num):
    """Tabix creates the index in the same folder as the input file."""
    print(f"Indexing {refdirs['ref']} with {label}:{num} intervals...")

    runtime = 0
    mem = 0
    idx_size_mb = 0
    if (options.tool == "tabix" or
        options.tool == "bedtools_sorted" or
        options.tool == "bedtools_tabix" or
        options.tool == "bedtk_sorted"):
            sort_rt, sort_mem = tool_call(f"sort -k1,1 -k2,2n -k3,3n {refdirs['ref'] / f'{label}.bed'} > {refdirs['idx'] / f'{label}.bed'}", options.logfile)
            runtime += sort_rt
            if sort_mem > mem:
                mem = sort_mem

            bgzip_rt, bgzip_mem = tool_call(f"bgzip -f {refdirs['idx'] / f'{label}.bed'} > {refdirs['idx'] / f'{label}.bed.gz'}", options.logfile)
            runtime += bgzip_rt
            if bgzip_mem > mem:
                mem = bgzip_mem

            # determine size of the index (in MB) - gzipped and tabixed
            bgzip_size = os.stat(refdirs['idx'] / f'{label}.bed.gz').st_size
            bgzip_size_mb = round(bgzip_size/(1024*1024), 5)
            idx_size_mb += bgzip_size_mb

            # create tabix index
            if options.tool == "tabix" or options.tool == "bedtools_tabix":
                tabix_rt, tabix_mem = tool_call(f"tabix -f -C -p bed {refdirs['idx'] / f'{label}.bed'}.gz", options.logfile)
                runtime += tabix_rt
                if tabix_mem > mem:
                    mem = tabix_mem

                csi_size = os.stat(refdirs['idx'] / f'{label}.bed.gz.csi').st_size
                csi_size_mb = round(csi_size/(1024*1024), 5)
                idx_size_mb += csi_size_mb


    elif options.tool == "giggle":
        sort_rt, sort_mem = tool_call(f"bash /giggle/scripts/sort_bed {refdirs['ref'] / f'{label}.bed'} {refdirs['idx']} 4", options.logfile)
        runtime += sort_rt
        if sort_mem > mem:
            mem = sort_mem

        giggle_rt, giggle_mem = tool_call(f"giggle index -i {refdirs['idx'] / f'{label}.bed.gz'} -o {refdirs['idx'] / f'{label}_index'} -f -s", options.logfile)
        runtime += giggle_rt
        if giggle_mem > mem:
            mem = giggle_mem

        indexpath = Path(options.datadir) / "bench" / options.benchname / options.tool
        """For some reason the giggle index is not created in ./giggle/idx/<index> but in ./giggle/<index> - so use this path"""
        giggle_size = os.stat(indexpath / f'{label}_index').st_size
        giggle_size_mb = round(giggle_size/(1024*1024), 5)
        idx_size_mb += giggle_size_mb

    elif options.tool == "igd":
        # copy the reference file to its own index directory
        idxindir = refdirs['idx'] / f'{label}_in'
        idxoutdir = refdirs['idx'] / f'{label}_out'
        idxindir.mkdir(parents=True, exist_ok=True)
        idxoutdir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(refdirs['ref'] / f'{label}.bed', idxindir / f'{label}.bed')

        igd_rt, igd_mem = tool_call(f"igd create {idxindir} {idxoutdir} {label}", options.logfile)
        runtime += igd_rt
        if igd_mem > mem:
            mem = igd_mem

        igd_size = os.stat(idxoutdir).st_size
        igd_size_mb = round(igd_size/(1024*1024), 5)
        idx_size_mb += igd_size_mb

    elif options.tool == "gia_sorted":
        sort_rt, sort_mem = tool_call(f"gia sort -i {refdirs['ref'] / f'{label}.bed'} -T bed4 -o {refdirs['idx'] / f'{label}.bed'}", options.logfile)
        runtime += sort_rt
        if sort_mem > mem:
            mem = sort_mem

    return runtime, mem, idx_size_mb

def query_call(options, label, num, reffiles, queryfile):
    tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)

    query_rt = 0
    query_mem = 0
    if options.tool == "tabix":
        query_rt, query_mem = tool_call(f"tabix {reffiles['idx']} -R {queryfile} > {tmpfile.name}", options.logfile)

    elif options.tool == "bedtools":
        query_rt, query_mem = tool_call(f"bedtools intersect -wa -a {reffiles['ref-unsrt']} -b {queryfile} > {tmpfile.name}", options.logfile)

    elif options.tool == "bedtools_sorted":
        # first sort the query file
        query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
        sort_rt, sort_mem = tool_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}", options.logfile)

        query_rt += sort_rt
        if sort_mem > query_mem:
            query_mem = sort_mem
        # query intervals
        bedtools_rt, bedtools_mem = tool_call(f"bedtools intersect -wa -a {reffiles['ref-srt']} -b {query_sorted.name} > {tmpfile.name}", options.logfile)
        query_rt += bedtools_rt
        if bedtools_mem > query_mem:
            query_mem = bedtools_mem

        query_sorted.close()

    elif options.tool == "bedtools_tabix":
        # first sort the query file
        query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
        sort_rt, sort_mem = tool_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}", options.logfile)

        query_rt += sort_rt
        if sort_mem > query_mem:
            query_mem = sort_mem

        bedtools_rt, bedtools_mem = tool_call(f"bedtools intersect -wa -a {reffiles['ref-srt']} -b {queryfile} > {tmpfile.name}", options.logfile)
        query_rt += bedtools_rt
        if bedtools_mem > query_mem:
            query_mem = bedtools_mem

    elif options.tool == "bedops":
        # first sort the query file
        query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
        sort_rt, sort_mem = tool_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}", options.logfile)
        query_rt += sort_rt
        if query_mem > query_mem:
            query_mem = query_mem

        bedops_rt = 0
        bedops_mem = 0
        if "basic" in str(queryfile):
            bedops_rt, bedops_mem = tool_call(f"bedops --element-of 1 {reffiles['ref-srt']} {query_sorted.name} > {tmpfile.name}", options.logfile)
        elif "complex" in str(queryfile):
            bedops_rt, bedops_mem = tool_call(f"bedmap --echo-map --multidelim '\n' {query_sorted.name} {reffiles['ref-srt']} > {tmpfile.name}", options.logfile)
        query_rt += bedops_rt
        if bedops_mem > query_mem:
            query_mem = bedops_mem

        query_sorted.close()

    elif options.tool == "bedmaps":
        # first sort the query file
        query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
        sort_rt, sort_mem = tool_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}", options.logfile)
        query_rt += sort_rt
        if sort_mem > query_mem:
            query_mem = sort_mem

        bedops_rt, bedops_mem = tool_call(f"bedmap --echo-map --multidelim '\n' {query_sorted.name} {reffiles['ref-srt']} > {tmpfile.name}", options.logfile)
        query_rt += bedops_rt
        if bedops_mem > query_mem:
            query_mem = bedops_mem

        query_sorted.close()

    elif options.tool == "giggle":
        query_sorted_dir = tempfile.TemporaryDirectory()
        sort_rt, sort_mem = tool_call(f" bash /giggle/scripts/sort_bed {queryfile} {query_sorted_dir.name} 4", options.logfile)
        query_rt += sort_rt
        if sort_mem > query_mem:
            query_mem = sort_mem

        indexpath = Path(options.datadir) / "bench" / options.benchname / options.tool
        """For some reason the giggle index is not created in ./giggle/idx/<index> but in ./giggle/<index> - so use this path"""
        giggle_rt, giggle_mem = tool_call(f"/giggle/bin/giggle search -i {indexpath / f'{label}_index'} -q {Path(query_sorted_dir.name) / f'{queryfile.name}.gz'} -v > {tmpfile.name}", options.logfile)
        query_rt += giggle_rt
        if giggle_mem > query_mem:
            query_mem = giggle_mem

        query_sorted_dir.cleanup()

    elif options.tool == "granges":
        # rename reference and query to .tsv
        ref_tsv = tempfile.NamedTemporaryFile(mode='w', delete=False)
        query_tsv = tempfile.NamedTemporaryFile(mode='w', delete=False)

        ref_tsv_name = ref_tsv.name + ".tsv"
        query_tsv_name = query_tsv.name + ".tsv"

        shutil.copy2(reffiles['ref-srt'], ref_tsv_name)
        shutil.copy2(queryfile, query_tsv_name)

        granges_rt, granges_mem = tool_call(f"granges filter --genome {reffiles['chromlens']} --left {ref_tsv_name} --right {query_tsv_name} > {tmpfile.name}", options.logfile)
        query_rt += granges_rt
        if granges_mem > query_mem:
            query_mem = granges_mem

    elif options.tool == "gia":
        query_rt, query_mem = tool_call(f"gia intersect -a {queryfile} -b {reffiles['ref-unsrt']} -t > {tmpfile.name}", options.logfile)

    elif options.tool == "gia_sorted":
        query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
        sort_rt, sort_mem = tool_call(f"gia sort -i {queryfile} -T bed4 -o {query_sorted.name}", options.logfile)
        query_rt += sort_rt
        if sort_mem > query_mem:
            query_mem = sort_mem

        gia_rt, gia_mem = tool_call(f"gia intersect --sorted -a {query_sorted.name} -b {reffiles['idx']} -t > {tmpfile.name}", options.logfile)
        query_rt += gia_rt
        if gia_mem > query_mem:
            query_mem = gia_mem

        query_sorted.close()

    elif options.tool == "bedtk":
        query_rt, query_mem = tool_call(f"bedtk isec {queryfile} {reffiles['ref-unsrt']} > {tmpfile.name}", options.logfile)

    elif options.tool == "bedtk_sorted":
        query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
        sort_rt, sort_mem = tool_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}", options.logfile)
        query_rt += sort_rt
        if sort_mem > query_mem:
            query_mem = sort_mem

        bedtk_rt, bedtk_mem = tool_call(f"bedtk isec {query_sorted.name} {reffiles['ref-srt']} > {tmpfile.name}", options.logfile)
        query_rt += bedtk_rt
        if bedtk_mem > query_mem:
            query_mem = bedtk_mem
        query_sorted.close()

    elif options.tool == "igd":
        tmpfile2 = tempfile.NamedTemporaryFile(mode='w', delete=False)
        idxpath = Path(options.datadir) / "bench" / options.benchname / options.tool / "idx"
        igd_rt, igd_mem = tool_call(f"igd search {idxpath / f'{label}_out' / f'{label}.igd'} -q {queryfile} -f > {tmpfile2.name}", options.logfile)
        query_rt += igd_rt
        if igd_mem > query_mem:
            query_mem = igd_mem

        # process the igd output to match the output of other tools (e.g., BED format)
        fh = open(tmpfile2.name)
        entries = []
        chrom = ""
        for line in fh:
            if line.startswith("Query"):
                chrom = line.split(",")[0].split()[1]
            elif line[0].isdigit():
                # extract start/end positons
                parts = line.split()
                start = parts[1].strip()
                end = parts[2].strip()
                entries.append(f"{chrom}\t{start}\t{end}\n")
        fh.close()

        fh = open(tmpfile.name, "w")
        for entry in entries:
            fh.write(entry)
        fh.close()

    elif options.tool == "ailist":
        tmpfile2 = tempfile.NamedTemporaryFile(mode='w', delete=False)
        ailist_rt, alilist_mem = tool_call(f"ailist {reffiles['ref-unsrt']} {queryfile} > {tmpfile2.name}", options.logfile)
        query_rt += ailist_rt
        if alilist_mem > query_mem:
            query_mem = alilist_mem

        # process the ailist output to match the output of other tools (e.g., BED format)
        # extract the lines that contain the overlaps (4th column contains the number of overlaps)
        fh = open(tmpfile2.name)
        for line in fh:
            if line.split()[3] != "0":
                tmpfile.write(line)
        fh.close()

    elif options.tool == "ucsc":
        query_rt, query_mem = tool_call(f"bedIntersect {reffiles['ref-unsrt']} {queryfile} {tmpfile.name}", options.logfile)

    tmpfile.close()



    return query_rt, query_mem, tmpfile
