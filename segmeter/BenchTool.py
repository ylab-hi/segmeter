import subprocess
import time
import os
import tempfile
from pathlib import Path
import shutil

import utility

class BenchTool:
    def __init__(self, options):
        self.options = options
        self.refdirs = self.get_refdirs() # get the reference directories
        self.querydirs = self.get_querydirs() # get the query directories

    def get_refdirs(self):
        """returns the input directories for the reference and query intervals"""
        refdirs = {}
        refdirs["ref"] = Path(self.options.outdir) / "sim" / self.options.format / "ref"
        refdirs["truth-basic"] = Path(self.options.outdir) / "sim" / self.options.format / "basic" / "truth"
        refdirs["truth-complex"] = Path(self.options.outdir) / "sim" / self.options.format / "complex" / "truth"

        # create directories for the index (if necessary - specified in idx_based_tools)
        if self.options.tool in self.options.idx_based_tools:
            refdirs["idx"] = Path(self.options.outdir) / "bench" / self.options.tool / "idx"
            refdirs["idx"].mkdir(parents=True, exist_ok=True) # need to be created (store the index files)

        return refdirs

    def get_querydirs(self):
        """returns the directories for the query intervals"""
        querydirs = {}
        querydirs["basic"] = {}
        for query in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
            querydirs["basic"][query] = Path(self.options.outdir) / "sim" / self.options.format / "basic" / "query" / query
        for query in ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]:
            querydirs["basic"][query] = Path(self.options.outdir) / "sim" / self.options.format / "basic" / "query" / query
        querydirs["complex"] = {}
        for query in ["mult"]:
            querydirs["complex"]["mult"] = Path(self.options.outdir) / "sim" / self.options.format / "complex" / "query" / query

        return querydirs

    def get_reffiles(self, label):
        reffiles = {}
        reffiles["ref-unsrt"] = self.refdirs["ref"] / f"{label}.bed"
        reffiles["ref-srt"] = self.refdirs["ref"] / f"{label}_sorted.bed"
        reffiles["ref"] = self.refdirs["ref"] / f"{label}.bed.gz"
        reffiles["truth-basic"] = self.refdirs["truth-basic"] / f"{label}.bed"
        reffiles["truth-complex"] = self.refdirs["truth-complex"] / f"{label}.bed"

        # create index if necessary
        if (self.options.tool == "tabix" or
            self.options.tool == "bedtools_sorted" or
            self.options.tool == "bedtools_tabix" or
            self.options.tool == "gia_sorted" or
            self.options.tool == "bedtk_sorted"):
                reffiles["idx"] = self.refdirs["idx"] / f"{label}.bed.gz"

        # if self.options.tool == "igd":


        # add genome length
        simpath = Path(self.options.outdir) / "sim" / self.options.format
        reffiles["chromlens"] =  simpath / f"{label}_chromlens.txt"

        return reffiles

    def get_queryfiles(self, label):
        """return the query files for the different query types according to the label"""
        queryfiles = {}
        queryfiles["basic"] = {}
        queryfiles["complex"] = {}
        for query in self.querydirs["basic"].keys():
            queryfiles["basic"][query] = {}
            for subset in range(10, 101, 10):
                queryfiles["basic"][query][subset] = self.querydirs["basic"][query] / f"{label}_{subset}p.bed"
        for query in self.querydirs["complex"].keys():
            queryfiles["complex"][query] = {}
            for bin in range(10, 101, 10):
                queryfiles["complex"][query][bin] = self.querydirs["complex"][query] / f"{label}_{bin}bin.bed"

        return queryfiles

    def init_stat(self, subset):
        """this function initializes the precision fields for the different query types"""
        fields = {}
        fields["basic"] = {}
        fields["basic"][subset] = {"TP": 0, "FP": 0, "TN": 0, "FN": 0, "negatives": []}
        fields["complex"] = {}
        fields["complex"][subset] = {"dist": 0}

        return fields

    def load_truth(self, truth_basic_file, truth_complex_file):
        truth = {}
        truth["basic"] = {}
        truth["complex"] = {}

        fh_basic = open(truth_basic_file)
        for line in fh_basic:
            fields = line.strip().split("\t")
            truth["basic"][(fields[0], fields[1], fields[2])] = ((fields[3], fields[4], fields[5]), fields[6])
        fh_basic.close()

        fh_complex = open(truth_complex_file)
        for line in fh_complex:
            fields = line.strip().split("\t")
            truth["complex"][(fields[0], fields[1], fields[2])] = fields[4] # only store number of records
        fh_complex.close()

        return truth

    def program_call(self, call):
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
        rss_value = utility.get_rss_from_stderr(stderr_output, rss_label)
        if rss_value:
            rss_value_mb = rss_value/(1000000)
            mem = rss_value_mb

        return runtime, mem


    def create_index(self, label, num):
        """Tabix creates the index in the same folder as the input file."""
        print(f"Indexing {self.refdirs['ref']} with {label}:{num} intervals...")

        runtime = 0
        mem = 0
        idx_size_mb = 0
        if (self.options.tool == "tabix" or
            self.options.tool == "bedtools_sorted" or
            self.options.tool == "bedtools_tabix" or
            self.options.tool == "gia_sorted" or
            self.options.tool == "bedtk_sorted"):

                sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {self.refdirs['ref'] / f'{label}.bed'} > {self.refdirs['idx'] / f'{label}.bed'}")
                runtime += sort_rt
                if sort_mem > mem:
                    mem = sort_mem

                bgzip_rt, bgzip_mem = self.program_call(f"bgzip -f {self.refdirs['idx'] / f'{label}.bed'} > {self.refdirs['idx'] / f'{label}.bed.gz'}")
                runtime += bgzip_rt
                if bgzip_mem > mem:
                    mem = bgzip_mem

                # determine size of the index (in MB) - gzipped and tabixed
                bgzip_size = os.stat(self.refdirs['idx'] / f'{label}.bed.gz').st_size
                bgzip_size_mb = round(bgzip_size/(1024*1024), 5)
                idx_size_mb += bgzip_size_mb

                # create tabix index
                if self.options.tool == "tabix" or self.options.tool == "bedtools_tabix":
                    tabix_rt, tabix_mem = self.program_call(f"tabix -f -C -p bed {self.refdirs['idx'] / f'{label}.bed'}.gz")
                    runtime += tabix_rt
                    if tabix_mem > mem:
                        mem = tabix_mem

                    csi_size = os.stat(self.refdirs['idx'] / f'{label}.bed.gz.csi').st_size
                    csi_size_mb = round(csi_size/(1024*1024), 5)
                    idx_size_mb += csi_size_mb


        elif self.options.tool == "giggle":
            sort_rt, sort_mem = self.program_call(f"bash /giggle/scripts/sort_bed {self.refdirs['ref'] / f'{label}.bed'} {self.refdirs['idx']} 4")
            runtime += sort_rt
            if sort_mem > mem:
                mem = sort_mem

            giggle_rt, giggle_mem = self.program_call(f"giggle index -i {self.refdirs['idx'] / f'{label}.bed.gz'} -o {self.refdirs['idx'] / f'{label}_index'} -f -s")
            runtime += giggle_rt
            if giggle_mem > mem:
                mem = giggle_mem

            indexpath = Path(self.options.outdir) / "bench" / self.options.tool
            """For some reason the giggle index is not created in ./giggle/idx/<index> but in ./giggle/<index> - so use this path"""
            giggle_size = os.stat(indexpath / f'{label}_index').st_size
            giggle_size_mb = round(giggle_size/(1024*1024), 5)
            idx_size_mb += giggle_size_mb

        elif self.options.tool == "igd":
            # copy the reference file to its own index directory
            idxindir = self.refdirs['idx'] / f'{label}_in'
            idxoutdir = self.refdirs['idx'] / f'{label}_out'
            idxindir.mkdir(parents=True, exist_ok=True)
            idxoutdir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(self.refdirs['ref'] / f'{label}.bed', idxindir / f'{label}.bed')

            igd_rt, igd_mem = self.program_call(f"igd create {idxindir} {idxoutdir} {label}")
            runtime += igd_rt
            if igd_mem > mem:
                mem = igd_mem

            igd_size = os.stat(idxoutdir).st_size
            igd_size_mb = round(igd_size/(1024*1024), 5)
            idx_size_mb += igd_size_mb


        return runtime, mem, idx_size_mb


    def query_interval(self, label, num, reffiles, queryfile):
        tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)

        query_rt = 0
        query_mem = 0
        if self.options.tool == "tabix":
            query_rt, query_mem = self.program_call(f"tabix {reffiles['idx']} -R {queryfile} > {tmpfile.name}")

        elif self.options.tool == "bedtools":
            query_rt, query_mem = self.program_call(f"bedtools intersect -wa -a {reffiles['ref-unsrt']} -b {queryfile} > {tmpfile.name}")

        elif self.options.tool == "bedtools_sorted":
            # first sort the query file
            query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}")

            query_rt += sort_rt
            if sort_mem > query_mem:
                query_mem = sort_mem
            # query intervals
            bedtools_rt, bedtools_mem = self.program_call(f"bedtools intersect -wa -a {reffiles['ref-srt']} -b {query_sorted.name} > {tmpfile.name}")
            query_rt += bedtools_rt
            if bedtools_mem > query_mem:
                query_mem = bedtools_mem

            query_sorted.close()

        elif self.options.tool == "bedtools_tabix":
            # first sort the query file
            query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}")

            query_rt += sort_rt
            if sort_mem > query_mem:
                query_mem = sort_mem

            bedtools_rt, bedtools_mem = self.program_call(f"bedtools intersect -wa -a {reffiles['ref-srt']} -b {queryfile} > {tmpfile.name}")
            query_rt += bedtools_rt
            if bedtools_mem > query_mem:
                query_mem = bedtools_mem

        elif self.options.tool == "bedops":
            # first sort the query file
            query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}")
            query_rt += sort_rt
            if query_mem > query_mem:
                query_mem = query_mem

            bedops_rt = 0
            bedops_mem = 0
            if "basic" in str(queryfile):
                bedops_rt, bedops_mem = self.program_call(f"bedops --element-of 1 {reffiles['ref-srt']} {query_sorted.name} > {tmpfile.name}")
            elif "complex" in str(queryfile):
                bedops_rt, bedops_mem = self.program_call(f"bedmap --echo-map --multidelim '\n' {query_sorted.name} {reffiles['ref-srt']} > {tmpfile.name}")
            query_rt += bedops_rt
            if bedops_mem > query_mem:
                query_mem = bedops_mem

            query_sorted.close()

        elif self.options.tool == "bedmaps":
            # first sort the query file
            query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}")
            query_rt += sort_rt
            if sort_mem > query_mem:
                query_mem = sort_mem

            bedops_rt, bedops_mem = self.program_call(f"bedmap --echo-map --multidelim '\n' {query_sorted.name} {reffiles['ref-srt']} > {tmpfile.name}")
            query_rt += bedops_rt
            if bedops_mem > query_mem:
                query_mem = bedops_mem

            query_sorted.close()

        elif self.options.tool == "giggle":
            query_sorted_dir = tempfile.TemporaryDirectory()
            sort_rt, sort_mem = self.program_call(f" bash /giggle/scripts/sort_bed {queryfile} {query_sorted_dir.name} 4")
            query_rt += sort_rt
            if sort_mem > query_mem:
                query_mem = sort_mem

            indexpath = Path(self.options.outdir) / "bench" / self.options.tool
            """For some reason the giggle index is not created in ./giggle/idx/<index> but in ./giggle/<index> - so use this path"""
            giggle_rt, giggle_mem = self.program_call(f"/giggle/bin/giggle search -i {indexpath / f'{label}_index'} -q {Path(query_sorted_dir.name) / f'{queryfile.name}.gz'} -v > {tmpfile.name}")
            query_rt += giggle_rt
            if giggle_mem > query_mem:
                query_mem = giggle_mem

            query_sorted_dir.cleanup()

        elif self.options.tool == "granges":
            # rename reference and query to .tsv
            ref_tsv = tempfile.NamedTemporaryFile(mode='w', delete=False)
            query_tsv = tempfile.NamedTemporaryFile(mode='w', delete=False)

            ref_tsv_name = ref_tsv.name + ".tsv"
            query_tsv_name = query_tsv.name + ".tsv"

            shutil.copy2(reffiles['ref-srt'], ref_tsv_name)
            shutil.copy2(queryfile, query_tsv_name)

            granges_rt, granges_mem = self.program_call(f"granges filter --genome {reffiles['chromlens']} --left {ref_tsv_name} --right {query_tsv_name} > {tmpfile.name}")
            query_rt += granges_rt
            if granges_mem > query_mem:
                query_mem = granges_mem

        elif self.options.tool == "gia":
            query_rt, query_mem = self.program_call(f"gia intersect -a {queryfile} -b {reffiles['ref-unsrt']} -t > {tmpfile.name}")

        elif self.options.tool == "gia_sorted":
            query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}")
            query_rt += sort_rt
            if sort_mem > query_mem:
                query_mem = sort_mem

            gia_rt, gia_mem = self.program_call(f"gia intersect -a {query_sorted.name} -b {reffiles['idx']} -t > {tmpfile.name}")
            query_rt += gia_rt
            if gia_mem > query_mem:
                query_mem = gia_mem

            query_sorted.close()

        elif self.options.tool == "bedtk":
            query_rt, query_mem = self.program_call(f"bedtk isec {queryfile} {reffiles['ref-unsrt']} > {tmpfile.name}")

        elif self.options.tool == "bedtk_sorted":
            query_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False)
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {queryfile} > {query_sorted.name}")
            query_rt += sort_rt
            if sort_mem > query_mem:
                query_mem = sort_mem

            bedtk_rt, bedtk_mem = self.program_call(f"bedtk isec {query_sorted.name} {reffiles['ref-srt']} > {tmpfile.name}")
            query_rt += bedtk_rt
            if bedtk_mem > query_mem:
                query_mem = bedtk_mem
            query_sorted.close()

        elif self.options.tool == "igd":
            tmpfile2 = tempfile.NamedTemporaryFile(mode='w', delete=False)
            igd_rt, igd_mem = self.program_call(f"igd search {self.refdirs['idx'] / f'{label}_out' / f'{label}.igd'} -q {queryfile} -f > {tmpfile2.name}")
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

        tmpfile.close()

        return query_rt, query_mem, tmpfile


    def query_intervals(self, label, num, subset):
        reffiles = self.get_reffiles(label) # get the reference files
        queryfiles = self.get_queryfiles(label)

        # load truths
        truth = self.load_truth(reffiles["truth-basic"], reffiles["truth-complex"])

        query_times = {}
        query_memory = {}
        query_precision = self.init_stat(subset)

        # print(f"Searching for overlaps in {subset}% of the queries in {num} intervals")

        last_message_len = 0
        for dtype in ["basic", "complex"]:
            query_times[dtype] = {}
            query_memory[dtype] = {}

            for qtype in queryfiles[dtype].keys():
                query_times[dtype][qtype] = {}
                query_memory[dtype][qtype] = {}

                # for subset in queryfiles[dtype][qtype]:
                # print(f"\rSearching for overlaps in {subset}% {dtype} '{qtype}' queries from {num} intervals {' '*20}", end="")
                current_message = f"\rSearching for overlaps in {subset}% {dtype} '{qtype}' queries from {num} intervals"
                clear_message = current_message + " " * (last_message_len - len(current_message))
                print(f"\r{clear_message}", end="")
                last_message_len = len(current_message)

                query_times[dtype][qtype][subset] = 0 # initialize the time
                query_memory[dtype][qtype][subset] = 0 # initialize the memory

                # determine the runtime and memory requirements
                query_rt, query_mem, query_result = self.query_interval(label, num, reffiles, queryfiles[dtype][qtype][subset])
                query_times[dtype][qtype][subset] += round(query_rt, 5)
                if query_mem > query_memory[dtype][qtype][subset]:
                    query_memory[dtype][qtype][subset] = query_mem

                # determine the precision of the tool
                precision = self.get_precision(queryfiles[dtype][qtype][subset], query_result, truth[dtype], dtype, qtype)
                if dtype == "basic":
                    query_precision[dtype][subset]["TP"] += precision["basic"]["TP"]
                    query_precision[dtype][subset]["FP"] += precision["basic"]["FP"]
                    query_precision[dtype][subset]["TN"] += precision["basic"]["TN"]
                    query_precision[dtype][subset]["FN"] += precision["basic"]["FN"]
                    query_precision[dtype][subset]["negatives"] += precision["basic"]["negatives"]
                elif dtype == "complex":
                    query_precision[dtype][subset]["dist"] += precision["complex"]["dist"]


        # final done message
        current_message = f"\rMatching all queries for {num} intervals using {self.options.tool} done!"
        clear_message = current_message + " " * (last_message_len - len(current_message))
        print(f"{clear_message}")
        last_message_len = len(current_message)

        return query_times, query_memory, query_precision

    def get_precision(self, queryfile, tmpfile, truth, dtype, qtype):
        """Check the file for the precision of the tool"""
        precision = {
            "basic": {"TP": 0, "FP": 0, "TN": 0, "FN": 0, "negatives": []},
            "complex": {"dist": 0}
        }

        if dtype == "basic":
            results = []
            # iterate/store the results from the queries
            fht = open(tmpfile.name)
            for line in fht:
                cols = line.strip().split("\t")
                result = tuple(cols[0:3])
                results.append(result)
            fht.close()

            # iterate/store the truth values from the queries
            fhq = open(queryfile)
            for line in fhq:
                cols = line.strip().split("\t")
                query = tuple(cols[0:3])
                truth_intvl = truth[query][0]
                truth_intvlid = truth[query][1]

                qgroup = utility.get_query_group("basic", qtype)
                if qgroup == "interval":
                    if truth_intvl in results:
                        precision["basic"]["TP"] += 1
                    else:
                        precision["basic"]["FN"] += 1
                        precision["basic"]["negatives"].append(f"{truth_intvlid}\tFN\n")

                elif qgroup == "gap":
                    if truth_intvl in results:
                        precision["basic"]["FP"] += 1
                        precision["basic"]["negatives"].append(f"{truth_intvlid}\tFP\n")
                    else:
                        precision["basic"]["TN"] += 1
            fhq.close()
        elif dtype == "complex":
            # if queryfile is empty, skip the precision calculation
            if os.stat(queryfile).st_size == 0:
                return precision

            results = []
            fht = open(tmpfile.name)
            lines = fht.readlines()
            results_entries_num = len(lines)
            fht.close()

            truth_entries_num = 0
            fht = open(queryfile)
            for line in fht:
                cols = line.strip().split("\t")
                query = tuple(cols[0:3])
                truthquery = truth[query]
                truth_entries_num += int(truthquery)
            fht.close()
            if results_entries_num != 0:
                precision["complex"]["dist"] = abs(results_entries_num - truth_entries_num)

        return precision
