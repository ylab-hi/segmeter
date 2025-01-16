import subprocess
import time
import os
import tempfile
from pathlib import Path

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
        refdirs["idx"] = Path(self.options.outdir) / "bench" / self.options.tool / "idx"
        refdirs["idx"].mkdir(parents=True, exist_ok=True) # need to be created (store the index files)
        refdirs["truth-basic"] = Path(self.options.outdir) / "sim" / self.options.format / "basic" / "truth"
        refdirs["truth-complex"] = Path(self.options.outdir) / "sim" / self.options.format / "complex" / "truth"

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
        reffiles["idx"] = self.refdirs["idx"] / f"{label}.bed.gz"
        reffiles["truth-basic"] = self.refdirs["truth-basic"] / f"{label}.bed"
        reffiles["truth-complex"] = self.refdirs["truth-complex"] / f"{label}.bed"

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

    def init_stat(self):
        """this function initializes the precision fields for the different query types"""
        fields = {}
        fields["basic"] = {}
        fields["complex"] = {}
        for subset in range(10,101,10):
            fields["basic"][subset] = {"TP": 0, "FP": 0, "TN": 0, "FN":0}
            fields["complex"][subset] = {"dist": 0}
        return fields

    def load_truth(self, truth_basic_file, truth_complex_file):
        truth = {}
        truth["basic"] = {}
        truth["complex"] = {}

        fh_basic = open(truth_basic_file)
        for line in fh_basic:
            fields = line.strip().split("\t")
            truth["basic"][(fields[0], fields[1], fields[2])] = (fields[3], fields[4], fields[5])
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
        if self.options.tool == "tabix":
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {self.refdirs['ref'] / f'{label}.bed'} > {self.refdirs['idx'] / f'{label}.bed'}")
            runtime += sort_rt
            if sort_mem > mem:
                mem = sort_mem

            bgzip_rt, bgzip_mem = self.program_call(f"bgzip -f {self.refdirs['idx'] / f'{label}.bed'} -o {self.refdirs['idx'] / f'{label}.bed'}.gz")
            runtime += bgzip_rt
            if bgzip_mem > mem:
                mem = bgzip_mem

            tabix_rt, tabix_mem = self.program_call(f"tabix -f -C -p bed {self.refdirs['idx'] / f'{label}.bed'}.gz")
            runtime += tabix_rt
            if tabix_mem > mem:
                mem = tabix_mem

            # determine size of the index (in MB)
            idx_size = os.stat(self.refdirs['idx'] / f'{label}.bed.gz.csi').st_size
            idx_size_mb = round(idx_size/(1024*1024), 5)

        elif self.options.tool == "bedtools_sorted":
            sort_rt, sort_mem = self.program_call(f"sort -k1,1 -k2,2n -k3,3n {self.refdirs['ref'] / f'{label}.bed'} > {self.refdirs['idx'] / f'{label}.bed'}")
            runtime = sort_rt
            mem = sort_mem

            idx_size = os.stat(self.refdirs['idx'] / f'{label}.bed').st_size
            idx_size_mb = round(idx_size/(1024*1024), 5)

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
            if query_mem > query_mem:
                query_mem = query_mem
            # query intervals
            query_rt, query_mem = self.program_call(f"bedtools intersect -wa -a {reffiles['ref-srt']} -b {query_sorted.name} > {tmpfile.name}")

            query_sorted.close()

        elif self.options.tool == "bedops":
            query_rt, query_mem = self.program_call(f"bedops --intersect {queryfile} {reffiles['ref-srt']} > {tmpfile.name}")


        tmpfile.close()

        return query_rt, query_mem, tmpfile


    def query_intervals(self, label, num):
        reffiles = self.get_reffiles(label) # get the reference files
        queryfiles = self.get_queryfiles(label)

        # load truths
        truth = self.load_truth(reffiles["truth-basic"], reffiles["truth-complex"])

        query_times = {}
        query_memory = {}
        query_precision = self.init_stat()

        for dtype in ["basic", "complex"]:
            query_times[dtype] = {}
            query_memory[dtype] = {}

            for qtype in queryfiles[dtype].keys():
                query_times[dtype][qtype] = {}
                query_memory[dtype][qtype] = {}

                for subset in queryfiles[dtype][qtype]:
                    print(f"\rSearching for overlaps in {subset}% of {num} {dtype} '{qtype}' queries...", end="")
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
                    elif dtype == "complex":
                        query_precision[dtype][subset]["dist"] += precision["complex"]["dist"]


                print("done!")

        # print(f"memory {query_memory}")
        return query_times, query_memory, query_precision

    def query_intervals_mem(self, reffiles, queryfile):
        # print("\t... measuring memory requirements...")
        rss_label = utility.get_time_rss_label()
        verbose = utility.get_time_verbose_flag()

        max_rss = 0
        result = subprocess.run(["/usr/bin/time", verbose] + self.query_call(reffiles, queryfile),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stderr_output = result.stderr
        rss_value = utility.get_rss_from_stderr(stderr_output, rss_label)
        if rss_value:
            rss_value_mb = rss_value/(1000000)
            if rss_value_mb > max_rss:
                max_rss = rss_value_mb
        return max_rss

    def get_precision(self, queryfile, tmpfile, truth, dtype, qtype):
        """Check the file for the precision of the tool"""
        precision = {
            "basic": {"TP": 0, "FP": 0, "TN": 0, "FN": 0},
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
                truthquery = truth[query]

                qgroup = utility.get_query_group("basic", qtype)
                if qgroup == "interval":
                    if truthquery in results:
                        precision["basic"]["TP"] += 1
                    else:
                        precision["basic"]["FN"] += 1
                elif qgroup == "gap":
                    if truthquery in results:
                        precision["basic"]["FP"] += 1
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
