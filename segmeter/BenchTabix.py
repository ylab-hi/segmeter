import subprocess
import time
import tempfile
from pathlib import Path

import utility

class BenchTabix:
    def __init__(self, options):
        self.options = options
        self.refdirs = self.get_refdirs() # get the reference directories
        self.querydirs = self.get_querydirs() # get the query directories

    def get_refdirs(self):
        """returns the input directories for the reference and query intervals"""
        refdirs = {}
        refdirs["ref"] = Path(self.options.outdir) / "sim" / self.options.format / "ref"
        refdirs["idx"] = Path(self.options.outdir) / "bench" / "tabix" / "idx"
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
        reffiles = {}
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
            fields["complex"][subset] = {"TP": 0, "FP": 0}
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

    def create_index(self, label, num):
        """Tabix creates the index in the same folder as the input file."""
        print(f"Indexing {self.refdirs['ref']} with {label}:{num} intervals...")
        start_time = time.time() # start the timer
        # sort bgzip and create index
        subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", f"{self.refdirs['ref'] / f'{label}.bed'}"],
            stdout=open(self.refdirs['idx'] / f'{label}.bed', 'w'))
        subprocess.run(["bgzip", "-f", f"{self.refdirs['idx'] / f'{label}.bed'}"])
        subprocess.run(["tabix", "-f", "-C", "-p", "bed", f"{self.refdirs['idx'] / f'{label}.bed'}.gz"])
        end_time = time.time() # end the timer
        duration = round(end_time - start_time, 5)
        return duration

    def create_index_mem(self, label, num):
        """Monitor the memory usage of the index creation"""
        print("\t... measuring memory requirements...")
        rss_label = utility.get_time_rss_label()
        verbose = utility.get_time_verbose_flag()
        result = subprocess.run(["/usr/bin/time", verbose, "tabix", "-f", "-C", "-p", "bed", f"{self.refdirs['idx'] / f'{label}.bed'}.gz"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stderr_output = result.stderr
        rss_value = utility.get_rss_from_stderr(stderr_output, rss_label)

        if rss_value:
            rss_value_mb = rss_value/(1000000)
            return rss_value_mb
        else:
            return -1

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
                    fh = open(queryfiles[dtype][qtype][subset])
                    query_times[dtype][qtype][subset] = 0 # initialize the time
                    query_memory[dtype][qtype][subset] = 0 # initialize the memory
                    for line in fh:
                        cols = line.strip().split("\t")
                        searchstr = f"{cols[0]}:{cols[1]}-{cols[2]}"
                        tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
                        start_time = time.time()
                        subprocess.run(["tabix", f"{reffiles['idx']}", f"{searchstr}"], stdout=tmpfile)
                        end_time = time.time()
                        query_times[dtype][qtype][subset] += round(end_time - start_time, 5)
                        tmpfile.close()

                        key = (cols[0], cols[1], cols[2])
                        precision = self.get_precision(tmpfile, truth[dtype][key], dtype, qtype)
                        query_precision[dtype][subset]["TP"] += precision["TP"]
                        query_precision[dtype][subset]["FP"] += precision["FP"]
                        if dtype == "basic":
                            query_precision[dtype][subset]["TN"] += precision["TN"]
                            query_precision[dtype][subset]["FN"] += precision["FN"]

                        # repeat the process for the memory measurement
                        query_memory[dtype][qtype][subset] = self.query_intervals_mem(
                            reffiles["idx"],
                            queryfiles[dtype][qtype][subset])

                print("done!")

        # print(f"memory {query_memory}")
        return query_times, query_memory, query_precision

    def query_intervals_mem(self, reffile, queryfile):
        # print("\t... measuring memory requirements...")
        rss_label = utility.get_time_rss_label()
        verbose = utility.get_time_verbose_flag()

        max_rss = 0
        # iterate through the query file
        with open(queryfile, "r") as fh:
            for line in fh:
                cols = line.strip().split("\t")
                searchstr = f"{cols[0]}:{cols[1]}-{cols[2]}"
                result = subprocess.run(["/usr/bin/time", verbose, "tabix", f"{reffile}", f"{searchstr}"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                stderr_output = result.stderr
                rss_value = utility.get_rss_from_stderr(stderr_output, rss_label)
                if rss_value:
                    rss_value_mb = rss_value/(1000000)
                    if rss_value_mb > max_rss:
                        max_rss = rss_value_mb
        return max_rss

    def get_precision(self, tmpfile, truth, dtype, qtype):
        """Check the file for the precision of the tool"""
        fho = open(tmpfile.name)
        found = False # flag to check if the interval could be found
        num_elements = 0
        for line in fho:
            splitted = line.strip().split("\t")
            if dtype == "basic":
                if tuple(splitted[0:3]) == truth:
                    found = True
            elif dtype == "complex":
                num_elements += 1
        fho.close()
        if dtype == "complex":
            if num_elements == int(truth):
                found = True

        precision = {"TP": 0, "FP": 0, "TN": 0, "FN": 0}
        if dtype == "basic":
            qgroup = utility.get_query_group("basic", qtype)
            if qgroup == "interval":
                if found:
                    precision["TP"] += 1
                else:
                    precision["FN"] += 1
            elif qgroup == "gap":
                if found:
                    precision["FP"] += 1
                else:
                    precision["TN"] += 1
        elif dtype == "complex":
            if found:
                precision["TP"] += 1
            else:
                precision["FP"] += 1

        return precision