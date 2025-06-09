# standard
import os
from pathlib import Path

# class
import utility
import calls
import shutil

class BenchTool:
    def __init__(self, options):
        self.options = options

        self.refdirs = self.get_refdirs() # get the reference directories
        if not self.options.simdata: # load the simulated data
            self.querydirs = self.get_querydirs() # get the query directories

    def get_refdirs(self):
        """returns the input directories for the reference and query intervals"""
        refdirs = {}
        if self.options.simdata:
            refdirs["ref"] = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format / "ref"
            refdirs["truth-basic"] = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format / "basic" / "truth"
            refdirs["truth-complex"] = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format / "complex" / "truth"

            # create directories for the index (if necessary - specified in idx_based_tools)
            if self.options.tool in self.options.idx_based_tools:
                refdirs["idx"] = Path(self.options.datadir) / "bench" / self.options.benchname / self.options.tool / "idx"
                refdirs["idx"].mkdir(parents=True, exist_ok=True) # need to be created (store the index files)
        else:
            refdirs["ref"] = Path(self.options.datadir) / "ref"
            refdirs["ref"].mkdir(parents=True, exist_ok=True) # need to be created (store the reference files)
            # copy target file to the reference directory (cp target.bed)
            target_file = Path(self.options.target)
            shutil.copy(target_file, refdirs["ref"] / "target.bed")

            #
            if self.options.tool in self.options.idx_based_tools:
                refdirs["idx"] = Path(self.options.datadir) / "bench" / self.options.benchname / self.options.tool / "idx"
                refdirs["idx"].mkdir(parents=True, exist_ok=True)

        return refdirs

    def get_querydirs(self):
        """returns the directories for the query intervals"""
        querydirs = {}
        querydirs["basic"] = {}
        for query in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
            querydirs["basic"][query] = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format / "basic" / "query" / query
        for query in ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]:
            querydirs["basic"][query] = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format / "basic" / "query" / query
        querydirs["complex"] = {}
        for query in ["mult"]:
            querydirs["complex"]["mult"] = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format / "complex" / "query" / query

        return querydirs

    def get_reffiles(self, label):
        reffiles = {}
        reffiles["ref-unsrt"] = self.refdirs["ref"] / f"{label}.bed"
        reffiles["ref-srt"] = self.refdirs["ref"] / f"{label}_sorted.bed"
        reffiles["ref"] = self.refdirs["ref"] / f"{label}.bed.gz"

        if self.options.simdata: # if simulated data, use the truth files
            reffiles["truth-basic"] = self.refdirs["truth-basic"] / f"{label}.bed"
            reffiles["truth-complex"] = self.refdirs["truth-complex"] / f"{label}.bed"

        # create index if necessary
        if (self.options.tool == "tabix" or
            self.options.tool == "bedtools_sorted" or
            self.options.tool == "bedtools_tabix" or
            self.options.tool == "bedtk_sorted"):
                reffiles["idx"] = self.refdirs["idx"] / f"{label}.bed.gz"

        elif self.options.tool == "gia_sorted":
            reffiles["idx"] = self.refdirs["idx"] / f"{label}.bed"

        # add genome length
        simpath = Path(self.options.datadir) / "sim" / self.options.simname / self.options.format
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


    def create_index(self, label, num):
        # runtime, mem, idx_size = calls.index_call(self.options, self.refdirs, label, num)
        return calls.index_call(self.options, self.refdirs, label, num)

    def query_interval_file(self, label, queryfile):
        """This function queries the intervals in the reference file with the query files"""
        reffiles = self.get_reffiles(label)
        query_rt, query_mem, query_result = calls.query_call(self.options, label, 0, reffiles, queryfile)

        return query_rt, query_mem, query_result

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
                query_rt, query_mem, query_result = calls.query_call(self.options, label, num, reffiles, queryfiles[dtype][qtype][subset])
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
            else:
                precision["complex"]["dist"] = truth_entries_num

        return precision
