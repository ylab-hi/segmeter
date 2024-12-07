from pip._internal import req
import os
from pathlib import Path
import time
import subprocess
import tempfile
import platform
import re

# Class
import utility

class BenchBase:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        query_time = {}
        query_precision = {}
        query_memory = {}

        if options.tool == "tabix":
            self.tool = BenchTabix(options)
            for i, (label, num) in enumerate(intvlnums.items()):
                outfile = Path(options.outdir) / "bench" / "tabix" / "simple" / f"{label}_idx_stats.txt"
                idx_time = self.tool.create_index(label, num)
                idx_mem = self.tool.create_index_mem(label, num)
                self.save_idx_stats(num, idx_time, idx_mem, outfile) # write stats for the index creation
                print(f"Detect overlaps for {num} intervals...({i+1} out of {len(intvlnums)})")
                query_time, query_memory, query_precision = self.tool.query_intervals(label, num)



            # index_time = self.tool.create_index()
            # index_time_file = Path(options.outdir) / 'bench' / f'{options.tool}_index_time.tsv'
            # self.save_index_time(index_time, index_time_file)
            # query_time, query_memory, query_precision = self.tool.query_intervals()


        elif options.tool == "bedtools":
            self.tool = BenchBEDTools(options, intvlnums)

        # query_time_file = Path(options.outdir) / "bench" / f"{options.tool}_query_time.tsv"
        # self.save_query_time(query_time, query_time_file)
        # query_precision_file = Path(options.outdir) / "bench" / f"{options.tool}_query_precision.tsv"
        # self.save_query_precision(query_precision, query_precision_file)
        # query_memory_file = Path(options.outdir) / "bench" / f"{options.tool}_query_memory.tsv"
        # self.save_query_memory(query_memory, query_memory_file)
        #

    def save_idx_stats(self, num, idx_time, idx_mem, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\ttime(s)\tmax_RSS(MB)\n")
        fh.write(f"{num}\t{idx_time}\t{idx_mem}\n")
        fh.close()

    def create_refdirs(self):
        indir = Path(self.options.outdir) / "sim" / "BED" / "simple" / "ref" # load reference intervals
        outdir = Path(self.options.outdir) / "bench" / "tabix" / "simple" / "ref"
        outdir.mkdir(parents=True, exist_ok=True)


    def save_query_time(self, query_time, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\tquery_type\tsubset\ttime\n")
        for key, value in query_time.items():
            for key2, value2 in value.items():
                for key3, value3 in value2.items():
                    fh.write(f"{key}\t{key2}\t{key3}%\t{value3}\n")
        fh.close()

    def save_query_memory(self, query_memory, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\tquery_type\tsubset\tmax_RSS\n")
        for key, value in query_memory.items():
            for key2, value2 in value.items():
                for key3, value3 in value2.items():
                    fh.write(f"{key}\t{key2}\t{key3}%\t{value3}\n")
        fh.close()


    def save_query_precision(self, query_precision, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\tsubset\tTP\tFP\tTN\tFN\tPrecision\tRecall\tF1\n")
        for key, value in query_precision.items():
            for key2, value2 in value.items():
                fh.write(f"{key}\t{key2}%\t{value2['TP']}\t{value2['FP']}\t{value2['TN']}\t{value2['FN']}\t")
                precision = value2['TP'] / (value2['TP'] + value2['FP'])
                recall = value2['TP'] / (value2['TP'] + value2['FN'])
                f1 = 2*(precision*recall) / (precision + recall)
                fh.write(f"{precision}\t{recall}\t{f1}\n")
        fh.close()

class BenchBEDTools:
    def __init__(self, options, intvlnums):
        self.options = options

class BenchTabix:
    def __init__(self, options):
        self.options = options
        self.refdirs = self.get_refdirs() # get the reference directories
        self.querydirs = self.get_querydirs() # get the query directories

    def get_refdirs(self):
        """returns the input directories for the reference and query intervals"""
        refdirs = {}
        refdirs[self.options.datatype] = {} # if its simple or clustered
        refdirs[self.options.datatype]["ref"] = Path(self.options.outdir) / "sim" / self.options.format / self.options.datatype / "ref"
        refdirs[self.options.datatype]["idx"] = Path(self.options.outdir) / "bench" / "tabix" / self.options.datatype / "idx"
        refdirs[self.options.datatype]["idx"].mkdir(parents=True, exist_ok=True) # need to be created (store the index files)
        refdirs[self.options.datatype]["truth"] = Path(self.options.outdir) / "sim" / self.options.format / self.options.datatype / "truth"
        return refdirs

    def get_querydirs(self):
        """returns the directories for the query intervals"""
        querydirs = {}
        if self.options.datatype == "basic":
            querydirs["basic"] = {}
            for query in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
                querydirs[query] = Path(self.options.outdir) / "sim" / self.options.format / self.options.datatype / "query" / query
            for query in ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]:
                querydirs["basic"][query] = Path(self.options.outdir) / "sim" / self.options.format / self.options.datatype / "query" / query
        elif self.options.datatype == "complex":
            querydirs["complex"] = {}

        return querydirs

    def get_reffiles(self, label):
        reffiles = {}
        if self.options.datatype == "basic":
            reffiles = {}
            reffiles["basic"] = {
                "ref": self.refdirs["ref"] / f"{label}.bed.gz",
                "truth": self.refdirs["truth"] / f"{label}.bed"
            }
        if self.options.datatype == "complex":
            reffiles = {}
        return reffiles

    def get_queryfiles(self, label):
        """return the query files for the different query types according to the label"""
        queryfiles = {}
        queryfiles[self.options.datatype] = {}
        if self.options.datatype == "basic":
            queryfiles["basic"] = {}
            for query in self.querydirs["basic"]["intvlqueries"].keys():
                queryfiles["basic"]["intvlqueries"][query] = {}
                for subset in range(10, 101, 10):
                    queryfiles["basic"]["intvlqueries"][query][subset] = self.querydirs[query] / f"{label}_{subset}p.bed"
            queryfiles["basic"]["gapqueries"] = {}
            for query in self.querydirs["basic"]["gapqueries"].keys():
                queryfiles["basic"]["gapqueries"][query] = {}
                for subset in range(10, 101, 10):
                    queryfiles["basic"]["gapqueries"][query][subset] = self.querydirs[query] / f"{label}_{subset}p.bed"
        elif self.options.datatype == "complex":
            queryfiles["complex"]["intvlqueries"] = {}

        return queryfiles


    def init_stat(self):
        """this function initializes the precision fields for the different query types"""
        fields = {}
        for subset in range(10,101,10):
            fields[subset] = {"TP": 0, "FP": 0, "TN": 0, "FN":0}
        return fields

    def load_truth(self, filename):
        truth = {}
        fh = open(filename)
        for line in fh:
            fields = line.strip().split("\t")
            truth[(fields[0], fields[1], fields[2])] = (fields[3], fields[4], fields[5])
        return truth

    def monitor_memory(self, reffile, filename):
        """This function monitors the memory usage of tabix queries (needs to be done separately)"""
        # Determine which option to use for /usr/bin/time based on the OS
        if platform.system() == "Darwin":  # macOS
            verbose = "-l"
            rss_label = "maximum resident set size"
        else:  # Linux
            verbose = "-v"
            rss_label = "Maximum resident set size (kbytes)"

        # Open the file and process each line
        with open(filename, "r") as fh:
            rss_value = None
            for line in fh:
                cols = line.strip().split("\t")
                searchstr = f"{cols[0]}:{cols[1]}-{cols[2]}"  # chr1:1000-2000

                verbose = "-v"
                if platform.system() == "Darwin":
                    verbose = "-l"

                result = subprocess.run(["/usr/bin/time", verbose, "tabix", f"{reffile}", f"{searchstr}"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


                # Extract RSS from the output
                stderr_output = result.stderr

                for line in stderr_output.split("\n"):
                    if rss_label in line:
                        # Extract the numerical value from the line
                        match = re.search(r"(\d+)", line)
                        if match:
                            rss_value = int(match.group(1))
                            break

            if rss_value:
                rss_value_mb = rss_value/(1000000)
                return rss_value_mb
            else:
                return -1


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

        truth = self.load_truth(reffiles["truth"]) # load the ground truth

        query_times = {}
        query_memory = {}
        query_precision = self.init_stat()

        for qtype in queryfiles[self.options.datatype].keys():
            query_times[qtype] = {}
            query_memory[qtype] = {}
            for subset in queryfiles[self.options.datatype][qtype]:
                print(f"\rSearching for overlaps in {subset}% of {num//10} '{qtype}' intervals...", end="")
                # query_times[qtype][subset] = 0.0
                fh = open(queryfiles[self.options.datatype][qtype][subset])
                for line in fh:
                    cols = line.strip().split("\t")
                    searchstr = f"{cols[0]}:{cols[1]}-{cols[2]}"
                    tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
                    start_time = time.time() # start taking the time
                    subprocess.run(["tabix", f"{reffiles['ref']}", f"{searchstr}"], stdout=tmpfile)
                    end_time = time.time()
                    query_times[qtype][subset] = round(end_time - start_time, 5)
                    tmpfile.close()

                    # check file (and determine precision)
                    fho = open(tmpfile.name)
                    key = (cols[0], cols[1], cols[2])
                    found = False # flag to check if the interval could be found
                    for line in fho:
                        splitted = line.strip().split("\t")
                        if tuple(splitted[0:3]) == truth[key]:
                            found = True
                    fho.close()
                    qgroup = utility.get_query_group(self.options.datatype, qtype) # intvl or gap
                    if qgroup == "interval":
                        if found:
                            query_precision[qtype][subset]["TP"] += 1
                        else:
                            query_precision[qtype][subset]["FN"] += 1
                    elif qgroup == "gap":
                        if found:
                            query_precision[qtype][subset]["FP"] += 1
                        else:
                            query_precision[qtype][subset]["TN"] += 1

                # memory measurement
                rss_value = self.monitor_memory(reffiles["ref"], queryfiles[self.options.datatype][qtype][subset])
                query_memory[qtype][subset] = rss_value


            print("done!")

        return query_times, query_memory, query_precision





        # for i, (label, num) in enumerate(self.intvlnums.items()):
        #     # load reference interval and query intervals
        #     reffiles = self.get_reffiles(refdirs, label)
        #     queryfiles = self.get_queryfiles(querydirs, label)

        #     print(f"Detect overlaps for {num} intervals...({i+1} out of {len(self.intvlnums)})")
        #     query_times[label] = {}
        #     query_precision[label] = self.init_stat()
        #     query_memory[label] = {}

        #     # load ground truth
        #     truth = self.load_truth(reffiles["truth"])
        #     # print(truth)
        #     #

        #     for qtype in queryfiles.keys():
        #         query_times[label][qtype] = {}
        #         query_memory[label][qtype] = {}

        #         for subset in queryfiles[qtype].keys():
        #             # counter that keeps track of the total number of intervals
        #             # total_intvls = utility.file_linecounter(queryfiles[qtype][subset])
        #             print(f"\rSearching for overlaps in {subset}% of {num//10} '{qtype}' intervals...", end="")
        #             searched_intvls = 0
        #             time_duration = 0.0
        #             memory_stats = []

        #             fh = open(queryfiles[qtype][subset])
        #             for line in fh:
        #                 cols = line.strip().split("\t")
        #                 searchstr = f"{cols[0]}:{cols[1]}-{cols[2]}" # chr1:1000-2000
        #                 tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
        #                 start_time = time.time() # start taking the time
        #                 subprocess.run(["tabix", f"{reffiles['ref']}", f"{searchstr}"], stdout=tmpfile)
        #                 end_time = time.time()
        #                 time_duration += round(end_time - start_time, 5)
        #                 tmpfile.close()
        #                 searched_intvls += 1

        #                 # check File
        #                 fho = open(tmpfile.name)
        #                 key = (cols[0], cols[1], cols[2])
        #                 found = False # flag to check if the interval could be found
        #                 for line in fho:
        #                     splitted = line.strip().split("\t")
        #                     if tuple(splitted[0:3]) == truth[key]:
        #                         found = True
        #                 binary = self.get_qtype_bin(qtype)
        #                 if binary == "pos":
        #                     if found:
        #                         query_precision[label][subset]["TP"] += 1
        #                     else:
        #                         query_precision[label][subset]["FN"] += 1
        #                 elif binary == "neg":
        #                     if found:
        #                         query_precision[label][subset]["FP"] += 1
        #                     else:
        #                         query_precision[label][subset]["TN"] += 1
        #                 fho.close()
        #             query_times[label][qtype][subset] = time_duration
        #             query_precision[label]

        #             # memory measurement
        #             rss_value = self.monitor_memory(reffiles["ref"], queryfiles[qtype][subset])
        #             query_memory[label][qtype][subset] = rss_value

        #         print("done!")
        # return query_times, query_memory, query_precision
