import os
from pathlib import Path
import time
import subprocess
import tempfile

# Class
import utility

class BenchBase:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        query_time = {}
        query_precision = {}
        if options.tool == "tabix":
            self.tool = BenchTabix(options, intvlnums)
            index_time = self.tool.create_index()
            index_time_file = Path(options.outdir) / 'bench' / f'{options.tool}_index_time.tsv'
            self.save_index_time(index_time, index_time_file)
            query_time, query_precision = self.tool.query_intervals()

        elif options.tool == "bedtools":
            self.tool = BenchBEDTools(options, intvlnums)

        query_time_file = Path(options.outdir) / "bench" / f"{options.tool}_query_time.tsv"
        self.save_query_time(query_time, query_time_file)
        query_precision_file = Path(options.outdir) / "bench" / f"{options.tool}_query_precision.tsv"
        self.save_query_precision(query_precision, query_precision_file)

    def save_index_time(self, index_time, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\ttime(s)\n")
        for key, value in index_time.items():
            fh.write(f"{key}\t{value}\n")
        fh.close()

    def save_query_time(self, query_time, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\tquery_type\tsubset\ttime\n")
        for key, value in query_time.items():
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
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        # create dirs (needed for later use)

    def get_refdirs(self):
        """returns the input directories for the reference and query intervals"""
        refdirs = {}
        refdirs["ref"] = Path(self.options.outdir) / "bench" / "tabix" / self.options.datatype / "ref"
        refdirs["truth"] = Path(self.options.outdir) / "sim" / "BED" / self.options.datatype / "truth"
        return refdirs

    def get_reffiles(self, refdirs, label):
        reffiles = {}
        if self.options.datatype == "simple":
            reffiles = {
                "ref": refdirs["ref"] / f"{label}.bed.gz",
                "truth": refdirs["truth"] / f"{label}.bed"
            }
        return reffiles

    def get_querydirs(self):
        """returns the directories for the query intervals"""
        querydirs = {}
        for query in ["perfect","5p-partial","3p-partial","enclosed","contained",
            "perfect-gap","left-adjacent-gap","right-adjacent-gap","mid-gap1","mid-gap2"]:
            querydirs[query] = Path(self.options.outdir) / "sim" / "BED" / self.options.datatype / "query" / query

        return querydirs

    def get_queryfiles(self, querydirs, label):
        """return the query files for the different query types according to the label"""
        queryfiles = {}
        if self.options.datatype == "simple":
            for query in querydirs.keys():
                queryfiles[query] = {}
                for prct in range(10,101,10):
                    queryfiles[query][prct] = querydirs[query] / f"{label}_{prct}p.bed"
        return queryfiles

    def get_qtype_bin(self, qtype):
        """returns the binary (pos, neg) or the qtype"""
        if qtype in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
            return "pos"
        elif qtype in ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]:
            return "neg"
        else:
            return None

    def init_stat(self):
        """this function initializes the precision fields for the different query types"""
        fields = {}
        for subset in range(10,101,10):
            fields[subset] = {"TP": 0, "FP": 0, "TN": 0, "FN":0}
        return fields

    def create_index(self):
        """Tabix creates the index in the same folder as the input file."""
        # create folder
        indir = Path(self.options.outdir) / "sim" / "BED" / "simple" / "ref" # load reference intervals
        outdir = Path(self.options.outdir) / "bench" / "tabix" / "simple" / "ref"
        outdir.mkdir(parents=True, exist_ok=True)

        creation_times = {}
        print("Sort the BED file, compress (bgzip) and create the index...")
        for label, num in self.intvlnums.items():
            print(f"Processing {indir / f'{label}.bed'} with {label}:{num} intervals...")

            start_time = time.time() # start timer

            # sort, bgzip and create index - change to subprocess
            subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", f"{indir / f'{label}.bed'}"], stdout=open(outdir / f'{label}.bed', 'w'))
            subprocess.run(["bgzip", "-f", f"{outdir / f'{label}.bed'}"])
            subprocess.run(["tabix", "-f", "-C", "-p", "bed", f"{outdir / f'{label}.bed'}.gz"])

            end_time = time.time() # end timer
            duration = round(end_time - start_time, 5)
            creation_times[label] = duration
        return creation_times

    def query_intervals(self):
        refdirs = self.get_refdirs()
        querydirs = self.get_querydirs()

        query_times = {}
        query_precision = {}
        for i, (label, num) in enumerate(self.intvlnums.items()):
            # load reference interval and query intervals
            reffiles = self.get_reffiles(refdirs, label)
            queryfiles = self.get_queryfiles(querydirs, label)

            print(f"Detect overlaps for {num} intervals...({i+1} out of {len(self.intvlnums)})")
            query_times[label] = {}
            query_precision[label] = self.init_stat()

            # load ground truth
            truth = self.load_truth(reffiles["truth"])
            # print(truth)
            #

            for qtype in queryfiles.keys():
                query_times[label][qtype] = {}

                for subset in queryfiles[qtype].keys():
                    # counter that keeps track of the total number of intervals
                    total_intvls = utility.file_linecounter(queryfiles[qtype][subset])
                    print(f"\rSearching for overlaps in {subset}% of {num} '{qtype}' intervals...", end="")
                    searched_intvls = 0
                    time_duration = 0.0

                    fh = open(queryfiles[qtype][subset])
                    for line in fh:
                        cols = line.strip().split("\t")
                        searchstr = f"{cols[0]}:{cols[1]}-{cols[2]}" # chr1:1000-2000
                        tmpfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
                        start_time = time.time() # start taking the time
                        subprocess.run(["tabix", f"{reffiles['ref']}", f"{searchstr}"], stdout=tmpfile)
                        end_time = time.time()
                        time_duration += round(end_time - start_time, 5)
                        tmpfile.close()
                        searched_intvls += 1

                        # check File
                        fho = open(tmpfile.name)
                        key = (cols[0], cols[1], cols[2])
                        found = False # flag to check if the interval could be found
                        for line in fho:
                            splitted = line.strip().split("\t")
                            if tuple(splitted[0:3]) == truth[key]:
                                found = True
                        binary = self.get_qtype_bin(qtype)
                        if binary == "pos":
                            if found:
                                query_precision[label][subset]["TP"] += 1
                            else:
                                query_precision[label][subset]["FN"] += 1
                        elif binary == "neg":
                            if found:
                                query_precision[label][subset]["FP"] += 1
                            else:
                                query_precision[label][subset]["TN"] += 1
                        fho.close()
                    query_times[label][qtype][subset] = time_duration
                    query_precision[label]


                print("done!")
        return query_times, query_precision

    def load_truth(self, filename):
        truth = {}
        fh = open(filename)
        for line in fh:
            fields = line.strip().split("\t")
            truth[(fields[0], fields[1], fields[2])] = (fields[3], fields[4], fields[5])
        return truth
