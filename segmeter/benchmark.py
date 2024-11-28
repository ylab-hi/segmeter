import os
from pathlib import Path
import time
import subprocess
import logging
import utility

class BenchBase:
    def __init__(self, options, intvlnums, logger):
        self.options = options
        self.intvlnums = intvlnums
        self.logger = logger

        query_time = {}
        if options.tool == "tabix":
            self.tool = BenchTabix(options, intvlnums)
            index_time = self.tool.create_index()
            index_time_file = Path(options.outdir) / 'bench' / f'{options.tool}_index_time.tsv'
            self.save_index_time(index_time, index_time_file)

            # query_time = self.tool.query_intervals()

        elif options.tool == "bedtools":
            self.tool = BenchBEDTools(options, intvlnums)

        query_time_file = Path(options.outdir) / 'bench' / f'{options.tool}_query_time.tsv'
        self.save_query_time(query_time, query_time_file)

    def save_index_time(self, index_time, filename):
        fh = open(filename, "w")
        fh.write(f"intvlnum\ttime(s)\n")
        for key, value in index_time.items():
            fh.write(f"{key}\t{value}\n")
        fh.close()

    def save_query_time(self, query_time, filename):
        fh = open(filename, "w")
        fh.write(f"intvlnum\tquery_type\tsubset\ttime\n")
        for key, value in query_time.items():
            for key2, value2 in value.items():
                for key3, value3 in value2.items():
                    fh.write(f"{key}\t{key2}\t{key3}%\t{value3}\n")
        fh.close()

class BenchBEDTools:
    def __init__(self, options, intvlnums):
        self.options = options

class BenchTabix:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        # create dirs (needed for later use)

        # create a logger
        logfile = Path(options.outdir) / 'tabix.log'
        self.logger = utility.setup_logger(logfile, options.loglevel, "[tabix]")

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
        for query in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
            querydirs[query] = Path(self.options.outdir) / "bench" / "tabix" / self.options.datatype / "query" / query
        return querydirs

    def get_queryfiles(self, )

    def get_outdirs(self):
        query = Path(self.options.outdir) / "bench" / "tabix" / self.options.datatype / "query"
        outdirs = {}
        outdirs["perfect"] = query / "perfect"
        outdirs["3p-partial"] = query / "3p-partial"
        # outdirs[""]

        # outdir_query = Path(self.options.outdir) / "bench" / "tabix" / "simple" / "query"

        # outdir_identical = outdir / "identical"
        # outdir_partial = outdir / "partial"
        # outdir_enclosed = outdir / "enclosed"
        # outdir_identical.mkdir(parents=True, exist_ok=True)
        # outdir_partial.mkdir(parents=True, exist_ok=True)
        # outdir_enclosed.mkdir(parents=True, exist_ok=True)
        #


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
        indirs = self.get_indirs()

        # indir_ref = Path(self.options.outdir) / "sim" / "BED" / "simple" / "ref"
        # indir_query = Path(self.options.outdir) / "sim" / "BED" / "simple" / "query"
        #

        # outdir = Path(self.options.outdir) / "bench" / "tabix" / "simple" / "query"
        # outdir_identical = outdir / "identical"
        # outdir_partial = outdir / "partial"
        # outdir_enclosed = outdir / "enclosed"
        # outdir_identical.mkdir(parents=True, exist_ok=True)
        # outdir_partial.mkdir(parents=True, exist_ok=True)
        # outdir_enclosed.mkdir(parents=True, exist_ok=True)

        query_times = {}
        for i, (label, num) in enumerate(self.intvlnums.items()):
            print(f"Detect overlaps for {num} intervals...({i+1} out of {len(self.intvlnums)})")
            infiles = {
               "identical": indir_ref / f"{label}.bed",
               "partial": indir_query / "partial" / f"{label}.bed",
               "enclosed": indir_query / "enclosed" / f"{label}.bed"
            }
            query_times[label] = {}

            # load ground truth
            truth = self.load_truth(truthdir / f"{label}.bed")
            print(truth)


            # subsets 10,20,30
            subsets = [num for num in range(10, 101, 10)]
            for qtype, infile in infiles.items():
                query_times[label][qtype] = {}
                fh = open(infile)
                total = 0 # total number of (query) intervals
                for line in fh:
                    total += 1
                fh.close()
                for subset in subsets: # iterate through percentages (10%, 20%,...)
                    print(f"\rSubsetting {num} '{qtype}' intervals (by {subset}%) and searching for overlaps...", end="")
                    num_searched = 0
                    duration = 0.0

                    fh = open(infile)
                    for line in fh:
                        splitted = line.strip().split("\t")
                        searchst = f"{splitted[0]}:{splitted[1]}-{splitted[2]}"
                        start_time = time.time()
                        subprocess.run(["tabix", f"{indexdir / f'{label}.bed.gz'}", f"{searchst}"], stdout=open(outdir / "tmp", "w"))
                        end_time = time.time()
                        duration += end_time - start_time
                        num_searched += 1

                        # check if the search result is correct
                        print("bla results")
                        fh = open(outdir / "tmp")
                        for line in fh:
                            print(line)
                        print("bla end")

                        progress = float(num_searched/total)
                        if progress >= float(subset/100):
                            end_time = time.time()
                            query_times[label][qtype][subset] = round(duration, 3)
                            break
                print("done!")

        return query_times

    def load_truth(self, filename):
        truth = {}
        fh = open(filename)
        for line in fh:
            fields = line.strip().split("\t")
            truth[(fields[0], fields[1], fields[2])] = (fields[3], fields[4], fields[5])
        return truth
