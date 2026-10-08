# Standard
from pathlib import Path
import shutil
import sys

# Class
from BenchTool import BenchTool
import calls

class BenchBase:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        self.validate()
        if options.tool == "bedtools_tabix":
            print("WARNING: bedtools_tabix is deprecated and will be removed in 0.15.0: it measures bedtools_sorted "
                  "plus a tabix index that bedtools cannot use; use tabix to measure the index.", file=sys.stderr)

        benchpath = Path(options.datadir) / "bench" / self.options.benchname / options.tool
        benchpath.mkdir(parents=True, exist_ok=True)

        # list of index-based tools
        self.options.idx_based_tools = [
            "tabix", "bedtools_sorted", "bedtools_tabix", "giggle",
            "bedtk_sorted", "igd", "bedops", "bedmaps"
        ]

        if not self.options.simdata:
            # check if query and target files are provided
            if not options.query or not options.target:
                raise ValueError("Query and target files must be provided for benchmarking")

            # check if query and target files exist
            if not Path(options.query).exists():
                raise FileNotFoundError(f"Query file {options.query} does not exist")
            if not Path(options.target).exists():
                raise FileNotFoundError(f"Target file {options.target} does not exist")

        # open log file; closed also when a tool call fails, so its stderr lands in the log
        options.logfile = open(benchpath / "log.txt", "w")
        try:
            self.run(benchpath)
        finally:
            options.logfile.close()

    def run(self, benchpath):
        options, intvlnums = self.options, self.intvlnums
        self.tool = BenchTool(options)

        if not self.options.simdata:
            # determine if index has to be created
            if options.tool in self.options.idx_based_tools:
                print(f"Create index for {options.tool}...")
                idx_time, idx_mem, idx_size = calls.index_call(options, self.tool.refdirs, "target")
                # save index stats
                (benchpath / "index_stats.txt").write_text(f"time(s)\tmax_RSS(MB)\tindex_size(MB)\n{idx_time}\t{idx_mem}\t{idx_size}\n")
            print(f"Query intervals using {options.tool} for provided query and target...")
            query_rt, query_mem, query_result = calls.query_call(options, "target", self.tool.get_reffiles("target"), Path(options.query)) # giggle needs a Path
            shutil.move(query_result.name, benchpath / "result.bed") # the overlaps found by the tool
            # save query stats
            (benchpath / "query_stats.txt").write_text(f"time(s)\tmax_RSS(MB)\n{query_rt}\t{query_mem}\n")
        else:
            # determine the subsets (to be used)
            subsets = self.parse_param_subset()
            for i, (label, num) in enumerate(intvlnums.items()):
                print(f"Detect overlaps using {options.tool} for {num} intervals...({i+1} out of {len(intvlnums)})")
                labelpath = benchpath / label
                labelpath.mkdir(parents=True, exist_ok=True)

                # if the tool is index-based, create index (and record stats)
                if options.tool in self.options.idx_based_tools:
                    outfile_idx = labelpath / f"{label}_idx_stats.txt"
                    idx_time, idx_mem, idx_size = calls.index_call(options, self.tool.refdirs, label)
                    self.save_idx_stats(num, idx_time, idx_mem, idx_size, outfile_idx)

                statspath = labelpath / "stats"
                precisionpath = labelpath / "precision"
                statspath.mkdir(parents=True, exist_ok=True)
                precisionpath.mkdir(parents=True, exist_ok=True)

                # parse queries - but separately for each subset (e.g, 10,20,30,40,...)
                for subset in subsets:
                    outfile_stats = statspath / f"{label}_query_stats_{subset}.txt"
                    query_time, query_memory, query_precision = self.tool.query_intervals(label, num, subset)
                    self.save_query_stats(num, query_time, query_memory, outfile_stats)

                    # save query precision stats
                    outfile_precision = precisionpath / f"{label}_query_precision_{subset}.txt"
                    outfile_negatives = precisionpath / f"{label}_query_precision_negatives_{subset}.txt"
                    self.save_query_prec_stats(num, query_precision, outfile_precision, outfile_negatives)


    def save_idx_stats(self, num, idx_time, idx_mem, idx_size, filename):
        Path(filename).write_text(f"intvlnum\ttime(s)\tmax_RSS(MB)\tindex_size(MB)\n{num}\t{idx_time}\t{idx_mem}\t{idx_size}\n")

    def save_query_stats(self, num, query_time, query_mem, filename):
        with open(filename, "w") as fh:
            fh.write("intvlnum\tdata_type\tquery_type\ttime\tmax_RSS(MB)\n")
            for key, value in query_time["basic"].items(): # save basic query stats
                for key2, value2 in value.items():
                    fh.write(f"{num}\tbasic\t{key}_{key2}%\t{value2}\t{query_mem['basic'][key][key2]}\n")
            for key, value in query_time["complex"].items(): # save complex query stats
                for key2, value2 in value.items():
                    fh.write(f"{num}\tcomplex\t{key}_{key2}bin\t{value2}\t{query_mem['complex'][key][key2]}\n")

    def save_query_prec_stats(self, num, query_precision, filename, filename_negatives):
        with open(filename, "w") as fh:
            fh.write("intvlnum\tsubset\tTP\tFP\tTN\tFN\tPrecision\tRecall\tF1\n")
            for key, value in query_precision["basic"].items():
                precision = 0
                recall = 0
                f1 = 0
                if value["TP"] > 0:
                    precision = value["TP"] / (value["TP"] + value["FP"])
                    recall = value["TP"] / (value["TP"] + value["FN"])
                    f1 = 2 * ((precision * recall) / (precision + recall))
                fh.write(f"{num}\t{key}%\t{value['TP']}\t{value['FP']}\t{value['TN']}\t{value['FN']}\t")
                fh.write(f"{precision}\t{recall}\t{f1}\n")
            fh.write("\nintvlnum\tbin\tdistance\n")
            for key, value in query_precision["complex"].items():
                fh.write(f"{num}\t{key}bin\t{value['dist']}\n")

        with open(filename_negatives, "w") as fh:
            fh.write("intvlnum\tsubset\tFP\n")
            for value in query_precision["basic"].values():
                fh.writelines(value["negatives"])

    def parse_param_subset(self):
        subset = []
        for sub in self.options.subset.split(","):
            if "-" in sub:
                start, end = sub.split("-")
                # add all values between start and end (by 10 increments)
                for i in range(int(start), int(end)+1, 10):
                    subset.append(i)
            else:
                subset.append(int(sub))
        return subset

    def validate(self):
        if not Path(self.options.datadir).exists():
            raise FileNotFoundError(f"Directory {self.options.datadir} does not exist")
        simpath = Path(self.options.datadir) / "sim" / self.options.simname
        if self.options.simdata and not simpath.exists():
            raise FileNotFoundError(f"Directory {simpath} does not exist")
        if self.options.tool is None:
            raise ValueError("Tool not specified")
