# Standard
from pathlib import Path

# Class
from BenchTool import BenchTool

class BenchBase:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        query_time = {}
        query_precision = {}
        query_memory = {}

        benchpath = Path(options.outdir) / "bench" / options.tool
        benchpath.mkdir(parents=True, exist_ok=True)

        self.tool = BenchTool(options)
        if options.tool == "tabix":
            for i, (label, num) in enumerate(intvlnums.items()):
                outfile = Path(options.outdir) / "bench" / "tabix" / f"{label}_idx_stats.txt"
                idx_time = self.tool.create_index(label, num)
                idx_mem = self.tool.create_index_mem(label, num)
                self.save_idx_stats(num, idx_time, idx_mem, outfile) # write stats for the index creation
                print(f"Detect overlaps for {num} intervals...({i+1} out of {len(intvlnums)})")
                query_time, query_memory, query_precision = self.tool.query_intervals(label, num)

                outfile_query = Path(options.outdir) / "bench" / "tabix" / f"{label}_query_stats.txt"
                self.save_query_stats(num, query_time, query_memory, outfile_query)

                outfile_precision = Path(options.outdir) / "bench" / "tabix" / f"{label}_query_precision.txt"
                self.save_query_prec_stats(num, query_precision, outfile_precision)

        elif options.tool == "bedtools":
            for i, (label, num) in enumerate(intvlnums.items()):
                print(f"Detect overlaps for {num} intervals...({i+1} out of {len(intvlnums)})")
                query_time, query_memory, query_precision = self.tool.query_intervals(label, num)

                outfile_query = Path(options.outdir) / "bench" / "bedtools" / f"{label}_query_stats.txt"
                self.save_query_stats(num, query_time, query_memory, outfile_query)

                outfile_precision = Path(options.outdir) / "bench" / "bedtools" / f"{label}_query_precision.txt"
                self.save_query_prec_stats(num, query_precision, outfile_precision)

            options.tool = "bedtools_sorted"
            for i, (label, num) in enumerate(intvlnums.items()):
                print()



    def save_idx_stats(self, num, idx_time, idx_mem, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\ttime(s)\tmax_RSS(MB)\n")
        fh.write(f"{num}\t{idx_time}\t{idx_mem}\n")
        fh.close()

    def save_query_stats(self, num, query_time, query_mem, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\tdata_type\tquery_type\ttime\tmax_RSS(MB)\n")
        for key, value in query_time["basic"].items(): # save basic query stats
            for key2, value2 in value.items():
                fh.write(f"{num}\tbasic\t{key}_{key2}%\t{value2}\t{query_mem['basic'][key][key2]}\n")
        for key, value in query_time["complex"].items(): # save complex query stats
            for key2, value2 in value.items():
                fh.write(f"{num}\tcomplex\t{key}_{key2}bin\t{value2}\t{query_mem['complex'][key][key2]}\n")

    def save_query_prec_stats(self, num, query_precision, filename):
        fh = open(filename, "w")
        fh.write("intvlnum\tsubset\tTP\tFP\tTN\tFN\tPrecision\tRecall\tF1\n")
        for key, value in query_precision["basic"].items():
            precision = value["TP"] / (value["TP"] + value["FP"])
            recall = value["TP"] / (value["TP"] + value["FN"])
            f1 = 2 * ((precision * recall) / (precision + recall))
            fh.write(f"{num}\t{key}%\t{value['TP']}\t{value['FP']}\t{value['TN']}\t{value['FN']}\t")
            fh.write(f"{precision}\t{recall}\t{f1}\n")
        fh.write("\nintvlnum\tbin\tdistance\n")
        for key, value in query_precision["complex"].items():
            fh.write(f"{num}\t{key}bin\t{value['dist']}\n")
        fh.close()
