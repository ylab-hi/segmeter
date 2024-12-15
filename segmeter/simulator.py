import random
from pathlib import Path
import subprocess

# class
import utility

class SimBase:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums

        if options.format == "BED":
            self.format = SimBED(options, self.intvlnums)

class SimBED:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums
        self.chroms = self.init_chroms()

    def init_chroms(self):
        """Initialize the chromosomes"""
        chroms = {}
        # consists of a list of chromosome that haven't exeeded a maximum length
        chroms["space-left"] = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        chroms["leftgap"] = {} # contains the end position of the last gap (left of interval)
        for chr in chroms["space-left"]:
            chroms["leftgap"][chr] = {}
        return chroms

    def select_chrom(self):
        """Randomly select a chromosome"""
        gs_start, gs_end = [int(x) for x in self.options.gapsize.split("-")] # get the gap size
        selection = random.choice(self.chroms["space-left"])
        if self.chroms["leftgap"][selection] == {}: # no interval has been placed on this chromosome (yet)
            # therefore simulate a (left) gap
            self.chroms["leftgap"][selection] = self.simulate_gap(1, random.randint(gs_start, gs_end))
            return selection
        else:
            if self.chroms["leftgap"][selection]["end"] < self.options.max_chromlen:
                return selection
            else:
                self.chroms["space-left"].remove(selection) # chromosome has reached maximum length
                allchroms = list(self.chroms["leftgap"].keys()) # get all chromosomes
                scaffolds = [chr for chr in allchroms if "SCF" in allchroms] # filter for SCF
                scaffold_name = f"SCF{len(scaffolds)+1}"
                self.chroms["space-left"].append(scaffold_name)
                self.chroms["leftgap"][scaffold_name] = self.simulate_gap(1, random.randint(gs_start, gs_end))

            return selection


    def init_leftgap(self):
        leftgap = {}
        gs_start, gs_end = [int(x) for x in self.options.gapsize.split("-")]
        for chr in self.chroms:
            leftgap[chr] = self.simulate_gap(1, random.randint(gs_start, gs_end))
        return leftgap

    def simulate_gap(self, start, end):
        gap = {}
        gap["start"] = start
        gap["end"] = end
        gap["mid"] = int((start + end) // 2)
        return gap

    def update_leftgap(self, chrom, start, end, mid):
        self.chroms["leftgap"][chrom]["start"] = start
        self.chroms["leftgap"][chrom]["end"] = end
        self.chroms["leftgap"][chrom]["mid"] = mid

    def det_rightmost_start(self): # needs to be set to 0
        """Simulate start position of first interval (aka right-most position)"""
        rightmost = {}
        for chr in self.chroms:
            rightmost[chr] = random.randint(100, 10000)
        return rightmost

    def create_datadirs(self, outdir):
        refdir = outdir / "ref"

        truthdirs = {}
        truthdirs["basic"] = outdir / "basic" / "truth"
        truthdirs["complex"] = outdir / "complex" / "truth"

        querydirs = {}
        querydirs["basic"] = {}
        for query_pos in ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]:
            querydirs["basic"][query_pos] = outdir / "basic" / "query" / query_pos
        for query_neg in ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]:
            querydirs["basic"][query_neg] = outdir / "basic" / "query" / query_neg
        querydirs["complex"] = {}
        for query_pos in ["mult"]:
            querydirs["complex"][query_pos] = outdir / "complex" / "query" / query_pos

        # create folder
        refdir.mkdir(parents=True, exist_ok=True)
        truthdirs["basic"].mkdir(parents=True, exist_ok=True)
        truthdirs["complex"].mkdir(parents=True, exist_ok=True)
        for dir in querydirs["basic"].keys():
            querydirs["basic"][dir].mkdir(parents=True, exist_ok=True)
        querydirs["complex"]["mult"].mkdir(parents=True, exist_ok=True)

        return refdir, truthdirs, querydirs

    def open_datafiles(self, label, refdir, truthdirs, querydirs):
        datafiles = {}
        datafiles["ref"] = open(refdir / f"{label}.bed", 'w')
        datafiles["truth-basic"] = open(truthdirs["basic"] / f"{label}.bed", 'w')
        datafiles["truth-complex"] = open(truthdirs["complex"] / f"{label}.bed", 'w')
        datafiles["queries-basic"] = {}
        for key in querydirs["basic"].keys():
            datafiles["queries-basic"][key] = open(querydirs["basic"][key] / f"{label}.bed", 'w')
        datafiles["queries-complex"] = {}
        for key in querydirs["complex"].keys():
            datafiles["queries-complex"][key] = open(querydirs["complex"][key] / f"{label}.bed", 'w')

        return datafiles

    def close_datafiles_basic(self, datafiles):
        """Close all (basic) datafiles"""
        datafiles["truth-basic"].close()
        for key in datafiles["queries-basic"].keys():
            datafiles["queries-basic"][key].close()

    def close_datafiles_complex(self, datafiles):
        """Close all (complex) datafiles"""
        datafiles["truth-complex"].close()
        for key in datafiles["queries-complex"].keys():
            datafiles["queries-complex"][key].close()

    def sort_datafiles(self, datatype, label, truthdirs, querydirs):
        """Sort some of the datafiles
        datatype -> [basic, complex]"""
        truth_in = truthdirs[datatype] / f"{label}.bed"
        truth_out = truthdirs[datatype] / f"{label}_sorted.bed"
        utility.sort_BED(truth_in, truth_out)

        for key in querydirs[datatype].keys():
            infile = querydirs[datatype][key] / f"{label}.bed"
            outfile = querydirs[datatype][key] / f"{label}_sorted.bed"
            utility.sort_BED(infile, outfile)

    def subset_basic_queryfiles(self, querydirs, label, num):
        for key in querydirs["basic"].keys():
            infile = querydirs["basic"][key] / f"{label}.bed"
            for subset in range(10,101,10):
                outfile = querydirs["basic"][key] / f"{label}_{subset}p.bed"
                with open(outfile, 'w') as subset_bed:
                    subprocess.run(["shuf", "-n", str(int(num * (subset / 100))), str(infile)], stdout=subset_bed)
                sorted = querydirs["basic"][key] / f"{label}_{subset}p_sorted.bed"
                utility.sort_BED(outfile, sorted)

    def sim_intervals(self):
        outpath = Path(self.options.outdir) / "sim" / "BED"
        refdir, truthdirs, querydirs = self.create_datadirs(outpath)
        is_start, is_end = [int(x) for x in self.options.intvlsize.split("-")]

        for label, num in self.intvlnums.items():
            print(f"Simulate intervals for {label}:{num}...")
            datafiles = self.open_datafiles(label, refdir, truthdirs, querydirs)

            for i in range(1, (int(num))+1):
                chrom = self.select_chrom()
                intvl = {} # create/simulate new interval
                intvl["id"] = f"intvl_{i}" # create an interval ID
                intvl["chrom"] = chrom
                intvl["start"] = self.chroms["leftgap"][chrom]["end"]+1
                intvl_size = random.randint(is_start, is_end)
                intvl["end"] = intvl["start"] + (intvl_size-1)
                intvl["mid"] = int((intvl["start"] + intvl["end"]) // 2) # determine middle of interval

                gs_start, gs_end = [int(x) for x in self.options.gapsize.split("-")]
                rightgap = self.simulate_gap(intvl["end"]+1, intvl["end"]+1+random.randint(gs_start,gs_end))

                self.sim_basic_queries(datafiles, intvl, rightgap)

            datafiles["ref"].close() # close the reference file
            self.close_datafiles_basic(datafiles) # close the basic datafiles
            utility.sort_BED(datafiles["ref"] / f"{label}.bed", datafiles["ref"] / f"{label}_sorted.bed") # sort the reference file
            self.sort_datafiles("basic", label, truthdirs, querydirs) # sort the truth and query files

            # basic queries should also be subsetted
            if self.options.datatype == "basic":
                self.subset_basic_queryfiles(querydirs, label, num)

            # create file for chrlens
            # fh_chromlens = open(outdir / f"{label}_chromlens.txt",'w')
            # for chr in self.chroms["leftgap"]:
            #     if self.chroms["leftgap"][chr] != {}:
            #         fh_chromlens.write(f"{chr}\t{self.chroms['leftgap'][chr]['end']}\n")
            #     else:
            #         fh_chromlens.write(f"{chr}\t0\n")
            # fh_chromlens.close()

    def sim_basic_queries(self, datafiles, intvl, rightgap):
        intvl_id = intvl["id"]
        chrom = intvl["chrom"]
        start_intvl = intvl["start"]
        end_intvl = intvl["end"]
        mid_intvl = intvl["mid"]

        # write reference and perfect query (is the same)
        datafiles["ref"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}\n")
        datafiles["queries-basic"]["perfect"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect\n")

        # write partial (5' and 3') overlaps
        start_partial_5p = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
        end_partial_5p = random.randint(start_intvl, mid_intvl)
        datafiles["queries-basic"]["5p-partial"].write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{intvl_id}_5p\n")
        start_partial_3p = random.randint(mid_intvl+1, end_intvl)
        end_partial_3p = random.randint(end_intvl+1, rightgap["mid"])
        datafiles["queries-basic"]["3p-partial"].write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{intvl_id}_3p\n")

        # write enclosing (containing) intervals (both with respect to query and reference)
        start_contained = random.randint(start_intvl, mid_intvl) # query enclosed by reference
        end_contained = random.randint(mid_intvl+1, end_intvl)
        datafiles["queries-basic"]["contained"].write(f"{chrom}\t{start_contained}\t{end_contained}\t{intvl_id}_contained\n")
        start_enclosed = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
        end_enclosed = random.randint(end_intvl+1, rightgap["mid"])
        datafiles["queries-basic"]["enclosed"].write(f"{chrom}\t{start_enclosed}\t{end_enclosed}\t{intvl_id}_enclosed\n")

        # add no overlaps (falls within gaps)
        datafiles["perfect-gap"].write(f"{chrom}\t{self.chroms['leftgap'][chrom]['start']}\t{self.chroms['leftgap'][chrom]['end']}\t{intvl_id}_perfect-gap\n")
        end_left_adjacent = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, self.chroms["leftgap"][chrom]["end"])
        datafiles["left-adjacent-gap"].write(f"{chrom}\t{self.chroms['leftgap'][chrom]['start']}\t{end_left_adjacent}\t{intvl_id}_left-adjacent\n")
        start_right_adjacent = random.randint(self.chroms["leftgap"][chrom]["start"], self.chroms["leftgap"][chrom]["mid"])
        datafiles["right-adjacent-gap"].write(f"{chrom}\t{start_right_adjacent}\t{self.chroms['leftgap'][chrom]['end']}\t{intvl_id}_right-adjacent\n")
        start_mid_gap1 = random.randint(self.chroms["leftgap"][chrom]["start"], self.chroms["leftgap"][chrom]["mid"])
        end_mid_gap1 = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, self.chroms["leftgap"][chrom]["end"])
        datafiles["mid-gap1"].write(f"{chrom}\t{start_mid_gap1}\t{end_mid_gap1}\t{intvl_id}_mid-gap1\n")
        start_mid_gap2 = random.randint(self.chroms["leftgap"][chrom]["start"], self.chroms["leftgap"][chrom]["mid"])
        end_mid_gap2 = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, self.chroms["leftgap"][chrom]["end"])
        datafiles["mid-gap2"].write(f"{chrom}\t{start_mid_gap2}\t{end_mid_gap2}\t{intvl_id}_mid-gap2\n")

        # also write entries to truth
        datafiles["truth"].write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_5p:{intvl_id}\n")
        datafiles["truth"].write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_3p:{intvl_id}\n")
        datafiles["truth"].write(f"{chrom}\t{start_contained}\t{end_contained}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_contained:{intvl_id}\n")
        datafiles["truth"].write(f"{chrom}\t{start_enclosed}\t{end_enclosed}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_enclosed:{intvl_id}\n")
        datafiles["truth"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect:{intvl_id}\n")

        self.chroms["leftgap"][chrom] = rightgap # make rightgap to leftgap (for next iteration)

    def simulate_complex(self):
        print()
