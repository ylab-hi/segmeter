import platform
import random
from pathlib import Path

# class
import utility

class SimBED:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums
        self.gs_start, self.gs_end = [int(x) for x in options.gapsize.split("-")] # min and max gap size

    def init_chroms(self):
        """Initialize the chromosomes"""
        chroms = {}
        # consists of a list of chromosome that haven't exeeded a maximum length
        chroms["all"] = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'] # main chromosomes
        firstchrom = random.choice(chroms["all"]) # randomly select a chromosome
        chroms["space-left"] = [firstchrom] # start with a single chromosome (e.g. chr1)
        chroms["all"].remove(firstchrom) # remove the chromosome from the list (e.g. chr1)
        chroms["leftgap"] = {} # contains the end position of the last gap (left of interval)
        chroms["intvl"] = {} # stores the number of intervals on each chromosome
        for chrom in chroms["space-left"]:
            chroms["leftgap"][chrom] = {}
            chroms["intvl"][chrom] = 0
        return chroms

    def select_chrom(self, chroms):
        """Randomly select a chromosome"""
        # check if at least one has less than 10 intervals
        if len(chroms["all"]) > 0: # only if there are (main) chromosomes left
            if all(num >= 10 for num in chroms["intvl"].values()): # no chromosome left with less than 10 intervals
                nextchrom = random.choice(chroms["all"]) # randomly select a chromosome
                chroms["all"].remove(nextchrom) # remove the chromosome from the list
                chroms["space-left"].append(nextchrom) # add a new chromosome
                chroms["leftgap"][nextchrom] = {} # reset the gap
                chroms["intvl"][nextchrom] = 0 # reset the interval counter

        selection = random.choice(chroms["space-left"])
        if chroms["leftgap"][selection] == {}: # no interval has been placed on this chromosome (yet)
            # therefore simulate a gap
            chroms["leftgap"][selection] = self.simulate_gap(1, random.randint(self.gs_start, self.gs_end))
            return selection
        elif chroms["leftgap"][selection]["end"] < self.options.max_chromlen:
            return selection
        else: # chromosome is full, continue on a new scaffold
            chroms["space-left"].remove(selection)
            scaffolds = [chrom for chrom in chroms["leftgap"] if "SCF" in chrom]
            scaffold_name = f"SCF{len(scaffolds)+1}"
            chroms["space-left"].append(scaffold_name)
            chroms["leftgap"][scaffold_name] = self.simulate_gap(1, random.randint(self.gs_start, self.gs_end))
            chroms["intvl"][scaffold_name] = 0
            return scaffold_name

    def simulate_gap(self, start, end):
        gap = {}
        gap["start"] = start
        gap["end"] = end
        gap["mid"] = int((start + end) // 2)
        return gap

    def update_intvl_counter(self, chroms, chrom):
        chroms["intvl"][chrom] += 1 # every chromosome is registered with 0 by select_chrom

    def create_datadirs(self, datadir):
        refdir = datadir / "ref"

        truthdirs = {}
        truthdirs["basic"] = datadir / "basic" / "truth"
        truthdirs["complex"] = datadir / "complex" / "truth"

        querydirs = {}
        querydirs["basic"] = {}
        for query in utility.BASIC_QUERIES:
            querydirs["basic"][query] = datadir / "basic" / "query" / query
        querydirs["complex"] = {}
        querydirs["complex"]["mult"] = datadir / "complex" / "query" / "mult"

        # create folder
        refdir.mkdir(parents=True, exist_ok=True)
        truthdirs["basic"].mkdir(parents=True, exist_ok=True)
        truthdirs["complex"].mkdir(parents=True, exist_ok=True)
        for querydir in querydirs["basic"].values():
            querydir.mkdir(parents=True, exist_ok=True)
        querydirs["complex"]["mult"].mkdir(parents=True, exist_ok=True)

        return refdir, truthdirs, querydirs

    def open_datafiles(self, label, refdir, truthdirs, querydirs):
        datafiles = {}
        datafiles["ref"] = open(refdir / f"{label}.bed", 'w')
        datafiles["truth-basic"] = open(truthdirs["basic"] / f"{label}.bed", 'w')
        datafiles["queries-basic"] = {}
        for key in querydirs["basic"].keys():
            datafiles["queries-basic"][key] = open(querydirs["basic"][key] / f"{label}.bed", 'w')

        return datafiles

    def close_datafiles_basic(self, datafiles):
        """Close all (basic) datafiles"""
        datafiles["truth-basic"].close()
        for key in datafiles["queries-basic"].keys():
            datafiles["queries-basic"][key].close()

    def sort_datafiles(self, label, truthdirs, querydirs):
        """Sort the basic truth and query files"""
        utility.sort_BED(truthdirs["basic"] / f"{label}.bed", truthdirs["basic"] / f"{label}_sorted.bed")
        for querydir in querydirs["basic"].values():
            utility.sort_BED(querydir / f"{label}.bed", querydir / f"{label}_sorted.bed")

    def subset_basic_queryfiles(self, querydirs, label, num):
        """Random samples of 10%, 20%, ..., 100% of the queries of every basic type, in random order (as `shuf -n` gave them before #15)"""
        for key in querydirs["basic"].keys():
            with open(querydirs["basic"][key] / f"{label}.bed") as fh:
                lines = fh.readlines() # one query per interval
            for subset in range(10,101,10):
                outfile = querydirs["basic"][key] / f"{label}_{subset}p.bed"
                with open(outfile, 'w') as subset_bed:
                    subset_bed.writelines(random.sample(lines, int(num * (subset / 100))))
                utility.sort_BED(outfile, querydirs["basic"][key] / f"{label}_{subset}p_sorted.bed")

    def sim_intervals(self):
        outpath = Path(self.options.datadir) / "sim" / self.options.simname / "BED"
        refdir, truthdirs, querydirs = self.create_datadirs(outpath)
        is_start, is_end = [int(x) for x in self.options.intvlsize.split("-")]

        # seed once, so that the same seed and parameters give the same data; the seed is recorded for unseeded runs (#15)
        seed = self.options.seed if self.options.seed is not None else random.randrange(2**32)
        random.seed(seed)
        with open(outpath / "parameters.txt", "w") as fh:
            fh.write(f"seed\t{seed}\npython\t{platform.python_version()}\n") # the draws of `random` are an implementation detail of CPython
            for key in ("intvlnums", "intvlsize", "gapsize", "max_chromlen", "max_span"):
                fh.write(f"{key}\t{getattr(self.options, key)}\n")

        for label, num in self.intvlnums.items():
            print(f"Simulate intervals for {label}:{num}...")
            datafiles = self.open_datafiles(label, refdir, truthdirs, querydirs)

            chroms = self.init_chroms()

            for i in range(1, (int(num))+1):
                chrom = self.select_chrom(chroms)
                intvl = {} # create/simulate new interval
                intvl["id"] = f"intvl_{i}" # create an interval ID
                intvl["chrom"] = chrom
                intvl["start"] = chroms["leftgap"][chrom]["end"]+1
                intvl_size = random.randint(is_start, is_end)
                intvl["end"] = intvl["start"] + (intvl_size-1)
                intvl["mid"] = int((intvl["start"] + intvl["end"]) // 2) # determine middle of interval

                # update counter
                self.update_intvl_counter(chroms, chrom)

                rightgap = self.simulate_gap(intvl["end"]+1, intvl["end"]+1+random.randint(self.gs_start, self.gs_end))

                self.sim_basic_queries(chroms, datafiles, intvl, rightgap)

            datafiles["ref"].close() # close the reference file
            self.close_datafiles_basic(datafiles) # close the basic datafiles
            utility.sort_BED(refdir / f"{label}.bed", refdir / f"{label}_sorted.bed") # sort the reference file
            self.sort_datafiles(label, truthdirs, querydirs) # sort the truth and query files

            # basic queries should also be subsetted
            self.subset_basic_queryfiles(querydirs, label, num)

            # create file for chrlens
            fh_chromlens = open(outpath / f"{label}_chromlens.txt",'w')
            for chrom in chroms["leftgap"]:
                if chroms["leftgap"][chrom] != {}:
                    fh_chromlens.write(f"{chrom}\t{chroms['leftgap'][chrom]['end']}\n")
                else:
                    fh_chromlens.write(f"{chrom}\t0\n")
            fh_chromlens.close()

            # create file for chrnums (number of intervals on each chromosome)
            fh_chrnums = open(outpath / f"{label}_chrnums.txt",'w')
            for chrom in chroms["intvl"]:
                fh_chrnums.write(f"{chrom}\t{chroms['intvl'][chrom]}\n")
            fh_chrnums.close()

            self.sim_complex_queries(refdir, truthdirs, querydirs, num, label)

    def sim_basic_queries(self, chroms, datafiles, intvl, rightgap):
        intvl_id = intvl["id"]
        chrom = intvl["chrom"]
        start_intvl = intvl["start"]
        end_intvl = intvl["end"]
        mid_intvl = intvl["mid"]

        # write reference and perfect query (is the same)
        datafiles["ref"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}\n")
        datafiles["queries-basic"]["perfect"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect\n")

        # write partial (5' and 3') overlaps
        start_partial_5p = random.randint(chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
        end_partial_5p = random.randint(start_intvl+1, mid_intvl)
        datafiles["queries-basic"]["5p-partial"].write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{intvl_id}_5p\n")
        start_partial_3p = random.randint(mid_intvl+1, end_intvl-1)
        end_partial_3p = random.randint(end_intvl+1, rightgap["mid"])
        datafiles["queries-basic"]["3p-partial"].write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{intvl_id}_3p\n")

        # write enclosing (containing) intervals (both with respect to query and reference)
        start_contained = random.randint(start_intvl, mid_intvl) # query enclosed by reference
        end_contained = random.randint(mid_intvl+1, end_intvl)
        datafiles["queries-basic"]["contained"].write(f"{chrom}\t{start_contained}\t{end_contained}\t{intvl_id}_contained\n")
        start_enclosed = random.randint(chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
        end_enclosed = random.randint(end_intvl+1, rightgap["mid"])
        datafiles["queries-basic"]["enclosed"].write(f"{chrom}\t{start_enclosed}\t{end_enclosed}\t{intvl_id}_enclosed\n")

        # add no overlaps (falls within gaps)
        datafiles["queries-basic"]["perfect-gap"].write(f"{chrom}\t{chroms['leftgap'][chrom]['start']}\t{chroms['leftgap'][chrom]['end']}\t{intvl_id}_perfect-gap\n")
        end_left_adjacent = random.randint(chroms["leftgap"][chrom]["mid"]+1, chroms["leftgap"][chrom]["end"])
        datafiles["queries-basic"]["left-adjacent-gap"].write(f"{chrom}\t{chroms['leftgap'][chrom]['start']}\t{end_left_adjacent}\t{intvl_id}_left-adjacent\n")
        start_right_adjacent = random.randint(chroms["leftgap"][chrom]["start"], chroms["leftgap"][chrom]["mid"])
        datafiles["queries-basic"]["right-adjacent-gap"].write(f"{chrom}\t{start_right_adjacent}\t{chroms['leftgap'][chrom]['end']}\t{intvl_id}_right-adjacent\n")
        start_mid_gap1 = random.randint(chroms["leftgap"][chrom]["start"], chroms["leftgap"][chrom]["mid"])
        end_mid_gap1 = random.randint(chroms["leftgap"][chrom]["mid"]+1, chroms["leftgap"][chrom]["end"])
        datafiles["queries-basic"]["mid-gap1"].write(f"{chrom}\t{start_mid_gap1}\t{end_mid_gap1}\t{intvl_id}_mid-gap1\n")
        start_mid_gap2 = random.randint(chroms["leftgap"][chrom]["start"], chroms["leftgap"][chrom]["mid"])
        end_mid_gap2 = random.randint(chroms["leftgap"][chrom]["mid"]+1, chroms["leftgap"][chrom]["end"])
        datafiles["queries-basic"]["mid-gap2"].write(f"{chrom}\t{start_mid_gap2}\t{end_mid_gap2}\t{intvl_id}_mid-gap2\n")

        # also write entries to truth
        datafiles["truth-basic"].write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_5p:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_3p:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_contained}\t{end_contained}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_contained:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_enclosed}\t{end_enclosed}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_enclosed:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect:{intvl_id}\n")

        datafiles["truth-basic"].write(f"{chrom}\t{chroms['leftgap'][chrom]['start']}\t{chroms['leftgap'][chrom]['end']}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect-gap:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{chroms['leftgap'][chrom]['start']}\t{end_left_adjacent}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_left-adjacent:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_right_adjacent}\t{chroms['leftgap'][chrom]['end']}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_right-adjacent:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_mid_gap1}\t{end_mid_gap1}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_mid-gap1:{intvl_id}\n")
        datafiles["truth-basic"].write(f"{chrom}\t{start_mid_gap2}\t{end_mid_gap2}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_mid-gap2:{intvl_id}\n")

        chroms["leftgap"][chrom] = rightgap # make rightgap to leftgap (for next iteration)

    def sim_complex_queries(self, refdir, truthdirs, querydirs, num, label):
        reffile = refdir / f"{label}_sorted.bed" # reference file
        truthfile = truthdirs["complex"] / f"{label}.bed"

        outpath = Path(self.options.datadir) / "sim" / self.options.simname / "BED"
        with open(outpath / f"{label}_chrnums.txt") as fh:
            maxchrms = max((int(line.split("\t")[1]) for line in fh), default=0) # maximum number of intervals on a chromosome
        if self.options.max_span: # queries do not cover more than max_span intervals
            maxchrms = min(maxchrms, self.options.max_span)
        bins = {}
        frac10 = int(maxchrms * 0.1)
        for i in range(1,11):
            start = frac10*(i-1)+1
            end = frac10*i
            bins[(start, end)] = open(querydirs["complex"]["mult"] / f"{label}_{i*10}bin.bed", 'w')

        fh_truth = open(truthfile, 'w')

        fh = open(reffile, 'r')
        curr_chr = ""
        intvls = []
        for line in fh:
            splitted = line.strip().split("\t")
            chrom = splitted[0]

            if curr_chr != "" and curr_chr != chrom:
                self.sim_overlaps(intvls, bins, fh_truth)
                intvls = []

            intvls.append(splitted)
            curr_chr = chrom

        if len(intvls) > 1: # simulate if still intervals left
            self.sim_overlaps(intvls, bins, fh_truth)

        fh.close()
        for i in bins.keys():
            bins[i].close()
        fh_truth.close()

    def sim_overlaps(self, intvls, bins, fh_truth):
        intvlnum = len(intvls)
        if intvlnum > 1:
            chrom = intvls[0][0]
            for i in range(2, min(intvlnum, self.options.max_span or intvlnum)+1):
                start_intvl = random.randint(0, intvlnum-i)
                end_intvl = start_intvl + i - 1

                query_start = intvls[start_intvl][1]
                query_end = intvls[end_intvl][2]

                for bnds in bins.keys():
                    if i >= bnds[0] and i <= bnds[1]:
                        bins[bnds].write(f"{chrom}\t{query_start}\t{query_end}\tmult_{i}\n")
                        break

                fh_truth.write(f"{chrom}\t{query_start}\t{query_end}\tmult_{i}\t{i}\n")
