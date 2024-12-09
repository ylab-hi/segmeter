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
        datadirs = {}
        datadirs["ref"] = outdir / "ref"
        datadirs["truth"] = outdir / "truth"

        if self.options.datatype == "basic":
            # positive queries
            queries_pos = ["perfect", "5p-partial", "3p-partial", "enclosed", "contained"]
            for query in queries_pos:
                datadirs[query] = outdir / "query" / query

            # negative queries
            queries_neg = ["perfect-gap", "left-adjacent-gap", "right-adjacent-gap", "mid-gap1", "mid-gap2"]
            for query in queries_neg:
                datadirs[query] = outdir / "query" / query

        elif self.options.datatype == "complex":
            # positive queries
            queries_pos = ["span-intvl"]
            for query in queries_pos:
                datadirs[query] = outdir / "query" / query

            queries_neg = ["span-gap"]
            for query in queries_neg:
                datadirs[query] = outdir / "query" / query

        for key in datadirs.keys():
            datadirs[key].mkdir(parents=True, exist_ok=True)

        return datadirs

    def open_datafiles(self, datadirs, label):
        datafiles = {}
        for key in datadirs.keys():
            datafiles[key] = open(datadirs[key] / f"{label}.bed", 'w')
        return datafiles

    def close_datafiles(self, datafiles):
        for key in datafiles.keys():
            datafiles[key].close()

    def sort_datafiles(self, datadirs, label):
        for key in datadirs.keys():
            infile = datadirs[key] / f"{label}.bed"
            sorted = datadirs[key] / f"{label}_sorted.bed"
            utility.sort_BED(infile, sorted)

    def subset_datafiles(self, datadirs, label, num):
        for key in datadirs.keys():
            if (key == "ref" or key == "truth"):
                continue
            infile = datadirs[key] / f"{label}.bed"
            for subset in range(10,101,10):
                outfile = datadirs[key] / f"{label}_{subset}p.bed"
                with open(outfile, 'w') as subset_bed:
                    subprocess.run(["shuf", "-n", str(int(num * (subset / 100))), str(infile)], stdout=subset_bed)
                sorted = datadirs[key] / f"{label}_{subset}p_sorted.bed"
                utility.sort_BED(outfile, sorted) # also sort the subset files

    def simulate(self):
        outdir = Path(self.options.outdir) / "sim" / "BED" / self.options.datatype # create base output folder
        datadirs = self.create_datadirs(outdir)
        is_start, is_end = [int(x) for x in self.options.intvlsize.split("-")]

        for label, num in self.intvlnums.items():
            print(f"Simulate intervals for {label}:{num}...")
            datafiles = self.open_datafiles(datadirs, label)

            for i in range(1, (int(num)//10)+1):
                intvl_id = f"intvl_{i}" # create an interval ID

                # chrom = random.choice(self.chroms)
                chrom = self.select_chrom()
                start_intvl = self.chroms["leftgap"][chrom]["end"]+1
                intvl_size = random.randint(is_start, is_end)
                end_intvl = start_intvl + (intvl_size-1)
                mid_intvl = int((start_intvl + end_intvl) // 2) # determine middle of interval

                gs_start, gs_end = [int(x) for x in self.options.gapsize.split("-")]
                rightgap = self.simulate_gap(end_intvl+1, end_intvl+1+random.randint(gs_start,gs_end))

                # write reference and perfect query (is the same)
                datafiles["ref"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}\n")
                datafiles["perfect"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect\n")

                # write partial (5' and 3') overlaps
                start_partial_5p = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
                end_partial_5p = random.randint(start_intvl, mid_intvl)
                datafiles["5p-partial"].write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{intvl_id}_5p\n")
                start_partial_3p = random.randint(mid_intvl+1, end_intvl)
                end_partial_3p = random.randint(end_intvl+1, rightgap["mid"])
                datafiles["3p-partial"].write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{intvl_id}_3p\n")

                # write enclosing (containing) intervals (both with respect to query and reference)
                start_contained = random.randint(start_intvl, mid_intvl) # query enclosed by reference
                end_contained = random.randint(mid_intvl+1, end_intvl)
                datafiles["contained"].write(f"{chrom}\t{start_contained}\t{end_contained}\t{intvl_id}_contained\n")
                start_enclosed = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
                end_enclosed = random.randint(end_intvl+1, rightgap["mid"])
                datafiles["enclosed"].write(f"{chrom}\t{start_enclosed}\t{end_enclosed}\t{intvl_id}_enclosed\n")

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

            self.close_datafiles(datafiles)
            self.sort_datafiles(datadirs, label) # sort the ref and query files
            self.subset_datafiles(datadirs, label, num)

            # create file for chrlens
            fh_chromlens = open(outdir / f"{label}_chromlens.txt",'w')
            for chr in self.chroms:
                fh_chromlens.write(f"{chr}\t{self.chroms['leftgap'][chr]['end']}\n")
            fh_chromlens.close()

    def simulate_basic(self, datafiles, num):
        for i in range(1, (int(num)//10)+1):
            intvl_id = f"intvl_{i}" # create an interval ID

            # chrom = random.choice(self.chroms)
            chrom = self.select_chrom()
            start_intvl = self.chroms["leftgap"][chrom]["end"]+1
            intvl_size = random.randint(is_start, is_end)
            end_intvl = start_intvl + (intvl_size-1)
            mid_intvl = int((start_intvl + end_intvl) // 2) # determine middle of interval

            gs_start, gs_end = [int(x) for x in self.options.gapsize.split("-")]
            rightgap = self.simulate_gap(end_intvl+1, end_intvl+1+random.randint(gs_start,gs_end))

            # write reference and perfect query (is the same)
            datafiles["ref"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}\n")
            datafiles["perfect"].write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_perfect\n")

            # write partial (5' and 3') overlaps
            start_partial_5p = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
            end_partial_5p = random.randint(start_intvl, mid_intvl)
            datafiles["5p-partial"].write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{intvl_id}_5p\n")
            start_partial_3p = random.randint(mid_intvl+1, end_intvl)
            end_partial_3p = random.randint(end_intvl+1, rightgap["mid"])
            datafiles["3p-partial"].write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{intvl_id}_3p\n")

            # write enclosing (containing) intervals (both with respect to query and reference)
            start_contained = random.randint(start_intvl, mid_intvl) # query enclosed by reference
            end_contained = random.randint(mid_intvl+1, end_intvl)
            datafiles["contained"].write(f"{chrom}\t{start_contained}\t{end_contained}\t{intvl_id}_contained\n")
            start_enclosed = random.randint(self.chroms["leftgap"][chrom]["mid"]+1, start_intvl-1)
            end_enclosed = random.randint(end_intvl+1, rightgap["mid"])
            datafiles["enclosed"].write(f"{chrom}\t{start_enclosed}\t{end_enclosed}\t{intvl_id}_enclosed\n")

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
