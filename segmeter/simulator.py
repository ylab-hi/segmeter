import os
import random
from pathlib import Path

class SimBase:
    def __init__(self, options):
        self.options = options
        self.intvlnums = self.det_intvlnums()

        if options.format == "BED":
            self.format = SimBED(options, self.intvlnums)

    def det_intvlnums(self):
        intnums = {}
        nummap = {'K': 1000, 'M': 1000000}
        for num in self.options.intvlnums.split(","):
            if num[-1] in nummap:
                intnums[num] = int(num[:-1]) * nummap[num[-1]]
            else:
                intnums[num] = int(num)
        return intnums

class SimBED:
    def __init__(self, options, intvlnums):
        self.options = options
        self.intvlnums = intvlnums
        self.chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        # self.rightmost = self.det_rightmost_start()
        self.leftgap = self.init_leftgap()

    def init_leftgap(self):
        leftgap = {}
        for chr in self.chroms:
            leftgap[chr] = {}
            leftgap[chr]["start"] = 1
            leftgap[chr]["end"] = random.randint(100, 10000)
            leftgap[chr]["mid"] = int((leftgap[chr]["start"] + leftgap[chr]["end"]) // 2)
        return leftgap

    def det_rightmost_start(self): # needs to be set to 0
        """Simulate start position of first interval (aka right-most position)"""
        rightmost = {}
        for chr in self.chroms:
            rightmost[chr] = random.randint(100, 10000)
        return rightmost

    def simulate_ref(self):
        outdir = Path(self.options.outdir) / "sim" / "BED" # create base output folder
        outdir_ref = outdir / "ref"
        outdir_partial = outdir / "partial" # create query with partial 5' and 3' overlaps
        outdir_enclosed = outdir / "enclosed" # create query with enclosed intervals

        outdir_ref.mkdir(parents=True, exist_ok=True)
        outdir_partial.mkdir(parents=True, exist_ok=True)
        outdir_enclosed.mkdir(parents=True, exist_ok=True)

        for label, num in self.intvlnums.items():
            print(label)
            simdata_ref = outdir_ref / f"{label}.bed" # reference intervals
            simdata_partial = outdir_partial / f"{label}.bed"
            simdata_enclosed = outdir_enclosed / f"{label}.bed"

            # open files
            fh_ref = open(simdata_ref, 'w')
            fh_partial = open(simdata_partial, 'w')
            fh_enclosed = open(simdata_enclosed, 'w')
            for i in range(1, int(num)+1):
                chrom = random.choice(self.chroms)
                start_intvl = self.leftgap[chrom]+1
                intvl_size = random.randint(100,50000)
                end_intvl = start_intvl + (intvl_size-1)
                mid_intvl = int((start_intvl + end_intvl) // 2) # determine middle of interval

                rightgap = (end_intvl+1, end_intvl+1+random.randint(100,10000))
                mid_rightgap = int((rightgap[0] + rightgap[1]) // 2)

                # write reference
                fh_ref.write(f"{chrom}\t{start}\t{end}\n")

                # write partial (5' and 3' overlaps)
                start_partial_5p = random.randint(self.leftgap[chrom]["mid"]+1, start_intvl-1)
                end_partial_5p = random.randint(start_intvl, mid_intvl)
                start_partial_3p = random.randint(mid_intvl+1, end_intvl)
                end_partial_3p = random.randint(end_intvl+1, mid_rightgap)



                self.rightmost[chrom] = end
            fh_ref.close()
            simdatafile_sorted = outdir_ref / f"{label}_sorted.bed"
            os.system(f"sort -k1,1 -k2,2n -k3,3n  {simdata_ref} > {simdatafile_sorted}")

    def simulate_intvldist(self):
        """simulate the interval distance"""
        intvldist = random.randint(100, 10000)
        return intvldist

    # def simulate_queries(self):
    #     outdir_ref = Path(self.options.outdir) / "sim" /
    #     for label, num in self.intvlnums.items():
    #         ref_file =
    #         fh_ref = open()
