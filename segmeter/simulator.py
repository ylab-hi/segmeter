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
            leftgap[chr] = self.simulate_gap(1, random.randint(100, 10000))
        return leftgap

    def simulate_gap(self, start, end):
        gap = {}
        gap["start"] = start
        gap["end"] = end
        gap["mid"] = int((start + end) // 2)
        return gap

    def update_leftgap(self, chrom, start, end, mid):
        self.leftgap[chrom]["start"] = start
        self.leftgap[chrom]["end"] = end
        self.leftgap[chrom]["mid"] = mid

    def det_rightmost_start(self): # needs to be set to 0
        """Simulate start position of first interval (aka right-most position)"""
        rightmost = {}
        for chr in self.chroms:
            rightmost[chr] = random.randint(100, 10000)
        return rightmost

    def simulate(self):
        outdir = Path(self.options.outdir) / "sim" / "BED" / "simple" # create base output folder
        outdir_ref = outdir / "ref"
        outdir_partial = outdir / "partial" # create query with partial 5' and 3' overlaps
        outdir_enclosed = outdir / "enclosed" # create query with enclosed intervals
        outdir_truth = outdir / "truth" # create truth file

        outdir_ref.mkdir(parents=True, exist_ok=True)
        outdir_partial.mkdir(parents=True, exist_ok=True)
        outdir_enclosed.mkdir(parents=True, exist_ok=True)
        outdir_truth.mkdir(parents=True, exist_ok=True)

        for label, num in self.intvlnums.items():
            print(f"Simulate intervals for {label}:{num}...")
            simdata_ref = outdir_ref / f"{label}.bed" # reference intervals
            simdata_partial = outdir_partial / f"{label}.bed"
            simdata_enclosed = outdir_enclosed / f"{label}.bed"
            simdata_truth = outdir_truth / f"{label}.bed"

            # open files
            fh_ref = open(simdata_ref, 'w')
            fh_partial = open(simdata_partial, 'w')
            fh_enclosed = open(simdata_enclosed, 'w')
            fh_truth = open(simdata_truth, 'w')

            for i in range(1, int(num)+1):
                intvl_id = f"intvl_{i}"

                chrom = random.choice(self.chroms)
                start_intvl = self.leftgap[chrom]["end"]+1
                intvl_size = random.randint(100,50000)
                end_intvl = start_intvl + (intvl_size-1)
                mid_intvl = int((start_intvl + end_intvl) // 2) # determine middle of interval

                rightgap = self.simulate_gap(end_intvl+1, end_intvl+1+random.randint(100,10000))

                # write reference
                fh_ref.write(f"{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}\n")

                # write partial (5' and 3' overlaps)
                start_partial_5p = random.randint(self.leftgap[chrom]["mid"]+1, start_intvl-1)
                end_partial_5p = random.randint(start_intvl, mid_intvl)
                start_partial_3p = random.randint(mid_intvl+1, end_intvl)
                end_partial_3p = random.randint(end_intvl+1, rightgap["mid"])
                fh_partial.write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{intvl_id}_5p\n")
                fh_partial.write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{intvl_id}_3p\n")

                fh_truth.write(f"{chrom}\t{start_partial_5p}\t{end_partial_5p}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_5p:{intvl_id}\n")
                fh_truth.write(f"{chrom}\t{start_partial_3p}\t{end_partial_3p}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_3p:{intvl_id}\n")

                # write enclosing intervals (both with respect to query and reference)
                start_enclosed_query = random.randint(start_intvl, mid_intvl) # query enclosed by reference
                end_enclosed_query = random.randint(mid_intvl+1, end_intvl)
                fh_enclosed.write(f"{chrom}\t{start_enclosed_query}\t{end_enclosed_query}\t{intvl_id}_contained\n")

                start_enclosed_ref = random.randint(self.leftgap[chrom]["mid"]+1, start_intvl-1)
                end_enclosed_ref = random.randint(end_intvl+1, rightgap["mid"])
                fh_enclosed.write(f"{chrom}\t{start_enclosed_ref}\t{end_enclosed_ref}\t{intvl_id}_enclosed\n")

                fh_truth.write(f"{chrom}\t{start_enclosed_query}\t{end_enclosed_query}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_contained:{intvl_id}\n")
                fh_truth.write(f"{chrom}\t{start_enclosed_ref}\t{end_enclosed_ref}\t{chrom}\t{start_intvl}\t{end_intvl}\t{intvl_id}_enclosed:{intvl_id}\n")

                self.leftgap[chrom] = rightgap # make rightgap to leftgap (for next iteration)


            fh_ref.close()
            fh_partial.close()
            fh_enclosed.close()
            fh_truth.close()

            # simdatafile_sorted = outdir_ref / f"{label}_sorted.bed"
            os.system(f"sort -k1,1 -k2,2n -k3,3n  {simdata_ref} > {outdir_ref / f'{label}_sorted.bed'}")
            os.system(f"sort -k1,1 -k2,2n -k3,3n  {simdata_partial} > {outdir_partial / f'{label}_sorted.bed'}")
            os.system(f"sort -k1,1 -k2,2n -k3,3n  {simdata_enclosed} > {outdir_enclosed / f'{label}_sorted.bed'}")
