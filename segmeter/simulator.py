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
        self.rightmost = self.det_rightmost_start()

    def det_rightmost_start(self): # needs to be set to 0
        """Simulate start position of first interval (aka right-most position)"""
        rightmost = {}
        for chr in self.chroms:
            rightmost[chr] = random.randint(100, 10000)
        return rightmost

    def simulate_basic(self):
        outdir = Path(self.options.outdir) / "sim" / "BED"
        outdir_ref = outdir / "ref"
        outdir_ref.mkdir(parents=True, exist_ok=True)

        for label, num in self.intvlnums.items():
            simdata_ref = outdir_ref / f"{label}.bed"
            fh_ref = open(simdata_ref, 'w')
            for i in range(1, int(num)+1):
                chrom = random.choice(self.chroms)
                # determine space between last interval (rightmost) - no overlaps
                intvldist = random.randint(100, 10000)
                start = self.rightmost[chrom] + intvldist
                intvl_size = random.randint(100, 50000)
                end = start + intvl_size

                self.rightmost[chrom] = end
                fh_ref.write(f"{chrom}\t{start}\t{end}\n")
            fh_ref.close()
            simdatafile_sorted = outdir_ref / f"{label}_sorted.bed"
            os.system(f"sort -k1,1 -k2,2n -k3,3n  {simdata_ref} > {simdatafile_sorted}")
