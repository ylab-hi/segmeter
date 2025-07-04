# standard
from pathlib import Path
import argparse

# class
from simulator import SimBase
from benchmark import BenchBase

def main():
    options = parse_arguments()

    output = Path(options.datadir) / options.modus
    output.mkdir(parents=True, exist_ok=True)
    intvlnums = det_intvlnums(options.intvlnums)

    if options.modus == "sim":
        bed = SimBase(options, intvlnums)
        bed.format.sim_intervals()
    elif options.modus == "bench":
        bench = BenchBase(options, intvlnums)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Benchmarking tool for interval files")
    parser.add_argument("modus", type=str, help="modus in benchmarking", choices=["sim", "bench"])
    parser.add_argument("-r", "--simdata", action="store_true", help="use simulated data for benchmarking. Note that this will use the simulated data (mode sim) as input for the benchmarking")
    parser.add_argument("-q", "--query", type=str, help="query file used for benchmarking (not used when benchmarking simulated data)")
    parser.add_argument("--target", type=str, help="target file used for benchmarking (not used when benchmarking simulated data)")
    parser.add_argument("-o", "--datadir", type=str, help="output folder for the benchmark/simulation results. Note this also serves as input folder for the benchmarking", required=True)
    parser.add_argument("-f", "--format", type=str, help="format of the files to benchmark", choices=["BED"], default="BED")
    parser.add_argument("-n", "--intvlnums", type=str, help="""Number of intervals to simulate (should be divisible by 10).""", default="10")
    parser.add_argument("-s", "--subset", type=str, help="subset of the intervals to use for benchmarking. Format should be either XX-YY or XX,YY-ZZ", default="10-100")
    parser.add_argument("-m", "--max_chromlen", type=int, help="maximum length of the simulated chromosomes", default=1000000000)
    parser.add_argument("-b", "--benchname", type=str, help="name of the benchmark, used for the output folder", default="bench_001")
    parser.add_argument("-c", "--simname", type=str, help="name of the simulation, used for the output folder. Only used in simulation mode", default="sim_001")
    parser.add_argument("-t", "--tool", type=str, help="tool to benchmark",
        choices=["tabix", "bedtools", "bedtools_sorted", "bedtools_tabix", "bedops", "bedmaps", "giggle", "granges", "gia", "bedtk", "bedtk_sorted", "igd", "ailist",
            "ucsc", "awk", "intervaltree"])
    parser.add_argument("-g", "--gapsize", type=str, help="random size of the gaps (min and max) between the intervals", default="100-5000")
    parser.add_argument("-i", "--intvlsize", type=str, help="random size (min and max) of the intervals", default="100-10000")

    args = parser.parse_args()

    return args


def det_intvlnums(intvlnums):
    intnums = {}
    nummap = {'K': 1000, 'M': 1000000}
    if "," in intvlnums:
        intvlnums = intvlnums.split(",")
        for num in intvlnums:
            if num[-1] in nummap:
                intnums[num] = int(num[:-1]) * nummap[num[-1]]
            else:
                intnums[num] = int(num)
    else:
        if intvlnums[-1] in nummap:
            intnums[intvlnums] = int(intvlnums[:-1]) * nummap[intvlnums[-1]]
        else:
            intnums[intvlnums] = int(intvlnums)
    return intnums


main()
