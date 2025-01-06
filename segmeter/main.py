from pathlib import Path
import argparse
import os
import time
import logging

from simulator import SimBase
from benchmark import BenchBase
import utility

def main():
    options = parse_arguments()

    output = Path(options.outdir) / options.modus
    output.mkdir(parents=True, exist_ok=True)
    intvlnums = det_intvlnums(options.intvlnums)

    if options.modus == "sim":
        bed = SimBase(options, intvlnums)
        bed.format.sim_intervals()
    elif options.modus == "bench":
        bench = BenchBase(options, intvlnums)

    # look up the files to benchmark
    # fh = open(options.list, "r")
    # for line in fh:
    #     inputfile = line.strip()

        # if options.tool == "tabix":
        #     tabix_call(inputfile, options.output)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Benchmarking tool for interval files")
    parser.add_argument("modus", type=str, help="modus in benchmarking", choices=["sim", "bench"])
    parser.add_argument("-f", "--format", type=str, help="format of the files to benchmark", choices=["BED"], default="BED")
    parser.add_argument("-n", "--intvlnums", type=str, help="""Number of intervals to simulate (should be divisible by 10).""", default="10")
    parser.add_argument("-c", "--clustfrac", type=str, help="fraction of the intervals (as specified in --intvlsnums) that are part of a cluster", default="0.1-0.2")
    parser.add_argument("-m", "--max_chromlen", type=int, help="maximum length of the simulated chromosomes", default=1000000000)
    parser.add_argument("-o", "--outdir", type=str, help="output folder for the benchmark/simulation results. ", required=True)
    parser.add_argument("-t", "--tool", type=str, help="tool to benchmark", choices=["tabix", "bedtools"])
    parser.add_argument("-g", "--gapsize", type=str, help="random size of the gaps (min and max) between the intervals", default="100-5000")
    parser.add_argument("-i", "--intvlsize", type=str, help="random size (min and max) of the intervals", default="100-10000")
    # parser.add_argument("-l", "--loglevel", type=str, help="Log level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO")
    # parser.add_argument("-d", "--datatype", type=str, help="datatype of the intervals", choices=["basic", "complex", "custom"], default="basic")

    args = parser.parse_args()

    # Get the directory of the script
    script_dir = Path(__file__).parent.resolve()

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


# def tabix_call(inputfile, outputdir):
#     startIndex = time.time() # start the time measurement
#     output = Path(outputdir) / "tabix"
#      # create output file
#     if not os.path.exists(output):
#         os.mkdir(output)
#     outputfile = output / f"{Path(inputfile).stem}.gz"
#     print(outputfile)
#     os.system(f"bgzip -c {inputfile} > {outputfile} && tabix -p bed {outputfile}")


main()
