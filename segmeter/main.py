from pathlib import Path
import argparse
import os
import time

from simulator import SimBase

def main():
    options = parse_arguments()
    print(options)

    output = Path(options.outdir) / options.modus
    output.mkdir(parents=True, exist_ok=True)

    if options.modus == "sim":
        bed = SimBase(options)
        bed.format.simulate_basic()



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
    parser.add_argument("-n", "--intvlnums", type=str, help="number of intervals to simulate", default="10")
    parser.add_argument("-o", "--outdir", type=str, help="output folder for the benchmark results", required=True)

    args = parser.parse_args()

    # Get the directory of the script
    script_dir = Path(__file__).parent.resolve()

    return args


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
