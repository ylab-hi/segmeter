import argparse
import os
import time
from pathlib import Path
import simulate

def main():
    options = parse_arguments()
    print(options)

    if(options.modus == "BED"):
        bed = simulate.BED(options)

    # create folder
    if not os.path.exists(options.output):
        os.mkdir(options.output)

    # look up the files to benchmark
    fh = open(options.list, "r")
    for line in fh:
        inputfile = line.strip()

        if options.tool == "tabix":
            tabix_call(inputfile, options.output)



def parse_arguments():
    parser = argparse.ArgumentParser(description="Benchmarking tool for interval files")
    parser.add_argument("modus", type=str, help="modus in benchmarking", choices=["sim", "bench"])
    parser.add_argument("-f", "--format", type=str, help="format of the files to benchmark", choices=["BED"])
    # parser.add_argument("-l", "--list", type=str, help="list containing the files to benchmark")
    parser.add_argument("-o", "--outputdir", type=str, help="output folder for the benchmark results")

    args = parser.parse_args()
    return args


def tabix_call(inputfile, outputdir):
    startIndex = time.time() # start the time measurement
    output = Path(outputdir) / "tabix"
     # create output file
    if not os.path.exists(output):
        os.mkdir(output)
    outputfile = output / f"{Path(inputfile).stem}.gz"
    print(outputfile)
    os.system(f"bgzip -c {inputfile} > {outputfile} && tabix -p bed {outputfile}")


main()
