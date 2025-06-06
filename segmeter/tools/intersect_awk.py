#!/usr/bin/env python3

import argparse
import subprocess
import shlex
import sys
import os

def main():
    options = parse_arguments()

    # Validate input files
    if not os.path.isfile(options.target):
        print(f"Error: Target file '{options.target}' not found", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(options.query):
        print(f"Error: Query file '{options.query}' not found", file=sys.stderr)
        sys.exit(1)

    # AWK script that mimics bedtools intersect -wa behavior
    # First pass: read target intervals and store by chromosome
    # Second pass: for each query interval, check for overlaps and print if found
    awk_script = r'''
    BEGIN {
        FS = OFS = "\t"
    }
    NR==FNR {
        # First file (target): store intervals by chromosome
        chr = $1
        start = $2
        end = $3
        key = start "," end
        if (!(chr in intervals)) {
            intervals[chr] = ""
        }
        if (intervals[chr] == "") {
            intervals[chr] = key
        } else {
            intervals[chr] = intervals[chr] ";" key
        }
        next
    }
    {
        # Second file (query): check for overlaps
        chr = $1
        qstart = $2
        qend = $3

        if (chr in intervals) {
            # Split stored intervals for this chromosome
            n = split(intervals[chr], interval_list, ";")
            for (i = 1; i <= n; i++) {
                split(interval_list[i], coords, ",")
                tstart = coords[1]
                tend = coords[2]

                # Check for overlap: intervals overlap if qend > tstart && qstart < tend
                if (qend > tstart && qstart < tend) {
                    print $0
                    break  # Print once per query interval (like -wa)
                }
            }
        }
    }
    '''

    # Construct and run the AWK command
    awk_command = f"awk '{awk_script}' {shlex.quote(options.target)} {shlex.quote(options.query)}"

    try:
        result = subprocess.run(awk_command, shell=True, check=True,
                              capture_output=options.quiet, text=True)
        if options.quiet and result.stdout:
            print(result.stdout, end='')
    except subprocess.CalledProcessError as e:
        print(f"Error running awk command: {e}", file=sys.stderr)
        if e.stderr:
            print(f"AWK error: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: awk command not found. Please ensure awk is installed and in PATH.", file=sys.stderr)
        sys.exit(1)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Find overlapping intervals between two BED files (replicates bedtools intersect -wa)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -t reference.bed -q query.bed
  %(prog)s --target genes.bed --query peaks.bed --quiet
        """
    )
    parser.add_argument("-t", "--target",
                       type=str, required=True,
                       help="Target/reference intervals BED file")
    parser.add_argument("-q", "--query",
                       type=str, required=True,
                       help="Query intervals BED file")
    parser.add_argument("--quiet",
                       action="store_true",
                       help="Suppress error messages from subprocess")
    parser.add_argument("--version",
                       action="version",
                       version="%(prog)s 1.0")

    return parser.parse_args()

if __name__ == "__main__":
    main()
