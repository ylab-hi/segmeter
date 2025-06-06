#!/usr/bin/env python3

import argparse
import sys
from intervaltree import IntervalTree, Interval

def main():
    options = parse_arguments()

    # read input files
    # print(f"Reading target intervals from {options.target}...", file=sys.stderr)
    ref_intvls = read_target_intervals(options.target)

    # print(f"Reading query intervals from {options.query}...", file=sys.stderr)
    query_intervals(ref_intvls, options)

def read_target_intervals(file_path):
    trees = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            if chrom not in trees:
                trees[chrom] = IntervalTree()
            trees[chrom][start:end] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'intvlid': parts[3] if len(parts) > 3 else None
            }
    return trees

def query_intervals(trees, options):
    ofh = open(options.output, 'w') # output file handle
    with open(options.query, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            intvlid = parts[3] if len(parts) > 3 else None
            if chrom in trees:
                overlaps = trees[chrom][start:end]
                for interval in overlaps:
                    t_data = interval.data
                    t_start, t_end = interval.begin, interval.end

                    if options.report == "target":
                        ofh.write(f"{t_data['chrom']}\t{t_start}\t{t_end}\t{t_data['intvlid']}\n")
                    elif options.report == "query":
                        ofh.write(f"{chrom}\t{start}\t{end}\t{intvlid}\n")

def parse_arguments():
    parser = argparse.ArgumentParser(
            description="Find overlaps between two BED files using Interval Tree",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
    Examples:
      %(prog)s -t target.bed -q query.bed
      %(prog)s -t target.bed -q query.bed --min-overlap 100 --format detailed
      %(prog)s -t target.bed -q query.bed --stats
            """)

    parser.add_argument("-t", "--target", type=str, required=True,
                        help="Target BED file with intervals to search in")
    parser.add_argument("-q", "--query", type=str, required=True,
                        help="Query BED file with intervals to find overlaps with")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output file to write overlaps to")
    parser.add_argument("-r", "--report", type=str, default="target",
                        choices=['target', 'query'],
                        help="Report overlaps from target or query intervals (default: target)")
    parser.add_argument("--format", choices=['bed', 'detailed'], default='bed',
                        help="Output format (default: bed)")
    parser.add_argument("--stats", action='store_true',
                        help="Print overlap statistics")
    parser.add_argument("--version", action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()
