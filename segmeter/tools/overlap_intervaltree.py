#!/usr/bin/env python3

import argparse
import sys
from intervaltree import IntervalTree, Interval

def main():
    options = parse_arguments()

    # read input files
    # print(f"Reading target intervals from {options.target}...", file=sys.stderr)
    target = read_target_intervals(options.target)

    # print(f"Reading query intervals from {options.query}...", file=sys.stderr)
    query_intervals(target, options.query)


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
            }
    return trees

def query_intervals(trees, query_file_path):
    with open(query_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chrom in trees:
                overlaps = trees[chrom][start:end]
                for interval in overlaps:
                    t_data = interval.data
                    t_start, t_end = interval.begin, interval.end

                    print(f"{chrom}\t{start}\t{end}")

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
    parser.add_argument("--min-overlap", type=int, default=1,
                        help="Minimum overlap length required (default: 1)")
    parser.add_argument("--format", choices=['bed', 'detailed'], default='bed',
                        help="Output format (default: bed)")
    parser.add_argument("--stats", action='store_true',
                        help="Print overlap statistics")
    parser.add_argument("--version", action='version', version='%(prog)s 1.0')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()
