#!/usr/bin/env python3
"""Overlaps of two BED files with pygenogrove: the target intervals are inserted into a BED grove (one index per
chromosome), every query interval is intersected, and the target interval is printed once per query that hits it,
like `bedtools intersect -wa`. A grove key is closed and 0-based, BED is half-open, so a query [start, end) is
searched as [start, end-1]; the entry-deriving insert does the same for the targets."""
import argparse

import pygenogrove as pg


def main():
    parser = argparse.ArgumentParser(description="Find overlaps between two BED files using pygenogrove")
    parser.add_argument("-t", "--target", required=True, help="Target BED file with intervals to search in")
    parser.add_argument("-q", "--query", required=True, help="Query BED file with intervals to find overlaps with")
    parser.add_argument("-o", "--output", required=True, help="Output file: the target intervals hit, one line per (query, target) pair")
    parser.add_argument("-k", "--order", type=int, default=64, help="order of the B+ tree (default: %(default)s)")
    options = parser.parse_args()

    grove = pg.BedGrove(options.order)
    for entry in pg.BedReader(options.target):
        grove.insert(entry.chrom, entry)

    with open(options.output, "w") as out:
        for query in pg.BedReader(options.query):
            # '*' matches every strand, as the other tools ignore the strand column
            for hit in grove.intersect(pg.GenomicCoordinate("*", query.start, query.end - 1), query.chrom):
                data = hit.data
                out.write(f"{data.chrom}\t{data.start}\t{data.end}\n")


if __name__ == "__main__":
    main()
