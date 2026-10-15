#!/usr/bin/env python3
"""Overlaps of two BED files with polars-bio: both files are scanned as tab-separated frames, tagged as 0-based
half-open (an untagged frame is taken as 1-based closed, which also counts touching intervals), and
`polars_bio.overlap(..., overlap_output="left")` returns every target row once per query that overlaps it, like
`bedtools intersect -wa`. The target rows are written as they were read. polars and DataFusion are pinned to one
thread, like every other benchmarked tool."""
import argparse
import os

os.environ.setdefault("POLARS_MAX_THREADS", "1")  # read by polars at import
import polars as pl  # noqa: E402
import polars_bio as pb  # noqa: E402

pb.set_option("datafusion.execution.target_partitions", "1")


def scan_bed(path):
    # the chromosome column stays text: polars would infer Int64 from a file whose first rows are numbered chromosomes
    frame = pl.scan_csv(path, separator="\t", has_header=False, comment_prefix="#",
                        schema_overrides={"column_1": pl.String, "column_2": pl.Int64, "column_3": pl.Int64})
    names = frame.collect_schema().names()
    frame = frame.rename(dict(zip(names[:3], ("chrom", "start", "end"))))
    frame.config_meta.set(coordinate_system_zero_based=True)
    return frame


def main():
    parser = argparse.ArgumentParser(description="Find overlaps between two BED files using polars-bio")
    parser.add_argument("-t", "--target", required=True, help="Target BED file with intervals to search in")
    parser.add_argument("-q", "--query", required=True, help="Query BED file with intervals to find overlaps with")
    parser.add_argument("-o", "--output", required=True, help="Output file: the target intervals hit, one line per (query, target) pair")
    options = parser.parse_args()

    hits = pb.overlap(scan_bed(options.target), scan_bed(options.query), overlap_output="left", output_type="polars.DataFrame")
    hits.write_csv(options.output, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
