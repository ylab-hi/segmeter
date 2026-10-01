"""Checks for the helpers in `utility`.
Run with `python3 tests/core/test_utility.py` (or pytest). No external tools needed."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
import utility


def test_chrom_sort_key():
    """The order genomap/granges expect in a genome file (numbers, arms, then X Y M, then the rest)."""
    names = ["chrX", "scaffold_1", "chr10", "chr2R", "chrM", "chr2", "chrY", "chr2L", "chr1", "chrMT", "contig"]
    assert sorted(names, key=utility.chrom_sort_key) == [
        "chr1", "chr2", "chr2L", "chr2R", "chr10", "chrX", "chrY", "chrM", "chrMT", "contig", "scaffold_1"]


if __name__ == "__main__":
    test_chrom_sort_key()
    print("ok")
