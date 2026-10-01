"""Checks for the helpers in `utility`.
Run with `python3 tests/core/test_utility.py` (or pytest). No external tools needed."""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
import utility


def test_chrom_sort_key():
    """The order genomap/granges expect in a genome file (numbers, arms, then X Y M Z W O, then the rest)."""
    names = ["chrX", "scaffold_1", "chr10", "chr2R", "chrM", "chr2", "chrY", "chr2L", "chr1", "chrMT", "contig",
             "chrW", "chrZ", "chrO", "Mt", "22", "X", "3"]
    assert sorted(names, key=utility.chrom_sort_key) == [
        "chr1", "chr2", "chr2L", "chr2R", "3", "chr10", "22", "chrX", "X", "chrY", "chrM", "chrMT", "Mt",
        "chrZ", "chrW", "chrO", "contig", "scaffold_1"]


def test_chrom_lengths():
    """Largest end per chromosome over several files; headers and blank lines are skipped, other junk is an error."""
    with tempfile.TemporaryDirectory() as tmp:
        a, b, bad = Path(tmp) / "a.bed", Path(tmp) / "b.bed", Path(tmp) / "bad.bed"
        a.write_text("track name=x\nbrowser position chr1\n# comment\nchr1\t10\t200\tx\n\nchr2\t5\t50\r\n")
        b.write_text("chr1\t0\t300\nchrX\t1\t2") # no trailing newline, chromosome only in b
        bad.write_text("chr1 10 200\n")
        assert utility.chrom_lengths([a, b]) == {"chr1": 300, "chr2": 50, "chrX": 2}
        try:
            utility.chrom_lengths([bad])
        except ValueError as e:
            assert "bad.bed" in str(e) and "chr1 10 200" in str(e)
        else:
            raise AssertionError("space-separated line accepted")


if __name__ == "__main__":
    test_chrom_sort_key()
    test_chrom_lengths()
    print("ok")
