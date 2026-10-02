"""Checks for the precision scoring.
Run with `python3 tests/core/test_scoring.py` (or pytest). No external tools needed."""
import atexit
import os
import sys
import tempfile
import time
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
from BenchTool import BenchTool


def write_tmp(lines):
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    fh.writelines(lines)
    fh.close()
    atexit.register(os.unlink, fh.name)
    return fh


def test_basic_scoring():
    n = 50000
    ref = [("chr1", str(i * 100), str(i * 100 + 50)) for i in range(n)]
    truth = {r: (r, f"intvl_{i}") for i, r in enumerate(ref)}
    queries = write_tmp("\t".join(r) + "\n" for r in ref)
    results = write_tmp("\t".join(r) + "\n" for r in ref[: n // 2]) # tool finds half

    bench = BenchTool.__new__(BenchTool)
    start = time.time()
    hit = bench.get_precision(queries.name, results, truth, "basic", "perfect")["basic"]
    gap = bench.get_precision(queries.name, results, truth, "basic", "mid-gap1")["basic"]
    elapsed = time.time() - start

    assert (hit["TP"], hit["FN"]) == (n // 2, n // 2)
    assert (gap["FP"], gap["TN"]) == (n // 2, n // 2)
    assert elapsed < 5, f"scoring took {elapsed:.1f}s, is it quadratic again?"


def test_complex_scoring():
    truth = {("chr1", "0", "500"): "3", ("chr1", "100", "900"): "5"}
    queries = write_tmp("\t".join(q) + "\tmult\n" for q in truth)
    results = write_tmp("chr1\t0\t50\n" for _ in range(6)) # 8 expected, 6 found

    bench = BenchTool.__new__(BenchTool)
    assert bench.get_precision(queries.name, results, truth, "complex", "mult")["complex"]["dist"] == 2


def test_simdata_querydirs():
    """bench -r crashed since v0.13.1 because querydirs were only set without --simdata (#21)."""
    options = types.SimpleNamespace(simdata=True, tool="bedtools", idx_based_tools=[],
                                    datadir="data", simname="sim_001", format="BED")
    bench = BenchTool(options)
    assert set(bench.querydirs) == {"basic", "complex"}
    assert bench.querydirs["basic"]["perfect"] == Path("data/sim/sim_001/BED/basic/query/perfect")
    assert bench.querydirs["complex"] == {"mult": Path("data/sim/sim_001/BED/complex/query/mult")}


if __name__ == "__main__":
    test_basic_scoring()
    test_complex_scoring()
    test_simdata_querydirs()
    print("ok")
