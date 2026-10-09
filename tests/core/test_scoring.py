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
    """The complex score of a bin: TP/FP/FN on the set of references the queries cover, and the distance as the summed
    multiplicity differences, so a missing and an extra pair do not cancel (#74). The references a query covers come from
    the sorted reference by bisect and must agree with the truth file's count."""
    ref = [("chr1", str(i * 100), str(i * 100 + 50)) for i in range(10)] + [("chr2", "0", "50")] # chr1: 0-50 .. 900-950
    q1, q2 = ("chr1", "0", "450"), ("chr1", "300", "850") # cover ref 0..4 (5) and ref 3..8 (6): 11 pairs, 9 distinct references
    queries = write_tmp("\t".join(q) + "\tmult\n" for q in (q1, q2))
    expected_pairs = ref[0:5] + ref[3:9]
    bench = BenchTool.__new__(BenchTool)
    truth = bench.load_truth(write_tmp([]).name, write_tmp(f"{chr(9).join(q)}\tmult\t{n}\n" for q, n in ((q1, 5), (q2, 6))).name,
                             write_tmp("\t".join(r) + "\n" for r in ref).name)["complex"]
    assert truth["records"] == {q1: "5", q2: "6"} and sorted(truth["ref"]) == ["chr1", "chr2"] and list(truth["ref"]["chr1"][0]) == [i * 100 for i in range(10)]

    exact = write_tmp("\t".join(r) + "\n" for r in expected_pairs)
    assert bench.get_precision(queries.name, exact, truth, "complex", "mult")["complex"] == {"TP": 9, "FP": 0, "FN": 0, "dist": 0}

    # ref 3 reported once instead of twice (one pair missing), ref 9 reported although no query covers it (an extra pair),
    # ref 0 reported twice instead of once (an extra pair), ref 8 not at all (one pair missing): the line count is unchanged
    # (11 lines), the set and the pairs are not
    wrong = write_tmp("\t".join(r) + "\n" for r in ref[0:5] + ref[4:8] + [ref[9], ref[0]])
    score = bench.get_precision(queries.name, wrong, truth, "complex", "mult")["complex"]
    assert score == {"TP": 8, "FP": 1, "FN": 1, "dist": 4}, score

    assert bench.get_precision(queries.name, write_tmp([]), truth, "complex", "mult")["complex"] == {"TP": 0, "FP": 0, "FN": 9, "dist": 11}
    assert bench.get_precision(write_tmp([]).name, exact, truth, "complex", "mult")["complex"] == {"TP": 0, "FP": 0, "FN": 0, "dist": 0} # an empty bin

    wrong_truth = write_tmp(f"{chr(9).join(q1)}\tmult\t4\n")
    truth_bad = bench.load_truth(write_tmp([]).name, wrong_truth.name, write_tmp("\t".join(r) + "\n" for r in ref).name)["complex"]
    try:
        bench.get_precision(write_tmp([f"{chr(9).join(q1)}\tmult\n"]).name, exact, truth_bad, "complex", "mult")
        raise AssertionError("a truth count that disagrees with the reference must raise")
    except ValueError as e:
        assert "covers 5 references, the truth file says 4" in str(e)

    commented = write_tmp(["# a header\n"] + ["\t".join(r) + "\n" for r in expected_pairs])
    assert bench.get_precision(queries.name, commented, truth, "complex", "mult")["complex"] == {"TP": 9, "FP": 0, "FN": 0, "dist": 0}


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
