"""Checks for the benchmark harness around the tool calls.
Run with `python3 tests/core/test_benchmark.py` (or pytest). No external tools needed."""
import sys
import tempfile
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
import calls
from BenchTool import BenchTool
from benchmark import BenchBase
from main import det_intvlnums
import utility


def raises(exception, func, *args):
    try:
        func(*args)
    except exception:
        return True
    return False


def test_det_intvlnums():
    assert det_intvlnums("10") == {"10": 10}
    assert det_intvlnums("1K") == {"1K": 1000}
    assert det_intvlnums("10,1K,3M") == {"10": 10, "1K": 1000, "3M": 3000000}
    assert raises(ValueError, det_intvlnums, "x")


def test_validate():
    with tempfile.TemporaryDirectory() as tmp:
        bench = BenchBase.__new__(BenchBase)
        bench.options = types.SimpleNamespace(datadir=str(Path(tmp) / "missing"), simname="sim_001", simdata=True, tool="bedtools")
        assert raises(FileNotFoundError, bench.validate) # no data directory
        bench.options.datadir = tmp
        assert raises(FileNotFoundError, bench.validate) # no simulation to benchmark
        bench.options.simdata = False # the simulation is not needed for a target/query pair
        bench.validate()
        bench.options.tool = None
        assert raises(ValueError, bench.validate)


def test_negatives_file():
    """Every FP/FN of a subset is written to the negatives file."""
    precision = {"basic": {10: {"TP": 1, "FP": 1, "TN": 0, "FN": 1, "negatives": ["intvl_1_5p:intvl_1\tFN\n", "intvl_2_mid-gap1:intvl_2\tFP\n"]}},
                 "complex": {10: {"TP": 3, "FP": 1, "FN": 0, "dist": 3}}}
    with tempfile.TemporaryDirectory() as tmp:
        stats, negatives = Path(tmp) / "precision.txt", Path(tmp) / "negatives.txt"
        BenchBase.__new__(BenchBase).save_query_prec_stats(100, precision, stats, negatives)
        assert negatives.read_text().splitlines()[1:] == ["intvl_1_5p:intvl_1\tFN", "intvl_2_mid-gap1:intvl_2\tFP"]
        assert stats.read_text().splitlines()[1] == "100\t10%\t1\t1\t0\t1\t0.5\t0.5\t0.5"
        assert stats.read_text().splitlines()[-2] == "intvlnum\tbin\tTP\tFP\tFN\tPrecision\tRecall\tF1\tdistance"
        assert stats.read_text().splitlines()[-1] == f"100\t10bin\t3\t1\t0\t0.75\t1.0\t{2 * 0.75 / 1.75}\t3"


def test_query_intervals():
    """Time and memory of the tool call are recorded per query type, and the scores are summed over the types.
    The tool is replaced by one that takes 1.2345678 s and 42 MB and reports the first of two intervals."""
    ref = [("chr1", "100", "200"), ("chr1", "300", "400")]
    with tempfile.TemporaryDirectory() as tmp:
        sim = Path(tmp) / "sim" / "sim_001" / "BED"
        bench = BenchTool(types.SimpleNamespace(simdata=True, tool="bedtools", idx_based_tools=[],
                                                datadir=tmp, simname="sim_001", format="BED"))
        for path in [*bench.refdirs.values(), *bench.querydirs["basic"].values(), bench.querydirs["complex"]["mult"]]:
            path.mkdir(parents=True)
        truth = []
        for n, query in enumerate(utility.BASIC_QUERIES): # each query type has its own coordinates
            rows = [("chr1", str(1000 * n + i), str(1000 * n + i + 1)) for i in range(2)]
            (bench.querydirs["basic"][query] / "L_10p.bed").write_text("".join("\t".join(row) + "\n" for row in rows))
            truth += ["\t".join(row + r + (f"intvl_{i}_{query}:intvl_{i}",)) + "\n" for i, (row, r) in enumerate(zip(rows, ref))]
        (sim / "basic" / "truth" / "L.bed").write_text("".join(truth))
        (sim / "ref" / "L_sorted.bed").write_text("".join("\t".join(r) + "\n" for r in ref)) # the complex score reads it (#74)
        (sim / "complex" / "truth" / "L.bed").write_text("chr1\t100\t400\tmult_2\t2\n")
        (bench.querydirs["complex"]["mult"] / "L_10bin.bed").write_text("chr1\t100\t400\tmult_2\n")
        result = Path(tmp) / "result.bed"
        def query_call(options, label, reffiles, queryfile): # returns the closed temporary file with the tool's output
            result.write_text("\t".join(ref[0]) + "\n")
            return 1.2345678, 42.0, types.SimpleNamespace(name=str(result))

        original = calls.query_call
        calls.query_call = query_call
        try:
            times, memory, precision = bench.query_intervals("L", 2, 10)
        finally:
            calls.query_call = original
        assert not result.exists() # the output is removed once scored (#2)

    assert list(times["basic"]) == list(utility.BASIC_QUERIES) and list(times["complex"]) == ["mult"]
    assert all(value == {10: 1.23457} for dtype in times.values() for value in dtype.values())
    assert all(value == {10: 42.0} for dtype in memory.values() for value in dtype.values())
    basic = precision["basic"][10] # five interval types and five gap types, one of two intervals reported
    assert (basic["TP"], basic["FN"], basic["FP"], basic["TN"]) == (5, 5, 5, 5)
    assert len(basic["negatives"]) == 10
    assert precision["complex"][10] == {"TP": 1, "FP": 0, "FN": 1, "dist": 1} # two references covered, one reported


if __name__ == "__main__":
    test_det_intvlnums()
    test_validate()
    test_negatives_file()
    test_query_intervals()
    print("ok")
