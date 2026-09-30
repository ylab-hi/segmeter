"""Checks for the simulator.
Run with `python3 tests/core/test_simulator.py` (or pytest). No external tools needed."""
import glob
import io
import sys
import tempfile
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
from simulator import SimBED


def test_max_span():
    intvls = [["chr1", str(i * 100), str(i * 100 + 50)] for i in range(200)]
    for max_span, expected in [(None, 199), (1000, 199), (50, 49)]: # spans 2..200 and 2..50
        sim = SimBED(types.SimpleNamespace(max_span=max_span), {})
        queries, truth = io.StringIO(), io.StringIO()
        sim.sim_overlaps(intvls, {(1, 200): queries}, truth)
        spans = [int(line.split("\t")[4]) for line in truth.getvalue().splitlines()]
        assert spans == list(range(2, expected + 2))
        assert len(queries.getvalue().splitlines()) == expected


def test_max_span_bins():
    """The decile bins are computed from the capped span, and every bin query has a truth entry."""
    intvls = 1005 # one chromosome with more intervals than the cap
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "sim" / "sim_001" / "BED"
        sim = SimBED(types.SimpleNamespace(max_span=1000, datadir=tmp, simname="sim_001"), {})
        refdir, truthdirs, querydirs = sim.create_datadirs(out)
        with open(refdir / "L_sorted.bed", "w") as fh:
            fh.writelines(f"chr1\t{i * 100}\t{i * 100 + 50}\tintvl_{i}\n" for i in range(intvls))
        (out / "L_chrnums.txt").write_text(f"chr1\t{intvls}\n")
        sim.sim_complex_queries(refdir, truthdirs, querydirs, intvls, "L")

        truth = {tuple(l.split("\t")[:3]) for l in open(truthdirs["complex"] / "L.bed")}
        assert len(truth) == 999 # spans 2..1000
        bins = {}
        for f in glob.glob(str(querydirs["complex"]["mult"] / "L_*bin.bed")):
            spans = [int(l.split("\t")[3].split("_")[1]) for l in open(f)]
            bins[int(Path(f).stem.split("_")[1].rstrip("bin"))] = spans
            assert all(tuple(l.split("\t")[:3]) in truth for l in open(f))
        assert sorted(bins) == list(range(10, 101, 10))
        assert (min(bins[100]), max(bins[100])) == (901, 1000)
        assert sum(len(s) for s in bins.values()) == 999 # every span is in exactly one bin


if __name__ == "__main__":
    test_max_span()
    test_max_span_bins()
    print("ok")
