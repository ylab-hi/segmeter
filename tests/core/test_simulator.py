"""Checks for the simulator.
Run with `python3 tests/core/test_simulator.py` (or pytest). No external tools needed."""
import contextlib
import glob
import io
import platform
import sys
import tempfile
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
from simulator import SimBED


def test_select_chrom_scaffold():
    """A full chromosome is replaced by a scaffold that starts with a gap; the scaffold name counts the
    existing scaffolds (#7 fixed a membership test that always gave SCF1)."""
    sim = SimBED(types.SimpleNamespace(gapsize="100-100", max_chromlen=1000), {})
    chroms = {"all": [], "space-left": ["chr1"], "intvl": {"chr1": 10},
              "leftgap": {"chr1": {"start": 1, "end": 1000, "mid": 500}}}
    assert sim.select_chrom(chroms) == "SCF1"
    assert chroms["space-left"] == ["SCF1"] and chroms["intvl"]["SCF1"] == 0
    assert chroms["leftgap"]["SCF1"] == {"start": 1, "end": 100, "mid": 50}
    chroms["leftgap"]["SCF1"]["end"] = 1000 # SCF1 is full too
    assert sim.select_chrom(chroms) == "SCF2"


def test_max_span():
    intvls = [["chr1", str(i * 100), str(i * 100 + 50)] for i in range(200)]
    for max_span, expected in [(None, 199), (1000, 199), (50, 49)]: # spans 2..200 and 2..50
        sim = SimBED(types.SimpleNamespace(gapsize="100-100", max_span=max_span), {})
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
        sim = SimBED(types.SimpleNamespace(gapsize="100-100", max_span=1000, datadir=tmp, simname="sim_001"), {})
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


def simulate(tmp, seed):
    """A full `segmeter sim -n 10` run; returns the simulation folder and every file in it, relative path -> bytes."""
    options = types.SimpleNamespace(datadir=tmp, simname="sim_001", intvlnums="10", intvlsize="100-10000", gapsize="100-5000",
                                    max_chromlen=1000000000, max_span=None, seed=seed)
    with contextlib.redirect_stdout(io.StringIO()): # "Simulate intervals for ..."
        SimBED(options, {"10": 10}).sim_intervals()
    out = Path(tmp) / "sim" / "sim_001" / "BED"
    return out, {str(p.relative_to(out)): p.read_bytes() for p in out.rglob("*") if p.is_file()}


def test_seed():
    """The same seed gives the same data, another seed other data, an unseeded run records the seed it drew, and a second run
    into the same folder appends its block to parameters.txt instead of replacing the first (#15)."""
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b, tempfile.TemporaryDirectory() as c:
        out, files_a = simulate(a, None)
        params = dict(line.split("\t") for line in (out / "parameters.txt").read_text().splitlines() if line)
        assert params["intvlnums"] == "10" and params["max_span"] == "None" and params["python"] == platform.python_version()
        seed = int(params["seed"])
        _, files_b = simulate(b, seed)
        _, files_c = simulate(c, seed + 1)
        assert files_a == files_b
        assert files_a["ref/10.bed"] != files_c["ref/10.bed"]
        assert len(files_a["basic/query/perfect/10_30p.bed"].splitlines()) == 3 # 30% of 10 queries
        assert len({l for f in files_a if f.startswith("basic/query/perfect/10_") for l in files_a[f].splitlines()}) == 10 # samples of the 10 queries
        simulate(a, seed + 1) # a later run into the same simname (e.g. another -n) keeps the first run's block
        blocks = (out / "parameters.txt").read_text().split("\n\n")
        assert [b.split("\n")[0] for b in blocks if b] == [f"seed\t{seed}", f"seed\t{seed + 1}"]


if __name__ == "__main__":
    test_select_chrom_scaffold()
    test_max_span()
    test_max_span_bins()
    test_seed()
    print("ok")
