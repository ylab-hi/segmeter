"""Every tool of `segmeter bench` reports exactly the reference intervals that overlap a query, in
simulated-data mode (`bench -r`) and in real-data mode (`bench --target --query`, #19/#33), and
reports one output line per (query, reference) pair, which the complex score counts (#34).
Run with `python3 tests/tools/test_query.py` (or pytest). Needs `/usr/bin/time`; a tool that is not
installed is skipped, so the full table only runs across the project containers."""
import collections
import contextlib
import importlib.util
import io
import os
import random
import shutil
import sys
import tempfile
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
import calls
import utility
from BenchTool import BenchTool
from benchmark import BenchBase

# (tool, requirement, runs index_call, query directory: bedops picks its branch from the query path)
# the index_call column mirrors idx_based_tools in benchmark.py
TOOLS = [
    ("tabix", "tabix", True, "query"),
    ("bedtools", "bedtools", False, "query"),
    ("bedtools_sorted", "bedtools", True, "query"),
    ("bedops", "bedops", True, "basic"),
    ("bedops", "bedops", True, "complex"),
    ("bedops", "bedops", True, "query"), # arbitrary target/query pair, empty before #7
    ("bedmaps", "bedmap", True, "query"),
    ("giggle", "/giggle/bin/giggle", True, "query"),
    ("granges", "granges+bedtools", False, "query"),
    ("gia", "gia", False, "query"),
    ("bedtk", "bedtk+bedtools", False, "query"),
    ("bedtk_sorted", "bedtk+bedtools", True, "query"),
    ("awk", "awk", False, "query"),
    ("intervaltree", "py:intervaltree", False, "query"),
    ("igd", "igd", True, "query"),
    ("ailist", "ailist", False, "query"),
    ("ucsc", "bedIntersect+bedtools", False, "query"),
]
READS_INDEX = {"tabix", "bedtools_sorted", "bedtk_sorted", "bedops", "bedmaps", "igd"} # query reads refdirs["idx"]
CHROMS = ["chr1", "chr2", "chr10", "chrX"]
LABEL, NUM = "L", 5000


def available(requirement):
    """'a+b': every binary (or path, or 'py:' module) must be present."""
    for req in requirement.split("+"):
        if req.startswith("py:"):
            if importlib.util.find_spec(req[3:]) is None:
                return False
        elif not (shutil.which(req) or Path(req).exists()):
            return False
    return True


def random_bed(n, prefix, seed):
    random.seed(seed)
    for i in range(n):
        start = random.randint(0, 1_000_000)
        yield random.choice(CHROMS), start, start + random.randint(1, 5000), f"{prefix}_{i}"


def write_bed(path, rows):
    with open(path, "w") as fh:
        fh.writelines("\t".join(map(str, row)) + "\n" for row in rows)


def make_data(datadir):
    """Reference and queries laid out like a simulation, so BenchTool resolves the paths."""
    ref, queries = list(random_bed(NUM, "intvl", 7)), list(random_bed(NUM // 10, "query", 11))
    sim = datadir / "sim" / "sim_001" / "BED"
    (sim / "ref").mkdir(parents=True)
    write_bed(sim / "ref" / f"{LABEL}.bed", ref)
    utility.sort_BED(sim / "ref" / f"{LABEL}.bed", sim / "ref" / f"{LABEL}_sorted.bed")
    # not in natural chromosome order, like the simulator's random order: granges needs it reordered (#36)
    (sim / f"{LABEL}_chromlens.txt").write_text("".join(f"{c}\t2000000\n" for c in reversed(CHROMS)))
    # ponytail: O(ref * query) oracle, fine for 5000 x 500
    hits = {q: [r for r in ref if r[0] == q[0] and q[1] < r[2] and r[1] < q[2]] for q in queries}
    for qdir in {row[3] for row in TOOLS}:
        (datadir / qdir).mkdir()
        write_bed(datadir / qdir / "Q.bed", queries)
    # like the simulated complex queries: every one has hits, and the path says "complex", which selects the bedmap branch of
    # bedops and the duplicates pass of bedtk, granges and ucsc (#69)
    write_bed(datadir / "complex" / "C.bed", [q for q in queries if hits[q]])
    pairs = collections.Counter((c, str(s), str(e)) for q in queries for c, s, e, _ in hits[q]) # reference -> number of queries hitting it
    touching = sum(r[0] == q[0] and (q[1] == r[2] or r[1] == q[2]) for q in queries if hits[q] for r in ref)
    return set(pairs), pairs, touching


def intervals(text):
    return {tuple(line.split("\t")[:3]) for line in text.splitlines() if line and not line.startswith("#")}


def run_tool(tool, indexed, datadir, queryfile, index=True):
    """One query call; `index=False` reuses the index of an earlier call."""
    options = types.SimpleNamespace(simdata=True, tool=tool, idx_based_tools=[tool] if indexed else [],
                                    datadir=str(datadir), simname="sim_001", format="BED",
                                    benchname="bench_001", logfile=io.StringIO())
    bench = BenchTool(options)
    if indexed and index:
        calls.index_call(options, bench.refdirs, LABEL)
    logged = options.logfile.tell()
    _, _, out = calls.query_call(options, LABEL, bench.get_reffiles(LABEL), queryfile)
    if tool in READS_INDEX: # the query must read what the index step wrote, not the simulator's sorted copy (#39)
        assert str(bench.refdirs["idx"]) in options.logfile.getvalue()[logged:], f"{tool}: query does not read its index"
    text = Path(out.name).read_text()
    Path(out.name).unlink() # the caller removes the tool's output, like BenchTool does
    # the set of references (basic score) and how often each is reported (complex score: one line per (query, reference) pair)
    return intervals(text), collections.Counter(tuple(line.split("\t")[:3]) for line in text.splitlines() if line and not line.startswith("#"))


def run_real(tool, datadir, target, queryfile):
    """`bench --target --query`: the whole BenchBase run, as main.py starts it; the overlaps land in result.bed."""
    (datadir / "real").mkdir(exist_ok=True)
    options = types.SimpleNamespace(simdata=False, tool=tool, target=str(target), query=str(queryfile),
                                    datadir=str(datadir / "real"), simname="sim_001", format="BED",
                                    benchname="bench_001")
    with contextlib.redirect_stderr(io.StringIO()): # the tools' stderr
        BenchBase(options, {})
    return intervals((datadir / "real" / "bench" / "bench_001" / tool / "result.bed").read_text())


def test_query_tools():
    with tempfile.TemporaryDirectory() as tmp:
        datadir = Path(tmp)
        expected, pairs, touching = make_data(datadir)
        assert 0 < len(expected) < sum(pairs.values()) and touching > 0 # the seeds give duplicates and a touching pair (giggle must not report it, #70); a reseed must keep both
        tempfile.tempdir = str(datadir / "tmp") # every temporary file of the queries lands here (#2)
        (datadir / "tmp").mkdir()
        failed = []
        try:
            for tool, requirement, indexed, qdir in TOOLS:
                if not available(requirement):
                    print(f"skip {tool}: {requirement} not installed")
                    continue
                runs = {"sim": run_tool(tool, indexed, datadir, datadir / qdir / "Q.bed")[0]}
                if qdir == "query": # bedops' basic/complex rows are simulated-data only
                    runs["real"] = run_real(tool, datadir, datadir / "sim" / "sim_001" / "BED" / "ref" / f"{LABEL}.bed",
                                            datadir / qdir / "Q.bed")
                for mode, got in runs.items():
                    print(f"{'ok' if got == expected else 'FAIL'} {tool} ({qdir}, {mode}): {len(got)} of {len(expected)} overlaps")
                    if got != expected:
                        failed.append(f"{tool} ({qdir}, {mode})")
                # complex score (#34): get_precision counts the output lines against the (query, reference) pairs,
                # so a tool must report a reference once per query that hits it (bedtk, granges and ucsc get their
                # duplicates back in query_call for complex query files, #69); the basic rows above run on the raw output
                if tool != "bedops" or qdir == "complex": # once per tool
                    _, counts = run_tool(tool, indexed, datadir, datadir / "complex" / "C.bed", index=False)
                    print(f"{'ok' if counts == pairs else 'FAIL'} {tool} ({qdir}, complex): {sum(counts.values())} lines for {sum(pairs.values())} pairs, "
                          f"{sum(counts[r] != pairs[r] for r in counts.keys() | pairs.keys())} references with a wrong count")
                    if counts != pairs:
                        failed.append(f"{tool} ({qdir}, complex)")
                leftover = os.listdir(tempfile.tempdir)
                assert not leftover, f"{tool} leaves temporary files behind: {leftover}"
        finally:
            tempfile.tempdir = None # also after a failure, so later tempfile calls do not use the deleted directory
        assert not failed, f"{failed} differ from the oracle"


if __name__ == "__main__":
    test_query_tools()
    print("ok")
