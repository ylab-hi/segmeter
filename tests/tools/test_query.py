"""Every tool of `segmeter bench` reports exactly the reference intervals that overlap a query.
Run with `python3 tests/tools/test_query.py` (or pytest). Needs `/usr/bin/time`; a tool that is not
installed is skipped, so the full table only runs across the project containers."""
import importlib.util
import io
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

# (tool, requirement, runs index_call, query directory: bedops picks its branch from the query path)
# the index_call column mirrors idx_based_tools in benchmark.py
TOOLS = [
    ("tabix", "tabix", True, "query"),
    ("bedtools", "bedtools", False, "query"),
    ("bedtools_sorted", "bedtools", True, "query"),
    ("bedtools_tabix", "bedtools+tabix", True, "query"),
    ("bedops", "bedops", True, "basic"),
    ("bedops", "bedops", True, "complex"),
    ("bedops", "bedops", True, "query"), # arbitrary target/query pair, empty before #7
    ("bedmaps", "bedmap", True, "query"),
    ("giggle", "/giggle/bin/giggle", True, "query"),
    ("granges", "granges", False, "query"),
    ("gia", "gia", False, "query"),
    ("gia_sorted", "gia", True, "query"),
    ("bedtk", "bedtk+bedtools", False, "query"),
    ("bedtk_sorted", "bedtk+bedtools", True, "query"),
    ("awk", "awk", False, "query"),
    ("intervaltree", "py:intervaltree", False, "query"),
    ("igd", "igd", True, "query"),
    ("ailist", "ailist", False, "query"),
    ("ucsc", "bedIntersect", False, "query"),
]
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
    (sim / f"{LABEL}_chromlens.txt").write_text("".join(f"{c}\t2000000\n" for c in CHROMS))
    for qdir in {row[3] for row in TOOLS}:
        (datadir / qdir).mkdir()
        write_bed(datadir / qdir / "Q.bed", queries)
    # ponytail: O(ref * query) oracle, fine for 5000 x 500
    return {(c, str(s), str(e)) for c, s, e, _ in ref
            if any(qc == c and qs < e and s < qe for qc, qs, qe, _ in queries)}


def intervals(text):
    return {tuple(line.split("\t")[:3]) for line in text.splitlines() if line and not line.startswith("#")}


def run_tool(tool, indexed, datadir, queryfile):
    options = types.SimpleNamespace(simdata=True, tool=tool, idx_based_tools=[tool] if indexed else [],
                                    datadir=str(datadir), simname="sim_001", format="BED",
                                    benchname="bench_001", logfile=io.StringIO())
    bench = BenchTool(options)
    if indexed:
        calls.index_call(options, bench.refdirs, LABEL, NUM)
    _, _, out = calls.query_call(options, LABEL, NUM, bench.get_reffiles(LABEL), queryfile)
    return intervals(Path(out.name).read_text())


def test_query_tools():
    with tempfile.TemporaryDirectory() as tmp:
        datadir = Path(tmp)
        expected = make_data(datadir)
        assert 0 < len(expected) < NUM
        failed = []
        for tool, requirement, indexed, qdir in TOOLS:
            if not available(requirement):
                print(f"skip {tool}: {requirement} not installed")
                continue
            got = run_tool(tool, indexed, datadir, datadir / qdir / "Q.bed")
            print(f"{'ok' if got == expected else 'FAIL'} {tool} ({qdir}): {len(got)} of {len(expected)} overlaps")
            if got != expected:
                failed.append(f"{tool} ({qdir})")
        assert not failed, f"{failed} report other overlaps than expected"


if __name__ == "__main__":
    test_query_tools()
    print("ok")
