"""Checks for the measurement in `tool_call`.
Run with `python3 tests/core/test_calls.py` (or pytest). Needs `/usr/bin/time`."""
import contextlib
import io
import shutil
import sys
import tempfile
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
import calls
import utility

TOOLS = ["tabix", "bedtools", "bedtools_sorted", "bedops", "bedmaps", "giggle", "granges", "gia",
         "bedtk", "bedtk_sorted", "igd", "ailist", "ucsc", "awk", "intervaltree"]


def test_missing_rss_line():
    """A missing RSS line in the /usr/bin/time output must give 0 MB, not -1/1024 (#7); a present one is
    converted with the unit of the platform, bytes on macOS and kilobytes on Linux (#10)."""
    assert utility.get_rss_from_stderr("whatever", "no such label") == -1
    real = utility.get_time_rss_label
    utility.get_time_rss_label = lambda: "no such label"
    try:
        _, mem = calls.tool_call("true", io.StringIO())
    finally:
        utility.get_time_rss_label = real
    assert mem == 0, f"mem={mem}"
    _, mem = calls.tool_call("true", io.StringIO())
    assert 0.01 < mem < 100, f"mem={mem} MB for `true`: RSS unit wrong? (bytes on macOS, kB on Linux, #10)"


def test_failed_call_raises():
    """A measured command that fails must raise instead of yielding a runtime, a memory value and an empty
    result (#9); the stderr of the command still goes to the log."""
    for call, code in [("false", 1), ("sh -c 'echo broken >&2; exit 3'", 3), ("segmeter_no_such_tool", 127)]:
        log = io.StringIO()
        try:
            calls.tool_call(call, log)
        except RuntimeError as error:
            message = str(error)
        else:
            raise AssertionError(f"{call!r} did not raise")
        assert f"exit code {code}" in message and call in message, message
        assert f"Executing: {call}" in log.getvalue()
        if code == 3:
            assert "broken" in log.getvalue() and "broken" in message # the command's stderr is logged and reported
    log = io.StringIO()
    calls.tool_call("sh -c 'echo fine >&2'", log) # stderr output alone is not a failure
    assert "fine" in log.getvalue()


def test_query_steps_counted():
    """Every measured step of a query adds its runtime and folds its memory in with max, whichever step is the
    peak (#11: bedops dropped the memory of its query sort). The tools are replaced by a stub, so no tool is needed."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for name in ["ref.bed", "ref_sorted.bed", "query.bed"]:
            (tmp / name).write_text("chr1\t100\t200\tintvl_1\n")
        (tmp / "chromlens.txt").write_text("chr1\t1000\n")
        (tmp / "idx").mkdir()
        reffiles = {"ref-unsrt": tmp / "ref.bed", "ref-srt": tmp / "ref_sorted.bed", "idx": tmp / "idx" / "L.bed",
                    "chromlens": tmp / "chromlens.txt"}
        for memory in [[100, 10], [10, 100]]: # the peak is the first, then the last step
            def tool_call(call, logfile): # 0.5 s and the next memory value; creates the redirected output, like a tool
                if ">" in call:
                    Path(call.rsplit(">", 1)[1].strip()).write_text("")
                steps.append(call)
                return 0.5, memory[min(len(steps), len(memory)) - 1]
            original = calls.tool_call, calls.subprocess.run
            calls.tool_call, calls.subprocess.run = tool_call, lambda *args, **kwargs: None # the unmeasured steps: bedtk's bedtools pass, sorted_genome's sort
            try:
                for tool in TOOLS:
                    steps = []
                    options = types.SimpleNamespace(tool=tool, datadir=str(tmp), benchname="b", logfile=io.StringIO())
                    rt, mem, out = calls.query_call(options, "L", reffiles, tmp / "query.bed")
                    Path(out.name).unlink()
                    assert steps, f"{tool}: no branch in query_call"
                    assert rt == 0.5 * len(steps), f"{tool}: {rt} s for {len(steps)} steps"
                    assert mem == max(memory[:len(steps)]), f"{tool}: {mem} MB, steps {memory[:len(steps)]}"
                    steps = [] # an empty query file (an empty complex bin) runs no tool and reports 0 s, 0 MB (#9)
                    (tmp / "empty.bed").write_text("")
                    rt, mem, out = calls.query_call(options, "L", reffiles, tmp / "empty.bed")
                    assert (steps, rt, mem, Path(out.name).read_text()) == ([], 0, 0, ""), f"{tool}: {steps}"
                    Path(out.name).unlink()
            finally:
                calls.tool_call, calls.subprocess.run = original


def test_index_steps_counted():
    """Every measured step of an index adds its runtime and folds its memory in with max, and the index size is
    the sum of the index files (#52). The tools are replaced by a stub that creates the files a tool would."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "ref").mkdir()
        (tmp / "ref" / "L.bed").write_text("chr1\t100\t200\tintvl_1\n")
        refdirs = {"ref": tmp / "ref", "idx": tmp / "idx"}
        for memory in [[100, 10, 10], [10, 100, 10], [10, 10, 100]]: # the peak in every position (tabix has three steps)
            def tool_call(call, logfile): # 0.5 s and the next memory value; writes 1 MB into every index file a tool would create
                steps.append(call)
                made = [Path(call.rsplit(">", 1)[1].strip())] if ">" in call else []
                if call.startswith("tabix "):
                    made.append(Path(call.split()[-1] + ".csi"))
                if call.startswith("giggle index"): # a directory of files, next to idx/ (see index_call)
                    made.append(tmp / "bench" / "b" / "giggle" / "L_index" / "cache.0.dat")
                if call.startswith("igd create"): # igd create <in> <out> <label>: a directory with <label>.igd
                    made.append(Path(call.split()[3]) / "L.igd")
                for path in made:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(b"x" * 2**20)
                return 0.5, memory[min(len(steps), len(memory)) - 1]
            original = calls.tool_call
            calls.tool_call = tool_call
            try:
                for tool in TOOLS:
                    steps = []
                    shutil.rmtree(tmp / "idx", ignore_errors=True)
                    (tmp / "idx").mkdir()
                    options = types.SimpleNamespace(tool=tool, datadir=str(tmp), benchname="b", logfile=io.StringIO())
                    with contextlib.redirect_stdout(io.StringIO()): # "Indexing ... with <tool>..."
                        rt, mem, size = calls.index_call(options, refdirs, "L")
                    assert rt == 0.5 * len(steps), f"{tool}: {rt} s for {len(steps)} steps"
                    assert mem == max(memory[:len(steps)], default=0), f"{tool}: {mem} MB, steps {memory[:len(steps)]}"
                    expected = {"tabix": 2.0, "giggle": 1.0, "igd": 1.0}.get(tool, 0) # 1 MB per index file
                    assert size == expected, f"{tool}: index size {size} MB, expected {expected}"
            finally:
                calls.tool_call = original
        shutil.rmtree(tmp / "idx")
        (tmp / "idx").mkdir()
        for missing in [tmp / "nope", tmp / "idx"]: # no index: a missing path, an empty directory
            assert raises(RuntimeError if missing.exists() else FileNotFoundError, calls.index_size, missing), missing


def raises(exception, func, *args):
    try:
        func(*args)
    except exception:
        return True
    return False


if __name__ == "__main__":
    test_missing_rss_line()
    test_failed_call_raises()
    test_query_steps_counted()
    test_index_steps_counted()
    print("ok")
