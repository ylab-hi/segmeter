"""Checks for the measurement in `tool_call`.
Run with `python3 tests/core/test_calls.py` (or pytest). Needs `/usr/bin/time`."""
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[2] / "segmeter"))
import calls
import utility


def test_missing_rss_line():
    """A missing RSS line in the /usr/bin/time output must give 0 MB, not -1/1024 (#7)."""
    assert utility.get_rss_from_stderr("whatever", "no such label") == -1
    real = utility.get_time_rss_label
    utility.get_time_rss_label = lambda: "no such label"
    try:
        _, mem = calls.tool_call("true", io.StringIO())
    finally:
        utility.get_time_rss_label = real
    assert mem == 0, f"mem={mem}"
    _, mem = calls.tool_call("true", io.StringIO())
    assert mem > 0, f"mem={mem} with the RSS line present"


if __name__ == "__main__":
    test_missing_rss_line()
    print("ok")
