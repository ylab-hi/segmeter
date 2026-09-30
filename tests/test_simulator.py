"""Checks for the simulator.
Run with `python3 tests/test_simulator.py` (or pytest). No external tools needed."""
import io
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "segmeter"))
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


if __name__ == "__main__":
    test_max_span()
    print("ok")
