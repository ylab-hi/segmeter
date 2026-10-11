"""Artifact collection must skip private indexes and preserve publication paths."""
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parents[2] / 'scripts'))
from collect_benchmark_artifacts import collect


def test_reports_only():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / 'results-head'
        source.mkdir()
        (source / 'versions.json').write_text(json.dumps({'sizes': '100K', 'repeats': 1}))
        (source / 'benchmark.json').write_text('[]')
        tool = source / 'bench/run-0/giggle'
        stats = tool / '100K/stats'
        stats.mkdir(parents=True)
        (stats / '100K_query_stats_100.txt').write_text('query stats')
        precision = tool / '100K/precision'
        precision.mkdir()
        (precision / '100K_query_precision_100.txt').write_text('scores')
        (tool / 'log.txt').write_text('commands')
        (tool / '100K/100K_idx_stats.txt').write_text('index stats')
        index = tool / '100K_index'
        index.mkdir()
        (index / 'private.txt').write_text('binary index')
        index.chmod(0)
        destination = root / 'artifacts'
        original_glob = Path.glob
        def guarded_glob(path, pattern):
            assert path.name in ('stats', 'precision'), f'Traversed tool files: {path}'
            return original_glob(path, pattern)
        try:
            with patch.object(Path, 'glob', guarded_glob):
                collect(source, destination)
        finally:
            index.chmod(0o700)
        staged = destination / 'results-head'
        assert (staged / 'benchmark.json').read_text() == '[]'
        assert (staged / stats.relative_to(source) / '100K_query_stats_100.txt').read_text() == 'query stats'
        assert (staged / precision.relative_to(source) / '100K_query_precision_100.txt').read_text() == 'scores'
        assert (staged / tool.relative_to(source) / 'log.txt').read_text() == 'commands'
        assert (staged / 'bench/run-0/giggle/100K/100K_idx_stats.txt').read_text() == 'index stats'
        assert not (staged / index.relative_to(source)).exists()
        collect(root / 'missing-base', destination)


if __name__ == '__main__':
    test_reports_only()
    print('ok')
