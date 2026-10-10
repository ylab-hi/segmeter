"""Exporter and correctness gate tests; no external tools needed."""
import json
import subprocess
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parents[2] / 'scripts'))
from continuous_benchmark import ContainerRuntime, export, read_results, sha256
from validate_publication import validate


def test_export():
    data = {('bedtools', '10K', 'basic', 'perfect_100%'): [1, 2, 3],
            ('gia', '10K', 'basic', 'perfect_100%'): [2, 4, 6],
            ('gia', '10K', 'index', 'build'): [3, 4, 5]}
    relative = export(data, True)
    assert [r['value'] for r in relative] == [1, 2]
    absolute = export(data, False)
    assert [r['value'] for r in absolute] == [2, 4, 4]
    assert json.loads(relative[1]['extra'])['samples_seconds'] == [2, 4, 6]


def test_gate():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        folder = root / 'bench/run-0/bedtools/10K'
        (folder / 'precision').mkdir(parents=True)
        (folder / 'stats').mkdir()
        precision = folder / 'precision/10K_query_precision_100.txt'
        valid = ('intvlnum\tsubset\tTP\tFP\tTN\tFN\tPrecision\tRecall\tF1\n'
                 '10000\t100%\t5\t0\t5\t0\t1\t1\t1\n\n'
                 'intvlnum\tbin\tTP\tFP\tFN\tPrecision\tRecall\tF1\tdistance\n'
                 '10000\t100bin\t5\t0\t0\t1\t1\t1\t0\n')
        precision.write_text(valid)
        stats = folder / 'stats/10K_query_stats_100.txt'
        stats.write_text('intvlnum\tdata_type\tquery_type\ttime\tmax_RSS(MB)\n' +
                         ''.join(f'10000\tbasic\tq{i}\t0.5\t1\n' for i in range(11)))
        assert len(read_results(root, ['bedtools'], '10K', 1, [100])) == 11
        for invalid in [valid.replace('\t1\t1\t1\t0\n', '\t1\t1\t1\t2\n'),
                        valid.replace('\t5\t0\t5\t0\t', '\t5\t1\t5\t0\t'),
                        valid.split('\n\n')[0]]:
            precision.write_text(invalid)
            try:
                read_results(root, ['bedtools'], '10K', 1, [100])
            except ValueError:
                pass
            else:
                raise AssertionError('Invalid results accepted')


def test_comparison():
    script = Path(__file__).parents[2] / 'scripts/compare_benchmarks.py'
    with tempfile.TemporaryDirectory() as tmp:
        base, head = Path(tmp) / 'base', Path(tmp) / 'head'
        for folder in (base, head):
            folder.mkdir()
            (folder / 'dataset.json').write_text('{"ref": "hash"}')
            (folder / 'benchmark.json').write_text('[{"name": "case", "value": 1}]')
        def compare():
            return subprocess.run([sys.executable, str(script), str(base), str(head)], capture_output=True).returncode
        assert compare() == 0
        (head / 'benchmark.json').write_text('[{"name": "case", "value": 1.3}]')
        assert compare() == 1
        (head / 'benchmark.json').write_text('[{"name": "other", "value": 1}]')
        assert compare() == 1
        (head / 'dataset.json').write_text('{"ref": "changed"}')
        assert compare() == 1


def test_container_runtimes():
    with tempfile.TemporaryDirectory(prefix='segmeter images ') as tmp:
        directory = Path(tmp).resolve()
        image = directory / 'others.sif'
        image.write_bytes(b'test image')
        source, output = directory / 'source', directory / 'results'
        for executable in ('singularity', 'apptainer'):
            runtime = ContainerRuntime(executable, 'unused', directory)
            command = runtime.command('others', 'python3', '/segmeter/segmeter/main.py',
                                      source=source, output=output)
            assert command == [executable, 'exec', '--cleanenv', '--containall', '--pwd', '/tmp',
                               '--bind', f'{source}:/segmeter:ro', '--bind', f'{output}:/data',
                               str(image), 'python3', '/segmeter/segmeter/main.py']
            with patch('continuous_benchmark.subprocess.check_output') as capture:
                assert runtime.identity('others') == {'path': str(image), 'sha256': sha256(image)}
                capture.assert_not_called()  # no Docker inspection for SIF files
            try:
                runtime.image('giggle')
            except FileNotFoundError:
                pass
            else:
                raise AssertionError('Missing SIF accepted')
        docker = ContainerRuntime('docker', 'segmeter-bench')
        command = docker.command('others', 'true', source=source, output=output)
        assert command[-2:] == ['segmeter-bench:others', 'true']
        assert command[command.index('-v') + 1] == f'{source}:/segmeter:ro'


def test_publication():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        versions = {'source_commit': 'a' * 40, 'repeats': 3, 'tier': 'full',
                    'sizes': '10,100,1K,10K,100K,1M', 'images': {'others': 'hash'},
                    'tools': {'bedtools': '2.31.1'}, 'source_hashes': {'main.py': 'hash'}}
        (root / 'segmeter-full-versions.json').write_text(json.dumps(versions))
        (root / 'segmeter-full-dataset.json').write_text('{"ref": "hash"}')
        result_file = root / 'segmeter-full-benchmark.json'
        result_file.write_text('[{"name": "case", "unit": "seconds", "value": 1}]')
        assert validate(root) == 'a' * 40
        for results in ([], [{'name': 'case', 'unit': 'seconds', 'value': float('nan')}],
                        [{'name': 'case', 'unit': 'ratio to bedtools', 'value': 1}]):
            result_file.write_text(json.dumps(results))
            try:
                validate(root)
            except ValueError:
                pass
            else:
                raise AssertionError('Invalid publication accepted')
        versions['sizes'] = '1K'
        result_file.write_text('[{"name": "case", "unit": "seconds", "value": 1}]')
        (root / 'segmeter-full-versions.json').write_text(json.dumps(versions))
        try:
            validate(root)
        except ValueError:
            pass
        else:
            raise AssertionError('Short validation run accepted for publication')
        versions.update(sizes='10,100,1K,10K,100K', dataset={'kind': 'zenodo'})
        (root / 'segmeter-full-versions.json').write_text(json.dumps(versions))
        assert validate(root) == 'a' * 40


if __name__ == '__main__':
    test_export()
    test_gate()
    test_comparison()
    test_container_runtimes()
    test_publication()
    print('ok')
