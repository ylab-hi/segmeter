"""Cluster preparation tests without Slurm, Singularity or network access."""
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / 'scripts/cluster'))
from prepare_zenodo import extract, verify


def test_archive():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        archive = root / 'data.tar.gz'
        def write_archive(name):
            with tarfile.open(archive, 'w:gz') as tar:
                member = tarfile.TarInfo(name)
                member.size = 4
                tar.addfile(member, io.BytesIO(b'data'))
        name = 'simdata/sim/sim_001/BED/ref/1K.bed'
        write_archive(name)
        extract(archive, root / 'extracted')
        assert (root / 'extracted' / name).read_bytes() == b'data'
        for unsafe in ('../outside', '/absolute'):
            write_archive(unsafe)
            try:
                extract(archive, root / 'bad')
            except ValueError:
                pass
            else:
                raise AssertionError('Unsafe archive accepted')
        archive.write_bytes(b'data')
        with patch('prepare_zenodo.CHECKSUM', '8d777f385d3dfec8815d20f7496026dc'):
            verify(archive)
            archive.write_bytes(b'changed')
            try:
                verify(archive)
            except ValueError:
                pass
            else:
                raise AssertionError('Checksum mismatch accepted')


def test_launcher():
    with tempfile.TemporaryDirectory(prefix='segmeter cluster ') as tmp:
        root = Path(tmp)
        bin_dir = root / 'bin'
        bin_dir.mkdir()
        # Stub only external services. Run the actual Bash launcher, argument
        # parser, cache handling and checkout snapshot code.
        program = '''import json, os, pathlib, sys
name = pathlib.Path(sys.argv[0]).name
with open(os.environ['CLUSTER_TEST_LOG'], 'a') as fh:
    fh.write(json.dumps([name, *sys.argv[1:]]) + '\\n')
if name == 'singularity':
    pathlib.Path(sys.argv[3]).write_bytes(b'fake SIF')
elif name == 'sbatch':
    print('4242')
elif name == 'flock':
    pass
elif name == 'python3' and sys.argv[1].endswith('prepare_zenodo.py'):
    root = pathlib.Path(sys.argv[3]) / 'simdata'
    root.mkdir(parents=True, exist_ok=True)
    (root / 'zenodo.json').write_text('{}')
else:
    os.execv(os.environ['CLUSTER_REAL_PYTHON'], [os.environ['CLUSTER_REAL_PYTHON'], *sys.argv[1:]])
'''
        for name in ('singularity', 'sbatch', 'flock', 'python3'):
            file = bin_dir / name
            file.write_text(f'#!{sys.executable}\n' + program)
            file.chmod(0o755)
        initialization = root / 'modules.sh'
        initialization.write_text('module() { printf "%s\\n" "$*" >> "$CLUSTER_MODULE_LOG"; }\n')
        log = root / 'commands.jsonl'
        env = {**os.environ, 'PATH': str(bin_dir) + ':' + os.environ['PATH'],
               'BASH_ENV': str(initialization), 'CLUSTER_TEST_LOG': str(log),
               'CLUSTER_MODULE_LOG': str(root / 'modules.log'), 'CLUSTER_REAL_PYTHON': sys.executable}
        command = ['bash', str(ROOT / 'scripts/cluster/run-benchmark.sh'),
                   '--module', 'singularity/test', '--image-tag', 'v0.14.1',
                   '--work-dir', str(root / 'work'), '--dataset', 'zenodo',
                   '--account', 'lab', '--partition', 'cpu', '--sizes', '1K']
        for _ in range(2):
            result = subprocess.run(command, env=env, check=True, capture_output=True, text=True)
            assert 'Submitted Slurm job 4242' in result.stdout
        calls = [json.loads(line) for line in log.read_text().splitlines()]
        assert len([c for c in calls if c[0] == 'singularity']) == 3  # cached on second launch
        jobs = [c for c in calls if c[0] == 'sbatch']
        assert len(jobs) == 2
        for job in jobs:
            assert '--account=lab' in job and '--partition=cpu' in job
            assert '--dataset-dir' in job and job[-2:] == ['--sizes', '1K']
        snapshots = list((root / 'work/results').glob('*/checkout'))
        assert len(snapshots) == 2
        assert all((p / 'source-commit.txt').is_file() and (p / 'segmeter/main.py').is_file() for p in snapshots)
        assert (root / 'modules.log').read_text().splitlines() == ['load singularity/test'] * 2


if __name__ == '__main__':
    test_archive()
    test_launcher()
    print('ok')
