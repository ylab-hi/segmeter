#!/usr/bin/env python3
"""Run reproducible container benchmarks and export github-action-benchmark JSON."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import time

TOOLS = {
    'others': ['bedtools', 'bedtools_sorted', 'tabix', 'bedops', 'bedmaps',
               'bedtk', 'bedtk_sorted', 'igd', 'ailist', 'ucsc', 'awk', 'intervaltree'],
    'giggle': ['giggle'], 'rust-tools': ['gia', 'granges'],
}

EXECUTABLES = {'bedtools_sorted': 'bedtools', 'bedmaps': 'bedmap',
               'bedtk_sorted': 'bedtk', 'ucsc': 'bedIntersect'}


def tool_versions(runtime, container):
    # Older published images have no embedded manifest. Record exact installed
    # binary fingerprints there, rather than claiming the current Dockerfile pins.
    script = '''import hashlib, importlib.metadata, json, pathlib, shutil, sys
p = pathlib.Path('/opt/segmeter-tool-versions.json')
if p.exists():
    print(p.read_text())
else:
    result = {}
    for tool, executable in json.loads(sys.argv[2]).items():
        if tool == 'intervaltree':
            result[tool] = {'version': importlib.metadata.version('intervaltree')}
        else:
            binary = shutil.which(executable)
            if binary is None:
                raise RuntimeError('Missing tool: ' + executable)
            result[tool] = {'binary': binary, 'sha256': hashlib.sha256(pathlib.Path(binary).read_bytes()).hexdigest()}
    print(json.dumps({sys.argv[1]: result}))
'''
    executables = {tool: EXECUTABLES.get(tool, tool) for tool in TOOLS[container]}
    return json.loads(runtime.capture(container, 'python3', '-c', script, container,
                                      json.dumps(executables)))[container]


def run(*args):
    subprocess.run([str(a) for a in args], check=True)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class ContainerRuntime:
    def __init__(self, runtime, images, sif_dir=None):
        self.runtime, self.images = runtime, images
        self.sif_dir = sif_dir

    def image(self, container):
        if self.runtime == 'docker':
            return f'{self.images}:{container}'
        image = self.sif_dir.resolve() / f'{container}.sif'
        if not image.is_file():
            raise FileNotFoundError(f'Missing {image}; pull the images before submitting the job')
        return str(image)

    def command(self, container, *command, source=None, output=None):
        image = self.image(container)
        if self.runtime == 'docker':
            args = ['docker', 'run', '--rm', '--platform', 'linux/amd64']
            if source is not None:
                args += ['-v', f'{source}:/segmeter:ro', '-v', f'{output}:/data']
        else:
            args = [self.runtime, 'exec', '--cleanenv', '--containall', '--pwd', '/tmp']
            if source is not None:
                args += ['--bind', f'{source}:/segmeter:ro', '--bind', f'{output}:/data']
        return [*args, image, *command]

    def capture(self, container, *command):
        return subprocess.check_output(self.command(container, *command), text=True).strip()

    def identity(self, container):
        image = self.image(container)
        if self.runtime == 'docker':
            data = subprocess.check_output(['docker', 'image', 'inspect', image], text=True)
            return {'id': json.loads(data)[0]['Id']}
        return {'path': image, 'sha256': sha256(image)}


def read_results(root, tools, sizes, repeats, subsets):
    """Fail closed on missing data or incorrect answers; retain individual samples."""
    samples = {}
    for repeat in range(repeats):
        for tool in tools:
            for size in sizes.split(','):
                folder = root / 'bench' / f'run-{repeat}' / tool / size
                for subset in subsets:
                    precision = folder / 'precision' / f'{size}_query_precision_{subset}.txt'
                    blocks = precision.read_text().strip().split('\n\n')
                    if len(blocks) != 2:
                        raise ValueError(f'Missing scores: {precision}')
                    for block in blocks:
                        rows = list(csv.DictReader(block.splitlines(), delimiter='\t'))
                        if len(rows) != 1:
                            raise ValueError(f'Missing scores: {precision}')
                        for row in rows:
                            expected_score = 1 if int(row['TP']) else 0
                            if (int(row['FP']) or int(row['FN']) or
                                    float(row['Precision']) != expected_score or
                                    float(row['Recall']) != expected_score or int(row.get('distance', 0)) != 0):
                                raise ValueError(f'Incorrect answers: {precision}: {row}')
                    stats = folder / 'stats' / f'{size}_query_stats_{subset}.txt'
                    with stats.open() as fh:
                        rows = list(csv.DictReader(fh, delimiter='\t'))
                    if len(rows) != 11:
                        raise ValueError(f'Incomplete statistics: {stats}')
                    for row in rows:
                        key = (tool, size, row['data_type'], row['query_type'])
                        samples.setdefault(key, []).append(float(row['time']))
                index = folder / f'{size}_idx_stats.txt'
                if index.exists():
                    with index.open() as fh:
                        row = next(csv.DictReader(fh, delimiter='\t'))
                    samples.setdefault((tool, size, 'index', 'build'), []).append(float(row['time(s)']))
    return samples


def export(samples, relative):
    result = []
    for key, values in sorted(samples.items()):
        tool, size, dtype, query = key
        value = statistics.median(values)
        if relative and dtype == 'index':
            continue  # baseline bedtools has no index
        if relative:
            baseline = statistics.median(samples[('bedtools', size, dtype, query)])
            if baseline <= 0:
                raise ValueError('Baseline time must be positive')
            value /= baseline
        result.append({'name': '/'.join(key), 'unit': 'ratio to bedtools' if relative else 'seconds',
                       'value': value, 'extra': json.dumps({'samples_seconds': values})})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('.'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tier', choices=['smoke', 'full'], default='smoke')
    parser.add_argument('--images', default='segmeter-bench')
    parser.add_argument('--build', action='store_true')
    parser.add_argument('--sizes', help='Override sizes for local validation')
    parser.add_argument('--runtime', choices=['docker', 'singularity', 'apptainer'], default='docker')
    parser.add_argument('--sif-dir', type=Path, help='Directory containing others.sif, giggle.sif and rust-tools.sif')
    parser.add_argument('--dataset-dir', type=Path, help='Use existing DATADIR/sim/sim_001/BED data instead of simulating')
    args = parser.parse_args()
    if args.runtime != 'docker' and (args.sif_dir is None or args.build):
        parser.error('Singularity/Apptainer requires --sif-dir and pre-pulled images; --build is Docker-only')
    source, output = args.source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sizes = args.sizes or ('10K,100K' if args.tier == 'smoke' else
                           ('10,100,1K,10K,100K' if args.dataset_dir else '10,100,1K,10K,100K,1M'))
    repeats = 1 if args.tier == 'smoke' else 3
    subsets = [100] if args.tier == 'smoke' else list(range(10, 101, 10))
    metadata = {'tier': args.tier, 'sizes': sizes, 'seed': 1729, 'max_span': 100, 'repeats': repeats, 'images': {}}
    metadata['dataset'] = {'kind': 'simulated'}
    if args.dataset_dir:
        data_dir = args.dataset_dir.resolve()
        provenance = data_dir / 'zenodo.json'
        metadata.update(seed=None, max_span=None)
        metadata['dataset'] = {'kind': 'zenodo' if provenance.exists() else 'external'}
        if provenance.exists():
            metadata['dataset']['provenance'] = json.loads(provenance.read_text())
        for size in sizes.split(','):
            ref = data_dir / 'sim/sim_001/BED/ref' / f'{size}.bed'
            if not ref.is_file():
                parser.error(f'Dataset does not contain size {size}: {ref}')
    runtime = ContainerRuntime(args.runtime, args.images, args.sif_dir)
    metadata['runtime'] = subprocess.check_output([args.runtime, '--version'], text=True).strip()
    metadata['host'] = {'hostname': platform.node(), 'platform': platform.platform(),
                        'machine': platform.machine(), 'slurm': {
                            k: os.environ[k] for k in ('SLURM_JOB_ID', 'SLURM_JOB_NODELIST',
                                                       'SLURM_CPUS_PER_TASK', 'SLURM_JOB_PARTITION') if k in os.environ}}
    if Path('/proc/cpuinfo').exists():
        (output / 'cpuinfo.txt').write_bytes(Path('/proc/cpuinfo').read_bytes())
    for container in TOOLS:
        image = f'{args.images}:{container}'
        if args.build:
            run('docker', 'build', '--platform', 'linux/amd64', '-t', image, '-f', source / 'containers' / container / 'Dockerfile', source)
        metadata['images'][container] = runtime.identity(container)
        metadata.setdefault('tools', {}).update(tool_versions(runtime, container))
        packages = {'others': ['bedtools', 'tabix', 'mawk'],
                    'giggle': ['tabix', 'libhts-dev'], 'rust-tools': ['bedtools']}[container]
        metadata.setdefault('packages', {})[container] = runtime.capture(container, 'dpkg-query', '-W', *packages)
        # Save this checkout's build recipes separately from the installed image identities.
        (output / f'{container}.Dockerfile').write_bytes((source / 'containers' / container / 'Dockerfile').read_bytes())
    commit_file = source / 'source-commit.txt'
    metadata['source_commit'] = (commit_file.read_text().strip() if commit_file.exists() else
                                 subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip())
    snapshot = output / 'source'
    shutil.copytree(source / 'segmeter', snapshot / 'segmeter', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    metadata['source_hashes'] = {str(p.relative_to(snapshot)): sha256(p)
                                 for p in sorted((snapshot / 'segmeter').rglob('*.py'))}
    metadata['tools']['awk'] = runtime.capture('others', 'dpkg-query', '-W', 'mawk')
    (output / 'versions.json').write_text(json.dumps(metadata, indent=2))

    def execute(container, *command):
        run(*runtime.command(container, 'env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1',
                             'MKL_NUM_THREADS=1', 'python3', '/segmeter/segmeter/main.py', *command,
                             source=snapshot, output=output))

    wall_samples = {}
    if args.dataset_dir:
        shutil.copytree(data_dir / 'sim', output / 'sim', ignore=shutil.ignore_patterns('.DS_Store'))
    else:
        started = time.monotonic()
        execute('others', 'sim', '-o', '/data', '-n', sizes, '--seed', '1729', '--max_span', '100')
        wall_samples[('segmeter', 'all', 'simulation', 'wall')] = [time.monotonic() - started]
    hashes = {str(p.relative_to(output)): sha256(p)
              for p in sorted((output / 'sim').rglob('*')) if p.is_file()}
    (output / 'dataset.json').write_text(json.dumps(hashes, indent=2))
    for repeat in range(repeats):
        for container, tools in TOOLS.items():
            for tool in tools:
                started = time.monotonic()
                execute(container, 'bench', '-r', '-o', '/data', '-n', sizes,
                        '-s', '100' if args.tier == 'smoke' else '10-100',
                        '-b', f'run-{repeat}', '-t', tool)
                wall_samples.setdefault((tool, 'all', 'end-to-end', 'wall'), []).append(time.monotonic() - started)
    samples = read_results(output, sum(TOOLS.values(), []), sizes, repeats, subsets)
    results = export(samples, args.tier == 'smoke') + export(wall_samples, False)
    (output / 'benchmark.json').write_text(json.dumps(results, indent=2))
    if args.tier == 'full':
        publish = output / 'publish'
        publish.mkdir()
        for name in ('benchmark', 'versions', 'dataset'):
            shutil.copyfile(output / f'{name}.json', publish / f'segmeter-full-{name}.json')


if __name__ == '__main__':
    main()
