"""Stage reports without traversing tool indexes, which may be private to root."""
import argparse
import json
from pathlib import Path
import shutil

from continuous_benchmark import TOOLS


def collect(source, destination):
    target = destination / source.name

    def copy(path):
        if path.is_file():
            out = target / path.relative_to(source)
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, out)

    for name in ('benchmark.json', 'versions.json', 'dataset.json'):
        copy(source / name)
    for container in TOOLS:
        copy(source / f'{container}.Dockerfile')
    versions = source / 'versions.json'
    if not versions.is_file():
        return  # useful partial reports are retained even after an early failure
    metadata = json.loads(versions.read_text())
    for repeat in range(metadata['repeats']):
        for tool in sum(TOOLS.values(), []):
            folder = source / 'bench' / f'run-{repeat}' / tool
            copy(folder / 'log.txt')
            for size in metadata['sizes'].split(','):
                label = folder / size
                copy(label / f'{size}_idx_stats.txt')
                for category in ('stats', 'precision'):
                    # Only these reporting directories are enumerated. Never
                    # walk the tool directory or its sibling *_index / idx/.
                    for report in (label / category).glob('*.txt'):
                        copy(report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('sources', type=Path, nargs='+')
    args = parser.parse_args()
    for source in args.sources:
        collect(source, args.destination)
