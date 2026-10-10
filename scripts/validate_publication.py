"""Validate cluster publication files and print their source commit for checkout."""
import json
import math
from pathlib import Path
import re
import sys


def validate(folder):
    versions = json.loads((folder / 'segmeter-full-versions.json').read_text())
    commit = versions['source_commit']
    if not re.fullmatch(r'[0-9a-f]{40}', commit):
        raise ValueError('Invalid source commit')
    if versions['repeats'] != 3 or versions.get('tier') != 'full':
        raise ValueError('Expected a full benchmark with three repeats')
    kind = versions.get('dataset', {}).get('kind', 'simulated')
    if kind not in ('simulated', 'zenodo'):
        raise ValueError('Unrecognized dataset for publication')
    expected_sizes = '10,100,1K,10K,100K' if kind == 'zenodo' else '10,100,1K,10K,100K,1M'
    if versions['sizes'] != expected_sizes:
        raise ValueError('Short validation runs cannot be published as the full benchmark')
    if not versions['images'] or not versions['tools'] or not versions['source_hashes']:
        raise ValueError('Missing provenance')
    if not json.loads((folder / 'segmeter-full-dataset.json').read_text()):
        raise ValueError('Missing dataset hashes')
    results = json.loads((folder / 'segmeter-full-benchmark.json').read_text())
    if not isinstance(results, list) or not results:
        raise ValueError('Missing benchmark results')
    names = set()
    for row in results:
        if not isinstance(row['name'], str) or not row['name'] or row['name'] in names:
            raise ValueError('Invalid or duplicate benchmark name')
        names.add(row['name'])
        if row['unit'] != 'seconds' or not math.isfinite(row['value']) or row['value'] < 0:
            raise ValueError('Invalid full benchmark value')
    return commit


if __name__ == '__main__':
    print(validate(Path(sys.argv[1])))
