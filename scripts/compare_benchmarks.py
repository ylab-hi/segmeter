"""Compare same-runner smoke results, failing at a 25% slowdown."""
import json
from pathlib import Path
import sys

base, head = map(Path, sys.argv[1:])
if json.loads((base / 'dataset.json').read_text()) != json.loads((head / 'dataset.json').read_text()):
    raise SystemExit('Datasets differ; performance comparison is invalid')
before = {r['name']: r['value'] for r in json.loads((base / 'benchmark.json').read_text())}
after = {r['name']: r['value'] for r in json.loads((head / 'benchmark.json').read_text())}
if before.keys() != after.keys():
    raise SystemExit('Benchmark cases differ')
failed = False
print('| Benchmark | Head / base |\n| --- | ---: |')
for name in sorted(before):
    if before[name] <= 0:
        raise SystemExit(f'Nonpositive baseline: {name}')
    ratio = after[name] / before[name]
    print(f'| {name} | {ratio:.3f} |')
    failed |= ratio > 1.25
sys.exit(1 if failed else 0)
