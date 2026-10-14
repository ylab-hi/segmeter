"""Compare same-runner medians, requiring a 25% slowdown beyond sample noise."""
import json
from pathlib import Path
import sys

base, head = map(Path, sys.argv[1:])
if json.loads((base / 'dataset.json').read_text()) != json.loads((head / 'dataset.json').read_text()):
    raise SystemExit('Datasets differ; performance comparison is invalid')
before = {r['name']: r for r in json.loads((base / 'benchmark.json').read_text())}
after = {r['name']: r for r in json.loads((head / 'benchmark.json').read_text())}
if before.keys() != after.keys():
    raise SystemExit('Benchmark cases differ')
failed = False
print('| Benchmark | Head / base | Result |\n| --- | ---: | --- |')
for name in sorted(before):
    if before[name]['value'] <= 0:
        raise SystemExit(f'Nonpositive baseline: {name}')
    ratio = after[name]['value'] / before[name]['value']
    def samples(record):
        return json.loads(record.get('extra', '{}')).get('comparison_samples', [record['value']])
    regression = ratio > 1.25 and min(samples(after[name])) > max(samples(before[name]))
    status = 'regression' if regression else ('noisy; ranges overlap' if ratio > 1.25 else 'pass')
    print(f'| {name} | {ratio:.3f} | {status} |')
    failed |= regression
sys.exit(1 if failed else 0)
