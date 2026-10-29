"""Validate complete README benchmark logs and retain every measured call."""
import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    p.add_argument('--csv', required=True, type=Path)
    p.add_argument('--evidence', type=Path, help='Compact checks and timings, without license details')
    a = p.parse_args()
    manifest = json.loads((a.directory / 'manifest.json').read_text())
    cases = ['ipolate', 'split', 'rangejoin']
    rows, evidence = [], []
    pattern = re.compile(r'^README_BENCH,(\w+),(\d+),(reference|ctools),(\d+),(\d+),\s*([\d.]+),(\d+)$', re.M)
    for n in manifest['sizes']:
        log = (a.directory / f'run_{n}.log').read_text()
        assert re.search(r'^README_DRIVER_RC=0$', log, re.M), f'{n}: driver failed or unfinished'
        assert re.search(rf'^README_BENCHMARK_COMPLETE N={n}$', log, re.M), f'{n}: incomplete run'
        checks = re.findall(r'^README_CHECK,(\w+),(\d+),PASS,(\d+)$', log, re.M)
        outputs = {case: 2 * n - n // 500 if case == 'rangejoin' else n for case in cases}
        assert sorted(checks) == sorted((case, str(n), str(outputs[case])) for case in cases), f'{n}: missing exact comparisons'
        seen = set()
        for case, size, method, threads, rep, seconds, output in pattern.findall(log):
            size, threads, rep, output = map(int, (size, threads, rep, output))
            seconds = float(seconds)
            assert case in cases and size == n and threads == manifest['threads']
            assert 1 <= rep <= manifest['reps'] and output == outputs[case]
            assert math.isfinite(seconds) and seconds > 0
            key = (case, method, rep)
            assert key not in seen, f'duplicate {n}: {key}'
            seen.add(key)
            rows.append((case, n, method, threads, rep, seconds, output))
        expected = {(case, method, rep) for case in cases for method in ['reference', 'ctools'] for rep in range(1, manifest['reps'] + 1)}
        assert seen == expected, f'{n}: missing timing records'
        evidence.extend(line for line in log.splitlines() if line.startswith(('README_CHECK,', 'README_BENCH,', 'README_BENCHMARK_COMPLETE ', 'README_DRIVER_RC=')))
    with a.csv.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['command', 'input_rows_per_dataset', 'implementation', 'threads', 'repetition', 'seconds', 'output_rows'])
        writer.writerows(rows)
    if a.evidence:
        a.evidence.write_text('\n'.join(evidence) + '\n')
    print('| Command | Input rows per dataset | Reference (s) | ctools (s) | Speedup |')
    print('|---|---:|---:|---:|---:|')
    for case in cases:
        for n in manifest['sizes']:
            times = {method: statistics.median(row[5] for row in rows if row[:3] == (case, n, method)) for method in ['reference', 'ctools']}
            print(f"| `{case}` | {n:,} | {times['reference']:.3f} | {times['ctools']:.3f} | {times['reference'] / times['ctools']:.1f}× |")
    print(f'Validated {len(rows)} timings and {len(cases) * len(manifest["sizes"])} exact full-output comparisons.')


if __name__ == '__main__':
    main()
