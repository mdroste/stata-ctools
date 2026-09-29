"""Collect completed transport logs into raw timings and paired summaries.

Usage: python3 validation/summarize_transport.py MATRIX_DIRECTORY OUTPUT_PREFIX
The matrix directory contains <case>-<variant>/transport.log. Warmups (r0) and
verification passes are recorded in raw CSV but excluded from summary medians.
"""
import csv
from pathlib import Path
import re
import statistics
import sys
from collections import defaultdict


def summarize(directory, output):
    rows = []
    for log in sorted(directory.glob('*/transport.log')):
        text = log.read_text()
        if '\nTRANSPORT_COMPLETE RC=0\n' not in text:
            continue
        for line in text.splitlines():
            if not line.startswith('TRANSPORT,'):
                continue
            _, label, mode, n, k, threads, load, store, verify = line.split(',')
            match = re.fullmatch(r'(.+)_([^_]+)_w(\d+)_h([01])_r(\d+)', label)
            if not match:
                raise ValueError(label)
            variant, shape, width, hints, rep = match.groups()
            # A directory may hold a single variant or several alternating
            # variants. Its suffix names the run, not necessarily this record.
            group = re.sub(r'-(?:baseline|candidate|control)$', '', log.parent.name)
            rows.append(dict(group=group, variant=variant, shape=shape,
                             width=int(width), hints=int(hints), mode=mode,
                             selected_rows=int(n), columns=int(k), threads=int(threads),
                             repetition=int(rep), verify=int(verify),
                             load_seconds=float(load), store_seconds=float(store)))
    if not rows:
        raise ValueError('No completed transport runs')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_suffix('.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    cases = defaultdict(lambda: defaultdict(list))
    keys = ['group', 'shape', 'width', 'hints', 'mode', 'selected_rows', 'columns', 'threads']
    for row in rows:
        if row['repetition'] and not row['verify']:
            cases[tuple(row[k] for k in keys)][row['variant']].append(row)
    summaries = []
    for key, variants in sorted(cases.items()):
        if not {'baseline','candidate'}.issubset(variants): continue
        result = dict(zip(keys,key))
        for variant in ('baseline','candidate'):
            result[variant+'_repetitions'] = len(variants[variant])
            for phase in ('load','store'):
                vals = [r[phase+'_seconds'] for r in variants[variant]]
                for statistic, fn in [('median',statistics.median),('min',min),('max',max)]:
                    result[f'{variant}_{phase}_{statistic}'] = fn(vals)
        for phase in ('load','store'):
            base = result[f'baseline_{phase}_median']
            new = result[f'candidate_{phase}_median']
            result[phase+'_speedup'] = base/new if new else ''
            baseline = {r['repetition']: r[phase+'_seconds'] for r in variants['baseline']}
            candidate = {r['repetition']: r[phase+'_seconds'] for r in variants['candidate']}
            paired = [baseline[rep]/candidate[rep] for rep in baseline.keys() & candidate.keys()
                      if candidate[rep] > 0]
            result[phase+'_paired_count'] = len(paired)
            for statistic, fn in [('median',statistics.median),('min',min),('max',max)]:
                result[f'{phase}_paired_speedup_{statistic}'] = fn(paired) if paired else ''
        summaries.append(result)
    if summaries:
        with output.with_name(output.name+'_summary').with_suffix('.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
            writer.writeheader(); writer.writerows(summaries)
    print(f'{len(rows)} observations; {len(summaries)} paired cases')


if __name__ == '__main__':
    summarize(Path(sys.argv[1]), Path(sys.argv[2]))
