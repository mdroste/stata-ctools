"""Summarize the single-binary numeric-scheduler experiment; require both complete blocks."""
import csv
import json
from pathlib import Path
import statistics

OUT = Path(__file__).resolve().parent/'numeric_sixm_toggle'
m = json.loads((OUT/'manifest.json').read_text())
raw, processes, summaries = [], [], []
for run in m['runs']:
    log = Path(run['log'])
    text = log.read_text() if log.exists() else ''
    if '\nTRANSPORT_COMPLETE RC=0\n' not in text:
        raise ValueError(f'Incomplete block: {run["name"]}')
    times, cleanups, orders = {}, {}, {}
    for line in text.splitlines():
        if line.startswith('TRANSPORT,'):
            _, label, mode, n, k, threads, load, store, verify = line.split(',')
            if label in times:
                raise ValueError(f'Duplicate timing: {label}')
            times[label] = dict(mode=mode, rows=int(n), columns=int(k), threads=int(threads),
                                load=float(load), store=float(store), verify=int(verify))
        elif line.startswith('TRANSPORT_CLEANUP,'):
            _, label, seconds = line.split(',')
            if label in cleanups:
                raise ValueError(f'Duplicate cleanup: {label}')
            cleanups[label] = float(seconds)
        elif line.startswith('TRANSPORT_ORDER_BYTES,'):
            _, label, size = line.split(',')
            if label in orders:
                raise ValueError(f'Duplicate order allocation: {label}')
            orders[label] = int(size)
    expected = {f'{v}_{c["name"]}_r{rep}' for v in ('columns','tiles') for c in m['cases'] for rep in range(m['repetitions']+2)}
    if any(set(d) != expected for d in (times, cleanups, orders)):
        raise ValueError(f'Missing or unexpected records: {run["name"]}')
    for case in m['cases']:
        for variant in ('columns', 'tiles'):
            rows = []
            for rep in range(m['repetitions']+2):
                label = f'{variant}_{case["name"]}_r{rep}'
                item = times[label]
                if (item['rows'], item['columns'], item['mode']) != (case['rows'], case['columns'], 'read'):
                    raise ValueError(f'Case metadata mismatch: {label}')
                if item['verify'] != int(rep == m['repetitions']+1):
                    raise ValueError(f'Verification marker mismatch: {label}')
                expected_bytes = item['rows']*4
                if orders[label] != expected_bytes:
                    raise ValueError(f'Unexpected identity-order allocation: {label}, {orders[label]} != {expected_bytes}')
                row = dict(block=run['block'], case=case['name'], variant=variant, repetition=rep,
                           **item, cleanup=cleanups[label], sort_order_bytes=orders[label])
                row['transport'] = row['load']+row['cleanup']
                raw.append(row)
                if rep and not item['verify']:
                    rows.append(row)
            processes.append(dict(block=run['block'], case=case['name'], variant=variant,
                                  **{phase:statistics.median(r[phase] for r in rows) for phase in ('load','cleanup','transport')}))
for case in m['cases']:
    for phase in ('load','cleanup','transport'):
        default = {x['block']:x[phase] for x in processes if x['case']==case['name'] and x['variant']=='columns'}
        noorder = {x['block']:x[phase] for x in processes if x['case']==case['name'] and x['variant']=='tiles'}
        ratios = [default[b]/noorder[b] for b in default if noorder[b] > 0]
        summaries.append(dict(case=case['name'], phase=phase, process_pairs=len(ratios),
                              columns_median=statistics.median(default.values()), tiles_median=statistics.median(noorder.values()),
                              paired_speedup_median=statistics.median(ratios) if ratios else '',
                              paired_speedup_min=min(ratios) if ratios else '', paired_speedup_max=max(ratios) if ratios else '',
                              columns_min=min(default.values()), columns_max=max(default.values()),
                              tiles_min=min(noorder.values()), tiles_max=max(noorder.values())))
for filename, rows in [('raw.csv',raw),('process.csv',processes),('summary.csv',summaries)]:
    with (OUT/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(f'{len(raw)} records; {len(processes)} process medians; {len(summaries)} phase comparisons')
