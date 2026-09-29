"""Refresh completed-process campaign tables; never launch Stata or interpret results."""
import csv
import importlib.util
import json
from pathlib import Path
import types
import sys

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT.parents[1] / 'validation/benchmark_transport_campaign.py'
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location('transport_campaign', SCRIPT)
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)
index, raw, phases, resources = [], [], [], []
COMPLETED_CAMPAIGNS = (
    'tiny_campaign', 'tiny_weighted_screen', 'numeric_screen', 'str2045_screen',
    'numeric_tile_sweep_full', 'store_screen', 'runtime_focus', 'no_sort_order_focused',
)
for name in COMPLETED_CAMPAIGNS:
    path = ROOT / name / 'manifest.json'
    manifest = json.loads(path.read_text())
    if not {'runs', 'cases', 'variants', 'repetitions'}.issubset(manifest):
        continue
    directory = path.parent
    completed = []
    failed = []
    for run in manifest['runs']:
        log = Path(run['log'])
        content = log.read_text() if log.exists() else ''
        if '\nTRANSPORT_COMPLETE RC=0\n' in content:
            completed.append(run)
        elif '\nTRANSPORT_COMPLETE RC=' in content:
            failed.append(run['name'])
    index.append(dict(campaign=directory.name, cases=len(manifest['cases']),
                      variants=len(manifest['variants']), expected_processes=len(manifest['runs']),
                      completed_processes=len(completed), failed_processes=len(failed),
                      lifecycle=manifest['lifecycle'], no_sort_order=manifest.get('no_sort_order', False),
                      manifest=str(path)))
    if len(completed) != len(manifest['runs']):
        raise ValueError(f'Whitelisted campaign is not complete: {directory.name}')
    campaign.summarize(types.SimpleNamespace(directory=directory, baseline=None, allow_incomplete=False))
    cases = {c['name']: c for c in manifest['cases']}
    for filename, collection in [('raw.csv', raw), ('summary.csv', phases), ('resources.csv', resources)]:
        source = directory / filename
        if not source.exists():
            continue
        for item in csv.DictReader(source.open()):
            row = dict(campaign=directory.name, **item)
            if 'case' in item:
                case = cases[item['case']]
                for key in ('rows', 'host_columns', 'shape', 'storage', 'width', 'actual_width',
                            'hints', 'layout', 'threads', 'mode', 'predicate', 'obs_range'):
                    row['case_' + key] = case.get(key, '')
            collection.append(row)
for name, rows in [('campaign_index.csv', index), ('campaign_raw.csv', raw),
                   ('campaign_phase_summary.csv', phases), ('campaign_resources.csv', resources)]:
    if not rows:
        continue
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (ROOT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(name, len(rows))
