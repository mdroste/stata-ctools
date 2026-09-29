"""Run frozen whole-command comparisons, exclusively through oldstata."""
from pathlib import Path
import csv, hashlib, json, shlex, shutil, statistics, subprocess, time

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
OUT = ROOT / 'commands_final'
DRIVER = REPO / 'validation/benchmark_transport_commands.do'
OUT.mkdir(exist_ok=True)
variants = {'baseline': ROOT / 'production_final_baseline_build', 'candidate': ROOT / 'production_final/build'}
manifest = dict(driver=str(DRIVER), driver_sha256=hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
                repetitions=7, blocks=2, threads=8, variants={}, runs=[])
for name, directory in variants.items():
    manifest['variants'][name] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in directory.iterdir() if p.suffix in ('.ado', '.plugin')}
for block in range(2):
    for variant in (('baseline', 'candidate') if block == 0 else ('candidate', 'baseline')):
        label = f'b{block}_{variant}'
        directory = OUT / label
        directory.mkdir(exist_ok=True)
        shutil.copy2(OUT / 'clock.plugin', directory / 'clock.plugin')
        arguments = [str(DRIVER), str(variants[variant]), str(directory), label, '7', '8']
        shell = 'oldstata -q -b do ' + shlex.join(arguments)
        manifest['runs'].append(dict(block=block, variant=variant, label=label,
                                    directory=str(directory), shell=shell))
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
raw = []
for run in manifest['runs']:
    directory = Path(run['directory'])
    log = directory / ('commands_' + run['label'] + '.log')
    if not log.exists() or '\nADAPTIVE_COMMANDS_COMPLETE RC=0\n' not in log.read_text():
        start = time.monotonic()
        print('COMMAND_START', run['label'], flush=True)
        with (directory / 'resources.txt').open('w') as output:
            result = subprocess.run(['/usr/bin/time', '-l', '/bin/zsh', '-lic', run['shell']],
                                    cwd=directory, stdout=output, stderr=subprocess.STDOUT)
        (directory / 'execution.json').write_text(json.dumps(dict(exit_code=result.returncode,
                     elapsed_monotonic=time.monotonic()-start), indent=2) + '\n')
        if result.returncode:
            raise SystemExit(f'Wrapper failed: {run["label"]}')
    text = log.read_text() if log.exists() else ''
    if '\nADAPTIVE_COMMANDS_COMPLETE RC=0\n' not in text:
        print(text[-4000:])
        raise SystemExit(f'Command validation failed: {run["label"]}')
    records = []
    for line in text.splitlines():
        if not line.startswith('COMMANDCSV,'):
            continue
        _, label, case, rows, columns, threads, rep, elapsed, load, store = line.split(',')
        assert label == run['label']
        records.append(dict(block=run['block'], variant=run['variant'], case=case,
                            rows=int(rows), columns=int(columns), threads=int(threads), rep=int(rep),
                            elapsed=float(elapsed), load=float(load), store=float(store)))
    assert len(records) == 32
    assert len({(r['case'], r['rep']) for r in records}) == 32
    assert all({r['rep'] for r in records if r['case'] == case} == set(range(8))
               for case in ('tiny_numeric', 'narrow_numeric', 'wide_numeric', 'str2045_full'))
    raw.extend(records)
    print('COMMAND_COMPLETE', run['label'], flush=True)

def save(name, records):
    with (OUT / name).open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)

process = []
for run in manifest['runs']:
    for case in ('tiny_numeric', 'narrow_numeric', 'wide_numeric', 'str2045_full'):
        rows = [r for r in raw if r['block'] == run['block'] and r['variant'] == run['variant']
                and r['case'] == case and r['rep'] > 0]
        process.append(dict(block=run['block'], variant=run['variant'], case=case,
                            **{phase: statistics.median(r[phase] for r in rows)
                               for phase in ('elapsed', 'load', 'store')}))
summary = []
for case in ('tiny_numeric', 'narrow_numeric', 'wide_numeric', 'str2045_full'):
    for phase in ('elapsed', 'load', 'store'):
        base = [r[phase] for r in process if r['case'] == case and r['variant'] == 'baseline']
        candidate = [r[phase] for r in process if r['case'] == case and r['variant'] == 'candidate']
        ratios = [b/c for b,c in zip(base,candidate) if c > 0]
        summary.append(dict(case=case, phase=phase, baseline_median=statistics.median(base),
                            candidate_median=statistics.median(candidate),
                            paired_speedup_median=statistics.median(ratios),
                            paired_speedup_min=min(ratios), paired_speedup_max=max(ratios)))
save('raw.csv', raw); save('process.csv', process); save('summary.csv', summary)
print('COMMANDS_VERIFIED', len(raw), 'records', flush=True)
