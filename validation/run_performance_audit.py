"""Freeze and run the 40-column audit in one isolated Stata process per label.

Use matching builds and B-C-C-B order for optimization claims; --summarize
validates every planned trial and reports medians without deleting outliers.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
import re
import shutil
import statistics
import subprocess
from pathlib import Path

from stata_runner import ROOT, run_do

CASES = ('export import sort merge sample bsample encode decode destring ipolate split '
         'rangejoin winsor rangestat binscatter reghdfe ivreghdfe ppmlhdfe psmatch2 qreg '
         'group_numeric group_mixed group_count group_label winsor_by reghdfe_strcluster '
         'sample_by bsample_by').split()


def summarize(out):
    records = []
    checks = {}
    for manifest_path in sorted(out.glob('*.manifest.json')):
        plan = json.loads(manifest_path.read_text())
        label = plan['label']
        log = (out / (label + '.log')).read_text()
        if f'\nAUDIT_COMPLETE label={label}\n' not in log or not re.search(
                r'^CTOOLS_LAUNCH=[0-9a-f]+ RC=0$', log, re.M):
            raise ValueError('incomplete run: ' + label)
        seen = set()
        checked = set()
        for line in log.splitlines():
            if line.startswith('AUDIT,'):
                _, lab, n, case, threads, rep, seconds, rc, rows, columns = line.split(',')
                key = (int(n), case, int(rep))
                if lab != label or key in seen or int(threads) != plan['threads']:
                    raise ValueError('invalid/duplicate record: ' + line)
                seen.add(key)
                elapsed = float(seconds)
                if int(rc) or not math.isfinite(elapsed) or elapsed <= 0:
                    raise ValueError('failed measurement: ' + line)
                records.append(dict(label=lab, observations=int(n), command=case,
                                    threads=int(threads), repetition=int(rep), seconds=elapsed,
                                    output_rows=int(rows), output_columns=int(columns)))
            elif line.startswith('AUDIT_CHECK,'):
                _, lab, n, case, threads, signature, groups = line.split(',')
                key = (int(n), case, int(threads))
                if lab != label or int(threads) != plan['threads'] or key in checked:
                    raise ValueError('invalid/duplicate output check: ' + line)
                checked.add(key)
                value = (signature, groups)
                if key in checks and checks[key] != value:
                    raise ValueError(f'output mismatch {key}: {checks[key]} != {value}')
                checks[key] = value
        expected = {(n, c, r) for n in plan['sizes'] for c in plan['cases']
                    for r in range(plan['reps'] + 1)}
        if seen != expected:
            raise ValueError(f'incomplete grid {label}: {seen ^ expected}')
        expected_checks = {(n, c, plan['threads']) for n in plan['sizes'] for c in plan['cases']}
        if checked != expected_checks:
            raise ValueError(f'incomplete output checks {label}: {checked ^ expected_checks}')
    if not records:
        raise ValueError('no measurements')
    with (out / 'raw.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    groups = {}
    for row in records:
        if row['repetition']:
            key = (row['label'][0], row['observations'], row['command'], row['threads'])
            groups.setdefault(key, []).append(row['seconds'])
    summary = [dict(build=k[0], observations=k[1], command=k[2], threads=k[3], trials=len(v),
                    median=statistics.median(v), minimum=min(v), maximum=max(v))
               for k, v in sorted(groups.items())]
    with (out / 'summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print(f"{row['build']} {row['observations']:>8} {row['command']:20} "
              f"t={row['threads']:2} median={row['median']:.6f} "
              f"range=[{row['minimum']:.6f}, {row['maximum']:.6f}]")
    print(f'PASS: {len(records)} calls, {len(checks)} matching output signatures')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--label', default='b1')
    parser.add_argument('--sizes', type=int, nargs='+', default=[2000000, 5000000])
    parser.add_argument('--threads', type=int, default=12)
    parser.add_argument('--reps', type=int, default=3)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=CASES)
    parser.add_argument('--summarize', action='store_true')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if args.summarize:
        summarize(out)
        return
    if (args.build is None or args.reps < 1 or args.threads < 1 or
            any(n < 1000 or n % 200 for n in args.sizes) or
            not re.fullmatch(r'[a-z][a-z0-9_]*', args.label)):
        parser.error('provide a build, positive threads/repetitions and row counts divisible by 200')
    manifest_path = out / (args.label + '.manifest.json')
    if manifest_path.exists():
        parser.error('label already exists; preserve prior measurements and use a new label')
    mac = platform.system() == 'Darwin'
    plugin = ('ctools_mac_arm.plugin' if platform.machine() == 'arm64' else
              'ctools_mac_x86.plugin') if mac else 'ctools_linux.plugin'
    frozen = out / (args.label + '_build')
    frozen.mkdir()
    for path in [args.build / plugin, *sorted(args.build.glob('*.ado'))]:
        shutil.copy2(path, frozen / path.name)
    harness = out / (args.label + '.do')
    shutil.copy2(ROOT / 'validation/benchmark_audit.do', harness)
    clock = out / 'clock.plugin'
    if not clock.exists():
        flags = ['-bundle', '-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE'] if mac else [
            '-shared', '-fPIC', '-DSYSTEM=STUNIX', '-D_POSIX_C_SOURCE=200809L']
        subprocess.run(['cc', '-std=c11', '-O2', '-I', str(ROOT / 'src'), *flags,
                        str(ROOT / 'validation/benchmark_clock.c'),
                        str(ROOT / 'src/stplugin.c'), '-o', str(clock)], check=True)
    paths = [*sorted(frozen.iterdir()), harness, clock]
    paths += [ROOT / 'validation/run_performance_audit.py', ROOT / 'Makefile']
    manifest = dict(label=args.label, sizes=args.sizes, threads=args.threads, reps=args.reps,
                    cases=args.cases, platform=platform.platform(), machine=platform.machine(),
                    seed=95183, columns=40,
                    files={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    run_do(harness, out / (args.label + '.log'),
           [frozen, out, args.label, ' '.join(map(str, args.sizes)), args.threads,
            args.reps, ' '.join(args.cases)], markers=[f'AUDIT_COMPLETE label={args.label}'])
    print('PASS ' + args.label)


if __name__ == '__main__':
    main()
