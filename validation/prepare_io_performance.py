"""Prepare/run an ABBA comparison of frozen cimport/cexport build directories.

Uses the repository Stata launcher (interactive login zsh's stata alias), a
monotonic plugin clock, and fresh completion nonces. Summaries exclude warmups.
CSV output bytes and decompressed XLSX parts must match between both builds;
Stata's cf checks imported values on the same input files.
With --run, generated data are disposable; logs, manifests and timing tables
remain in --output. Use --keep-data to retain data for later --summarize-only.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import zipfile
from stata_runner import benchmark_data, launch_driver, stata_command

ROOT = Path(__file__).resolve().parents[1]
LABELS = ('b1', 'c1', 'c2', 'b2')


def summarize(out, data=None):
    data = out if data is None else data
    (out/'comparisons.txt').unlink(missing_ok=True)
    measurements = {}
    records = []
    for label in LABELS:
        text = (out / (label + '.log')).read_text()
        if text.splitlines().count('IO_BENCHMARK_RC=0') != 1:
            raise ValueError('incomplete benchmark: ' + label)
        for line in text.splitlines():
            if not line.startswith('PERF,'):
                continue
            _, tag, case, op, threads, rep, seconds = line.split(',')
            if tag != label:
                raise ValueError('unexpected benchmark label')
            record = (label, case, op, int(threads), int(rep), float(seconds))
            records.append(record)
            if int(rep):
                key = (case, op, int(threads))
                measurements.setdefault(key, {}).setdefault(label[0], []).append(float(seconds))
    with (out/'timings.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(('run','case','operation','threads','rep','seconds'))
        writer.writerows(records)
    summary = []
    for (case, op, threads), samples in measurements.items():
        if not samples.get('b') or len(samples['b']) != len(samples.get('c', [])):
            raise ValueError('unbalanced samples')
        b, c = (statistics.median(samples[tag]) for tag in ('b','c'))
        summary.append(dict(case=case, operation=op, threads=threads,
                            baseline_seconds=b, candidate_seconds=c, speedup=b/c,
                            baseline_min=min(samples['b']), baseline_max=max(samples['b']),
                            candidate_min=min(samples['c']), candidate_max=max(samples['c']),
                            samples_per_build=len(samples['b'])))
    if not summary:
        raise ValueError('no measurements')
    with (out/'summary.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader(); writer.writerows(summary)
    checked = 0
    for case, op, threads in measurements:
        if not op.startswith('export_'):
            continue
        ext = 'csv' if op == 'export_csv' else 'xlsx'
        expected = data/f'b1_{case}_{threads}.{ext}'
        for label in ('c1','c2','b2'):
            actual = data/f'{label}_{case}_{threads}.{ext}'
            if ext == 'csv':
                if hashlib.sha256(expected.read_bytes()).digest() != hashlib.sha256(actual.read_bytes()).digest():
                    raise ValueError('CSV bytes differ: ' + str(actual))
            else:
                with zipfile.ZipFile(expected) as b, zipfile.ZipFile(actual) as c:
                    if b.namelist() != c.namelist() or c.testzip() is not None:
                        raise ValueError('invalid workbook: ' + str(actual))
                    for name in b.namelist():
                        if b.read(name) != c.read(name):
                            raise ValueError(f'XLSX part differs: {actual}: {name}')
            checked += 1
    result = f'PASS {checked} exact export comparisons; imported values checked by Stata cf'
    (out/'comparisons.txt').write_text(result+'\n')
    print(result)
    for row in summary:
        print(f"{row['case']:8} {row['operation']:12} t={row['threads']:2} "
              f"{row['baseline_seconds']:.6f} -> {row['candidate_seconds']:.6f} s "
              f"({row['speedup']:.2f}x)")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--candidate',type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--cases',nargs='+',choices=['numeric','mixed','wide','strings','small','large'],
                   default=['numeric','mixed','wide','strings'])
    p.add_argument('--reps',type=int,default=3)
    p.add_argument('--run',action='store_true')
    p.add_argument('--summarize-only',action='store_true')
    p.add_argument('--keep-data',action='store_true',
                   help='retain generated data after --run for debugging or --summarize-only')
    a = p.parse_args(); out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    if a.summarize_only:
        summarize(out);return
    if a.baseline is None or a.candidate is None or a.reps < 1:
        p.error('provide --baseline, --candidate, and positive --reps')
    (out/'comparisons.txt').unlink(missing_ok=True)
    mac=platform.system()=='Darwin'
    plugin=('ctools_mac_arm.plugin' if platform.machine()=='arm64' else 'ctools_mac_x86.plugin') if mac else 'ctools_linux.plugin'
    builds={'b':a.baseline.resolve(),'c':a.candidate.resolve()}
    manifest={}
    for tag, build in builds.items():
        manifest[tag]={str(f):hashlib.sha256(f.read_bytes()).hexdigest()
                       for f in [build/plugin,*sorted(build.glob('*.ado'))]}
    flags=['-std=c11','-O2','-I',str(ROOT/'src')]
    flags+=['-bundle','-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE'] if mac else ['-shared','-fPIC','-DSYSTEM=STUNIX','-D_POSIX_C_SOURCE=200809L']
    subprocess.run(['cc',*flags,str(ROOT/'validation/benchmark_clock.c'),str(ROOT/'src/stplugin.c'),'-o',str(out/'clock.plugin')],check=True)
    shutil.copy2(ROOT/'validation/benchmark_io_performance.do',out/'benchmark.do')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Prepared-only drivers keep stable paths; executed runs own a scratch directory.
    with benchmark_data(out, keep=not a.run or a.keep_data) as data:
        for label in LABELS:
            driver=out/(label+'.do')
            driver.write_text(f'''clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{builds[label[0]]}" {label} "{out}" "{' '.join(a.cases)}" {a.reps} "{data}"
local rc=_rc
di "IO_BENCHMARK_RC=`rc'"
log close
exit, clear
''')
            if a.run:
                for path, digest in manifest[label[0]].items():
                    if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:
                        raise ValueError('frozen build changed: '+path)
                print('Running '+label,flush=True)
                launch_driver(driver,out/(label+'.log'))
            else:
                print(shlex.join(stata_command(driver)))
        if a.run:
            summarize(out, data)


if __name__=='__main__':
    main()
