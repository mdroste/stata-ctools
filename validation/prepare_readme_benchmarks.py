"""Prepare native/SSC benchmarks; print stata commands, never launch Stata.

Pass a frozen, compiled source tree to --snapshot to isolate concurrent edits.
The manifest hashes sources, wrappers, references and the actual loaded plugin.
"""
from stata_runner import stata_command
import shlex
import argparse
import hashlib
import json
import platform
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot', required=True, type=Path)
    p.add_argument('--reference', required=True, type=Path)
    p.add_argument('--reference-dependency', action='append', type=Path, default=[],
                   help='Additional reference source to hash, e.g. rangestat.ado (repeatable)')
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--sizes', nargs='+', type=int, default=[1000000, 5000000, 20000000])
    p.add_argument('--threads', type=int, default=12)
    p.add_argument('--reps', type=int, default=4)
    a = p.parse_args()
    if a.threads < 1 or a.reps < 2 or a.reps % 2 or any(n < 1000 or n % 1000 for n in a.sizes):
        p.error('Use positive threads, an even repeat count >=2, and sizes divisible by 1000.')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    snapshot = a.snapshot.resolve()
    build = snapshot / 'build'
    refs = a.reference.resolve()
    mac = platform.system() == 'Darwin'
    plugin = ('ctools_mac_arm.plugin' if platform.machine() == 'arm64' else 'ctools_mac_x86.plugin') if mac else 'ctools_linux.plugin'
    if not (build / plugin).is_file():
        p.error('Build the frozen snapshot first.')
    flags = ['-std=c11', '-O2', '-I', str(snapshot / 'src')]
    flags += ['-bundle', '-DSYSTEM=APPLEMAC', '-D_DARWIN_C_SOURCE', '-mmacosx-version-min=11.0'] if mac else ['-shared', '-fPIC', '-DSYSTEM=STUNIX', '-D_POSIX_C_SOURCE=200809L']
    subprocess.run(['cc', *flags, str(ROOT / 'validation/benchmark_clock.c'), str(snapshot / 'src/stplugin.c'), '-o', str(out / 'clock.plugin')], check=True)
    shutil.copy2(ROOT / 'validation/benchmark_readme_commands.do', out / 'benchmark.do')
    files = sorted(f for f in (snapshot / 'src').rglob('*') if f.is_file())
    files += sorted(build.glob('*.ado')) + [build / plugin, snapshot / 'Makefile']
    files += sorted(refs.glob('*.ado')) + [out / 'clock.plugin', out / 'benchmark.do']
    files += [f.resolve() for f in a.reference_dependency]
    files += [ROOT / 'validation' / name for name in
              ['benchmark_clock.c', 'prepare_readme_benchmarks.py', 'summarize_readme_benchmarks.py']]
    for path in ['/Applications/Stata/ado/base/i/ipolate.ado', '/Applications/Stata/ado/base/s/split.ado']:
        if Path(path).exists():
            files.append(Path(path))
    hashes = {str(f): hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    source_identity = '\n'.join(f'{f.relative_to(snapshot)} {hashes[str(f)]}'
                                for f in sorted((snapshot / 'src').rglob('*')) if f.is_file())
    manifest = dict(platform=platform.platform(), machine=platform.machine(), sizes=a.sizes,
                    threads=a.threads, reps=a.reps, seed=95183,
                    source_tree_sha256=hashlib.sha256(source_identity.encode()).hexdigest(), files=hashes)
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    for n in a.sizes:
        driver = out / f'run_{n}.do'
        driver.write_text(f'''clear all
set more off
log using "{out}/run_{n}.log", text replace
capture noisily do "{out}/benchmark.do" "{build}" "{refs}" "{out}" {n} {a.threads} {a.reps}
local rc=_rc
di "README_DRIVER_RC=`rc'"
log close
exit, clear
''')
        print(shlex.join(stata_command(driver)))


if __name__ == '__main__':
    main()
