"""Prepare an ABBA benchmark of two frozen ctools build directories.

Compiles a benchmark-only monotonic timer and generates separate Stata drivers.
Use --run to execute all four sessions and clean generated datasets on success
or failure. Without --run, print the commands for manual execution instead.
"""
from stata_runner import benchmark_data, launch_driver, stata_command
from summarize_command_performance import CASES, summarize
import shlex
import argparse
import hashlib
import json
import platform
from pathlib import Path
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--baseline',required=True,type=Path)
    p.add_argument('--candidate',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--cases',nargs='+',default=[])
    p.add_argument('--run',action='store_true')
    p.add_argument('--keep-data',action='store_true',
                   help='retain generated comparison datasets after --run')
    a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    mac=platform.system()=='Darwin'
    plugin='ctools_mac_arm.plugin' if mac and platform.machine()=='arm64' else 'ctools_mac_x86.plugin' if mac else 'ctools_linux.plugin'
    manifest={}
    for tag,build in [('b',a.baseline.resolve()),('c',a.candidate.resolve())]:
        files=[build/plugin]+[build/(command+'.ado') for command in ('cipolate','csplit','crangejoin')]
        manifest[tag]={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    flags=['-std=c11','-O2','-I',str(ROOT/'src')]
    flags+=['-bundle','-DSYSTEM=APPLEMAC','-D_DARWIN_C_SOURCE'] if mac else ['-shared','-fPIC','-DSYSTEM=STUNIX','-D_POSIX_C_SOURCE=200809L']
    subprocess.run(['cc',*flags,str(ROOT/'validation/benchmark_clock.c'),str(ROOT/'src/stplugin.c'),'-o',str(out/'clock.plugin')],check=True)
    shutil.copy2(ROOT/'validation/benchmark_command_performance.do',out/'benchmark.do')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Use separate processes: loading two plugins containing static OpenMP
    # runtimes into one Stata session is not a supported comparison method.
    with benchmark_data(out, keep=not a.run or a.keep_data) as data:
        for label in ('b1','c1','c2','b2'):
            build=a.baseline.resolve() if label.startswith('b') else a.candidate.resolve()
            driver=out/(label+'.do')
            driver.write_text(f'''clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{build}" {label} "{out}" "{' '.join(a.cases)}" "{data}"
local rc=_rc
di "PERFORMANCE_DRIVER_RC=`rc'"
log close
exit, clear
''')
            if a.run:
                for path, digest in manifest[label[0]].items():
                    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
                        raise ValueError('Frozen build changed: '+path)
                print('Running '+label,flush=True)
                launch_driver(driver,out/(label+'.log'))
            else:
                print(shlex.join(stata_command(driver)))
        if a.run:
            summarize(out, out/'timings.csv', a.cases or CASES)

if __name__=='__main__':
    main()
