"""Prepare an ABBA benchmark of two frozen ctools build directories.

Compiles a benchmark-only monotonic timer and generates separate Stata drivers.
Invoke the printed stata command yourself; this script never launches Stata.
"""
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
    for label in ('b1','c1','c2','b2'):
        build=a.baseline.resolve() if label.startswith('b') else a.candidate.resolve()
        driver=out/(label+'.do')
        driver.write_text(f'''clear all
set more off
log using "{out}/{label}.log", text replace
capture noisily do "{out}/benchmark.do" "{build}" {label} "{out}" "{' '.join(a.cases)}"
local rc=_rc
di "PERFORMANCE_DRIVER_RC=`rc'"
log close
exit, clear
''')
        print(f"/bin/zsh -lic 'stata -q -b do \"{driver}\"'")

if __name__=='__main__':
    main()
