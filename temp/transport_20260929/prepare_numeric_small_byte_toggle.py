"""Prepare short-wide numeric cases against the existing scheduler-toggle binary.
No build or Stata invocation. Provenance remains the earlier toggle snapshot.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import shlex
import shutil
import sys

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO/'validation'))
spec = importlib.util.spec_from_file_location('campaign', REPO/'validation/benchmark_transport_campaign.py')
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)
ORIGIN = ROOT/'numeric_scheduler_toggle'
OUT = ROOT/'numeric_small_byte_toggle'
REPS = 11
m = json.loads((ORIGIN/'manifest.json').read_text())
assert m['compiled'] and campaign.sha(Path(m['plugin'])) == m['plugin_sha256']
assert campaign.sha(ORIGIN/'benchmark.c') == m['harness_sha256']
if (OUT/'manifest.json').exists():
    raise ValueError('Already prepared; preserve the recorded provenance')
OUT.mkdir(exist_ok=True)
for name in ('benchmark.c','_ctools_strw.ado'):
    shutil.copy2(ORIGIN/name, OUT/name)
cases = [dict(name=f'n{n}_k{k}_byte', rows=n, columns=k, storage='byte', mode='read')
         for n,k in [(15625,128),(50000,128),(100000,128),(4096,512),(50000,256)]]
cases.append(dict(name='n50000_k128_double', rows=50000, columns=128, storage='double', mode='read'))
cases = [campaign.validate_case(c) for c in cases]
runs=[]
for block in range(2):
    directory = OUT/f'b{block}'
    directory.mkdir(exist_ok=True)
    lines = ['clear all', 'set more off', 'set linesize 255',
             f'log using "{directory}/transport.log", text replace', f'adopath ++ "{OUT}"',
             f'program transport_io, plugin using("{m["plugin"]}")',
             'program run_campaign', 'version 16']
    for case in cases:
        expanded=[]
        for line in campaign.case_driver(case,'columns',REPS).splitlines():
            if line.startswith('plugin call '):
                rep = int(line.split('"')[1].rsplit('_r',1)[1])
                for flag in ((0,1) if (rep+block)%2==0 else (1,0)):
                    call = line if not flag else line.replace('"columns_','"tiles_',1)
                    expanded.append(call+f' "{flag}"')
            else:
                expanded.append(line)
        casefile=directory/(case['name']+'.do')
        casefile.write_text('\n'.join(expanded)+'\n')
        lines.append(f'do "{casefile}"')
    lines += ['end','capture noisily run_campaign','local rc = _rc',
              'di "TRANSPORT_COMPLETE RC=`rc\'"','log close','exit, clear']
    driver=directory/'run.do'
    driver.write_text('\n'.join(lines)+'\n')
    resources=directory/'resources.txt'
    shell='oldstata -q -b do '+shlex.quote(str(driver))
    command=shlex.join(['/usr/bin/time','-l','/bin/zsh','-lic',shell])+' 2>'+shlex.quote(str(resources))
    runs.append(dict(name=f'b{block}',block=block,variant='runtime_toggle',cases=[c['name'] for c in cases],
                     driver=str(driver),log=str(directory/'transport.log'),resources=str(resources),command=command))
summary=ROOT/'summarize_numeric_small_byte_toggle.py'
summary.write_text((ROOT/'summarize_numeric_scheduler_toggle.py').read_text().replace("/'numeric_scheduler_toggle'","/'numeric_small_byte_toggle'"))
m.update(schema='single_plugin_numeric_scheduler_short_wide_v1',cases=cases,runs=runs,repetitions=REPS,
         summarize_script=str(summary),reuse_manifest=str(ORIGIN/'manifest.json'),
         provenance_note='Reuses the original scheduler-toggle source, harness and binary without rebuilding. This is not the final production source hash. K>=128 cases exercise the same numeric gate decision as the corrected production threshold.')
(OUT/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
(OUT/'run-plan.tsv').write_text('run\tcommand\n'+''.join(r['name']+'\t'+r['command']+'\n' for r in runs))
print(f'Prepared {len(cases)} cases, 2 processes, 11 timed pairs each; existing plugin '+m['plugin_sha256'])
print('No build or Stata launch performed.')
