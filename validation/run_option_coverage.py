"""Run every explicit option comparison, retaining per-case failures and coverage."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tempfile
from option_cases import CASES
from stata_runner import ROOT, launch_driver

def quote(s):
    return '`"' + s + '"\''

def source(cases, log):
    lines=[f'cd "{ROOT}"', 'clear all', 'set more off', 'set linesize 255',
           'capture log close _all', f'log using "{log}", text replace',
           'capture noisily {', 'quietly do validation/validate_setup.do', 'quietly do validation/coverage_helpers.do',
           'local failed_cases=0', 'di "REFERENCE_STATA_VERSION=`c(stata_version)\' OS=`c(os)\'"']
    for cmd in ['qreg','reghdfe','ivreghdfe','ppmlhdfe','psmatch2','binscatter','rangestat','rangejoin','winsor2']:
        lines.append('which '+cmd)
    for index,c in cases:
        lines += [f'forvalues variant=1/3 {{', 'quietly coverage_fixture `variant\'']
        for line in c.get('setup','').splitlines(): lines.append('quietly '+line)
        opts=''
        for key in ('sortby','predict','scalars','required','ref_after'):
            if c.get(key): opts+=' '+key.replace('ref_after','refafter')+'('+quote(c[key])+')'
        lines += ['local before=$TESTS_FAILED', 'global FAILURE_COUNT=0',
                  f'capture noisily coverage_{c["oracle"]}, name("option-{index} variant-`variant\'") '+
                  'reference('+c['reference']+') candidate('+c['candidate']+')'+opts,
                  'local case_rc=_rc',
                  f'if `case_rc\' test_fail "option-{index}" "harness error r(`case_rc\')"',
                  'local bad=($TESTS_FAILED>`before\')', 'if `bad\' print_failures',
                  'local failed_cases=`failed_cases\'+`bad\'',
                  f'di "OPTION_CASE={index} VARIANT=`variant\' FAILED=`bad\' RC=`case_rc\'"', '}']
    lines += ['di "OPTION_COVERAGE_PASSED=$TESTS_PASSED FAILED=$TESTS_FAILED SKIPPED=$TESTS_SKIPPED"',
              'local rc=cond(`failed_cases\',9,0)',
              'di "OPTION_COVERAGE_COMPLETE"','}', 'local harness_rc=_rc', "if `harness_rc' local rc=`harness_rc'", 'log close','exit, clear']
    return '\n'.join(lines)+'\n'

def validate_cases(text, cases):
    results={}
    for i,v,f,rc in re.findall(r'^OPTION_CASE=(\d+) VARIANT=(\d+) FAILED=(\d+) RC=(\d+)$',text,re.M):
        if f not in ('0','1'): raise ValueError('invalid failure flag')
        key=(int(i),int(v))
        if key in results: raise ValueError('duplicate case result')
        results[key]=dict(failed=int(f),rc=int(rc))
    expected={(i,v) for i,_ in cases for v in (1,2,3)}
    if set(results)!=expected: raise ValueError(f'incomplete option cases: missing {expected-set(results)}')
    return results

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--commands',nargs='+')
    p.add_argument('--output',type=Path,default=ROOT/'temp/option_coverage')
    a=p.parse_args()
    cases=[(i,c) for i,c in enumerate(CASES) if not a.commands or c['command'] in a.commands]
    if not cases: raise SystemExit('No cases selected')
    a.output.mkdir(parents=True,exist_ok=True)
    log=a.output/'results.log'
    inputs=list((ROOT/'build').glob('*.ado'))+list((ROOT/'build').glob('*mac_arm.plugin'))
    hashes=lambda: {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    before_hashes=hashes()
    launch_error=None
    with tempfile.TemporaryDirectory(prefix='ctools-options-') as t:
        driver=Path(t)/'options.do';driver.write_text(source(cases,log))
        try: launch_driver(driver,log)
        except ValueError as error: launch_error=str(error) # Retain failed case details, but never turn launcher failure into success.
    text=log.read_text(errors='replace')
    results=validate_cases(text,cases)
    report=[]
    for i,c in cases:
        for variant in (1,2,3): report.append(dict(index=i,variant=variant,**c,**results[(i,variant)]))
    after_hashes=hashes()
    (a.output/'provenance.json').write_text(json.dumps(dict(before=before_hashes,after=after_hashes),indent=2)+'\n')
    (a.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    failures=[r for r in report if r['failed'] or r['rc']]
    for r in failures: print(f'FAIL {r["name"]} variant={r["variant"]} rc={r["rc"]}')
    print(f'{len(report)-len(failures)}/{len(report)} option/data cases passed; {len(failures)} failed; log: {log}')
    if before_hashes!=after_hashes: raise SystemExit("Build inputs changed during validation; rerun on a stable build")
    if failures or launch_error: raise SystemExit(1)
if __name__=='__main__': main()
