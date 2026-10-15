"""Mutation checks for option coverage accounting and complete result reporting."""
import json
from pathlib import Path
import shutil
import tempfile
from unittest.mock import patch
import option_coverage as coverage
from option_cases import CASES
from run_option_coverage import validate_cases, source


def reject(fn):
    try:
        fn()
    except ValueError:
        return
    raise AssertionError('invalid coverage accepted')


def main():
    rows=coverage.records()
    assert rows and not [r for r in rows if r['status']=='gap']
    # Transaction wrappers must not hide their implementation's public syntax.
    assert {'cencode','cdecode','cmerge','cwinsor','crangestat','creghdfe'} <= {r['command'] for r in rows}
    assert any(r['command']=='cencode' and r['option']=='replace' for r in rows)
    assert (coverage.ROOT/'docs/OPTION_COVERAGE.md').read_text()==coverage.report(rows)
    cases=[(0,CASES[0]),(1,CASES[1])]
    complete='\n'.join(f'OPTION_CASE={i} VARIANT={v} FAILED=0 RC=0' for i,_ in cases for v in (1,2,3))
    assert len(validate_cases(complete,cases))==6
    reject(lambda: validate_cases(complete.split('\n',1)[1],cases))
    reject(lambda: validate_cases(complete+'\n'+complete,cases))
    reject(lambda: validate_cases(complete+'\nOPTION_CASE=99 VARIANT=1 FAILED=0 RC=0',cases))
    reject(lambda: validate_cases(complete.replace('FAILED=0','FAILED=2',1),cases))
    failed=validate_cases(complete.replace('FAILED=0 RC=0','FAILED=1 RC=198',1),cases)
    assert failed[(0,1)]==dict(failed=1,rc=198)
    driver=source(cases,Path('/tmp/option coverage.log'))
    assert driver.rstrip().endswith('exit, clear') and 'OPTION_COVERAGE_COMPLETE' in driver
    assert coverage.exercised(dict(options=['generate'],candidate='cencode s, gen(out)'), 'generate','GENerate')
    assert not coverage.exercised(dict(options=['generate'],candidate='cencode s, agent(out)'), 'generate','GENerate')
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        shutil.copytree(coverage.ROOT/'package/ado',root/'package/ado')
        shutil.copytree(coverage.ROOT/'validation',root/'validation',ignore=shutil.ignore_patterns('*.log','__pycache__'))
        path=root/'package/ado/csort.ado'
        path.write_text(path.read_text().replace('syntax ', 'syntax ',1)+'\nprogram define csort\nsyntax, NEWcoverageoption\nend\n')
        with patch.object(coverage,'ROOT',root):
            gaps=[r for r in coverage.records() if r['status']=='gap']
            assert [(r['command'],r['option']) for r in gaps]==[('csort','newcoverageoption')]
            evidence=root/'validation/option_evidence.json'
            data=json.loads(evidence.read_text())
            data['ctools']['version']['evidence'][0]['anchor']='ABSENT_TEST_ASSERTION_1234'
            evidence.write_text(json.dumps(data))
            reject(coverage.records)
    print('PASS option inventory mutations, stale evidence, incomplete/duplicate results, and failed-case retention')


if __name__=='__main__':
    main()
