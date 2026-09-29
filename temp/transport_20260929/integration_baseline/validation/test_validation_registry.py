"""Regression checks for registry coverage and the shared Stata launch contract."""
import json
from pathlib import Path
import tempfile
from unittest.mock import patch
import suite_registry as registry
import stata_runner as launcher


def reject(fn):
    try: fn()
    except ValueError: return
    raise AssertionError('invalid input accepted')


def main():
    registry.check_registry()
    assert 'memory_safety' in registry.master_names()
    for marker in ('', '\nCTOOLS_LAUNCH=stale RC=0\n', '\nCTOOLS_LAUNCH=now RC=198\n'):
        reject(lambda: launcher.validate_completion(marker, 'now'))
    launcher.validate_completion('\nCTOOLS_LAUNCH=now RC=0\n','now')
    with tempfile.TemporaryDirectory(prefix='ctools registry spaces ') as tmp:
        root=Path(tmp);(root/'validation').mkdir()
        reg=root/'suites.json';master=root/'master.do'
        reg.write_text(json.dumps({'schema':1,'suites':[]}))
        with patch.multiple(registry, ROOT=root, REGISTRY=reg, MASTER=master):
            master.write_text(registry.master_source());registry.check_registry()
            (root/'validation/validate_unregistered.do').write_text('')
            reject(registry.check_registry)
        driver=root/'driver with spaces.do';log=root/'out.log'
        driver.write_text(f'log using "{log}", text replace\nlocal rc = 0\nlog close\nexit, clear\n')
        launched=[]
        def fake_run(argv, **kwargs):
            import shlex,re
            assert argv[:2]==['/bin/zsh','-lic']
            command=shlex.split(argv[2]);assert command[:4]==['stata','-q','-b','do']
            path=Path(command[4]);launched.append(path)
            text=path.read_text();assert text.rstrip().endswith('exit, clear')
            nonce=re.search('CTOOLS_LAUNCH=([a-f0-9]+)',text)[1]
            log.write_text(f'\nCTOOLS_LAUNCH={nonce} RC=0\n')
        with patch.object(launcher.subprocess,'run',fake_run): launcher.launch_driver(driver,log)
        assert launched and not launched[0].exists()
        driver.write_text('local rc = 0\n')
        reject(lambda: launcher.launch_driver(driver,log))
    print('PASS registry completeness, stale/missing/error logs, path quoting, and driver cleanup')


if __name__=='__main__':main()
