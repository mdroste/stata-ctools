"""Regression checks for registry coverage and the shared Stata launch contract."""
import json
from pathlib import Path
import platform
import shlex
import tempfile
from unittest.mock import patch
import suite_registry as registry
import stata_runner as launcher
import prepare_io_performance as benchmark
import prepare_command_performance as command_benchmark
import prepare_encode_performance as encode_benchmark
import prepare_group_performance as group_benchmark


def reject(fn):
    try: fn()
    except ValueError: return
    raise AssertionError('invalid input accepted')


def check_io_cleanup(root):
    for fail in (False, True):
        try:
            with launcher.io_workspace() as workspace:
                assert (workspace/'validation/io_test_helpers.do').is_file()
                assert (workspace/'build').is_dir()
                (workspace/'temp/io_allformats/native/output.csv').write_text('fixture')
                if fail:
                    raise ValueError('injected test failure')
        except ValueError:
            assert fail
        assert not workspace.exists()

    build = root/'build'
    build.mkdir()
    plugin = ('ctools_mac_arm.plugin' if platform.machine() == 'arm64' else
              'ctools_mac_x86.plugin') if platform.system() == 'Darwin' else 'ctools_linux.plugin'
    (build/plugin).write_bytes(b'fake frozen build')
    for mode in ('success', 'launch_failure', 'mismatch', 'keep'):
        out = root/mode
        out.mkdir()
        (out/'unrelated.csv').write_text('preserve me')
        (out/'comparisons.txt').write_text('stale pass')
        workspaces = []

        def fake_launch(driver, log):
            command = next(line for line in driver.read_text().splitlines()
                           if line.startswith('capture noisily do '))
            data = Path(shlex.split(command)[-1])
            workspaces.append(data)
            label = driver.stem
            log.write_text(f'IO_BENCHMARK_RC=0\nPERF,{label},small,export_csv,1,1,0.1\n')
            (data/f'{label}_small_1.csv').write_text(
                'different' if mode == 'mismatch' and label == 'c1' else 'same bytes')
            if mode == 'launch_failure':
                raise ValueError('injected Stata failure')

        argv = ['prepare_io_performance.py', '--baseline', str(build), '--candidate',
                str(build), '--output', str(out), '--cases', 'small', '--reps', '1', '--run']
        if mode == 'keep':
            argv.append('--keep-data')
        with patch('sys.argv', argv), patch.object(benchmark.subprocess, 'run'), \
                patch.object(benchmark, 'launch_driver', fake_launch):
            if mode in ('launch_failure', 'mismatch'):
                reject(benchmark.main)
                assert not (out/'comparisons.txt').exists()
            else:
                benchmark.main()
                assert (out/'comparisons.txt').read_text().startswith('PASS 3 ')
                assert (out/'timings.csv').is_file() and (out/'summary.csv').is_file()
        assert (out/'unrelated.csv').read_text() == 'preserve me'
        assert (out/'b1.log').is_file()
        assert workspaces
        if mode == 'keep':
            assert (out/'c1_small_1.csv').is_file()
            benchmark.summarize(out)
        else:
            assert all(not path.exists() for path in workspaces)


def check_dataset_cleanup(root):
    build = root/'build'
    for command in ('cipolate', 'csplit', 'crangejoin'):
        (build/(command+'.ado')).write_text('fake frozen ado')
    for name, runner in (('command', command_benchmark), ('encode', encode_benchmark),
                         ('group', group_benchmark)):
        for mode in ('success', 'launch_failure', 'summary_failure', 'keep', 'prepare'):
            if mode == 'prepare' and name != 'command':
                continue
            out = root/f'{name}_{mode}'
            out.mkdir()
            (out/'unrelated.dta').write_bytes(b'preserve me')
            workspaces = []

            def fake_launch(driver, log):
                command = next(line for line in driver.read_text().splitlines()
                               if line.startswith('capture noisily do '))
                data = Path(shlex.split(command)[-1])
                workspaces.append(data)
                (data/'expected_fixture.dta').write_bytes(b'generated comparison data')
                label = driver.stem
                if name == 'command':
                    lines = [f'PERFORMANCE_COMPLETE label={label}', 'PERFORMANCE_DRIVER_RC=0']
                    lines += [f'PERF,{label},ip_sorted,{threads},{rep},0.1,1000'
                              for threads in (1,4,12) for rep in range(1,7)]
                else:
                    lines = [f'{name.upper()}_BENCHMARK_RC=0',
                             f'PERF,{label},fixture,ctools,1000,1,1,0.1',
                             f'PERF,{label},fixture,reference,1000,0,1,0.2']
                log.write_text('\n'+'\n'.join(lines)+'\n')
                if mode == 'launch_failure':
                    raise ValueError('injected Stata failure')

            argv = [runner.__file__, '--baseline', str(build), '--candidate', str(build),
                    '--output', str(out)]
            if name == 'command':
                argv += ['--cases', 'ip_sorted']
                if mode != 'prepare':
                    argv.append('--run')
            if mode == 'keep':
                argv.append('--keep-data')
            with patch('sys.argv', argv), patch.object(runner.subprocess, 'run'), \
                    patch.object(runner, 'launch_driver', fake_launch):
                if mode == 'summary_failure':
                    with patch.object(runner, 'summarize', side_effect=ValueError('bad summary')):
                        reject(runner.main)
                elif mode == 'launch_failure':
                    reject(runner.main)
                else:
                    runner.main()
            assert (out/'unrelated.dta').read_bytes() == b'preserve me'
            if mode == 'prepare':
                assert not workspaces and (out/'b1.do').is_file()
                assert not list(out.glob('data-*'))
                continue
            assert workspaces and (out/'b1.log').is_file()
            if mode == 'keep':
                assert (out/'expected_fixture.dta').is_file()
            else:
                assert all(not path.exists() for path in workspaces)
            if mode in ('success', 'keep'):
                assert (out/'timings.csv').is_file()
                if name != 'command':
                    assert (out/'summary.csv').is_file()


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
            entry = dict(name='unregistered', path='validation/validate_unregistered.do',
                         runner='master', category='correctness', marker='unregistered',
                         fixtures=[], dependencies=[])
            reg.write_text(json.dumps({'schema':1,'suites':[entry]}))
            master.write_text(registry.master_source())
            reject(registry.check_registry)  # Registered, but no completion assignment.
            (root/entry['path']).write_text('global CTOOLS_COMPONENT_COMPLETE "unregistered"\n')
            registry.check_registry()
            master.write_text('stale generated list')
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
        check_io_cleanup(root)
        check_dataset_cleanup(root)
    print('PASS registry completeness, stale/missing/error logs, path quoting, and driver cleanup')
    print('PASS I/O cleanup on success/failure, retained logs, and benchmark keep-data')
    print('PASS command/encode/group dataset cleanup, failure paths, and manual preparation')


if __name__=='__main__':main()
