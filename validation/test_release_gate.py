"""Reject incomplete Stata logs and incomplete installation directories."""
import importlib.util
from pathlib import Path
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load('audit_runner', ROOT / 'validation/run_stata_audit.py')
packaging = load('staging', ROOT / 'scripts/stage_package.py')


def must_reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError('invalid release input was accepted')


def logs():
    nonce = 'current-run'
    log = '\n' + '\n'.join(
        f'CTOOLS_COMPONENT={name} RC=0 COMPLETE={name}' for name in runner.COMPONENTS)
    log += '\nCTOOLS_FULL_PASSED=2000 FAILED=0 SKIPPED=48 EXPECTED_SKIPPED=48\n'
    log += f'CTOOLS_DRIVER={nonce} RC=0\n'
    assert runner.validate_log(log, nonce) == (2000, 48)
    for bad in (
        log.replace(nonce, 'stale-run'),
        log.replace('COMPLETE=cpplmhdfe', 'COMPLETE='),
        log.replace('cmerge RC=0', 'cmerge RC=679'),
        log.replace('FAILED=0', 'FAILED=1'),
        log.replace('SKIPPED=48 EXPECTED', 'SKIPPED=49 EXPECTED'),
        log.replace('PASSED=2000', 'PASSED=0'),
        log.replace('CTOOLS_FULL_', 'MISSING_FULL_'),
        log.replace(f'CTOOLS_DRIVER={nonce} RC=0', f'CTOOLS_DRIVER={nonce} RC=1'),
    ):
        must_reject(lambda: runner.validate_log(bad, nonce))
    print('PASS complete-run, failure, stale-log, abort, and unexpected-skip gates')


def packages():
    with tempfile.TemporaryDirectory(prefix='ctools-package-test-') as tmp:
        directory = Path(tmp)
        source = directory / 'source'
        source.mkdir()
        # Real source manifests and helpers; dummy binaries test staging only.
        for entry in (ROOT / 'build').iterdir():
            if entry.is_file() and (entry.suffix in ('.ado', '.sthlp', '.pkg')
                                   or entry.name in ('stata.toc', 'LICENSE', 'THIRD_PARTY_NOTICES')):
                shutil.copy2(entry, source / entry.name)
        for target in packaging.PLATFORMS:
            (source / f'ctools_{target}.plugin').write_bytes(b'test fixture, not a plugin')
        for targets in (['mac_arm'], list(packaging.PLATFORMS)):
            staged = packaging.stage(source, directory / ('stage-' + targets[0] + str(len(targets))),
                                     targets, 'test-only')
            assert not packaging.check_package(staged)
            assert (staged / 'cpplmhdfe_p.ado').is_file()
            assert {p.name for p in staged.glob('*.plugin')} == {f'ctools_{p}.plugin' for p in targets}
            (staged / 'cpplmhdfe_p.ado').unlink()
            assert packaging.check_package(staged)
        (source / 'ctools_linux.plugin').unlink()
        must_reject(lambda: packaging.stage(source, directory / 'bad', list(packaging.PLATFORMS), 'test-only'))
        assert not (directory / 'bad').exists()
        (source / 'cpplmhdfe_p.ado').unlink()
        must_reject(lambda: packaging.stage(source, directory / 'bad-helper', ['mac_arm'], 'test-only'))
    print('PASS single/all-platform manifests, helper completeness, and atomic staging')


def inventory():
    public = {p.stem for p in (ROOT / 'build').glob('*.sthlp')} - {'ctools'}
    registered = set(runner.COMPONENTS) - {'audit_p1', 'audit_p2', 'sep22', 'sep24'}
    assert public == registered, (public - registered, registered - public)
    manifest = (ROOT / 'build/ctools.pkg').read_text().splitlines()
    for name in public:
        assert (ROOT / 'src' / name).is_dir(), name
        for suffix in ('.ado', '.sthlp'):
            assert (ROOT / 'build' / (name + suffix)).is_file(), name
            assert f'f {name}{suffix}' in manifest, name
    print(f'PASS {len(public)} public-command help/source/package/test inventory')


if __name__ == '__main__':
    logs()
    packages()
    inventory()
