"""Check install ZIP contents and reject incomplete or corrupted packages."""
import json
from pathlib import Path
import re
import sys
import tempfile
from zipfile import ZipFile

from validation import ROOT
from scripts.archive_package import archive_package
from scripts.package_contract import check_package
from scripts.stage_package import PLATFORMS, stage


def reject(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError('invalid release archive was accepted')


def main():
    workflow = (ROOT / '.github/workflows/build.yml').read_text()
    package_block = workflow.split('\n  package:\n')[1]
    package_condition = re.search(r'    if: >-\n(.*?)\n    runs-on:', package_block, re.S).group(1)
    build_jobs = ('build-linux', 'build-macos-arm', 'build-macos-intel', 'build-windows')
    for event in ('push', 'release', 'workflow_dispatch', 'pull_request'):
        for signing in ('success', 'failure', 'skipped', 'cancelled'):
            for failed_build in (None, *build_jobs):
                expression = package_condition.replace('always()', 'True').replace('!cancelled()', 'True')
                for job in build_jobs:
                    expression = expression.replace(f'needs.{job}.result', repr('failure' if job == failed_build else 'success'))
                expression = expression.replace('needs.sign-windows.result', repr(signing))
                expression = expression.replace('github.event_name', repr(event))
                expression = ' '.join(expression.replace('&&', 'and').replace('||', 'or').split())
                allowed = failed_build is None and (signing == 'success' or (event == 'pull_request' and signing == 'skipped'))
                assert eval(expression, {'__builtins__': {}}, {}) == allowed
    block = workflow.split('\n  release:\n')[-1]
    condition = re.search(r'    if: >-\n(.*?)\n    needs:', block, re.S).group(1)
    for event, ref, publish in (
        ('push', 'refs/heads/main', True),
        ('push', 'refs/heads/dev', False),
        ('pull_request', 'refs/heads/main', False),
        ('workflow_dispatch', 'refs/heads/main', False),
        ('release', 'refs/tags/v1.0.0', True),
    ):
        for status in ('success', 'failure', 'skipped', 'cancelled'):
            expression = condition.replace('always()', 'True').replace('!cancelled()', 'True')
            for key, value in {
                'needs.package.result': status,
                'github.event_name': event, 'github.ref': ref,
            }.items():
                expression = expression.replace(key, repr(value))
            expression = ' '.join(expression.replace('&&', 'and').replace('||', 'or').split())
            assert eval(expression, {'__builtins__': {}}, {}) == (publish and status == 'success')
    assert 'PLUGIN_PUBLICATION: ${{ needs.commit-plugins.result }}' in block
    assert 'test "$PLUGIN_PUBLICATION" = success' in block
    with tempfile.TemporaryDirectory(prefix='ctools-zip-') as tmp:
        root = Path(tmp)
        binaries = root / 'binaries'
        binaries.mkdir()
        # Packaging fixtures deliberately do not pretend to be executable plugins.
        for platform in PLATFORMS:
            (binaries / f'ctools_{platform}.plugin').write_bytes(platform.encode())
        package = stage(ROOT / 'build', root / 'package', list(PLATFORMS),
                        'archive-test-revision', binaries)
        (package / 'scratch.log').write_text('must not ship')
        archive = archive_package(package, root / 'dist')
        version = json.loads((ROOT / 'config/release.json').read_text())['version']
        assert archive.name == f'ctools-{version}.zip'
        info = json.loads((package / 'BUILD_INFO.json').read_text())
        assert info['version'] == version
        assert info['release_date'] == json.loads((ROOT / 'config/release.json').read_text())['date']
        with ZipFile(archive) as zipped:
            assert set(zipped.namelist()) == set(info['sha256']) | {'BUILD_INFO.json'}
            assert all('/' not in name for name in zipped.namelist())
            assert zipped.testzip() is None
            zipped.extractall(root / 'installed')
        assert not check_package(root / 'installed', integrity=True)
        assert json.loads((root / 'installed/BUILD_INFO.json').read_text())['revision'] == 'archive-test-revision'
        assert {p.name for p in (root / 'installed').glob('*.plugin')} == {
            f'ctools_{platform}.plugin' for platform in PLATFORMS}

        single = stage(ROOT / 'build', root / 'single', ['linux'], 'test', binaries)
        reject(lambda: archive_package(single, root / 'bad'))
        wrong_release = root / 'release.json'
        wrong_release.write_text(json.dumps({'version': '999.0.0'}))
        reject(lambda: archive_package(package, root / 'bad', wrong_release))
        wrong_release.write_text(json.dumps({'version': '../escape'}))
        reject(lambda: archive_package(package, root / 'bad', wrong_release))
        original_archive = archive.read_bytes()
        (package / 'ctools_windows.plugin').write_bytes(b'corrupt')
        reject(lambda: archive_package(package, root / 'dist'))
        assert archive.read_bytes() == original_archive
        (package / 'ctools_linux.plugin').unlink()
        reject(lambda: archive_package(package, root / 'bad'))
        assert not (root / 'bad').exists()
    print('PASS flat all-platform ZIP, extracted integrity, version, scope, and corruption checks')


if __name__ == '__main__':
    main()
