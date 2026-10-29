"""Exercise publication decisions, release checksums, and cross-shell cleanup."""
import fnmatch
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
from unittest.mock import patch
from types import SimpleNamespace
from validation.support.stata_runner import stata_command

from validation import ROOT
import sys
from scripts.package_contract import check_package
from scripts.stage_package import stage


def workflow_contract():
    source = (ROOT/'.github/workflows/build.yml').read_text()
    for job in ('commit-plugins', 'release'):
        block = source.split('\n  '+job+':\n')[-1]
        condition = re.search(r'    if: (.*?)\n    needs:', block, re.S).group(1)
        condition = condition.removeprefix('>-').strip()
        assert 'stata-validation' in block.split('    runs-on:')[0]
        for status in ('success', 'failure', 'skipped', 'cancelled'):
            expression = condition.replace('always()', 'True').replace('!cancelled()', 'True')
            expression = expression.replace('needs.package.result', repr('success'))
            expression = expression.replace('needs.stata-validation.result', repr(status))
            expression = expression.replace('github.event_name', repr('push' if job == 'commit-plugins' else 'release'))
            expression = expression.replace('github.ref', repr('refs/heads/main'))
            expression = ' '.join(expression.replace('&&', 'and').replace('||', 'or').split())
            assert eval(expression, {'__builtins__': {}}, {}) == (status == 'success'), (job, status)
    paths = re.findall(r"    paths:\n((?:      - .*\n)+)", source)
    assert len(paths) == 2
    for block in paths:
        patterns = re.findall(r"'([^']+)'", block)
        for name in ('build/new.ado', 'build/cimport.sthlp', 'build/ctools.pkg'):
            assert any(fnmatch.fnmatchcase(name, pattern) for pattern in patterns)
    assert 'platform: [mac_arm, mac_x86, windows, linux]' in source
    assert 'name: ctools_${{ matrix.platform }}' in source
    with patch.dict(os.environ, {'CTOOLS_STATA_EXECUTABLE': '/test path/Stata'}):
        assert stata_command('driver.do')[0] == '/test path/Stata'
        assert stata_command('driver.do')[-1] == str(Path('driver.do').resolve())
    with patch('validation.support.stata_runner.os', SimpleNamespace(name='nt', environ={'CTOOLS_STATA_EXECUTABLE': 'C:/Program Files/Stata/StataMP.exe'})):
        command = stata_command('driver.do')
        assert command[0] == 'C:/Program Files/Stata/StataMP.exe'
        assert command[1:3] == ['/e', 'do']
    print('PASS publication rejects skipped/failed/cancelled gates; canonical sources trigger CI')


def packages():
    with tempfile.TemporaryDirectory(prefix='ctools-release-integrity-') as tmp:
        directory = stage(ROOT/'build', Path(tmp)/'package', ['mac_arm'], 'test-only')
        assert not check_package(directory, integrity=True)
        info_path = directory/'BUILD_INFO.json'
        original = info_path.read_bytes()
        for name in ('ctools_mac_arm.plugin', 'cppmlhdfe_p.ado'):
            payload_path = directory/name
            payload = payload_path.read_bytes()
            payload_path.write_bytes(b'corrupt or truncated payload')
            assert any('checksum mismatch' in e for e in check_package(directory, integrity=True))
            payload_path.write_bytes(payload)
        for mutation in ('missing', 'scope', 'inventory', 'digest', 'manifest', 'invalid'):
            info = json.loads(original)
            if mutation == 'missing': info_path.unlink()
            elif mutation == 'invalid': info_path.write_text('null')
            else:
                if mutation == 'scope': info['platforms'] = ['linux']
                if mutation == 'inventory': info['sha256'].pop('cppmlhdfe_p.ado')
                if mutation == 'digest': info['sha256']['ctools.pkg'] = '0' * 64
                if mutation == 'manifest': (directory/'ctools.pkg').write_text((directory/'ctools.pkg').read_text()+'d changed\n')
                info_path.write_text(json.dumps(info))
            assert check_package(directory, integrity=True), mutation
            if mutation == 'manifest':
                p = directory/'ctools.pkg'; p.write_text(p.read_text().removesuffix('d changed\n'))
            info_path.write_bytes(original)
        assert not check_package(directory, integrity=True)
    print('PASS checksum corruption, inventory, scope, manifest, and missing build-info checks')


def cleanup():
    with tempfile.TemporaryDirectory(prefix='ctools-clean-') as tmp:
        directory = Path(tmp)
        build, source = directory/'build', directory/'src'
        build.mkdir(); source.mkdir(); (build/'obj').mkdir()
        artifacts = [build/'ctools.plugin', build/'ctools_windows.plugin', build/'x.o', source/'x.o', build/'obj/x.o']
        for artifact in artifacts: artifact.touch()
        result = subprocess.run(['make', '-s', '-f', str(ROOT/'Makefile'), 'clean', 'DETECTED_OS=Windows',
                                 f'BUILD_DIR={build}', f'SRC_DIR={source}', 'LIBOMP_PREFIX=/nonexistent'], cwd=ROOT, check=True)
        assert all(not p.exists() for p in artifacts)
        # A matching directory cannot be unlinked: genuine cleanup errors propagate.
        (build/'bad.o').mkdir()
        result = subprocess.run(['python3', 'scripts/clean_build.py', str(build), str(source)], cwd=ROOT,
                                capture_output=True)
        assert result.returncode != 0
    print('PASS MSYS/POSIX recipe cleanup, generic plugin, and error propagation')


if __name__ == '__main__':
    workflow_contract(); packages(); cleanup()
