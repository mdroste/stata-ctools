"""Exercise signing failure paths without a certificate or Apple network access."""

import base64
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from validation import ROOT


MOCK = r'''#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys

command = Path(sys.argv[0]).name
args = sys.argv[1:]
scenario = os.environ['SIGNING_SCENARIO']
with open(os.environ['SIGNING_JOURNAL'], 'a') as log:
    log.write(json.dumps([command, *args]) + '\n')
if command == 'security':
    if args[0] == 'create-keychain':
        Path(args[-1]).touch()
    if args[0] == 'import' and scenario == 'import-fails':
        sys.exit(1)
    if args[0] == 'delete-keychain':
        Path(args[-1]).unlink()
    if args[0] == 'list-keychains' and '-s' not in args:
        print('    "/Users/runner/Library/Keychains/login.keychain-db"')
    if args[0] == 'find-identity' and scenario != 'missing-identity':
        team = 'WRONGTEAM0' if scenario == 'wrong-team' else 'ABCDE12345'
        print('  1) ' + 'A' * 40 + f' "Developer ID Application: Test ({team})"')
if command == 'codesign':
    if '--sign' in args and scenario == 'sign-fails':
        sys.exit(1)
    if '--verify' in args and scenario == 'invalid-signature':
        sys.exit(1)
    if '--display' in args and scenario != 'missing-timestamp':
        print('Timestamp=Nov 2, 2026 at 12:00:00', file=sys.stderr)
    if '--check-notarization' in args and scenario == 'missing-ticket':
        sys.exit(1)
if command == 'ditto':
    Path(args[-1]).write_bytes(b'archive')
if command == 'xcrun':
    status = 'Invalid' if scenario == 'rejected' else 'Accepted'
    if args[1] == 'store-credentials' and scenario == 'invalid-credentials':
        sys.exit(1)
    if args[1] == 'submit':
        if scenario == 'timeout':
            status = 'In Progress'
        print(json.dumps({'id': '12345678-1234-1234-1234-123456789abc', 'status': status}))
        if scenario == 'timeout':
            sys.exit(1)
    if args[1] == 'log':
        if scenario == 'log-fails':
            sys.exit(1)
        Path(args[-1]).write_text(json.dumps({'status': status, 'issues': None}))
'''


def main():
    workflow = (ROOT / '.github/workflows/build.yml').read_text()
    signing = workflow.split('\n  sign-macos:\n')[1].split('\n  build-windows:\n')[0]
    condition = re.search(r'    if: >-\n(.*?)\n    needs:', signing, re.S).group(1)
    for event in ('push', 'pull_request', 'workflow_dispatch', 'release'):
        for ref in (
            'refs/heads/main', 'refs/heads/dev', 'refs/heads/develop',
            'refs/heads/untrusted', 'refs/tags/v1.2.3', 'refs/tags/other',
            'refs/pull/1/merge',
        ):
            expression = condition.replace('github.event_name', repr(event))
            expression = expression.replace('github.ref', repr(ref))
            expression = ' '.join(expression.replace('&&', 'and').replace('||', 'or').split())
            allowed = event != 'pull_request' and (
                ref in ('refs/heads/main', 'refs/heads/dev', 'refs/heads/develop')
                or (event == 'release' and ref.startswith('refs/tags/v'))
            )
            assert eval(expression, {'__builtins__': {}, 'startsWith': str.startswith}) == allowed
    assert 'environment: macos-signing' in signing
    assert 'pattern: macos-*-unsigned' in signing
    assert 'pattern: ctools_*' in workflow.split('\n  package:\n')[1]
    print('PASS macOS signing: trusted refs and artifact isolation')
    script = ROOT / 'scripts/sign_macos.sh'
    subprocess.run(['bash', '-n', str(script)], check=True)
    with tempfile.TemporaryDirectory(prefix='ctools-signing-test-') as tmp:
        root = Path(tmp)
        commands = root / 'bin'
        commands.mkdir()
        for name in ('security', 'codesign', 'xcrun', 'ditto', 'sleep'):
            path = commands / name
            path.write_text(MOCK)
            path.chmod(0o755)
        scenarios = (
            'accepted', 'missing-config', 'import-fails', 'missing-identity',
            'wrong-team', 'sign-fails', 'invalid-signature', 'missing-timestamp',
            'invalid-credentials', 'rejected', 'timeout', 'log-fails', 'missing-ticket',
        )
        for scenario in scenarios:
            work = root / scenario
            work.mkdir()
            plugins = work / 'plugins'
            plugins.mkdir()
            for arch in ('arm', 'x86'):
                (plugins / f'ctools_mac_{arch}.plugin').write_bytes(b'fixture')
            journal = work / 'journal.jsonl'
            summary = work / 'summary'
            env = dict(
                os.environ,
                PATH=str(commands) + os.pathsep + os.environ['PATH'],
                RUNNER_TEMP=str(work),
                GITHUB_STEP_SUMMARY=str(summary),
                APPLE_CERTIFICATE_P12_BASE64=base64.b64encode(b'fixture').decode(),
                APPLE_CERTIFICATE_PASSWORD='fixture-export-password',
                APPLE_TEAM_ID='ABCDE12345',
                APPLE_ID='fixture@example.invalid',
                APPLE_APP_SPECIFIC_PASSWORD='fixture-notary-password',
                SIGNING_SCENARIO=scenario,
                SIGNING_JOURNAL=str(journal),
            )
            if scenario == 'missing-config':
                del env['APPLE_CERTIFICATE_P12_BASE64']
            result = subprocess.run(
                ['bash', str(script), str(plugins), str(work / 'results')],
                env=env, text=True, capture_output=True,
            )
            assert (result.returncode == 0) == (scenario == 'accepted'), (
                scenario, result.stdout, result.stderr,
            )
            assert not list(work.glob('ctools-signing.*')), scenario
            assert summary.exists() == (scenario == 'accepted'), scenario
            calls = [json.loads(line) for line in journal.read_text().splitlines()] if journal.exists() else []
            if scenario != 'missing-config':
                assert any(c[:2] == ['security', 'delete-keychain'] for c in calls), scenario
                # The temporary keychain joins the search list before signing,
                # and the user's original list is restored on every exit.
                login = '/Users/runner/Library/Keychains/login.keychain-db'
                searches = [c for c in calls if c[:4] == ['security', 'list-keychains', '-d', 'user'] and '-s' in c]
                assert searches and searches[0][6:] == [login], (scenario, searches)
                assert Path(searches[0][5]).name == 'signing.keychain-db', (scenario, searches)
                assert searches[-1][5:] == [login], (scenario, searches)
                first_search = calls.index(searches[0])
                assert all(i > first_search for i, c in enumerate(calls) if c[0] == 'codesign'), scenario
            if scenario == 'accepted':
                signs = [c for c in calls if c[0] == 'codesign' and '--sign' in c]
                tickets = [c for c in calls if c[0] == 'codesign' and '--check-notarization' in c]
                assert len(signs) == len(tickets) == 2
                for call in signs:
                    assert '--timestamp' in call and '--keychain' in call
                assert {Path(c[-1]).name for c in tickets} == {
                    'ctools_mac_arm.plugin', 'ctools_mac_x86.plugin',
                }
            if scenario in ('rejected', 'timeout'):
                assert (work / 'results/notarization-log.json').exists()
                assert not any('--check-notarization' in c for c in calls)
            if scenario in ('import-fails', 'wrong-team', 'sign-fails', 'invalid-signature', 'missing-timestamp'):
                assert not any(c[:3] == ['xcrun', 'notarytool', 'submit'] for c in calls)
            print(f'PASS macOS signing: {scenario}')


if __name__ == '__main__':
    main()
