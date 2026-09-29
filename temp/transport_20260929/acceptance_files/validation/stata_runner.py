"""The single Stata launch boundary: interactive zsh alias, fresh completion nonce.

Stata batch exit status alone does not report do-file failure. Each driver must
save _rc before logging/cleanup; a nonce rejects stale or interrupted output.
"""
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import uuid

ROOT = Path(__file__).resolve().parents[1]


def stata_command(driver):
    return ['/bin/zsh', '-lic', 'stata -q -b do ' + shlex.quote(str(Path(driver).resolve()))]


def validate_completion(text, nonce):
    matches = re.findall(r'^CTOOLS_LAUNCH=' + re.escape(nonce) + r' RC=(\d+)$', text, re.M)
    if matches != ['0']:
        raise ValueError('Stata failed or did not finish this driver: ' + repr(matches))


def launch_driver(driver, log, *, rc_local='rc', cwd=ROOT):
    source = Path(driver).read_text()
    if not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*', rc_local):
        raise ValueError('invalid return-code local')
    if not source.rstrip().endswith('exit, clear'):
        raise ValueError('batch driver must end with exit, clear')
    if '\nlog close\n' not in source:
        raise ValueError('batch driver must close its log after reporting _rc')
    nonce = uuid.uuid4().hex
    source = source.replace('\nlog close\n', f'\ndi "CTOOLS_LAUNCH={nonce} RC=`{rc_local}\'"\nlog close\n')
    log = Path(log)
    log.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(prefix='ctools-stata-') as tmp:
        executable = Path(tmp) / 'driver.do'
        executable.write_text(source)
        subprocess.run(stata_command(executable), cwd=cwd, check=True)
    text = log.read_text(errors='replace') if log.exists() else ''
    try:
        validate_completion(text, nonce)
    except ValueError:
        print(text[-20000:])
        raise
    return text


def run_do(path, log, args=(), markers=()):
    # Stata compound quotes accommodate spaces and ordinary embedded quotes.
    def quote(value):
        value = str(value)
        if '\n' in value or '\r' in value or '`' in value:
            raise ValueError('unsupported Stata path characters')
        return '`"' + value + '"\''
    source = f'''clear all
set more off
set linesize 255
capture log close _all
log using {quote(log)}, text replace
capture noisily do {quote(Path(path).resolve())} {' '.join(quote(a) for a in args)}
local rc = _rc
log close
exit, clear
'''
    with tempfile.TemporaryDirectory(prefix='ctools-do-') as tmp:
        driver = Path(tmp) / 'driver.do'
        driver.write_text(source)
        text = launch_driver(driver, log)
    for marker in markers:
        if marker not in text.splitlines():
            raise ValueError(f'missing completion marker: {marker}')
    return text
