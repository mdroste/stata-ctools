"""The single Stata launch boundary: configured executable or local zsh alias.

Stata batch exit status alone does not report do-file failure. Each driver must
save _rc before logging/cleanup; a nonce rejects stale or interrupted output.
"""
from pathlib import Path
from contextlib import contextmanager
import re
import os
import shlex
import subprocess
import tempfile
import uuid

from validation import ROOT
@contextmanager
def benchmark_data(out, *, keep=False):
    """Share generated data across benchmark sessions, then remove it even on failure."""
    if keep:
        yield Path(out)
    else:
        with tempfile.TemporaryDirectory(prefix='data-', dir=out) as tmp:
            yield Path(tmp)


@contextmanager
def io_workspace():
    """Give legacy relative-path I/O tests disposable data, retaining external logs."""
    with tempfile.TemporaryDirectory(prefix='ctools-io-') as tmp:
        workspace = Path(tmp)
        for name in ('build', 'validation'):
            (workspace / name).symlink_to(ROOT / name, target_is_directory=True)
        # cio_export_test shares these output paths across the option suites.
        for name in ('native', 'result'):
            (workspace / 'temp/io_allformats' / name).mkdir(parents=True)
        yield workspace


def stata_command(driver):
    executable = os.environ.get('CTOOLS_STATA_EXECUTABLE')
    driver = str(Path(driver).resolve())
    if executable:
        # A path, not a shell command: spaces need no shell quoting.
        batch = ['/e', 'do'] if os.name == 'nt' else ['-q', '-b', 'do']
        return [executable, *batch, driver]
    if os.name == 'nt':
        raise ValueError('Set CTOOLS_STATA_EXECUTABLE to the licensed Stata executable')
    # Retain the existing local alias workflow; CI configures an executable.
    return ['/bin/zsh', '-lic', 'stata -q -b do ' + shlex.quote(driver)]


def validate_completion(text, nonce):
    matches = re.findall(r'^CTOOLS_LAUNCH=' + re.escape(nonce) + r' RC=(\d+)$', text, re.M)
    if matches != ['0']:
        raise ValueError('Stata failed or did not finish this driver: ' + repr(matches))


def quote(value):
    """Quote a Stata argument, including paths with spaces or ordinary quotes."""
    value = str(value)
    if '\n' in value or '\r' in value or '`' in value:
        raise ValueError('unsupported Stata path characters')
    return '`"' + value + '"\''


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
        # Stata opens driver.log before reading the do-file. Start in scratch
        # space, then restore the requested working directory for suite paths.
        executable.write_text(f'cd {quote(Path(cwd).resolve())}\n' + source)
        subprocess.run(stata_command(executable), cwd=tmp, check=True)
    text = log.read_text(errors='replace') if log.exists() else ''
    try:
        validate_completion(text, nonce)
    except ValueError:
        print(text[-20000:])
        raise
    return text


def run_do(path, log, args=(), markers=()):
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
