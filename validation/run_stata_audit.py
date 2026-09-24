"""Run the complete offline correctness gate using the required oldstata wrapper."""
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import uuid

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / 'validation/audit_ci.log'
COMPONENTS = 'csort cmerge cimport cexport creghdfe cqreg civreghdfe cdecode cencode cdestring csample cbsample cbinscatter cpsmatch crangestat cwinsor cpplmhdfe audit_p1 audit_p2 sep22 sep24'.split()


def validate_log(text, nonce):
    if f'\nCTOOLS_DRIVER={nonce} RC=0\n' not in text:
        raise ValueError('Stata driver failed or did not finish')
    for name in COMPONENTS:
        if f'\nCTOOLS_COMPONENT={name} RC=0 COMPLETE={name}\n' not in text:
            raise ValueError(f'incomplete component: {name}')
    result = re.search(r'^CTOOLS_FULL_PASSED=(\d+) FAILED=(\d+) SKIPPED=(\d+) EXPECTED_SKIPPED=(\d+)$', text, re.M)
    if not result:
        raise ValueError('missing full-suite summary')
    passed, failed, skipped, expected = map(int, result.groups())
    if passed == 0 or failed or skipped != expected:
        raise ValueError('failed checks, unexpected skips, or empty suite')
    return passed, skipped


def main():
    subprocess.run(['python3', 'scripts/fetch_validation_data.py', '--check-only'], cwd=ROOT, check=True)
    LOG.unlink(missing_ok=True)
    nonce = uuid.uuid4().hex
    with tempfile.TemporaryDirectory(prefix='ctools-stata-') as tmp:
        driver = Path(tmp) / 'audit.do'
        driver.write_text(f'''clear all
set more off
set linesize 255
capture log close _all
log using "{LOG}", text replace
local failed = 0
foreach cmd in reghdfe ivreghdfe ppmlhdfe psmatch2 binscatter rangestat winsor2 gstats {{
    capture which `cmd'
    if _rc {{
        di as error "Required validation reference missing: `cmd'"
        local failed = 1
    }}
}}
if !`failed' {{
    capture noisily do validation/validate_all.do
    local failed = _rc
}}
di "CTOOLS_DRIVER={nonce} RC=`failed'"
log close
exit, clear
''')
        subprocess.run(['/bin/zsh', '-lic', 'oldstata -q -b do ' + shlex.quote(str(driver))], cwd=ROOT, check=True)
    text = LOG.read_text()
    try:
        passed, skipped = validate_log(text, nonce)
    except ValueError as error:
        print(text[-20000:])
        raise SystemExit(f'{error}; see {LOG}')
    print(f'Complete Stata suite: {passed} passed; {skipped} explicitly documented method comparisons excluded. Log: {LOG}')

if __name__ == '__main__': main()
