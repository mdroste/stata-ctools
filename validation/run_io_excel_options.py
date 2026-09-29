"""Run native Excel/DBF option comparisons through the stata shell alias."""
from stata_runner import launch_driver
import re
import subprocess
from pathlib import Path
from io_excel_format_fixtures import fixtures
from io_excel_formula_fixtures import formula_fixtures, compare_formula_outputs
ROOT = Path(__file__).resolve().parents[1]

def main():
    directory = ROOT / 'temp/io_excel_options'
    directory.mkdir(parents=True, exist_ok=True)
    for prefix in ('native', 'result'):
        for path in directory.glob(f'formula_{prefix}_*.xlsx'):
            path.unlink()
    fixtures(directory)
    formula_fixtures(directory)
    for ext in ('xls', 'xlsx'):
        (directory/f'{ext}_invalid.{ext}').write_bytes(b'invalid file header\n')
        (directory/f'{ext}_truncated.{ext}').write_bytes((directory/f'display_formats.{ext}').read_bytes()[:16])
    log = directory / 'regressions.log'
    log.unlink(missing_ok=True)
    driver = directory / 'driver.do'
    driver.write_text(f'''clear all
capture log close _all
log using "{log}", text replace
capture noisily do "{ROOT}/validation/validate_io_excel_options.do" "{directory}"
local rc = _rc
di "IO_OPTIONS_DRIVER_RC=`rc'"
log close
exit, clear
''')
    launch_driver(driver, log)
    text = log.read_text(errors='replace')
    result = re.findall(r'IO_OPTIONS_DRIVER_RC=([0-9]+)', text)
    summaries = re.findall(r'IO_OPTIONS_PASSED=([0-9]+) FAILED=([0-9]+)', text)
    if not result or result[-1] != '0' or not summaries or summaries[-1][1] != '0':
        raise SystemExit('FAIL Excel option regressions; see '+str(log))
    formula_checks = compare_formula_outputs(directory)
    print('PASS '+summaries[-1][0]+' exact option comparisons; '+str(log))
    print(f'PASS {formula_checks} formula/calculation-chain structural comparisons')

if __name__ == '__main__':
    main()
