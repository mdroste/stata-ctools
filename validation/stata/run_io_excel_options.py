"""Run native Excel/DBF option comparisons through the stata shell alias."""
from validation.support.stata_runner import io_workspace, launch_driver
import re
from pathlib import Path
from validation.fixtures.io_excel_format_fixtures import fixtures
from validation.fixtures.io_excel_formula_fixtures import formula_fixtures, compare_formula_outputs
from validation import ROOT
def main():
    log = ROOT / 'temp/io_excel_options/regressions.log'
    log.parent.mkdir(parents=True, exist_ok=True)
    with io_workspace() as workspace:
        directory = workspace / 'temp/io_excel_options'
        directory.mkdir(parents=True, exist_ok=True)
        fixtures(directory)
        formula_fixtures(directory)
        for ext in ('xls', 'xlsx'):
            (directory/f'{ext}_invalid.{ext}').write_bytes(b'invalid file header\n')
            (directory/f'{ext}_truncated.{ext}').write_bytes((directory/f'display_formats.{ext}').read_bytes()[:16])
        log.unlink(missing_ok=True)
        driver = directory / 'driver.do'
        driver.write_text(f'''clear all
capture log close _all
log using "{log}", text replace
capture noisily do "{ROOT}/validation/stata/validate_io_excel_options.do" "{directory}"
local rc = _rc
di "IO_OPTIONS_DRIVER_RC=`rc'"
log close
exit, clear
''')
        launch_driver(driver, log, cwd=workspace)
        formula_checks = compare_formula_outputs(directory)
    text = log.read_text(errors='replace')
    result = re.findall(r'IO_OPTIONS_DRIVER_RC=([0-9]+)', text)
    summaries = re.findall(r'IO_OPTIONS_PASSED=([0-9]+) FAILED=([0-9]+)', text)
    if not result or result[-1] != '0' or not summaries or summaries[-1][1] != '0':
        raise SystemExit('FAIL Excel option regressions; see '+str(log))
    print('PASS '+summaries[-1][0]+' exact option comparisons; '+str(log))
    print(f'PASS {formula_checks} formula/calculation-chain structural comparisons')

if __name__ == '__main__':
    main()
