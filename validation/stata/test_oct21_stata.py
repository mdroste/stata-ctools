"""Licensed integration regressions for the October 21 release audit."""
from pathlib import Path
import os
import tempfile
from validation.fixtures.oct21_fixtures import fixtures
from validation.support.stata_runner import ROOT, run_do


def main():
    log = ROOT/'temp/oct21_stata.log'
    log.parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='ctools-oct21-stata-') as tmp:
        directory = Path(tmp)
        fixtures(directory)
        run_do(ROOT/'validation/stata/validate_oct21.do', log,
               [directory, os.environ.get('CTOOLS_TEST_PLUGIN_DIR', '')], ['CTOOLS_OCT21_COMPLETE'])
    print('PASS October 21 Stata regressions')


if __name__ == '__main__':
    main()
