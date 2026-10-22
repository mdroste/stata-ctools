"""Licensed Stata regression checks for the October 22 audit fixes."""
from pathlib import Path
import tempfile
from oct22_fixtures import create
from stata_runner import run_do

ROOT=Path(__file__).resolve().parents[1]

def main():
    with tempfile.TemporaryDirectory(prefix='ctools-oct22-stata-') as temporary:
        folder=Path(temporary)
        create(folder)
        log=ROOT/'temp/registry_logs/oct22_stata.log'
        log.parent.mkdir(parents=True,exist_ok=True)
        run_do(ROOT/'validation/validate_oct22.do', log, args=[folder], markers=['OCT22_REGRESSIONS_COMPLETE'])
        print('OCT22_REGRESSIONS_COMPLETE')

if __name__=='__main__': main()
