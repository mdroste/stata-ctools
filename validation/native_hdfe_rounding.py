"""Row-space conjugate gradient at the rounding of a column far from zero."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_rounding_native():
    source = (ROOT / 'validation/native/hdfe_rounding.c').read_text(encoding='utf-8')
    with tempfile.TemporaryDirectory(prefix='ctools-hdfe-rounding-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'hdfe_rounding_' + str(parallel), source, openmp=parallel)
