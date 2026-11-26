"""Weak FE connections must not be certified as converged by a zero CG step."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_numerics_native():
    source = (ROOT / 'validation/native/hdfe_numerics.c').read_text(encoding='utf-8')
    with tempfile.TemporaryDirectory(prefix='ctools-hdfe-numerics-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'hdfe_numerics_' + str(parallel), source, openmp=parallel)
