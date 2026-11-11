"""Cancellation tests, including OpenMP host-finished-first and nested teams."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_cancel_native():
    source = (ROOT / 'validation/native/hdfe_cancel.c').read_text()
    with tempfile.TemporaryDirectory(prefix='ctools-hdfe-cancel-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'hdfe_cancel_' + str(parallel), source, openmp=parallel)
