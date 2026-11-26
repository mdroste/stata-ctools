"""Coefficient-space fixed-effect projection against an independent reference."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_coef_native():
    source = (ROOT / 'validation/native/hdfe_coef.c').read_text(encoding='utf-8')
    source = source.replace(
        '#include "newcommands_mock.h"',
        '#include "' + (ROOT / 'validation/native/newcommands_mock.h').as_posix() + '"',
    )
    with tempfile.TemporaryDirectory(prefix='ctools-hdfe-coef-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'hdfe_coef_' + str(parallel), source, openmp=parallel)
