"""creghdfe's column preparation, bit for bit against the passes it replaces."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_creghdfe_columns_native():
    source = (ROOT / 'validation/native/creghdfe_columns.c').read_text(encoding='utf-8')
    source = source.replace(
        '#include "newcommands_mock.h"',
        '#include "' + (ROOT / 'validation/native/newcommands_mock.h').as_posix() + '"',
    )
    with tempfile.TemporaryDirectory(prefix='ctools-creghdfe-columns-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'creghdfe_columns_' + str(parallel), source, openmp=parallel)
