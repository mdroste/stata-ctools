"""Analytic HDFE slope projections and detailed metadata contracts."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_paths_native():
    source=(ROOT/'validation/native/hdfe_paths.c').read_text(encoding='utf-8')
    source=source.replace('#include "newcommands_mock.h"',
        '#include "'+(ROOT/'validation/native/newcommands_mock.h').as_posix()+'"')
    with tempfile.TemporaryDirectory(prefix='ctools-hdfe-paths-') as tmp:
        for parallel in (False,True):
            run_test(Path(tmp),'hdfe_paths_'+str(parallel),source,openmp=parallel)
