"""Grouped/parallel correctness and publication/allocator failure coverage."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_newcommands_paths_native():
    source=(ROOT/'validation/native/newcommands_paths.c').read_text(encoding='utf-8')
    source=source.replace('#include "newcommands_mock.h"',
        '#include "'+(ROOT/'validation/native/newcommands_mock.h').as_posix()+'"')
    with tempfile.TemporaryDirectory(prefix='ctools-newcommands-paths-') as tmp:
        for parallel in (False,True):
            run_test(Path(tmp),'newcommands_paths_'+str(parallel),source,
                     modules=('ctools_parse.c',),openmp=parallel)
