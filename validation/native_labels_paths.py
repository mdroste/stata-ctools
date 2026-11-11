"""Full encoder/decoder protocol and failure-path native coverage."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_labels_paths_native():
    source=(ROOT/'validation/native/labels_paths.c').read_text(encoding='utf-8')
    source=source.replace('#include "newcommands_mock.h"',
        '#include "'+(ROOT/'validation/native/newcommands_mock.h').as_posix()+'"')
    with tempfile.TemporaryDirectory(prefix='ctools-labels-paths-') as tmp:
        for parallel in (False,True):
            run_test(Path(tmp),'labels_paths_'+str(parallel),source,openmp=parallel)
