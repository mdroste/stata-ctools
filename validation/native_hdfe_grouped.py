"""One-pass projections on rows grouped by a fixed effect, bit for bit against the general passes."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_grouped_native():
    source = (ROOT / 'validation/native/hdfe_grouped.c').read_text(encoding='utf-8')
    with tempfile.TemporaryDirectory(prefix='ctools-hdfe-grouped-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'hdfe_grouped_' + str(parallel), source, openmp=parallel)
