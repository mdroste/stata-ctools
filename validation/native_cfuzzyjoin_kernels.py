"""cfuzzyjoin's string kernels (Levenshtein, soundex, grams, profiles) against references."""
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import run_test


def test_cfuzzyjoin_kernels_native():
    source = (ROOT / "validation/native/cfuzzyjoin_kernels.c").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="ctools-cfuzzyjoin-kernels-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "cfuzzyjoin_kernels_" + str(parallel), source, openmp=parallel)
