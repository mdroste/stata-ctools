"""cjoinby key index and output plan against a brute-force join."""
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import run_test


def test_cjoinby_paths_native():
    source = (ROOT / "validation/native/cjoinby_paths.c").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="ctools-cjoinby-paths-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "cjoinby_paths_" + str(parallel), source, openmp=parallel)
