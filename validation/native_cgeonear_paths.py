"""cgeonear's k-d tree search against brute force, and its distance kernels."""
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import run_test


def test_cgeonear_paths_native():
    source = (ROOT / "validation/native/cgeonear_paths.c").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="ctools-cgeonear-paths-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "cgeonear_paths_" + str(parallel), source, openmp=parallel)
