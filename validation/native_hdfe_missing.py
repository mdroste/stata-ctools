"""Parallel missing screening and stable in-place compaction contracts."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_missing_native():
    source = (ROOT / "validation/native/hdfe_missing.c").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="ctools-hdfe-missing-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "hdfe_missing_" + str(parallel), source, openmp=parallel)
