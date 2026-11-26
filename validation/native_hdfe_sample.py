"""Parallel FE remap, factor renumbering, components and nesting contracts."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_hdfe_sample_parallel_native():
    source = (ROOT / "validation/native/hdfe_sample_parallel.c").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="ctools-hdfe-sample-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "hdfe_sample_parallel_" + str(parallel), source, openmp=parallel)
