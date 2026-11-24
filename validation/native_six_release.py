"""Native failure regressions complementing the six-command Stata suite."""
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import harness, run_test


def test_six_release_native():
    source = (ROOT / "validation/native/cbsample_idpass.c").read_text(encoding="utf-8")
    source = source.replace(
        '#include "newcommands_mock.h"',
        '#include "' + (ROOT / "validation/native/newcommands_mock.h").as_posix() + '"',
    )
    with tempfile.TemporaryDirectory(prefix="ctools-six-release-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "bootstrap_idpass_" + str(parallel), source,
                     modules=("ctools_parse.c",), openmp=parallel)
            run_test(Path(tmp), "range_release_" + str(parallel),
                     harness("memory_rangestat"), openmp=parallel)
