"""Native tests for distinct counting and standalone fixed-effect projection."""

import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import run_test


def test_count_hdfe_native():
    source = (ROOT / "validation/native/count_hdfe.c").read_text()
    source = source.replace('#include "alloc.h"', f'#include "{ROOT}/validation/native/alloc.h"')
    with tempfile.TemporaryDirectory(prefix="ctools-count-hdfe-") as tmp:
        run_test(Path(tmp), "count_hdfe", source)
