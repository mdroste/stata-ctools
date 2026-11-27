"""HDFE sample construction that shares its steps' work, against a reference."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test

# The second set of thresholds takes the paths the first (the file's own)
# keeps short: whole-data probes, serial ranges of levels, tables per block,
# and singleton removal that counts rather than lists from a few single rows
# on, with a table per block for every factor it or the PPML selection counts.
OTHER = (
    "#define CTOOLS_HDFE_BLOCK_ROWS 64\n"
    "#define HDFE_SINGLETON_PROBE_ROWS 100\n"
    "#define HDFE_LISTED_MIN_DROPPED 3\n"
    "#define HDFE_LISTED_MIN_SEEDS 2\n"
    "#define HDFE_ORDER_RANGE_LEVELS 1000000\n"
    "#define REMAP_SCATTER_CODES 32\n"
    "#define HDFE_COUNT_FAR_SLOTS 64\n"
    "#define HDFE_COUNT_BLOCK_ENTRIES 4096\n"
    "#define HDFE_SEEDED_TABLE_ENTRIES 1048576\n"
    "#define PPML_COUNT_TABLE_LEVELS 1048576\n"
)


def test_hdfe_sample_ordered_native():
    source = (ROOT / "validation/native/hdfe_sample_ordered.c").read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix="ctools-hdfe-ordered-") as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), "hdfe_sample_ordered_" + str(parallel), source, openmp=parallel)
        run_test(Path(tmp), "hdfe_sample_ordered_other", OTHER + source, openmp=True)
