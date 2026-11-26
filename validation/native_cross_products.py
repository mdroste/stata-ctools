"""Cross products, compensated sums and sandwich meats against sequential references."""
import tempfile
from pathlib import Path
from validation import ROOT
from validation.native import run_test


def test_cross_products_native():
    source = (ROOT / 'validation/native/cross_products.c').read_text(encoding='utf-8')
    source = source.replace(
        '#include "newcommands_mock.h"',
        '#include "' + (ROOT / 'validation/native/newcommands_mock.h').as_posix() + '"',
    )
    # The clustered meat has two paths, chosen by how many cluster sums there
    # are: the default takes the designs with few clusters in one pass per
    # regressor, a budget of a few sums sends them all through the sorted
    # blocks. Both must give the bits of the definition.
    sorted_path = '#define CTOOLS_VCE_CLUSTER_SUMS 8\n' + source
    with tempfile.TemporaryDirectory(prefix='ctools-cross-') as tmp:
        for parallel in (False, True):
            run_test(Path(tmp), 'cross_products_' + str(parallel), source, openmp=parallel)
            run_test(Path(tmp), 'cross_products_sorted_' + str(parallel), sorted_path, openmp=parallel)
