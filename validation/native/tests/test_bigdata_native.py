"""Regression checks for fused sort stores: order, ownership, and checked SPI.

Run with python3 -m validation.native.tests.test_bigdata_native. Uses the existing native
test compiler helper and tests both serial and OpenMP/thread-pool execution.
"""
from validation.native.support.runner import harness
from pathlib import Path
import tempfile

from validation.native.support.runner import SPI, run_test


SOURCE = SPI + harness("bigdata_native_1.c")


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-bigdata-native-") as tmp:
        for omp in (False, True):
            run_test(Path(tmp), "sorted_store_" + str(omp), SOURCE,
                     ["ctools_types.c", "ctools_arena.c", "ctools_threads.c"], openmp=omp)


if __name__ == "__main__":
    main()
