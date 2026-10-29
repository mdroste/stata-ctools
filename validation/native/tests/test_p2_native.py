"""P2 hash, encoding, SPI fault-injection, and compiler-selection regressions."""
from validation.native.support.runner import harness
import os
from pathlib import Path
import platform
import subprocess
import tempfile

from validation import ROOT
from validation.native.support.runner import SPI, run_test






def main():
    with tempfile.TemporaryDirectory(prefix="ctools-p2-native-") as tmp:
        directory = Path(tmp)
        run_test(directory, "signed_labels", harness("p2_native_1.c"), ["ctools_arena.c"])
        run_test(directory, "encoding_bom", harness("p2_native_2.c"))
        run_test(directory, "spi_failures", SPI + harness("p2_native_3.c"), ["ctools_types.c", "ctools_arena.c", "ctools_threads.c"], openmp=True)

        run_test(directory, "csv_store_failures", SPI + harness("p2_native_4.c"), ["ctools_types.c", "cimport/cimport_parse.c"])
        run_test(directory, "xlsx_store_failures", SPI + harness("p2_native_5.c"))
        run_test(directory, "matrix_store_failures", SPI + harness("p2_native_6.c"))


if __name__ == "__main__":
    main()
