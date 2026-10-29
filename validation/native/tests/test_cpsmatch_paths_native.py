"""Scan-oracle parity, sentinels, radius gaps and allocation-failure ownership."""

import tempfile
from pathlib import Path

from validation.native.support.runner import ROOT, harness, run_test


def main():
    source = (ROOT / "validation/native/memory_safety/alloc.h").read_text()
    source += "\n" + harness("cpsmatch_paths.c")
    with tempfile.TemporaryDirectory(prefix="ctools-cpsmatch-paths-") as tmp:
        for omp in (False, True):
            run_test(Path(tmp), "cpsmatch_paths_" + str(omp), source, openmp=omp)
        oom = (ROOT / "validation/native/memory_safety/alloc.h").read_text()
        oom += "\n" + harness("cpsmatch_paths_oom.c")
        run_test(Path(tmp), "cpsmatch_paths_oom", oom, openmp=True)


if __name__ == "__main__":
    main()
