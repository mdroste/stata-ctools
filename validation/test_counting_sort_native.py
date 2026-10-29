"""Counting-sort eligibility, thread limits, and stable LSD fallback.

Serial and OpenMP builds cover threshold boundaries, shuffled incoming orders,
multiple keys, signed zero, system/extended missings, range and magnitude
boundaries, and fractional/outlier values outside the bounded sample.
"""

import os
import platform
import shlex
import subprocess
import tempfile
from pathlib import Path

from native.support.compiler import compile_native
from native.support.runner import ROOT, harness, openmp_prefix


def main():
    with tempfile.TemporaryDirectory(prefix="ctools-counting-sort-") as temporary:
        directory = Path(temporary)
        source = directory / "test.c"
        source.write_text(harness("counting_sort_native.c"))
        source_root = Path(os.environ.get("CTOOLS_TEST_SOURCE", ROOT / "src"))
        for omp in (False, True):
            binary = directory / ("test-omp" if omp else "test-serial")
            flags = [
                "-std=c11",
                "-O1",
                "-g",
                "-fno-fast-math",
                "-fno-omit-frame-pointer",
                "-pthread",
                "-I",
                str(source_root),
            ]
            sanitizers = os.environ.get("CTOOLS_SANITIZERS", "undefined")
            if sanitizers:
                flags += ["-fsanitize=" + sanitizers]
            if platform.system() == "Darwin":
                flags += ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE", "-Wl,-dead_strip"]
                if omp:
                    prefix = str(openmp_prefix())
                    flags += [
                        "-Xpreprocessor",
                        "-fopenmp",
                        "-I",
                        prefix + "/include",
                        prefix + "/lib/libomp.a",
                    ]
            else:
                flags += ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE"]
                if omp:
                    flags += ["-fopenmp"]
            compiler = shlex.split(os.environ.get("CC", "clang"))
            libraries = [] if platform.system() == "Darwin" else ["-lm"]
            compile_native(
                [
                    *compiler,
                    *flags,
                    str(source),
                    str(source_root / "ctools_arena.c"),
                    *libraries,
                    "-o",
                    str(binary),
                ],
                check=True,
            )
            subprocess.run(
                [str(binary)],
                check=True,
                timeout=300,
                env=dict(
                    os.environ,
                    ASAN_OPTIONS="detect_leaks=0:abort_on_error=1",
                    UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1",
                ),
            )
            print(f"PASS counting_sort openmp={omp}", flush=True)


if __name__ == "__main__":
    main()
