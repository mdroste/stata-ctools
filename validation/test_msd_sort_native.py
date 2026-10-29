"""MSD task scheduling, mixed keys, and long strings against a stable reference.

Includes the randomized LSD fixtures without editing them, then adds task/team
boundaries, long shared prefixes, and a depth-limited merge case. The source
override permits validating an isolated prototype before applying its patch.
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
    with tempfile.TemporaryDirectory(prefix="ctools-msd-sort-") as temporary:
        directory = Path(temporary)
        source = directory / "test.c"
        source.write_text(harness("msd_sort_native.c"))
        (directory / "sort_keys_native.c").write_text(harness("sort_keys_native.c"))
        source_root = Path(os.environ.get("CTOOLS_TEST_SOURCE", ROOT / "src"))
        implementation = Path(
            os.environ.get("CTOOLS_MSD_SOURCE", source_root / "ctools_sort_radix_msd.c")
        )
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
                    str(implementation),
                    *libraries,
                    "-o",
                    str(binary),
                ],
                check=True,
            )
            environments = ({}, {"OMP_THREAD_LIMIT": "2", "OMP_DYNAMIC": "TRUE"}) if omp else ({},)
            for constrained in environments:
                subprocess.run(
                    [str(binary)],
                    check=True,
                    timeout=300,
                    env=dict(
                        os.environ,
                        ASAN_OPTIONS="detect_leaks=0:abort_on_error=1",
                        UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1",
                        **constrained,
                    ),
                )
                print(f"PASS msd_sort openmp={omp} constrained={bool(constrained)}", flush=True)


if __name__ == "__main__":
    main()
