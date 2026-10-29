"""Check csort's ordered-input fast path and explicit-engine dispatch contract.

The adapter is tested independently of the sorting engines: exact numeric and
byte-string order, ties, extended missings, signed zeros, late inversions and
engine failures, in serial and OpenMP builds. Large ranges check prefix and
worker boundaries, ties and thread limits. Real-Stata payload/transport coverage is in
validate_csort_identity.do. Runs without Stata.
"""
from native.support.compiler import compile_native
from native.support.runner import ROOT, harness, openmp_prefix
import os
from pathlib import Path
import platform
import shlex
import subprocess
import tempfile


def main():
    source_dir = Path(os.environ.get("CTOOLS_TEST_SOURCE", ROOT / "src"))
    with tempfile.TemporaryDirectory(prefix="ctools-csort-ordered-") as tmp:
        directory = Path(tmp)
        source = directory / "test.c"
        source.write_text(harness("csort_ordered_native.c"))
        for omp in (False, True):
            binary = directory / ("test-" + str(omp))
            flags = ["-std=c11", "-O1", "-g", "-fno-fast-math", "-fno-omit-frame-pointer",
                     "-ffunction-sections", "-fdata-sections", "-pthread", "-I", str(source_dir)]
            sanitizers = os.environ.get("CTOOLS_SANITIZERS", "undefined")
            if sanitizers:
                flags += ["-fsanitize=" + sanitizers]
            if platform.system() == "Darwin":
                flags += ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE", "-Wl,-dead_strip",
                          "-Wl,-undefined,dynamic_lookup"]
                if omp:
                    prefix = str(openmp_prefix())
                    flags += ["-Xpreprocessor", "-fopenmp", "-I", prefix + "/include",
                              prefix + "/lib/libomp.a"]
            else:
                flags += ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE", "-Wl,--gc-sections"]
                if omp:
                    flags += ["-fopenmp"]
            compiler = shlex.split(os.environ.get("CC", "clang"))
            libraries = [] if platform.system() == "Darwin" else ["-lm"]
            compile_native([*compiler, *flags, str(source), *libraries, "-o", str(binary)],
                           check=True)
            for runtime_limit in ((8, 2) if omp else (1,)):
                subprocess.run([str(binary)], check=True, timeout=120,
                               env=dict(os.environ, ASAN_OPTIONS="detect_leaks=0:abort_on_error=1",
                                        UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1",
                                        OMP_THREAD_LIMIT=str(runtime_limit)))
                print(f"PASS csort ordered-input openmp={omp} limit={runtime_limit}", flush=True)


if __name__ == "__main__":
    main()
