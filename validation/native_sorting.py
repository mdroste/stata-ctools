"""Native checks; execution and compiler setup are shared in validation.native.

Functions retain their independent fixtures, failure sweeps, and runtime variants.
"""

import argparse
import os
import platform
import shlex
import subprocess
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import compile_native, harness, openmp_prefix


def test_counting_sort_native():
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


def test_csort_ordered_native():
    source_dir = Path(os.environ.get("CTOOLS_TEST_SOURCE", ROOT / "src"))
    with tempfile.TemporaryDirectory(prefix="ctools-csort-ordered-") as tmp:
        directory = Path(tmp)
        source = directory / "test.c"
        source.write_text(harness("csort_ordered_native.c"))
        for omp in (False, True):
            binary = directory / ("test-" + str(omp))
            flags = [
                "-std=c11",
                "-O1",
                "-g",
                "-fno-fast-math",
                "-fno-omit-frame-pointer",
                "-ffunction-sections",
                "-fdata-sections",
                "-pthread",
                "-I",
                str(source_dir),
            ]
            sanitizers = os.environ.get("CTOOLS_SANITIZERS", "undefined")
            if sanitizers:
                flags += ["-fsanitize=" + sanitizers]
            if platform.system() == "Darwin":
                flags += [
                    "-DSYSTEM=APPLEMAC",
                    "-D_DARWIN_C_SOURCE",
                    "-Wl,-dead_strip",
                    "-Wl,-undefined,dynamic_lookup",
                ]
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
                flags += ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE", "-Wl,--gc-sections"]
                if omp:
                    flags += ["-fopenmp"]
            compiler = shlex.split(os.environ.get("CC", "clang"))
            libraries = [] if platform.system() == "Darwin" else ["-lm"]
            compile_native(
                [*compiler, *flags, str(source), *libraries, "-o", str(binary)], check=True
            )
            for runtime_limit in (8, 2) if omp else (1,):
                subprocess.run(
                    [str(binary)],
                    check=True,
                    timeout=120,
                    env=dict(
                        os.environ,
                        ASAN_OPTIONS="detect_leaks=0:abort_on_error=1",
                        UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1",
                        OMP_THREAD_LIMIT=str(runtime_limit),
                    ),
                )
                print(f"PASS csort ordered-input openmp={omp} limit={runtime_limit}", flush=True)


def test_msd_sort_native():
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


def test_order_parallel():
    SOURCE = harness("order_parallel_source.c")
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=ROOT)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="ctools-order-") as tmp:
        cfile = Path(tmp) / "test.c"
        binary = Path(tmp) / "test"
        cfile.write_text(SOURCE)
        flags = [
            "-std=c11",
            "-O2",
            "-g",
            "-fsanitize=undefined",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-Wall",
            "-Wextra",
            "-Werror",
            "-I",
            str(args.source_root / "src"),
        ]
        if platform.system() == "Darwin":
            prefix = openmp_prefix()
            flags += [
                "-DSYSTEM=APPLEMAC",
                "-D_DARWIN_C_SOURCE",
                "-Xpreprocessor",
                "-fopenmp",
                "-I",
                str(prefix / "include"),
            ]
            libs = [str(prefix / "lib/libomp.a")]
        else:
            flags += ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE", "-fopenmp"]
            libs = []
        compile_native(
            [
                os.environ.get("CC", "clang"),
                *flags,
                str(cfile),
                str(args.source_root / "src/ctools_order.c"),
                *libs,
                "-lm",
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary)], check=True, timeout=120)
        subprocess.run(
            [str(binary)], env={**os.environ, "OMP_THREAD_LIMIT": "2"}, check=True, timeout=120
        )


def test_p1_native():
    CC = os.environ.get("CC", "clang")
    COMMON = harness("p1_native_common.c")

    def run_test(directory, name, source, *, openmp=False, sanitize=True):
        cfile = directory / f"{name}.c"
        binary = directory / name
        cfile.write_text(source)
        flags = ["-std=c11", "-O3", "-fno-fast-math", "-ffp-contract=off", "-I", str(ROOT / "src")]
        if sanitize:
            flags += [
                "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
                "-fno-omit-frame-pointer",
            ]
        libraries = []
        if platform.system() == "Darwin":
            flags += [
                "-DSYSTEM=APPLEMAC",
                "-D_DARWIN_C_SOURCE",
                "-Wl,-dead_strip",
                "-Wl,-undefined,dynamic_lookup",
            ]
            if openmp:
                prefix = os.environ.get("LIBOMP_PREFIX")
                if not prefix:
                    prefix = subprocess.check_output(
                        ["brew", "--prefix", "libomp"], text=True
                    ).strip()
                flags += [
                    "-Xpreprocessor",
                    "-fopenmp",
                    "-I",
                    prefix + "/include",
                    "-L",
                    prefix + "/lib",
                    "-lomp",
                ]
        else:
            flags += [
                "-DSYSTEM=STUNIX",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ]
            libraries = ["-lm"]
            if openmp:
                flags += ["-fopenmp"]
        compile_native([CC, *flags, str(cfile), *libraries, "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True, timeout=60)
        print(f"PASS {name}", flush=True)

    with tempfile.TemporaryDirectory(prefix="ctools-p1-native-") as tmp:
        directory = Path(tmp)
        run_test(directory, "compensated_sum", harness("p1_native_2.c"))
        for module, allocators in [
            ("ctools_sort_radix_lsd.c", ["radix_context", "string_context"]),
            ("ctools_sort_radix_msd.c", ["msd_context"]),
        ]:
            body = COMMON + f'\n#include "{module}"\nint main(void) {{\n'
            for allocator in allocators:
                body += f"\n    calls = 0; fail_at = 0;\n    {allocator}_t *ctx_{allocator} = {allocator}_alloc(4);\n    assert(ctx_{allocator});\n    int count_{allocator} = calls;\n    {allocator}_free(ctx_{allocator});\n    assert(live == 0);\n    for (int failure = 1; failure <= count_{allocator}; failure++) {{\n        calls = 0; fail_at = failure;\n        ctx_{allocator} = {allocator}_alloc(4);\n        assert(!ctx_{allocator});\n        assert(live == 0);\n    }}\n".replace(
                    f"{allocator}_t",
                    {
                        "radix_context": "radix_sort_context_t",
                        "string_context": "string_sort_context_t",
                        "msd_context": "msd_sort_context_t",
                    }[allocator],
                )
            run_test(directory, "test_" + module[:-2], body + "return 0; }\n")
        run_test(
            directory, "permutation_allocations", COMMON + harness("p1_native_3.c"), openmp=True
        )


def test_sort_keys_native():
    with tempfile.TemporaryDirectory(prefix="ctools-sort-keys-native-") as tmp:
        directory = Path(tmp)
        source = directory / "test.c"
        source.write_text(harness("sort_keys_native.c"))
        for omp in (False, True):
            binary = directory / ("test-" + str(omp))
            flags = [
                "-std=c11",
                "-O1",
                "-g",
                "-fno-fast-math",
                "-fno-omit-frame-pointer",
                "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
                "-pthread",
                "-I",
                str(ROOT / "src"),
            ]
            if platform.system() == "Darwin":
                flags += [
                    "-DSYSTEM=APPLEMAC",
                    "-D_DARWIN_C_SOURCE",
                    "-Wl,-dead_strip",
                    "-Wl,-undefined,dynamic_lookup",
                ]
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
                flags += [
                    "-DSYSTEM=STUNIX",
                    "-D_GNU_SOURCE",
                    "-ffunction-sections",
                    "-fdata-sections",
                    "-Wl,--gc-sections",
                    "-lm",
                ]
                if omp:
                    flags += ["-fopenmp"]
            compile_native(
                [
                    os.environ.get("CC", "clang"),
                    *flags,
                    str(source),
                    str(ROOT / "src/ctools_sort_radix_lsd.c"),
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
            print(f"PASS sort_keys openmp={omp}", flush=True)
