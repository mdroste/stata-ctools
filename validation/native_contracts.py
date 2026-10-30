"""Native checks; execution and compiler setup are shared in validation.native.

Functions retain their independent fixtures, failure sweeps, and runtime variants.
"""

import argparse
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path

from validation import ROOT
from validation.native import (
    CC,
    SPI,
    compile_native,
    harness,
    run_test,
    run_transfer,
    sort_source,
    unit,
)


def test_bigdata_native():
    SOURCE = SPI + harness("bigdata_native_1.c")
    with tempfile.TemporaryDirectory(prefix="ctools-bigdata-native-") as tmp:
        for omp in (False, True):
            run_test(
                Path(tmp),
                "sorted_store_" + str(omp),
                SOURCE,
                ["ctools_types.c", "ctools_arena.c", "ctools_threads.c"],
                openmp=omp,
            )


def test_contracts_native():
    SOURCE = harness("contracts_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-contracts-") as tmp:
        c = Path(tmp) / "test.c"
        exe = Path(tmp) / "test"
        c.write_text(SOURCE)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
            "-I" + str(ROOT / "src"),
            "-pthread",
            "-ffunction-sections",
            "-fdata-sections",
        ]
        flags += (
            ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
            if platform.system() == "Darwin"
            else ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE"]
        )
        flags += ["-Wl,-dead_strip"] if platform.system() == "Darwin" else ["-Wl,--gc-sections"]
        compile_native(
            [
                os.environ.get("CC", "clang"),
                *flags,
                str(c),
                *[
                    str(ROOT / "src" / p)
                    for p in [
                        "ctools_parse.c",
                        "ctools_ols.c",
                        "ctools_runtime.c",
                        "ctools_threads.c",
                        "ctools_arena.c",
                    ]
                ],
                "-lm",
                "-o",
                str(exe),
            ],
            check=True,
        )
        subprocess.run([str(exe)], check=True)


def test_cppmlhdfe_native():
    SOURCE = harness("cppmlhdfe_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-ppml-native-") as tmp:
        src = Path(tmp) / "test.c"
        binary = Path(tmp) / "test"
        src.write_text(SOURCE)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-fno-fast-math",
            "-fno-omit-frame-pointer",
            "-fsanitize=address,undefined",
            "-I",
            str(ROOT / "src"),
        ]
        if platform.system() == "Darwin":
            flags += ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip"]
        else:
            flags += [
                "-DSYSTEM=STUNIX",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ]
        compile_native(
            [
                os.environ.get(
                    "CC", "/usr/bin/clang" if platform.system() == "Darwin" else "clang"
                ),
                *flags,
                str(src),
                "-lm",
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary)], check=True, timeout=60)
        print(
            "PASS: native regressor/FE separation, overlap, tiny positive outcomes, heterogeneous slopes, simplex, LSMR (ASan/UBSan)"
        )


def test_cpsmatch_paths_native():
    source = (ROOT / "validation/native/alloc.h").read_text()
    source += "\n" + harness("cpsmatch_paths.c")
    with tempfile.TemporaryDirectory(prefix="ctools-cpsmatch-paths-") as tmp:
        for omp in (False, True):
            run_test(Path(tmp), "cpsmatch_paths_" + str(omp), source, openmp=omp)
        oom = (ROOT / "validation/native/alloc.h").read_text()
        oom += "\n" + harness("cpsmatch_paths_oom.c")
        run_test(Path(tmp), "cpsmatch_paths_oom", oom, openmp=True)


def test_cqreg_native(argv=()):

    def build(directory, output, *, baseline=False, fallback=False):
        flags = [
            "-O3",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-funroll-loops",
            "-Wno-deprecated-declarations",
            "-I" + str(directory),
            "-I" + str(ROOT / "src"),
            "-I" + str(ROOT / "src/cqreg"),
        ]
        link = ["-lm"]
        if platform.system() == "Darwin":
            flags += ["-DSYSTEM=APPLEMAC"]
            link += ["-framework", "Accelerate"]
        else:
            flags += ["-DSYSTEM=STUNIX"]
        if fallback:
            flags += ["-DUSE_BLAS=0"]
        if baseline:
            flags += ["-DCQREG_BASELINE"]
        if os.environ.get("CTOOLS_SANITIZERS"):
            flags += [
                "-fsanitize=" + os.environ["CTOOLS_SANITIZERS"],
                "-fno-sanitize-recover=all",
                "-g",
            ]
        prefix = os.environ.get("CTOOLS_TEST_OPENMP_PREFIX")
        if prefix:
            flags += ["-Xpreprocessor", "-fopenmp", "-I" + prefix + "/include"]
            link += [prefix + "/lib/libomp.a"]
        elif platform.system() != "Darwin":
            flags += ["-fopenmp"]
        sources = [unit(output.parent, "cqreg_native")]
        counted_types = output.with_suffix(".c")
        counted_types.write_text(
            f"""#include <stdlib.h>\nextern int cqreg_record_alloc(void **, size_t, size_t);\n#define posix_memalign cqreg_record_alloc\n#include "{directory / 'cqreg_types.c'}"\n"""
        )
        sources += [counted_types]
        sources += [directory / ("cqreg_" + s + ".c") for s in ("linalg", "blas", "sparsity")]
        compile_native(
            [os.environ.get("CC", "cc"), *flags, *map(str, sources), *link, "-o", str(output)],
            check=True,
            timeout=60,
        )

    def run(binary, args=()):
        env = dict(os.environ, OMP_NUM_THREADS="4")
        result = subprocess.run(
            [str(binary), *map(str, args)],
            env=env,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        print(result.stdout, end="", flush=True)
        return [line.split(",") for line in result.stdout.splitlines()]

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--rows", type=int, default=20000)
    parser.add_argument("--columns", type=int, default=6)
    parser.add_argument("--baseline-dir", type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args(argv)
    if not 100 <= args.rows <= 200000 or not 2 <= args.columns <= 32:
        parser.error("use 100..200000 rows and 2..32 columns to keep runs bounded")
    if not 1 <= args.repetitions <= 5:
        parser.error("use 1..5 repetitions")
    with tempfile.TemporaryDirectory(prefix="cqreg-native-") as tmp:
        tmp = Path(tmp)
        current = tmp / "current"
        build(ROOT / "src/cqreg", current)
        run(current)
        fallback = tmp / "fallback"
        build(ROOT / "src/cqreg", fallback, fallback=True)
        run(fallback)
        if args.benchmark:
            if args.baseline_dir:
                baseline = tmp / "baseline"
                build(args.baseline_dir.resolve(), baseline, baseline=True)
            new_times, old_times = ({}, {})
            for _ in range(args.repetitions):
                current_results = run(current, (args.rows, args.columns))
                if args.baseline_dir:
                    old_results = run(baseline, (args.rows, args.columns))
                    assert len(current_results) == len(old_results)
                    for new, old in zip(current_results, old_results):
                        assert new[:3] == old[:3]
                        if new[0] == "workspace":
                            assert int(new[4]) < int(old[4])
                            print(
                                f"workspace allocation saved: {int(old[4]) - int(new[4]):,} bytes"
                            )
                            continue
                        assert new[4:] == old[4:], (new, old)
                        new_times.setdefault(new[0], []).append(float(new[3]))
                        old_times.setdefault(old[0], []).append(float(old[3]))
            for name in new_times:
                before = statistics.median(old_times[name])
                after = statistics.median(new_times[name])
                print(f"{name}: {before:.6f}s -> {after:.6f}s ({before / after:.2f}x)")


def test_encode_labels_native():
    with tempfile.TemporaryDirectory(prefix="ctools-encode-labels-") as tmp:
        run_test(Path(tmp), "encode_label_file", harness("encode_label_file.c"), ["ctools_arena.c"])


def test_group_windows_native():
    with tempfile.TemporaryDirectory(prefix="ctools-group-moments-") as tmp:
        run_test(
            Path(tmp),
            "group_moments",
            harness("group_moments.c"),
            ["crangestat/crangestat_kernels.c"],
        )


def test_hdfe_failures_native():
    FIXTURES = ROOT / "validation/native"
    apple = platform.system() == "Darwin"
    env = dict(os.environ, ASAN_OPTIONS="detect_leaks=0", UBSAN_OPTIONS="halt_on_error=1")
    with tempfile.TemporaryDirectory(prefix="ctools-hdfe-failures-") as tmp:

        def build(name):
            binary = Path(tmp) / name
            flags = (
                ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip"]
                if apple
                else [
                    "-DSYSTEM=OPUNIX",
                    "-ffunction-sections",
                    "-fdata-sections",
                    "-Wl,--gc-sections",
                ]
            )
            compile_native(
                [
                    os.environ.get("CC", "/usr/bin/clang" if apple else "cc"),
                    "-std=gnu11",
                    "-O1",
                    "-g",
                    "-fno-fast-math",
                    "-ffp-contract=off",
                    "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "address,undefined"),
                    "-I",
                    str(ROOT / "src"),
                    "-I",
                    str(FIXTURES),
                    *flags,
                    str(unit(Path(tmp), "memory_" + name)),
                    "-lm",
                    "-pthread",
                    "-o",
                    str(binary),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=ROOT,
            )
            return binary

        def run(binary, *args):
            result = subprocess.run(
                [str(binary), *map(str, args)],
                check=True,
                capture_output=True,
                text=True,
                env=env,
                timeout=30,
            )
            return int(re.search("allocations=(\\d+)", result.stdout)[1])

        iv = build("iv")
        cases = 0
        for variant in [(0, 0, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]:
            count = run(iv, 0, *variant)
            for failure in range(1, count + 1):
                run(iv, failure, *variant)
                cases += 1
            print(
                f"PASS IV {variant}: {count} allocation failures, estimates and recovery",
                flush=True,
            )
        for save in (1, 2):
            run(iv, 0, 0, 0, 0, save)
            cases += 1
        for missing in (
            "__civreghdfe_has_cluster2",
            "__civreghdfe_tolerance",
            "__civreghdfe_center",
        ):
            run(iv, 0, 0, 0, 0, 0, missing)
            cases += 1
        print("PASS IV: failed scalar saves and unset configuration scalars are errors", flush=True)
        bins = build("bins_hdfe")
        for weighted in (0, 1):
            for first in (0, 1, 2):
                run(bins, 0, weighted, first)
            count = run(bins, 0, weighted)
            for failure in range(1, count + 1):
                run(bins, failure, weighted)
                cases += 1
            for maxiter in (0, 1):
                run(bins, 0, weighted, 0, maxiter)
        print(
            f"PASS batch HDFE: dummy-OLS parity, weighted orthogonality, column order, iteration limits; {cases} fault injections",
            flush=True,
        )


def test_match_winsor_optimizations_native():
    with tempfile.TemporaryDirectory(prefix="ctools-match-winsor-") as tmp:
        for omp in (False, True):
            run_test(
                Path(tmp),
                "match_winsor_" + str(omp),
                harness("match_winsor_optimizations.c"),
                openmp=omp,
            )


def test_memory_safety_native():
    FIXTURES = ROOT / "validation/native"
    SYSTEM = "APPLEMAC" if platform.system() == "Darwin" else "OPUNIX"
    CC = os.environ.get("CC", "/usr/bin/clang" if SYSTEM == "APPLEMAC" else "cc")
    SAN = os.environ.get("CTOOLS_SANITIZERS", "address,undefined,float-cast-overflow")
    ENV = dict(
        os.environ,
        UBSAN_OPTIONS="halt_on_error=1",
        ASAN_OPTIONS="detect_leaks=0",
        OMP_NUM_THREADS="4",
        OMP_THREAD_LIMIT="4",
    )
    with tempfile.TemporaryDirectory(prefix="ctools-memory-safety-") as tmp:

        def build(name, source=None, flags=()):
            exe = Path(tmp) / name
            omp = os.environ.get("CTOOLS_TEST_OPENMP_PREFIX")
            parallel = (
                ["-Xpreprocessor", "-fopenmp", "-I" + omp + "/include", omp + "/lib/libomp.a"]
                if omp
                else ["-fopenmp"] if SYSTEM != "APPLEMAC" else []
            )
            command = [
                CC,
                "-std=gnu11",
                "-O1",
                "-g",
                f"-DSYSTEM={SYSTEM}",
                "-I",
                str(ROOT / "src"),
                "-I",
                str(FIXTURES),
                f"-fsanitize={SAN}",
                "-ffunction-sections",
                "-fdata-sections",
                *flags,
                *parallel,
                str(unit(Path(tmp), "memory_" + (source or name))),
                "-o",
                str(exe),
                "-Wl,-dead_strip" if SYSTEM == "APPLEMAC" else "-Wl,--gc-sections",
                "-lm",
                "-pthread",
            ]
            compile_native(command, check=True, cwd=ROOT, capture_output=True, text=True)
            return exe

        def run(exe, *args):
            result = subprocess.run(
                [str(exe), *map(str, args)],
                cwd=ROOT,
                env=ENV,
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            m = re.search("rc=(\\d+) (?:allocations|calls)=(\\d+) alive=(\\d+)", result.stdout)
            if m:
                rc, calls, alive = map(int, m.groups())
                assert alive == 0, result.stdout
                return (rc, calls)
            return None

        cases = 0
        for name in ("regress", "ppml", "sampling", "bootstrap", "rangestat"):
            exe = build(
                name,
                "sampling" if name == "bootstrap" else name,
                ["-DBOOTSTRAP"] if name == "bootstrap" else [],
            )
            variants = (
                [(0, 0, 0), (1, 1, 0), (1, 2, 0), (1, 1, 2), (0, 0, 1), (1, 1, 3)]
                if name == "regress"
                else [(0, 0, 0), (1, 1, 0)] if name == "ppml" else [()]
            )
            for variant in variants:
                rc, count = run(exe, 0, *variant)
                assert rc == 0, (name, variant, rc)
                for failure in range(1, count + 1):
                    rc, _ = run(exe, failure, *variant)
                    assert rc in (0, 920, 430), (name, variant, failure, rc)
                    cases += 1
                if name != "ppml":
                    for store in (1, 6, 12):
                        rc, _ = run(exe, 0, *variant, store)
                        assert rc == 459, (name, variant, store, rc)
                        cases += 1
            if name in ("regress", "ppml"):
                for scalar in (1, 8, 20) if name == "regress" else (1, 2):
                    assert run(exe, 0, 1, 1, 0, 0, scalar)[0] == 459
                    cases += 1
                for matrix in (1, 2, 4) if name == "regress" else (1, 2):
                    assert run(exe, 0, 1, 1, 0, 0, 0, matrix)[0] == 459
                    cases += 1
            print(f"{name}: allocation and output failures passed", flush=True)
        for name in ("remap", "timsort", "threads"):
            run(build(name))
            print(f"{name}: passed", flush=True)
        exe = build("styles")
        for failure in range(5):
            run(exe, failure)
        print(f"styles: all parser/growth failures passed; {cases} command injections", flush=True)


def test_newcommands_native():
    SOURCE = harness("newcommands_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-new-native-") as tmp:
        cfile = Path(tmp) / "test.c"
        binary = Path(tmp) / "test"
        cfile.write_text(SOURCE)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
            "-fno-omit-frame-pointer",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-Wno-unknown-pragmas",
            "-I",
            str(ROOT / "src"),
        ]
        flags += (
            ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
            if platform.system() == "Darwin"
            else ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE"]
        )
        compile_native(
            [
                os.environ.get("CC", "clang"),
                *flags,
                str(cfile),
                str(ROOT / "src/ctools_parse.c"),
                "-lm",
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary)], check=True, timeout=60)


def test_p2_native():
    with tempfile.TemporaryDirectory(prefix="ctools-p2-native-") as tmp:
        directory = Path(tmp)
        run_test(directory, "signed_labels", harness("p2_native_1.c"), ["ctools_arena.c"])
        run_test(directory, "encoding_bom", harness("p2_native_2.c"))
        run_test(
            directory,
            "spi_failures",
            SPI + harness("p2_native_3.c"),
            ["ctools_types.c", "ctools_arena.c", "ctools_threads.c"],
            openmp=True,
        )
        run_test(
            directory,
            "csv_store_failures",
            SPI + harness("p2_native_4.c"),
            ["ctools_types.c", "cimport/cimport_parse.c"],
        )
        run_test(directory, "xlsx_store_failures", SPI + harness("p2_native_5.c"))
        run_test(directory, "matrix_store_failures", SPI + harness("p2_native_6.c"))


def test_refactor_boundaries_native():
    subprocess.run(["python3", str(ROOT / "scripts/sync_architecture.py"), "--check"], check=True)
    headers = [
        "ctools_data.h",
        "ctools_data_io.h",
        "ctools_sort.h",
        "ctools_convert.h",
        "ctools_alloc.h",
        "ctools_platform.h",
        "ctools_rng.h",
        "ctools_group.h",
        "ctools_linalg.h",
        "ctools_text.h",
        "formats/ctools_output.h",
        "formats/excel_address.h",
        "hdfe/hdfe_types.h",
        "hdfe/hdfe_solver.h",
        "hdfe/hdfe_remap.h",
    ]
    with tempfile.TemporaryDirectory(prefix="ctools-boundaries-") as tmp:
        directory = Path(tmp)
        for header in headers:
            source = directory / "header.c"
            source.write_text(f'#include "{header}"\nint main(void) {{ return 0; }}\n')
            system = "APPLEMAC" if platform.system() == "Darwin" else "STUNIX"
            subprocess.run(
                [
                    CC,
                    "-std=c11",
                    "-fsyntax-only",
                    "-DSYSTEM=" + system,
                    "-I",
                    str(ROOT / "src"),
                    str(source),
                ],
                check=True,
            )
        print("PASS standalone public headers", flush=True)
        for name in ("creghdfe/creghdfe_regress.c", "cppmlhdfe/cppmlhdfe_irls.c"):
            macros = subprocess.check_output(
                [
                    CC,
                    "-std=c11",
                    "-E",
                    "-dM",
                    "-DSD_FASTMODE",
                    "-DSYSTEM=" + system,
                    "-I",
                    str(ROOT / "src"),
                    str(ROOT / "src" / name),
                ],
                text=True,
            )
            for macro, callback in (("SF_vdata", "safevdata"), ("SF_vstore", "safestore")):
                definition = next(
                    (
                        line
                        for line in macros.splitlines()
                        if line.startswith("#define " + macro + "(")
                    )
                )
                assert callback in definition, (name, definition)
        print("PASS estimator checked-SPI include order", flush=True)
        run_test(
            directory,
            "refactor_boundaries",
            harness("refactor_boundaries.c"),
            ["ctools_group.c", "ctools_platform.c", "ctools_runtime.c", "formats/excel_address.c"],
        )


def test_rowstage_native():
    run_transfer("rowstage_native")


def test_sep22_native():
    SOURCE = harness("sep22_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-sep22-native-") as tmp:
        cfile = Path(tmp) / "test.c"
        binary = Path(tmp) / "test"
        cfile.write_text(SOURCE)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-fno-fast-math",
            "-fno-omit-frame-pointer",
            "-fsanitize=address,undefined",
            "-I",
            str(ROOT / "src"),
        ]
        if platform.system() == "Darwin":
            flags += ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"]
        else:
            flags += [
                "-DSYSTEM=STUNIX",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ]
        compile_native(
            [
                os.environ.get(
                    "CC", "/usr/bin/clang" if platform.system() == "Darwin" else "clang"
                ),
                *flags,
                str(cfile),
                "-lm",
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary)], check=True, timeout=60)
        print(
            "PASS: full-double clusters; original-row selection; singleton chains; every selection allocation failure (ASan/UBSan)"
        )


def test_sep24_native():
    ALLOC = harness("sep24_native_alloc.c")
    VCE = ALLOC + harness("sep24_native_5.c")
    ARENA = ALLOC + harness("sep24_native_6.c")
    POOL = harness("sep24_native_pool.c")
    ATOMIC = harness("sep24_native_atomic.c")
    KEY_GUARD = ALLOC + harness("sep24_native_8.c")
    RANK_SOLVE = ALLOC + harness("sep24_native_9.c")

    def compile_run(tmp, name, source, flags=(), env=None, args=()):
        cfile = tmp / (name + ".c")
        binary = tmp / name
        cfile.write_text(source)
        base = [
            "-std=c11",
            "-O1",
            "-g",
            "-fno-fast-math",
            "-fno-omit-frame-pointer",
            "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
            "-I",
            str(ROOT / "src"),
        ]
        if platform.system() == "Darwin":
            base += [
                "-DSYSTEM=APPLEMAC",
                "-D_DARWIN_C_SOURCE",
                "-Wl,-dead_strip",
                "-Wl,-undefined,dynamic_lookup",
            ]
        else:
            base += [
                "-DSYSTEM=STUNIX",
                "-D_POSIX_C_SOURCE=200809L",
                "-ffunction-sections",
                "-Wl,--gc-sections",
            ]
        compile_native(
            [
                os.environ.get(
                    "CC", "/usr/bin/clang" if platform.system() == "Darwin" else "clang"
                ),
                *base,
                str(cfile),
                *flags,
                "-pthread",
                "-lm",
                "-o",
                str(binary),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        run_env = dict(
            os.environ,
            ASAN_OPTIONS="detect_leaks=0:symbolize=0:abort_on_error=1",
            UBSAN_OPTIONS="halt_on_error=1:symbolize=0",
            **env or {},
        )
        subprocess.run(
            [str(binary), *map(str, args)],
            check=True,
            env=run_env,
            timeout=45,
            capture_output=True,
            text=True,
        )
        print("PASS", name, flush=True)

    with tempfile.TemporaryDirectory(prefix="ctools-sep24-native-") as directory:
        tmp = Path(directory)
        for name, source in [
            ("key_type_guard", KEY_GUARD),
            ("rank_solve", RANK_SOLVE),
            ("vce_faults", VCE),
            ("arena_reset", ARENA),
            ("pool_create_failure", POOL),
            ("atomic_output", ATOMIC),
        ]:
            compile_run(
                tmp,
                name,
                source,
                args=[tmp / "literal threads(2) path"] if name == "atomic_output" else [],
            )
        for alg, num, string in [
            (
                "sample",
                "sample_sort_numeric_impl(order,keys,128,4)",
                "sample_sort_string_impl(order,strings,128,4,2)",
            ),
            (
                "merge",
                "parallel_merge_sort_numeric(order,keys,128,4)",
                "parallel_merge_sort_string(order,strings,128,4,2)",
            ),
        ]:
            compile_run(tmp, alg + "_faults", sort_source(alg, numeric=num, string=string))
        for alg in ["sample", "merge", "radix_lsd", "counting", "ips4o"]:
            source = sort_source(alg)
            extra = [str(ROOT / "src/ctools_arena.c")] if alg in ("counting", "ips4o") else []
            compile_run(tmp, alg + "_noomp", source, extra)
            prefix = os.environ.get("CTOOLS_TEST_OPENMP_PREFIX")
            if prefix:
                flags = [
                    "-Xpreprocessor",
                    "-fopenmp",
                    "-I",
                    str(Path(prefix) / "include"),
                    str(Path(prefix) / "lib/libomp.a"),
                ]
                for limit in ["1", "2"]:
                    compile_run(
                        tmp,
                        alg + "_team" + limit,
                        source,
                        flags + extra,
                        {"OMP_THREAD_LIMIT": limit, "OMP_DYNAMIC": "TRUE"},
                    )


def test_spi_matrix_native():
    SOURCE = harness("spi_matrix_native_source.c")
    PUBLICATION = harness("spi_matrix_native_publication.c")
    system = "APPLEMAC" if platform.system() == "Darwin" else "OPUNIX"
    cc = os.environ.get("CC", "/usr/bin/clang" if system == "APPLEMAC" else "cc")
    with tempfile.TemporaryDirectory(prefix="ctools-spi-matrix-") as tmp:
        source = Path(tmp) / "matrix.c"
        source.write_text(SOURCE)
        for mode in ("SD_FASTMODE", "SD_SAFEMODE"):
            binary = Path(tmp) / mode
            compile_native(
                [
                    cc,
                    "-std=c11",
                    "-g",
                    "-O1",
                    "-fsanitize=address,undefined",
                    "-fno-omit-frame-pointer",
                    *(["-DSD_FASTMODE"] if mode == "SD_FASTMODE" else []),
                    f"-DSYSTEM={system}",
                    "-I",
                    str(ROOT / "src"),
                    str(source),
                    str(ROOT / "src/ctools_spi.c"),
                    "-o",
                    str(binary),
                ],
                check=True,
            )
            subprocess.run([str(binary)], check=True, timeout=30)
            print(f"PASS {mode}: invalid matrix writes, injected error, and recovery")
        source.write_text(PUBLICATION)
        binary = Path(tmp) / "publication"
        link = (
            ["-Wl,-dead_strip", "-Wl,-undefined,dynamic_lookup"]
            if system == "APPLEMAC"
            else ["-Wl,--gc-sections", "-lm"]
        )
        compile_native(
            [
                cc,
                "-std=gnu11",
                "-g",
                "-O1",
                "-fsanitize=address,undefined",
                "-ffunction-sections",
                "-fdata-sections",
                "-DSD_FASTMODE",
                f"-DSYSTEM={system}",
                "-I",
                str(ROOT / "src"),
                str(source),
                str(ROOT / "src/ctools_spi.c"),
                *link,
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary)], check=True, timeout=30)
        print(
            "PASS IV/quantile publication: allocation failure, 16 write failures, shapes, recovery"
        )


def test_split_tokens_native():
    SOURCE = harness("newcommands_native_source.c")
    EXTRA = harness("split_tokens_native_extra.c")
    with tempfile.TemporaryDirectory(prefix="ctools-split-tokens-") as tmp:
        cfile = Path(tmp) / "test.c"
        binary = Path(tmp) / "test"
        cfile.write_text(
            "#define CTOOLS_HARNESS_MAIN original_main\n"
            + SOURCE
            + "\n#undef CTOOLS_HARNESS_MAIN\n"
            + EXTRA
        )
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-fsanitize=" + os.environ.get("CTOOLS_SANITIZERS", "undefined"),
            "-fno-omit-frame-pointer",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-Wno-unknown-pragmas",
            "-I",
            str(ROOT / "src"),
        ]
        flags += (
            ["-DSYSTEM=APPLEMAC", "-D_DARWIN_C_SOURCE"]
            if platform.system() == "Darwin"
            else ["-DSYSTEM=STUNIX", "-D_GNU_SOURCE"]
        )
        compile_native(
            [
                os.environ.get("CC", "clang"),
                *flags,
                str(cfile),
                str(ROOT / "src/ctools_parse.c"),
                "-lm",
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run([str(binary)], check=True, timeout=60)


def test_thread_policy_native():
    SYSTEM = platform.system()
    CAP = 16
    SCRUB = ("CTOOLS_NUM_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "XDG_CONFIG_HOME")

    def build(tmp):
        source = tmp / "thread_policy.c"
        source.write_text(harness("thread_policy.c"))
        binary = tmp / ("thread_policy.exe" if SYSTEM == "Windows" else "thread_policy")
        default = "/usr/bin/clang" if SYSTEM == "Darwin" else "cc"
        command = [
            *shlex.split(os.environ.get("CC", default)),
            "-std=gnu11",
            "-O1",
            "-g",
            "-I",
            str(ROOT / "src"),
        ]
        sanitizers = os.environ.get("CTOOLS_SANITIZERS", "undefined")
        if sanitizers:
            command += ["-fno-omit-frame-pointer", "-fsanitize=" + sanitizers]
        command += {"Darwin": ["-DSYSTEM=APPLEMAC"], "Windows": ["-DSYSTEM=STWIN32"]}.get(
            SYSTEM, ["-DSYSTEM=STUNIX"]
        )
        command += [str(source), "-o", str(binary)]
        if SYSTEM != "Windows":
            command += ["-pthread", "-lm"]
        compile_native(command, check=True, capture_output=True, text=True)
        return binary

    def environment(home, **values):
        env = {k: v for k, v in os.environ.items() if k not in SCRUB}
        env.update(
            HOME=str(home),
            APPDATA=str(home / "AppData"),
            ASAN_OPTIONS="detect_leaks=0:abort_on_error=1",
            UBSAN_OPTIONS="halt_on_error=1",
        )
        env.update(values)
        return env

    def settings_path(home, env):
        if SYSTEM == "Windows":
            return str(home / "AppData" / "ctools" / "ctools.conf")
        if SYSTEM == "Darwin":
            return f"{home}/Library/Application Support/ctools/ctools.conf"
        base = env.get("XDG_CONFIG_HOME") or f"{home}/.config"
        return f"{base}/ctools/ctools.conf"

    def report(binary, env, prefix=()):
        out = subprocess.run(
            [*prefix, str(binary), "report"],
            env=env,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        ).stdout
        fields = dict(re.findall("(\\w+)=(.*?)(?=\\s+\\w+=|\\s*$)", out))
        return {k: v if k == "source" else int(v) for k, v in fields.items()}

    def check(name, got, **expected):
        wrong = {k: (got[k], v) for k, v in expected.items() if got[k] != v}
        assert not wrong, f"{name}: {wrong} in {got}"
        assert got["max"] == got["threads"] and got["clamp"] == got["limit"], (name, got)
        print("PASS", name, flush=True)

    with tempfile.TemporaryDirectory(prefix="ctools-thread-policy-") as directory:
        tmp = Path(directory)
        binary = build(tmp)
        homes = [tmp / "home"]
        if SYSTEM == "Linux":
            homes.append(tmp / "xdg-home")
        for home in homes:
            home.mkdir()
            (home / "AppData").mkdir()
            env = environment(home)
            if home.name == "xdg-home":
                env["XDG_CONFIG_HOME"] = str(home / "xdg")
            work = home / "work"
            work.mkdir()
            subprocess.run(
                [str(binary), "unit", str(work), settings_path(home, env)],
                env=env,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            print("PASS unit", home.name, flush=True)
        home = tmp / "cases"
        home.mkdir()
        (home / "AppData").mkdir()
        base = report(binary, environment(home))
        usable = base["available"]
        limit = base["limit"]
        assert 1 <= usable <= base["detected"] and limit == min(usable, 256), base
        check("automatic", base, threads=min(usable, CAP), source="automatic", saved=0)
        check(
            "omp list",
            report(binary, environment(home, OMP_NUM_THREADS="3,2")),
            threads=min(3, limit),
            source="OMP_NUM_THREADS",
            omp=3,
        )
        check(
            "ctools over omp",
            report(binary, environment(home, OMP_NUM_THREADS="3", CTOOLS_NUM_THREADS="2")),
            threads=min(2, limit),
            source="CTOOLS_NUM_THREADS",
        )
        check(
            "explicit above usable",
            report(binary, environment(home, CTOOLS_NUM_THREADS="100000")),
            threads=limit,
            source="CTOOLS_NUM_THREADS",
        )
        check(
            "thread limit",
            report(binary, environment(home, OMP_THREAD_LIMIT="1", CTOOLS_NUM_THREADS="8")),
            threads=1,
            limit=1,
            omp_limit=1,
        )
        check(
            "invalid ignored",
            report(binary, environment(home, CTOOLS_NUM_THREADS="abc", OMP_NUM_THREADS="0")),
            threads=min(usable, CAP),
            source="automatic",
            env=-1,
            omp=-1,
        )
        conf = Path(settings_path(home, environment(home)))
        conf.parent.mkdir(parents=True)
        conf.write_text("threads=5\n")
        check(
            "saved over omp",
            report(binary, environment(home, OMP_NUM_THREADS="3")),
            threads=min(5, limit),
            source="saved setting",
            saved=5,
        )
        check(
            "ctools over saved",
            report(binary, environment(home, CTOOLS_NUM_THREADS="2")),
            threads=min(2, limit),
            source="CTOOLS_NUM_THREADS",
            saved=5,
        )
        conf.write_text("threads=lots\n")
        check(
            "invalid saved ignored",
            report(binary, environment(home)),
            threads=min(usable, CAP),
            source="automatic",
            saved=-1,
        )
        conf.unlink()
        taskset = shutil.which("taskset") if SYSTEM == "Linux" else None
        if taskset:
            check(
                "affinity",
                report(binary, environment(home), [taskset, "-c", "0"]),
                available=1,
                limit=1,
                threads=1,
            )
            check(
                "affinity explicit",
                report(binary, environment(home, CTOOLS_NUM_THREADS="8"), [taskset, "-c", "0"]),
                available=1,
                threads=1,
                source="CTOOLS_NUM_THREADS",
            )


def test_transport_adaptive_native():
    run_transfer("transport_adaptive_native_1", passive=True)


def test_transport_native():
    run_transfer("transport_native_1")


def test_transport_scheduling():
    run_transfer("transport_scheduling_1")


def test_transport_store_native():
    run_transfer("transport_store_native_1")


def test_winsor_transport_native():
    allocator = ROOT / "validation/native/alloc.h"
    source = allocator.read_text() + "\n" + harness("winsor_transport.c")
    with tempfile.TemporaryDirectory(prefix="ctools-winsor-transport-") as tmp:
        run_test(Path(tmp), "winsor_transport", source, openmp=True)
