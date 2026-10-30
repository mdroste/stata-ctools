"""Compile native harnesses with the production modules they exercise.

White-box harnesses deliberately include implementation files to interpose SPI
and allocation failures. Keep extracted companion modules in that same test
translation unit so every allocation remains covered by fault injection.
Ordinary plugin builds compile these modules separately.
"""

import os
import platform
import re
import subprocess
import tempfile
from pathlib import Path

from validation import ROOT

SRC = ROOT / "src"
COMPANIONS = {
    "cimport/cimport_impl.c": ["cimport/cimport_cache.c"],
    "civreghdfe/civreghdfe_impl.c": [
        "civreghdfe/civreghdfe_result.c",
        "civreghdfe/civreghdfe_input.c",
    ],
    "creghdfe/creghdfe_regress.c": ["creghdfe/creghdfe_result.c", "creghdfe/creghdfe_input.c"],
    "crangestat/crangestat_impl.c": ["crangestat/crangestat_kernels.c"],
    "cpsmatch/cpsmatch_impl.c": [
        "cpsmatch/cpsmatch_kernels.c",
        "cpsmatch/cpsmatch_options.c",
        "cpsmatch/cpsmatch_stats.c",
        "cpsmatch/cpsmatch_radius.c",
    ],
    "ctools_data_io.c": ["ctools_transport.c", "ctools_data_plan.c", "ctools_data_store.c"],
    "ctools_types.c": ["ctools_platform.c"],
    "csort/csort_stream.c": ["ctools_platform.c"],
    "csample/csample_impl.c": ["ctools_group.c"],
    "cbsample/cbsample_impl.c": ["ctools_group.c"],
    "cexport/cexport_parse.c": ["formats/ctools_output.c"],
    "cexport/cexport_xlsx.c": ["ctools_text.c", "formats/excel_address.c"],
    "formats/xlsx_xml.c": ["formats/excel_address.c"],
    "cppmlhdfe/cppmlhdfe_irls.c": [
        "cppmlhdfe/cppmlhdfe_result.c",
        "cppmlhdfe/cppmlhdfe_input.c",
        "cppmlhdfe/cppmlhdfe_options.c",
        "cppmlhdfe/cppmlhdfe_solvers.c",
        "cppmlhdfe/cppmlhdfe_simplex.c",
        "cppmlhdfe/cppmlhdfe_relu.c",
        "cppmlhdfe/cppmlhdfe_moments.c",
    ],
    "io/cio.c": [
        "io/cio_context.c",
        "io/cio_readstat.c",
        "io/cio_dbf.c",
        "io/cio_metadata.c",
        "io/cio_excel.c",
        "io/cio_xls.c",
        "io/cio_xlsx.c",
        "io/cio_xport.c",
        "io/cio_xpf.c",
        "io/cio_shp.c",
        "formats/excel_address.c",
    ],
}


def included_sources(path, includes, seen=None):
    seen = set() if seen is None else seen
    path = path.resolve()
    if path in seen or not path.is_file():
        return seen
    seen.add(path)
    if path.name.endswith(".suite.c"):
        return seen  # Selected sections declare their own companion includes.
    for name in re.findall(r'^\s*#include\s+"([^"\n]+\.c)"', path.read_text(), re.M):
        for folder in [path.parent, *includes]:
            candidate = folder / name
            if candidate.is_file():
                included_sources(candidate, includes, seen)
                break
    return seen


def compile_native(command, **kwargs):
    provided = kwargs.pop("provided", ())
    command = list(map(str, command))
    source_root = Path(os.environ.get("CTOOLS_TEST_SOURCE", SRC)).resolve()
    includes = []
    for i, part in enumerate(command):
        if part == "-I":
            includes.append(Path(command[i + 1]))
        elif part.startswith("-I"):
            includes.append(Path(part[2:]))
    includes.append(source_root)
    sources = [
        (i, Path(p)) for i, p in enumerate(command) if p.endswith(".c") and Path(p).is_file()
    ]
    all_present = set().union(*(included_sources(p, includes) for _, p in sources))
    all_present.update((source_root / p).resolve() for p in provided)
    with tempfile.TemporaryDirectory(prefix="ctools-native-modules-") as tmp:
        for index, path in sources:
            present = included_sources(path, includes)
            additions = {}
            for primary, companions in COMPANIONS.items():
                if (source_root / primary).resolve() in present:
                    for companion in companions:
                        candidate = (source_root / companion).resolve()
                        if candidate.is_file() and candidate not in all_present:
                            additions.setdefault((source_root / primary).resolve(), []).append(
                                candidate
                            )
                            all_present.add(candidate)
            if additions:
                wrapper = Path(tmp) / f"unit-{index}.c"

                def include(match):
                    for folder in [path.parent, *includes]:
                        candidate = (folder / match[1]).resolve()
                        if candidate.is_file():
                            # Inject companions before the harness undefines its
                            # allocator/SPI macros, not at the end of the file.
                            return "\n".join(
                                f'#include "{p}"'
                                for p in [candidate, *additions.pop(candidate, [])]
                            )
                    return match[0]

                source = re.sub(
                    r'^\s*#include\s+"([^"\n]+)"', include, path.read_text(), flags=re.M
                )
                tail = [p for group in additions.values() for p in group]
                wrapper.write_text(
                    f'#line 1 "{path.resolve()}"\n'
                    + source
                    + "\n"
                    + "\n".join(f'#include "{p}"' for p in tail)
                    + "\n"
                )
                command[index] = str(wrapper)
        return subprocess.run(command, **kwargs)


"""Reusable native compilation, execution, and SPI fixture support."""


CC = os.environ.get("CC", "clang")


def openmp_prefix():
    """Configured static runtime, otherwise the host's Homebrew installation."""
    prefix = os.environ.get("CTOOLS_LIBOMP_PREFIX") or os.environ.get("LIBOMP_PREFIX")
    if not prefix:
        prefix = subprocess.check_output(["brew", "--prefix", "libomp"], text=True).strip()
    return Path(prefix)


def run_test(directory, name, source, modules=(), openmp=False):
    cfile = directory / f"{name}.c"
    binary = directory / name
    cfile.write_text(source)
    flags = [
        "-std=c11",
        "-O2",
        "-DSD_FASTMODE",
        "-fno-fast-math",
        "-ffp-contract=off",
        "-fsanitize=undefined",
        "-fno-omit-frame-pointer",
        "-pthread",
    ]
    for path in [ROOT / "src", ROOT / "src/cimport"]:
        flags += ["-I", str(path)]
    if platform.system() == "Darwin":
        flags += [
            "-DSYSTEM=APPLEMAC",
            "-D_DARWIN_C_SOURCE",
            "-Wl,-dead_strip",
            "-Wl,-undefined,dynamic_lookup",
        ]
        if openmp:
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
        ]
        if openmp:
            flags += ["-fopenmp"]
    libraries = [] if platform.system() == "Darwin" else ["-lm"]
    compile_native(
        [
            CC,
            *flags,
            str(cfile),
            *[str(ROOT / "src" / m) for m in modules],
            *libraries,
            "-o",
            str(binary),
        ],
        check=True,
    )
    subprocess.run([str(binary)], check=True, timeout=60)
    print(f"PASS {name}", flush=True)


def harness(name):
    """Include a named C case; the C preprocessor preserves isolated translation units."""
    name = name.removesuffix(".c")
    shared = {
        "p1_native_common": "sort_alloc.h",
        "sep24_native_alloc": "fault_alloc.h",
        "sort_keys_native": "sort_reference.h",
    }
    if name in shared:
        return f'#include "{ROOT / "validation/native" / shared[name]}"\n'
    matches = [
        p
        for p in (ROOT / "validation/native").glob("*.suite.c")
        if f"#if defined(CTOOLS_CASE_{name})" in p.read_text()
    ]
    if len(matches) != 1:
        raise ValueError(f"Native case {name}: expected one source, found {matches}")
    return f'#define CTOOLS_CASE_{name}\n#include "{matches[0]}"\n#undef CTOOLS_CASE_{name}\n'


def unit(directory, name):
    path = directory / (name + ".c")
    path.write_text(harness(name))
    return path


def sort_source(algorithm, numeric=None, string=None):
    definitions = f'#define SORT_IMPLEMENTATION "ctools_sort_{algorithm}.c"\n'
    if numeric is not None:
        definitions += f"#define NUMERIC_CALL {numeric}\n#define STRING_CALL {string}\n"
        return harness("sep24_native_alloc") + definitions + harness("sep24_native_7")
    definitions += f"#define SORT_ORDER ctools_sort_{algorithm}_order_only\n"
    definitions += f'#define STRING_VARIANTS {1 if algorithm == "counting" else 2}\n'
    return definitions + harness("sep24_native_sort")


SPI = f'#include "{ROOT / "validation/native/spi_mock.h"}"\n'


def run_transfer(name, *, passive=False):
    """Identical SPI transport build matrix; only adaptive scheduling needs passive waits."""
    source_dir = Path(os.environ.get("CTOOLS_TEST_SOURCE", ROOT / "src"))
    sanitizers = os.environ.get("CTOOLS_SANITIZERS", "undefined")
    with tempfile.TemporaryDirectory(prefix="ctools-transfer-") as temporary:
        directory = Path(temporary)
        source = directory / "test.c"
        source.write_text(SPI + harness(name))
        for omp in (False, True):
            binary = directory / ("test-" + str(omp))
            flags = [
                "-std=c11",
                "-O2",
                "-DSD_FASTMODE",
                "-fno-fast-math",
                "-fsanitize=" + sanitizers,
                "-fno-omit-frame-pointer",
                "-pthread",
                "-I" + str(source_dir),
            ]
            if platform.system() == "Darwin":
                flags += [
                    "-DSYSTEM=APPLEMAC",
                    "-D_DARWIN_C_SOURCE",
                    "-Wl,-dead_strip",
                    "-Wl,-undefined,dynamic_lookup",
                ]
                if omp:
                    prefix = os.environ.get("LIBOMP_PREFIX", "/opt/homebrew/opt/libomp")
                    flags += [
                        "-Xpreprocessor",
                        "-fopenmp",
                        "-I" + prefix + "/include",
                        prefix + "/lib/libomp.a",
                    ]
            else:
                flags += [
                    "-DSYSTEM=OPUNIX",
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
                    *[
                        str(source_dir / f)
                        for f in ("ctools_types.c", "ctools_arena.c", "ctools_threads.c")
                    ],
                    "-o",
                    str(binary),
                ],
                check=True,
            )
            env = dict(os.environ, UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1")
            if passive:
                env.update(KMP_BLOCKTIME="0", OMP_WAIT_POLICY="PASSIVE")
            subprocess.run([str(binary)], check=True, timeout=120, env=env)
            print(
                f"PASS {name} sanitizers={sanitizers} openmp={omp} source={source_dir}", flush=True
            )


def run_import_scan(name):
    """Run the CSV delimiter and bounded-scan fixtures with their common sanitizer recipe."""
    with tempfile.TemporaryDirectory(prefix="ctools-cimport-") as temporary:
        directory = Path(temporary)
        source = unit(directory, name)
        binary = directory / "test"
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-DSD_FASTMODE",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
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
                "-liconv",
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
        modules = [
            "ctools_types.c",
            "ctools_arena.c",
            "ctools_threads.c",
            "cimport/cimport_parse.c",
            "cimport/cimport_mmap.c",
            "cimport/cimport_encoding.c",
        ]
        compile_native(
            [
                os.environ.get("CC", "cc"),
                *flags,
                str(source),
                *[str(ROOT / "src" / m) for m in modules],
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run(
            [str(binary)],
            check=True,
            timeout=60,
            env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"},
        )
        print(f"PASS {name} (ASan/UBSan)")
