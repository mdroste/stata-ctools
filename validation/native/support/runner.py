"""Reusable native compilation, execution, and SPI fixture support."""

from .compiler import compile_native
import os
from pathlib import Path
import platform
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
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
    return (ROOT / "validation/native/harnesses" / name).read_text()


SPI = (Path(__file__).with_name("spi_mock.h")).read_text()
