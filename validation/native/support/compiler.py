"""Compile native harnesses with the production modules they exercise.

White-box harnesses deliberately include implementation files to interpose SPI
and allocation failures. Keep extracted companion modules in that same test
translation unit so every allocation remains covered by fault injection.
Ordinary plugin builds compile these modules separately.
"""

from pathlib import Path
import os
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
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
    for name in re.findall(r'^\s*#include\s+"([^"\n]+\.c)"', path.read_text(), re.M):
        for folder in [path.parent, *includes]:
            candidate = folder / name
            if candidate.is_file():
                included_sources(candidate, includes, seen)
                break
    return seen


def compile_native(command, **kwargs):
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
