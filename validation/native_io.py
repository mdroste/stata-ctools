"""Native checks; execution and compiler setup are shared in validation.native.

Functions retain their independent fixtures, failure sweeps, and runtime variants.
"""

import os
import platform
import struct
import subprocess
import tempfile
import zipfile
from pathlib import Path

from validation import ROOT
from validation.fixtures import codec_fixtures, shapefiles
from validation.fixtures import excel_fixtures as fixtures
from validation.native import compile_native, harness, run_import_scan, run_test, unit


def test_cimport_options_native():
    run_import_scan("cimport_options_native_source")


def test_cimport_stream_native():
    run_import_scan("cimport_stream_native")


def test_io_formats_native():
    SOURCE = harness("io_formats_native_source.c")
    vendor = ROOT / "src/vendor"
    with tempfile.TemporaryDirectory(prefix="ctools-binary-asan-") as tmp:
        directory = Path(tmp)
        shapefiles(directory)
        codec_fixtures(directory)
        fixtures(directory)
        cfile = directory / "test.c"
        cfile.write_text(SOURCE)
        sources = [
            *sorted((vendor / "readstat").rglob("*.c")),
            *sorted((vendor / "zlib").glob("*.c")),
            *sorted((vendor / "libxls/src").glob("*.c")),
            *sorted((ROOT / "src/vendor/miniz").glob("*.c")),
            *sorted((ROOT / "src/vendor/libdeflate").rglob("*.c")),
            ROOT / "src/cexport/cexport_parse.c",
            *[
                ROOT / "src" / name
                for name in [
                    "cimport/cimport_xlsx.c",
                    "formats/xlsx_xml.c",
                    "formats/xlsx_zip.c",
                    "ctools_types.c",
                    "ctools_arena.c",
                    "ctools_threads.c",
                ]
            ],
        ]
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-DHAVE_ZLIB=1",
            "-DSD_FASTMODE",
            "-D_GNU_SOURCE",
            "-D_DARWIN_C_SOURCE",
            "-pthread",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
            "-I",
            str(ROOT / "src"),
            "-I",
            str(vendor / "readstat"),
            "-I",
            str(vendor / "zlib"),
            "-I",
            str(vendor / "libxls/include"),
            "-I",
            str(vendor / "libxls/include/libxls"),
        ]
        if platform.system() == "Darwin":
            flags += ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip", "-liconv"]
        else:
            flags += [
                "-DSYSTEM=STUNIX",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
                "-lm",
            ]
        binary = directory / "test"
        compile_native(
            [os.environ.get("CC", "cc"), *flags, str(cfile), *map(str, sources), "-o", str(binary)],
            check=True,
            provided=("formats/excel_address.c",),
        )
        subprocess.run(
            [str(binary), str(directory)],
            check=True,
            timeout=60,
            env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"},
        )
        print("PASS binary codecs, truncation, and cache cleanup (ASan/UBSan)")


def test_io_parity_native():
    SOURCE = harness("io_parity_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-io-parity-") as directory:
        run_test(Path(directory), "io_parity", SOURCE, modules=("ctools_types.c",))


def test_io_performance_native():
    CRC = harness("io_perf_crc.c")
    XML = harness("io_perf_xml.c")
    CSV = harness("io_perf_csv.c")
    mac = platform.system() == "Darwin"
    vendor = ROOT / "src/vendor"
    with tempfile.TemporaryDirectory(prefix="ctools-io-performance-") as tmp:
        directory = Path(tmp)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-DSD_FASTMODE",
            "-fno-fast-math",
            "-ffp-contract=off",
            "-fno-strict-aliasing",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
            "-pthread",
            "-I",
            str(ROOT / "src"),
        ]
        flags += (
            [
                "-DSYSTEM=APPLEMAC",
                "-D_DARWIN_C_SOURCE",
                "-Wl,-dead_strip",
                "-Wl,-undefined,dynamic_lookup",
            ]
            if mac
            else [
                "-DSYSTEM=STUNIX",
                "-D_GNU_SOURCE",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ]
        )
        crc_sources = [
            vendor / "libdeflate/lib/crc32.c",
            vendor / "miniz/miniz.c",
            vendor / "libdeflate/lib/arm/cpu_features.c",
            vendor / "libdeflate/lib/x86/cpu_features.c",
        ]
        xml_sources = [
            *(
                ROOT / "src" / p
                for p in (
                    "ctools_types.c",
                    "ctools_arena.c",
                    "ctools_threads.c",
                    "formats/xlsx_xml.c",
                    "formats/xlsx_zip.c",
                )
            ),
            *sorted((vendor / "miniz").glob("*.c")),
            *sorted((vendor / "libdeflate").rglob("*.c")),
        ]
        csv_sources = [
            ROOT / "src" / p
            for p in (
                "ctools_types.c",
                "ctools_arena.c",
                "ctools_threads.c",
                "cimport/cimport_parse.c",
                "cimport/cimport_mmap.c",
                "cimport/cimport_encoding.c",
            )
        ]
        for name, source, modules in [
            ("crc", CRC, crc_sources),
            ("xml", XML, xml_sources),
            ("csv", CSV, csv_sources),
        ]:
            cfile = directory / (name + ".c")
            cfile.write_text(source)
            binary = directory / name
            compile_native(
                [
                    os.environ.get("CC", "cc"),
                    *flags,
                    str(cfile),
                    *map(str, modules),
                    *(["-liconv"] if mac else ["-lm"]),
                    "-o",
                    str(binary),
                ],
                check=True,
            )
            subprocess.run(
                [str(binary)], check=True, env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"}
            )
            print("PASS I/O performance " + name + " (ASan/UBSan)", flush=True)


def test_oct21_native():

    def xlsx_fixtures(folder):
        for name, method in [
            ("valid.xlsx", zipfile.ZIP_DEFLATED),
            ("stored.xlsx", zipfile.ZIP_STORED),
        ]:
            with zipfile.ZipFile(
                folder / name, "w", compression=method, compresslevel=0
            ) as archive:
                archive.writestr("sheet.xml", b"<v>123</v>")
        original = (folder / "valid.xlsx").read_bytes()
        data = bytearray(original)
        data[data.index(b"<v>123") + 3] ^= 8
        (folder / "bitflip.xlsx").write_bytes(data)
        for name, field in [("crc.xlsx", 16), ("length.xlsx", 24)]:
            data = bytearray(original)
            offset = data.index(b"PK\x01\x02") + field
            struct.pack_into("<I", data, offset, struct.unpack_from("<I", data, offset)[0] + 1)
            (folder / name).write_bytes(data)

    with tempfile.TemporaryDirectory(prefix="ctools-oct21-") as tmp:
        directory = Path(tmp)
        xlsx_fixtures(directory)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-DSD_FASTMODE",
            "-D_GNU_SOURCE",
            "-D_DARWIN_C_SOURCE",
            "-pthread",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
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
                "-lm",
            ]
        sources = [
            unit(directory, "oct21_xlsx"),
            *[
                ROOT / "src" / p
                for p in [
                    "formats/xlsx_xml.c",
                    "ctools_types.c",
                    "ctools_arena.c",
                    "ctools_threads.c",
                ]
            ],
            *sorted((ROOT / "src/vendor/miniz").glob("*.c")),
            *sorted((ROOT / "src/vendor/libdeflate").rglob("*.c")),
        ]
        binary = directory / "xlsx"
        compile_native(
            [os.environ.get("CC", "cc"), *flags, *map(str, sources), "-o", str(binary)], check=True
        )
        subprocess.run(
            [str(binary), str(directory)],
            check=True,
            env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"},
        )
        output = directory / "output"
        compile_native(
            [
                os.environ.get("CC", "cc"),
                *flags,
                str(unit(directory, "oct21_output")),
                "-o",
                str(output),
            ],
            check=True,
        )
        subprocess.run(
            [str(output), str(directory)],
            check=True,
            env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"},
        )


def test_oct22_native():
    with tempfile.TemporaryDirectory(prefix="ctools-oct22-") as tmp:
        directory = Path(tmp)
        flags = [
            "-std=c11",
            "-O1",
            "-g",
            "-DSD_FASTMODE",
            "-D_GNU_SOURCE",
            "-D_DARWIN_C_SOURCE",
            "-pthread",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
            "-I",
            str(ROOT / "src"),
        ]
        if platform.system() == "Darwin":
            flags += ["-DSYSTEM=APPLEMAC", "-Wl,-dead_strip"]
            libraries = []
        else:
            flags += [
                "-DSYSTEM=STUNIX",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ]
            libraries = ["-lm"]
        for name in ["xml", "export", "match_oom", "range"]:
            sources = [unit(directory, "oct22_" + name)]
            if name == "xml":
                sources += [
                    ROOT / "src" / p
                    for p in [
                        "formats/xlsx_xml.c",
                        "ctools_types.c",
                        "ctools_arena.c",
                        "ctools_threads.c",
                    ]
                ]
                sources += sorted((ROOT / "src/vendor/miniz").glob("*.c")) + sorted(
                    (ROOT / "src/vendor/libdeflate").rglob("*.c")
                )
            if name == "export":
                sources += [ROOT / "src/ctools_hash.c", ROOT / "src/ctools_arena.c"]
            binary = directory / name
            compile_native(
                [
                    os.environ.get("CC", "cc"),
                    *flags,
                    *map(str, sources),
                    *libraries,
                    "-o",
                    str(binary),
                ],
                check=True,
            )
            subprocess.run(
                [str(binary)], check=True, env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"}
            )


def test_xls_edit_native():
    SOURCE = harness("xls_edit_native_source.c")
    vendor = ROOT / "src/vendor"
    with tempfile.TemporaryDirectory(prefix="ctools-xls-edit-") as tmp:
        directory = Path(tmp)
        source = directory / "test.c"
        source.write_text(SOURCE)
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
            "-I",
            str(ROOT / "src"),
            "-I",
            str(vendor / "libxls/include"),
            "-I",
            str(vendor / "libxls/include/libxls"),
            "-D_DARWIN_C_SOURCE",
            "-D_GNU_SOURCE",
        ]
        if platform.system() == "Darwin":
            flags += [
                "-DSYSTEM=APPLEMAC",
                "-Wl,-dead_strip",
                "-Wl,-undefined,dynamic_lookup",
                "-liconv",
            ]
        else:
            flags += [
                "-DSYSTEM=STUNIX",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
                "-lm",
            ]
        compile_native(
            [
                os.environ.get("CC", "cc"),
                *flags,
                str(source),
                *map(str, sorted((vendor / "libxls/src").glob("*.c"))),
                "-o",
                str(binary),
            ],
            check=True,
        )
        subprocess.run(
            [str(binary), str(directory)],
            check=True,
            timeout=60,
            env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"},
        )
        print(
            "PASS XLS long SST continuations, formulas, preserved OLE streams, record bounds (ASan/UBSan)"
        )


def test_xlsx_edit_native():
    SOURCE = harness("xlsx_edit_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-xlsx-edit-") as tmp:
        directory = Path(tmp)
        source = directory / "test.c"
        source.write_text(SOURCE)
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
        compile_native(
            [
                os.environ.get("CC", "cc"),
                *flags,
                str(source),
                str(ROOT / "src/formats/xlsx_xml.c"),
                "-o",
                str(binary),
            ],
            check=True,
            provided=("formats/excel_address.c",),
        )
        subprocess.run(
            [str(binary)],
            check=True,
            timeout=60,
            env={**os.environ, "UBSAN_OPTIONS": "halt_on_error=1"},
        )
        print(
            "PASS XLSX edits preserve styles/formulas/XML, namespaces, truncated tags (ASan/UBSan)"
        )


def test_xlsx_strl_native():
    SOURCE = harness("xlsx_strl_native_source.c")
    with tempfile.TemporaryDirectory(prefix="ctools-xlsx-strl-") as tmp:
        directory = Path(tmp)
        source = directory / "test.c"
        source.write_text(SOURCE)
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
        modules = ["ctools_data_io.c", "ctools_types.c", "ctools_arena.c", "ctools_threads.c"]
        compile_native(
            [
                os.environ.get("CC", "cc"),
                *flags,
                str(source),
                *[str(ROOT / "src" / name) for name in modules],
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
        print(
            "PASS XLSX long text, filtering, SPI failures, binary rejection, UTF-16 limit (ASan/UBSan)"
        )
