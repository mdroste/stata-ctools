# Bundled libraries

ctools compiles the libraries below into its plugin. See
[third-party notices](../../docs/licenses/THIRD_PARTY_NOTICES) for versions and
licenses, and [build requirements](../../docs/development/PLATFORMS.md) for
system libraries needed on each platform.

## ReadStat and compression

[ReadStat 1.1.9](https://github.com/WizardMac/ReadStat/tree/v1.1.9) (MIT) provides
SAS, XPORT and SPSS support. The repository includes the C parsers generated
by Ragel, so Ragel is not needed to build ctools.

ReadStat uses the system `iconv` functions on Unix and a Windows code-page
adapter on Windows. `readstat_iconv.h` selects the implementation and prevents
repeated header definitions. SAS string conversion preserves spaces before a
NUL terminator, matching Stata. Records padded with spaces use ReadStat's
original trimming behavior.

Two changes in `sas/readstat_xport_write.c` match Stata's XPORT writer. A
format that is a bare name, with no width, is copied as the format name, since
the generated format parser rejects two-character names such as `YN`. A
version 5 member name is checked only for its eight-byte length, because Stata
names the member after the file name as written.

`zlib/zlib.h` maps the zlib functions used by ReadStat's SPSS `.zsav` reader to
miniz, the MIT-licensed library in `miniz/` also used for XLSX files. The header
provides this translation; ctools does not include the zlib library itself.

[libdeflate 1.25](https://github.com/ebiggers/libdeflate/tree/v1.25) (MIT)
provides decompression and CRC-32 checksums. XLSX export uses its fast CRC-32
implementation, choosing CPU-specific instructions where available and the
portable slice-by-eight algorithm otherwise. The CRC sources and tables are
unchanged from the v1.25 archive, whose SHA-256 checksum is
`d11473c1ad4c57d874695e8026865e38b47116bbcb872bfc622ec8f37a86017d`.

## Excel files

[libxls 1.6.3](https://github.com/libxls/libxls/tree/v1.6.3) (BSD 2 clause)
provides five C source files and their public headers. It uses `iconv` on Unix
and the same Windows adapter as ReadStat. Integer type names have `XLS` prefixes
to avoid collisions with `windows.h`. The locale and endian headers have
`xls_` prefixes so the compiler can distinguish them from system headers.

One local change, marked `ctools:` in `src/xls.c`: `xls_preparseWorkSheet` and
`xls_parseWorkSheet` skip a BOF..EOF substream nested in a worksheet (an
embedded chart), so its cached series values are not read as worksheet cells
and its EOF does not end the sheet. Native Stata reads such sheets the same way.

## Character encoding detection

The ICU 67.1 adaptation (Unicode license) includes byte maps, frequency tables
and C implementations of its default encoding detectors. See the
[ICU guide](icu_charset/README.md) for source files, checksums and instructions
for regenerating the tables. Building and using ctools requires neither Java
nor the ICU libraries.

## DuckDB C API header

`duckdb/duckdb.h` is the C API header of [DuckDB 2.0.1](https://github.com/duckdb/duckdb/tree/v2.0.1) (MIT), from the release's `libduckdb-src.zip` (SHA-256 `22627d9f8fb1c49c3374db8b8c4d8e2ffd186720057ba58aad1a7aaaf78c98f4`), unchanged. `csql` uses its types only: it defines `DUCKDB_STATIC_BUILD` so that no function is declared `dllimport`, and resolves the functions from the library it loads at run time. DuckDB itself is not compiled into the plugin.

## Compiler settings

`scripts/build_plugin.py` applies library-specific compiler options only to the
relevant sources. Only ReadStat uses `HAVE_ZLIB` and the `vendor/zlib` include
path; only libxls uses its configuration and public include directory.
Suppressed warnings apply only to bundled libraries. ctools' own code uses
`-Isrc` and explicit module paths, so library headers cannot replace system
headers accidentally. The renamed libxls headers also support the C tests.
