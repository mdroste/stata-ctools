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

## Character encoding detection

The ICU 67.1 adaptation (Unicode license) includes byte maps, frequency tables
and C implementations of its default encoding detectors. See the
[ICU guide](icu_charset/README.md) for source files, checksums and instructions
for regenerating the tables. Building and using ctools requires neither Java
nor the ICU libraries.

## Compiler settings

`scripts/build_plugin.py` applies library-specific compiler options only to the
relevant sources. Only ReadStat uses `HAVE_ZLIB` and the `vendor/zlib` include
path; only libxls uses its configuration and public include directory.
Suppressed warnings apply only to bundled libraries. ctools' own code uses
`-Isrc` and explicit module paths, so library headers cannot replace system
headers accidentally. The renamed libxls headers also support the C tests.
