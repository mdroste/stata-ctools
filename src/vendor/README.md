# Pinned C codecs

* ReadStat 1.1.9 (MIT), https://github.com/WizardMac/ReadStat/tree/v1.1.9
  Only library sources for SAS, XPORT, and SPSS are included. Generated Ragel C
  parsers are checked in; Ragel is not a build dependency.
* `zlib/zlib.h` is not zlib: it is a small shim that maps the zlib names
  ReadStat's SPSS .zsav code uses onto miniz (src/vendor/miniz, MIT), which
  ctools already compiles for XLSX I/O. No zlib sources are vendored.

Everything is compiled into the plugin; users need no additional dynamic libraries.
ReadStat uses system iconv on Unix and the Windows code-page adapter on Windows.
Local changes to upstream sources are documented here as they are introduced.

* libdeflate 1.25 (MIT), https://github.com/ebiggers/libdeflate/tree/v1.25
  The decompressor and CRC-32 implementation are included. XLSX export uses
  its accelerated CRC-32, including runtime CPU dispatch and the portable
  slice-by-eight fallback. The CRC sources and tables are unchanged upstream
  files, added from the v1.25 archive (SHA-256
  `d11473c1ad4c57d874695e8026865e38b47116bbcb872bfc622ec8f37a86017d`).

ReadStat changes: iconv includes route through `readstat_iconv.h`, which selects
the local code-page adapter on Windows. The header also has an include guard for repeated charset includes.

SAS string conversion preserves spaces before a NUL terminator, matching Stata;
space-padded records retain the upstream trimming behavior.

* libxls 1.6.3 (BSD 2 clause), https://github.com/libxls/libxls/tree/v1.6.3
  The five library C sources and public headers are included. Config selects
  iconv; Windows conversion uses the same system adapter as ReadStat. The
  integer typedefs have XLS prefixes to avoid collisions with windows.h.

libxls locale/endian headers are renamed with xls_ prefixes so the repository
include-directory discovery does not shadow system locale.h or endian.h.

* ICU 67.1 charset recognition (Unicode license). Statistical byte maps and
  frequency tables plus a portable C adaptation of the default recognizers;
  Java and ICU libraries are not runtime or build dependencies. Provenance,
  source checksums and regeneration instructions are in `icu_charset/README.md`.

## Unified ownership and compilation

Compression now lives alongside the codecs: `miniz/` and `libdeflate/` retain
the exact previously vendored sources and licenses. This relocation does not
upgrade upstream packages. See the root `THIRD_PARTY_NOTICES` for their versions
and licenses. ICU input provenance is in `icu_charset/` and `tools/generation/`.

`scripts/build_plugin.py` scopes warning suppression to vendor sources. Only
ReadStat sees `HAVE_ZLIB` and `vendor/zlib`; only libxls sees its config and public
include directory. First-party modules use `-Isrc` and explicit module paths.
The renamed libxls locale/endian headers remain for compatibility with native
test builds, but no recursive global include path can shadow system headers.
