# Import/export parity audit

The requested target is every Stata `import`/`export` format with file parsing,
serialization, and I/O executed in C. There is no runtime fallback to native
`import` or `export`. Full feature parity is **not yet established**. The file
format families exposed by the installed dispatcher have C engines. Newer
StataNow Parquet import remains unimplemented. FRED, Haver and HaverDirect are
excluded from the requested scope; the remaining compatibility work is below.

## Format coverage

This inventory follows the installed Stata dispatcher and help files, with
newer documented formats identified separately.

| Native subcommand | cimport | cexport | Verified scope and remaining work |
|---|---|---|---|
| `delimited` | C reader | C writer | Multiple/literal/collapsed delimiters, ranges, explicit names, long strings, quote limits, statistical encoding detection, returned charset/delimiters and Unicode numeric locales; broader encoding/error coverage remains |
| `excel` — XLSX | C ZIP/XML reader | C ZIP/XML writer | Names, types, dates, ranges, column selection, allstring(format), long text, describe and workbook editing; complex workbook features need expansion |
| `excel` — XLS | C BIFF/OLE reader | C BIFF8/OLE writer | Same tested interface, BIFF8 workbook editing, styles and SST continuation; complex workbook features need expansion |
| `dbase` | C III/IV reader | C writer | Numeric/string/date fields, special missing values, names, datafmt, dataset date characteristics; broader encoding/error-option coverage remains |
| `sas` | C SAS7BDAT reader with SAS7BCAT catalogs | C SAS7BDAT writer | Reader compared with native, including selection and catalogs; native SAS export is unavailable in this installation, so writer parity is unverified |
| `parquet` — newer StataNow | Not implemented | — | Not exposed by the installed dispatcher; still required for the broader all-format target |
| `sasxport5` | C XPORT reader | C XPORT writer | Multiple members, member selection, formats.xpf labels, renamed collisions, missing tags; describe results for multiple members; companion-file/error coverage needs expansion |
| `sasxport8` | C XPORT reader | C XPORT writer | Values, types, formats, case, long strings, SAS syntax companion; companion validation and full metadata/error coverage remain |
| `spss` | C SAV/ZSAV reader | C SAV writer | Labels, missing values, dates, encodings, selection, long strings; user-missing definitions/ranges and standard formats; additional mappings need broader fixtures |
| `shp` | C ESRI reader | C ESRI writer | All supported shape codes, multipart/null records, binary headers, coordinate metadata and SHX; larger/degenerate files and error behavior need expansion |
| `fred` | Not implemented | Not a native export subcommand | Excluded by user instruction |
| `haver`, `haverdirect` | Not implemented on Windows; native-equivalent platform rejection on macOS | Not native export subcommands | Excluded by user instruction |

`jdbc`, `odbc`, `infile`, `infix`, `outfile`, and `use` are separate commands,
not subcommands of `import`/`export`. Version aliases and platform behavior also
need testing before an all-format parity claim.

## Implemented file-format behavior

ReadStat 1.1.9 and libxls 1.6.3 are pinned and compiled into the plugin;
ReadStat's .zsav compression uses the vendored miniz through a zlib-compatible
shim. Source provenance and local changes are in
[`src/io/vendor/README.md`](../src/io/vendor/README.md); distribution licenses
are bundled in `ctools-codecs-LICENSE.txt` and `THIRD_PARTY_NOTICES`. Unix uses system iconv; Windows uses
system code pages through a C adapter, without an iconv DLL.

The shared C cache transports variable names, storage widths, labels, formats,
value labels and dataset characteristics. Text longer than 2,045 bytes and SHP
binary headers are transferred in bounded hex chunks and stored through Mata's
strL interface. The ado wrappers create variables, apply metadata, mark samples
and evaluate Stata display expressions where needed. They never invoke native
`import` or `export`.

Statistical readers preserve format-specific missing/date semantics rather
than imposing one common conversion. For example, XPORT5 retains tagged
missing codes while native SAS/XPORT8 imports collapse them; SPSS seconds use
an epoch in 1582 and Stata stores milliseconds. XPORT5's tested `formats.xpf`
companion round-trips value labels and member selection reads the selected
dataset from multi-member libraries.

Excel supports literal sheet names, explicit/implicit column lists, firstrow,
cell ranges, allstring with an optional numeric format, and workbook describe
without replacing data in memory. XLSX and BIFF8 XLS exports can add, modify
or replace worksheets; modify preserves cells outside the written rectangle,
and keepcellfmt retains styles of overwritten cells. XLSX editing preserves
untouched ZIP entries; XLS editing retains unrelated OLE streams. Both retain
untouched formula records/cells, styles and the workbook date system. Date
formats retain native display metadata; allstring dates/times are formatted in
C, including native padding. Daily-date values follow native rounding through
milliseconds. Like native, one date/datetime decision is made per column (any
time-of-day cell makes it a datetime) and every numeric cell converts to that
unit; mixed display formats fall back to %td/%tc. Built-in number formats take
native's display formats when a whole column shares one (2, 7 and 8 %14.2f;
3, 37 and 38 %10.0gc; 4, 39 and 40 %14.2fc; 9 %4.2f; 10 %6.4f; 11 %10.2e; 48
%10.1e; the built-in id decides even if the workbook redefines it); custom
number formats and mixed columns use %10.0g. TRUE/FALSE cells take no part in
that decision and import as missing in date/datetime columns. Allstring
date/time text is the stored Stata value shown through the cell's format,
truncated to whole seconds (serial 44197.1 is 02:23:59, as in native). XLS
text formula results, TRUE/FALSE and error cells import as strings, 1/0 and
missing. XLSX parts are resolved through the package relationships (chart,
dialog and macro sheets are refused by name with r(601) and listed without a
range); `<dimension>` is only a sizing hint, and namespace-prefixed markup,
implicit cell/row positions and either attribute quote are accepted. Shared
strings have no length limit, keep whitespace-only runs and omit phonetic
(`rPh`) text; rich inline strings keep their first run, as native does.
Explicit and implicit selectors permit the native extra blank column at the
right edge of the data area.

Deliberate Excel import differences, each observed natively (StataNow 19.5):
without `sheet()` the first *worksheet* is imported, whereas native imports
the first tab and so returns an empty dataset (rc 0) when a workbook begins
with a chart, dialog or macro sheet; a workbook with no worksheet is an empty
dataset in both. `describe` names chart/dialog/macro sheets (native lists them
unnamed with range A0:A-1). Native data losses that are not reproduced:
characters beyond U+FFFF become `cp & 0xFFFF` (U+1F600 → U+F600), cell text is
cut to 32,766 bytes (possibly mid-UTF-8), cells or rows lacking `r` after
explicitly positioned ones are dropped, duplicated or shifted, and worksheets
reached through `./`/`../` targets or stored outside the workbook folder are
unnamed and empty. Native also rejects with r(603) packages whose parts lack
SpreadsheetML content types and shared strings with an empty `<rPh/>`; the C
reader accepts both.
XLS output uses a BIFF8 shared-string table with CONTINUE records. XLSX edits
remove the calculation-chain part and its declarations, request recalculation,
and move shared-formula anchors to the first surviving cell when overwritten.
Relative A1 references are translated; native leaves whole-row/column ranges
unchanged. Array-formula anchors follow native's observed overwrite behavior.
Excel output is compared by reimported content and metadata; formula nodes,
cached values and calculation-chain declarations receive separate structural
comparisons. ZIP/OLE bytes and timestamps are not compared.

SHP imports reproduce Stata's binary `rec_header`, separator rows, shape_order,
coordinate storage, sorting and dataset characteristics. Native Stata discards
measure arrays for shape types 28 and 31 and then fails their export with 3300;
the C path matches that observed behavior in the regression fixtures.

Statistical, DBF and SHP readers clear data before parsing, following native's
failure behavior once replacement is authorized by `clear` or empty/unchanged
data. Excel validates its container before replacing data. Short and malformed
header failures compare native codes and the complete dataset left in memory.
XPORT8 companion loading follows native's partial-data behavior described below.
Writers stage output in sibling files and publish after a successful write.
Companion files are published individually: the whole group is not an atomic
transaction if a later publication fails.

## Remaining compatibility gaps

Newer StataNow supports [Parquet import](https://www.stata.com/manuals/dimportparquet.pdf),
including column selection, describe, rowrange and favormemory. It is absent
from this installation and has no C engine here. Implementing and validating
that engine remains part of the requested all-format target.

Delimited import implements delimiter sets, literal sequences, collapse,
explicit variable lists, UTF-8/16/32 and named code-page conversion, strL and
favorstrfixed, range-dependent inference and strict maxquotedrows enforcement.
Numeric locale parsing uses C profiles generated from public OpenJDK 17 CLDR
locale facts; no Java executes at runtime. All 1,016 installed locale names are
recognized. Automatic encoding detection adapts all default ICU 67.1
recognizers to C, including statistical language models, confidence ties and
the 8,000-byte sample. Native charset/delimiter results, literal versus shortcut
delimiter names, and empty-file charset validation are compared. The complete
Java charset alias inventory on every platform, malformed inputs, implicit
system locales outside this Mac's US locale, and exact returned/error results
need broader comparisons.
Delimited export supports long strings; display-format and file-state coverage
still needs expansion.

Excel editing is implemented for XLSX and BIFF8 XLS (Excel 1997/2003). BIFF5
editing is not supported; that format predates the documented native XLS scope.
Shared/array formulas intersecting edited rectangles and calculation chains
have differential structural fixtures. More complex formulas and templates,
legacy extended-ASCII workbook locale behavior,
and all error paths need further native comparisons. Accepting a locale option
alone does not establish equivalent Unicode conversion behavior.

SAS/SPSS format mappings, catalogs and catalog encodings need broader fixtures.
SPSS defined numeric/string missing values and ranges are converted to native
Stata missing values. XPORT5 describe reports the final member's N/k/size, and
uses lowercase member names. XPORT5 permits replacement of saved, unchanged
data without clear; other file readers reject data in memory. XPORT8 companion
handling follows the installed native command's behavior, which does not
attach labels and returns 4 after reading a nonempty companion. The C path
leaves the companion in memory with its declared widths, double numerics and
original names, matching native before compression/case conversion. DBF
exports accept native version(3), and default to version IV. DBF numeric fields
reproduce native's %20.0g text and its 3-, 5- and 10-column byte, int and long
fields; such a field widens by one column only when values that native
truncates are exported (byte -100 to -127, int from -10000, long from
-1000000000; native writes "-10", "-1000", "-100000000"). Companion-file
publication is not an atomic group transaction.

Three long-string SAS/ZSAV cases expose invalid native string display formats
(for example, %2104s) that Stata itself refuses to assign through its public
metadata interface. The C reader preserves the data in strL with valid %9s.
The native SPSS Z format assigns %8s to a numeric variable, which the public
metadata interface also rejects; the C reader uses %10.0g. Tests identify these
four exceptions as **KNOWN GAP**, normalize only the invalid reference format,
and do not count them as exact format-parity passes.

FRED/Haver engines are absent and explicitly excluded by the user. No provider
parity claim is made.

## Validation

The statistical/binary suite passes 183 checks: 179 exact comparisons and four
explicit invalid-native-format exceptions. Its state tests compare saved,
modified, empty and zero-observation datasets for every file family, plus
XPORT8 companion errors and the data/metadata left in memory. Forty additional
comparisons cover malformed/truncated headers, empty versus modified caller
data, unsupported charsets, and XPORT5 failure after saved data. SAS/SPSS/XPORT8
and DBF imports also compare native returned N/k. The delimited
option suite passes 227 exact comparisons, including returned results and
automatic regional encodings; the Excel/dBase option suite passes 196. Excel
fixtures cover built-in and custom dates/times, fractional serials, both file formats,
numeric display metadata and allstring output. Fifteen additional structural
comparisons check shared/array formula edits, rectangular overwrites, complete
group removal, sheet replacement/addition and calculation-chain cleanup.
Exclusive describe syntax and malformed-container describe failures are compared.
It covers values (hexadecimal doubles and binary strings), names, storage
types, formats, variable/data labels, label definitions, dataset characteristics,
sorting, workbook descriptions and native dataset state on corrupted-input errors.
Native commands create reference fixtures only. Independent ESRI fixtures and
C-generated compressed/long-string statistical fixtures complement those
references. Unicode values and metadata are compared in SPSS/XPORT8, and long
Unicode strings are round-tripped through both Excel writers and readers.

The four dedicated suites total 622 exact comparisons, 15 additional formula
structure comparisons, and four explicit exceptions. The original exact
CSV/XLSX suite passes 32 cases, including
byte-for-byte CSV export of 20,000 randomized rows and sensitive final-digit rounding cases.
It also covers native quoting (only delimiters and quotes), default date/time
formatting and labeled-value fallbacks, whitespace-only and padded fields,
leading-zero and small-exponent parsing, a multi-chunk parallel import, rows
longer than the export size sample, and `shell` after `clear all`.
The larger component suites pass 338 import checks and 213 export checks.
The C boundary and six existing P2 native tests pass under UBSan; binary/Excel
codec truncation/cache-cleanup tests run under ASan/UBSan. Charset recognition
also exercises short/random buffers and the 8,000-byte sample boundary under
ASan/UBSan. The dedicated XLSX long-text loader passes ASan/UBSan checks for
filtering, read failures, binary
rejection, UTF-16 limits and cleanup.

These results are from Apple Silicon Stata on this Mac. Intel macOS and Windows
cross-compilation/dependency checks do not establish runtime parity on those
platforms. Linux runtime/build validation is outstanding; Docker's daemon is
unavailable in this environment. Release dependency checks require macOS 11.0
or earlier and system libraries only, and Windows system DLLs only.

From the repository root:

```sh
python3 validation/test_io_parity_native.py
python3 validation/test_io_formats_native.py
python3 validation/test_xlsx_strl_native.py
python3 validation/test_cimport_options_native.py
python3 validation/test_xlsx_edit_native.py
python3 validation/test_xls_edit_native.py
python3 validation/run_io_delimited_options.py
python3 validation/run_io_excel_options.py
python3 validation/run_io_parity.py
python3 validation/run_io_formats.py
```

The Stata runners invoke the `stata` shell alias through an interactive login
zsh and check the final return-code marker. Fixtures/logs live in
`temp/io_parity/`, `temp/io_allformats/`, `temp/io_delimited_options/` and
`temp/io_excel_options/`; `--plugin-path` accepts a separately
built plugin. The original larger suites are useful supplementary checks, but
their positional/tolerance comparisons do not independently certify exact parity.
