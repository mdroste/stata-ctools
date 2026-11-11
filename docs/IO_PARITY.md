# Import and export formats

`cimport` and `cexport` read and write files entirely in C, and they don't need Java. The one exception is the online sources FRED and Haver, which `cimport` hands to Stata's own `import` commands. This page lists the supported formats, explains how each one is handled, and notes where ctools behaves differently from Stata.

## Supported formats

| Format | `cimport` | `cexport` |
|---|---|---|
| Delimited text (`delimited`) | Yes | Yes |
| Excel XLSX (`excel`) | Yes | Yes, including edits to existing workbooks |
| Excel XLS, Excel 97–2003 (`excel`) | Yes | Yes, including edits to existing workbooks |
| dBase (`dbase`) | dBase III and IV | dBase IV, or III with `version(3)` |
| SAS (`sas`) | .sas7bdat, with .sas7bcat catalogs | .sas7bdat |
| SAS XPORT v5 (`sasxport5`) | Yes | Yes |
| SAS XPORT v8 (`sasxport8`) | Yes | Yes |
| SPSS (`spss`) | .sav and .zsav | .sav |
| ESRI shapefiles (`shp`) | Yes | Yes |
| Apache Parquet (`parquet`) | Yes | No |
| FRED (`fred`), Haver (`haver`, `haverdirect`) | Through Stata's `import` | No |

Commands like `odbc`, `jdbc`, `infile`, and `infix` are separate from `import` and have no ctools versions. Stata has no `export parquet`.

Metadata preservation depends on the file format and options. Statistical formats can carry variable and value labels and display formats; delimited text cannot preserve a Stata dataset's full metadata. Text longer than 2,045 bytes is stored as strL on import, and so are shapefile binary headers.

## Delimited text

`cimport delimited` supports:

- `delimiters()` with several characters, a multi-character string, or `collapse`
- explicit variable names, row and column ranges, `favorstrfixed`, and `maxquotedrows()`
- UTF-8, UTF-16, UTF-32, and every charset name and alias that Java gives `import delimited` (single-byte charsets are decoded from Java's tables, so each byte reads as natively; on Windows, Big5-HKSCS, x-EUC-TW, x-ISCII91, ISO-2022-JP-2 and the IBM EBCDIC CJK pages give r(198))
- automatic encoding detection, using the same ICU detectors as Stata
- `parselocale()` with 1,016 locale names and their conventions for digits, decimal and grouping separators, and signs

Variable names, quoting, `NaN` and `Infinity` values, `d` exponents, overlapping `numericcols()` and `stringcols()`, and the variables that `rowrange()` keeps follow `import delimited`, as does a quote still open at the end of the file under `bindquotes(strict)` (its line ends are kept). `cimport` differs in the following case:

- **Forced numeric text.** With `numericcols()`, text that is not a number becomes missing as in Stata. A column of text only gets Stata's type, which follows the column's string width. When text and numbers mix, `import delimited` can store the same file as different types from one run to the next; `cimport` stores the type the numbers need.

Display formats follow `import delimited`: each numeric variable takes the format of the type that a sample of the file's lines needs (the first 513 lines, then the groups of 512 lines that start at lines 1,026, 2,050, 4,098 and so on, and a last group that is not full).

A column with non-integer values is stored as `float`, with or without `asfloat`, unless one of its integer values (every value beyond 2^53 is one) lies outside the `long` range; then it is `double`. `cimport` reads every row for this choice. `import delimited` starts from the type of its sample and widens it line by line, so its choice depends on where values lie: a column that is `float` in the sample stays `float` even when a later value is outside the `long` range, and a column that is `long` in the sample becomes `double` at a later non-integer value.

`cexport delimited` handles strL variables. With `datafmt`, a comma format such as `%12.0gc` is written without thousands separators and a decimal-comma format such as `%9,2f` with a decimal point, as `export delimited` writes them.

## Excel

### Import

`cimport excel` supports `sheet()`, `cellrange()`, `firstrow`, column lists (including renames like `identifier=A`), `allstring` with an optional numeric format, and `describe`, which reports on the workbook without touching the data in memory.

As with `import excel`, the workbook kind comes from the file's content, so `.xlsm`, `.xltx` or any other name reads. `cexport excel` keeps any file name, adding `.xlsx` only when there is no extension; a name ending in `.xls` gets a binary workbook. Native abbreviations work: `sh()`, `all`, `date()`, `miss()`; delimited `decimals()`, `groups()`, `q`, `novar`; XPORT5 `val()`.

Dates and numbers are handled as in `import excel`:

- Both the 1900 and 1904 date systems work.
- Each column is a date or a datetime column. A single cell with a time of day makes the whole column datetime.
- Built-in date/time formats receive their Stata display formats; a custom format displays as `%td` (a date), `%tc` (a date with a time, or an elapsed time such as `[h]:mm`) or `%tcHH:MM:SS` (a time of day). Columns whose cells mix formats use `%td` or `%tc`, or `%tcHH:MM:SS` when every date or time cell is a time. Cells without a style take the workbook's first cell style, and XML parts stored as UTF-16 are read.
- In text columns and header labels, date and time cells read as their formatted text; a value Stata cannot show as a date (beyond year 9999) reads as Stata shows it, in the format's width.
- TRUE/FALSE cells import as 1/0, error cells as missing, and text formula results as strings.

`cimport excel` differs from `import excel` as follows:

- **Default sheet.** Without `sheet()`, it reads the first worksheet. Stata reads the first tab, so a workbook that starts with a chart, dialog, or macro sheet gives you an empty dataset.
- **Non-data sheets.** Chart, dialog, and macro sheets can't be imported (error 601), and `describe` lists them by name.
- **Additional supported text and worksheet layouts.** It keeps emoji and other characters outside Unicode's Basic Multilingual Plane, cell text longer than 32,766 bytes, cells and rows that don't record their own position, and worksheets stored outside the workbook's usual folder.
- **Additional supported XLSX files.** It reads XLSX files that Stata rejects with error 603, such as files missing SpreadsheetML content types or with an empty phonetic run.
- **Worksheets too large for memory.** An XLSX worksheet whose used range would need more than about half of physical memory to import (at roughly 26 bytes a cell, counted from column A) fails at once with error 909, and the data in memory is kept. A single value at `XFD1048576` is such a worksheet. More columns than `maxvar` fail with error 900, as in Stata, also without replacing the data.
- **Package parts too large for memory.** A part of an XLSX file that would decompress to more than half of physical memory, or to more than `c(max_memory)` when that is set, fails with error 909 before anything is decompressed. A part declaring more data than its compressed size can hold fails with error 603, as in Stata.

### Export

`cexport excel` creates a new workbook or, with `sheet(..., modify)` or `sheet(..., replace)`, adds, updates, or replaces one worksheet in an existing workbook. Everything else in the workbook is kept as it was, including other sheets, formulas, styles, and the date system.

- `modify` leaves cells outside the range you write untouched, and `keepcellfmt` keeps the styles of the cells it overwrites, except that `%t` variables still receive their formats.
- If you overwrite part of a shared formula, the cells that remain keep their formula.
- After an XLSX edit, Excel recalculates formulas the next time the workbook is opened.
- XLS output uses the Excel 97–2003 (BIFF8) format. Older BIFF5 workbooks can't be edited.
- Charts and objects embedded in an edited worksheet are kept. A chart sheet can't be the target of `sheet()` (error 601).
- `sheet()`, `sheetmodify` and `sheetreplace` follow `export excel`: the separate options take precedence over a mode in `sheet()`, the first of `modify` and `replace` in `sheet()` applies, and without `sheet()` a mode applies to the first worksheet. A worksheet name Excel refuses (over 31 characters, containing `: \ / ? * [ ]`, or beginning or ending with an apostrophe) is error 9901, and with `sheet()` an existing file that is not a workbook is error 603, as in `export excel`.
- Cell formats follow `export excel`: `m/d/yy` for dates, `m/d/yy h:mm` for datetimes, `h:mm:ss AM/PM` for `%tc` formats that begin with an hour, minute or am/pm code; `0.00` for `%9.2f`, `#,##0.00` for `%9.2fc`, `#,##0` for any `%gc`, `0.00E+00` for `%9.2e`, and `0` for every other numeric format (so `%9.0g` shows 1234.57 as 1235). Blank missing cells keep the column's format; under a numeric `missing()` and under `keepcellfmt`, only `%t` variables receive formats.

Some Excel output differs from `export excel`:

- Dates: `cexport` writes `%td`, `%-td` and other period dates as Excel dates, `%tb` dates through their business calendar, and `%tg` values as plain numbers. `export excel` writes a `%-td` value (and an unconverted `%tC` value in .xls) as the raw Stata number, `%tb` as empty text, and `%tg` as a datetime.
- Fractional seconds: for a `%tc` format that shows seconds with `s` (such as `%tcHH:MM:SS.sss`), `export excel` writes empty text. `cexport` writes the datetime, formatted `yyyy-mm-dd hh:mm:ss.000` in .xlsx and `m/d/yy h:mm` in .xls. A `%-td` or `%-tc` value, which `cexport` writes as a date, has a date format; `export excel` writes the raw number formatted `0`.
- Long text: an Excel cell holds at most 32,767 characters. `cexport excel` refuses longer strL text with error 109; `export excel` truncates it to 32,766 characters.

## Parquet

`cimport parquet` is a C reader for Apache Parquet files with the syntax of StataNow 19.5's `import parquet`: a column list, `describe` (with `r(N_rows)` and `r(n_columns)` and the same listing), `clear`, `rowrange()`, and `favormemory`, which is accepted and changes nothing. It doesn't need StataNow, Java, or the Arrow library. Stata's command is built on Arrow, so `cimport` reproduces how Arrow presents each Parquet type and how Stata then stores it. Its tests check that both commands create the same variables, types, display formats, variable and value labels, and values, and print the same messages, for every fixture listed below.

The format support:

- **Footer and schema.** The Thrift compact-protocol footer, all physical types, logical types and the older converted types. Nested columns (groups, lists, maps, repeated fields) are listed by `describe` as `struct`, `list` or `map` and skipped with a note, as in Stata.
- **Pages.** All row groups, data pages version 1 and 2 (compressed or not), dictionary pages, and definition levels in the RLE/bit-packed hybrid or the old BIT_PACKED encoding. Arrow reads BIT_PACKED levels least significant bit first, though the format specifies the opposite order; `cimport` reads them as Arrow does.
- **Encodings.** PLAIN, PLAIN_DICTIONARY, RLE_DICTIONARY, RLE (booleans), DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, and BYTE_STREAM_SPLIT.
- **Compression.** Snappy and LZ4 (the Hadoop framing of the LZ4 codec, a bare block as Arrow accepts it, and LZ4_RAW) use decoders written for ctools; Gzip uses the bundled libdeflate; Zstandard and Brotli use decoders written for ctools from RFC 8878 and RFC 7932 (the Brotli static dictionary is embedded). A damaged Brotli page fails with error 691 and Arrow's message, leaving the data in memory unchanged. Neither command reads LZO; both fail with error 691.
- **Arrow metadata.** Files written by Arrow-based tools (pandas, pyarrow, polars) store their Arrow schema in the footer. As in Stata, a dictionary (categorical) column over text becomes a labeled numeric variable, large string and binary types are listed as `large_utf8` and `large_binary`, and durations are skipped as unsupported. A stored schema that can't be read fails with error 3900, as in Stata.

How values are stored follows `import parquet`, including its unusual choices:

- Integer types start from `byte` (8-bit and boolean), `int` (16-bit) or `long` (32-bit) and widen by one step when the largest value *in the row range* exceeds the type's maximum. A value below the type's minimum, such as -128 in an int8 column, is stored as missing. A 64-bit integer column is `long`, unless a value in the row range doesn't fit in a `long`; then every value is stored as text, the value's 64 bits in lowercase hexadecimal (`str16`).
- Float columns are `float`, or `double` when a value exceeds `c(maxfloat)`. NaN and infinities are missing. Doubles beyond Stata's largest value are stored unchanged.
- Dates are `long` with `%td`. Times of day are `double` milliseconds with no date format. Timestamps, including INT96, are `double` milliseconds with `%tc`; finer units are truncated toward zero, and the UTC flag is ignored.
- Decimals are text, as wide as the precision plus two, in Arrow's notation (`0.05`, `1E-10`).
- Strings are as wide as the longest value *in the whole column*, whatever the row range; past 2,045 bytes they are strL. Binary, fixed-length and enum columns are always strL.
- Dictionary columns are `byte`, `int` or `long` by the size of the column's dictionary (the dictionaries of all row groups, unified in order, including unused entries). Their value labels, named after the variable, modify any label of that name already in memory, and the display format is widened to the longest entry.
- Variable names keep letters, digits and underscores, drop leading digits, and are cut to 32 characters; an empty name becomes `v#`, and if any two names would coincide, every variable is named `v0`, `v1`, .... A column whose name is a Stata reserved word, such as `in` or `byte`, fails with error 3598, as in Stata.
- An empty file (no rows) fails with error 3900, as does a `rowrange()` whose end is not an integer.

Only the selected columns are decoded, and only the row groups the row range needs (string and dictionary columns are measured over every row group, as Stata measures them). Columns are decoded in parallel. On a 2-million-row, 5-column test file, `cimport parquet` took about a third of the time `import parquet` took.

Option `fast` skips the copy that restores the data in memory after an error, and `snapshot` makes it with `preserve`, as for the other formats. `import parquet` itself always restores.

## FRED and Haver

`cimport fred`, `cimport haver`, and `cimport haverdirect` (with the abbreviations `hav`, `have`, `haverd`, and so on) pass every argument unchanged to Stata's `import fred`, `import haver`, and `import haverdirect`. FRED is a web API that needs an API key, and Haver needs Haver Analytics' DLX library, which exists only on Windows; neither is a file format that a C reader could replace. Results, messages, stored results and the data left after an error are those of Stata's commands. Options `snapshot` and `fast` don't apply. Because the tests may not use the network, they only compare the errors both commands give before contacting FRED or opening a Haver database: missing series, conflicting options, data in memory, a missing series-list file, a missing database, and unknown options.

## dBase

`cexport dbase` writes dBase IV by default and dBase III with `version(3)`. As with `export dbase`, the dataset and the exported variables are limited to 255 variables and the exported variables to an observation width (`c(width)`) of 4,000 bytes; the record written can be longer (4,001 bytes for 200 doubles). Numeric field widths match `export dbase`, with one exception. Stata cuts off the last digit of very negative integers: bytes from -100 to -127, ints from -10,000, and longs from -1,000,000,000. `cexport` widens those fields by one column so the values are written in full.

## SAS, SPSS, and SAS XPORT

- **SAS:** reads .sas7bdat files with `bcat()` catalogs and `encoding()`, and accepts names and `if`/`in` before `using`, read as `import sas` reads them (exact, case-sensitive names; a hyphen separates names; unknown names ignored; r(111) without a match). The names and an `in` range are applied while reading.
- **SPSS:** reads .sav and .zsav files with `encoding()`, and accepts names and `if`/`in` before `using`, read as `import spss` reads them (the same rules as SAS). The names and an `in` range are applied while reading. User-defined missing values and ranges become Stata missing values.
- **XPORT v5:** supports `member()` for files with several datasets, `novallabels`, and `describe`. Value labels are written to and read from a companion `formats.xpf` file.
- **XPORT v8:** like `import sasxport8`, a nonempty companion file isn't applied. The import stops with error 4 and leaves the companion's contents in memory.

Missing values and dates follow Stata's own handling for each format. For example, XPORT v5 imports keep SAS special missing values and SAS and XPORT v8 imports don't, just as with Stata's commands. SPSS and XPORT v8 exports write `%tC` values without leap seconds, and monthly, quarterly, weekly and half-yearly dates as the first day of the period with a matching date format, as `export spss` and `export sasxport8` do. The conversion uses one temporary variable for each such variable, so an export of a dataset at `c(maxvar)` variables fails with error 900 where Stata's commands succeed.

As with Stata's commands, SPSS numeric formats are the display width plus three, COMMA for comma formats and E for `%e`. An XPORT v5 member is named after the file name, upper case, up to its first period and eight bytes, as written; XPORT v8 members are named `DATASET`. XPORT v5 names are upper case and at most eight bytes, and only value labels that variables use are written. A variable name that collides takes a numeric suffix (`ABCDEFG2`); a value-label name loses trailing digits (`lab1` is written as `LAB`), one that collides takes letters counting with A–J as digits (`a2` and `a3` become `A` and `AA`), and a SAS format name such as `date` becomes `USER`. Without `rename`, any other change is error 110 with the suggested names, variables first, as in `export sasxport5`.

ctools uses different display formats in the following cases:

- Some long-string SAS and SPSS variables come with formats like `%2104s`, which Stata can't assign. `cimport` stores them as strL with `%9s`.
- For numeric SPSS variables in Z format, Stata assigns the string format `%8s`. `cimport` uses `%10.0g`.
- For business-calendar (`%tb`) variables, `export spss` writes a string format on a numeric variable, so its file reads back as garbled strings, and `export sasxport8` writes a `$8` format. `cexport` writes the business-day numbers as plain numbers.

## Shapefiles

`cimport shp` produces the same layout as `import shp`: the binary `rec_header`, separator rows, `shape_order`, coordinates, sort order, and dataset characteristics. As in Stata, measure (M) values for shape types 28 and 31 are dropped on import, and exporting those shapes fails with error 3300.

## When an import or export fails

As in Stata, SAS, SPSS, XPORT, dBase, and shapefile imports clear the data in memory before reading the file, so a failed import leaves an empty dataset unless you specify `snapshot`. Excel imports check the file before replacing anything. See [Compatibility](COMPATIBILITY.md#when-a-command-fails) for `snapshot` and `fast`.

Exports write to temporary files in the destination directory and move them to the requested filenames after writing succeeds. As with Stata's commands, `replace` turns a destination that is a symbolic link into a regular file and leaves the link's target unchanged, while `sheet(..., modify)` writes through the link to its target. A write-protected destination is not replaced, with the error code Stata's command gives (603 for delimited and Excel, 602 for SPSS, XPORT v8, dBase and shapefiles, 608 for XPORT v5); unlike `export spss` and `export sasxport8`, `cexport` leaves that file in place. Exports containing several files, such as shapefiles and XPORT files with value labels, retain the originals until all new files are in place and restore them if an error occurs. A system crash can leave only some files replaced. If restoration fails, ctools retains the backups and reports their locations. On filesystems that cannot rename or link files without replacing an existing destination, ctools creates a new destination and copies into it. That filename is visible during copying and is removed if copying fails.

## Testing

Tests compare with Stata's `import` and `export` where the reference command is available (for Parquet, the `s2_parquet` case and the files written by `validation/support/parquet_fixtures.py`), checking values, variable names, storage types, display formats, labels, characteristics, sort order, and the data left in memory after an error. SAS export is tested by reading the generated files. The Stata installation used for testing does not provide SAS export. List the discovered cases with `python3 validation/run.py --list`; see the [testing guide](../validation/README.md).

## Third-party code

SAS and SPSS support builds on [ReadStat](https://github.com/WizardMac/ReadStat) 1.1.9, and XLS import uses [libxls](https://github.com/libxls/libxls) 1.6.3. Both are compiled into the plugin. The Parquet reader and its Snappy, LZ4 and Zstandard decoders are ctools' own; Gzip uses the bundled libdeflate. See [THIRD_PARTY_NOTICES](licenses/THIRD_PARTY_NOTICES) for licenses.
