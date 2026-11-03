# Import and export formats

`cimport` and `cexport` read and write files entirely in C. They never call Stata's `import` or `export`, and they don't need Java. This page lists the supported formats, explains how each one is handled, and notes where ctools behaves differently from Stata.

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

StataNow's `import parquet`, `import fred`, and `import haver` aren't supported. Commands like `odbc`, `jdbc`, `infile`, and `infix` are separate from `import` and have no ctools versions.

Metadata preservation depends on the file format and options. Statistical formats can carry variable and value labels and display formats; delimited text cannot preserve a Stata dataset's full metadata. Text longer than 2,045 bytes is stored as strL on import, and so are shapefile binary headers.

## Delimited text

`cimport delimited` supports:

- `delimiters()` with several characters, a multi-character string, or `collapse`
- explicit variable names, row and column ranges, `favorstrfixed`, and `maxquotedrows()`
- UTF-8, UTF-16, UTF-32, and named code pages
- automatic encoding detection, using the same ICU detectors as Stata
- `parselocale()` with 1,016 names from the pinned OpenJDK locale tables, so numbers like `1.234,5` parse correctly

Variable names, quoting, `NaN` and `Infinity` values, `d` exponents, overlapping `numericcols()` and `stringcols()`, and the variables that `rowrange()` keeps follow `import delimited`. `cimport` differs in a few rare cases:

- **Display formats.** `import delimited` chooses each numeric variable's display format from the type its first 513 lines and a few of its last lines need. `cimport` uses the first 513 lines only, so a `long`-range value only near the end of a file displays as `%8.0g` where Stata shows `%12.0g`.
- **Unclosed quotes.** Under `bindquotes(strict)`, a quote still open at the end of the file keeps the file's last newline in Stata's value; `cimport` drops it.
- **Forced numeric text.** With `numericcols()`, text that is not a number becomes missing as in Stata, but the storage type can differ.
- **Very wide rows.** A row with more than 32,767 fields fails with error 902; `import delimited` reads it when `maxvar` allows.

A column with non-integer values is stored as `float`, with or without `asfloat`, unless one of its integer values (every value beyond 2^53 is one) lies outside the `long` range; then it is `double`, as with `import delimited`. `cimport` reads every row for this choice. In files of more than about 262,144 rows, `import delimited` can miss such a value in a later part of the file and store the column as `float`.

`cexport delimited` handles strL variables. With `datafmt`, a comma format such as `%12.0gc` is written without thousands separators and a decimal-comma format such as `%9,2f` with a decimal point, as `export delimited` writes them.

## Excel

### Import

`cimport excel` supports `sheet()`, `cellrange()`, `firstrow`, column lists (including renames like `identifier=A`), `allstring` with an optional numeric format, and `describe`, which reports on the workbook without touching the data in memory.

Dates and numbers are handled as in `import excel`:

- Both the 1900 and 1904 date systems work.
- Each column is a date or a datetime column. A single cell with a time of day makes the whole column datetime.
- Recognized built-in and custom date/time formats receive corresponding Stata display formats. Other custom formats are not fully mapped; mixed columns follow the inferred column type.
- TRUE/FALSE cells import as 1/0, error cells as missing, and text formula results as strings.

`cimport excel` differs from `import excel` in a few ways, all deliberate:

- **Default sheet.** Without `sheet()`, it reads the first worksheet. Stata reads the first tab, so a workbook that starts with a chart, dialog, or macro sheet gives you an empty dataset.
- **Non-data sheets.** Chart, dialog, and macro sheets can't be imported (error 601), and `describe` lists them by name.
- **Content Stata loses.** It keeps emoji and other characters outside Unicode's Basic Multilingual Plane, cell text longer than 32,766 bytes, cells and rows that don't record their own position, and worksheets stored outside the workbook's usual folder.
- **Files Stata rejects.** It reads XLSX files that Stata rejects with error 603, such as files missing SpreadsheetML content types or with an empty phonetic run.
- **Worksheets too large for memory.** An XLSX worksheet whose used range would need more than about half of physical memory to import (at roughly 26 bytes a cell, counted from column A) fails at once with error 909, and the data in memory is kept. A single value at `XFD1048576` is such a worksheet. More columns than `maxvar` fail with error 900, as in Stata, also without replacing the data.
- **Package parts too large for memory.** A part of an XLSX file that would decompress to more than half of physical memory, or to more than `c(max_memory)` when that is set, fails with error 909 before anything is decompressed. A part declaring more data than its compressed size can hold fails with error 603, as in Stata.
- **Time formats and header cells.** A datetime column whose time-formatted cells (such as `h:mm:ss` or `mm:ss.0`) mix with General numbers displays as `%tc`; Stata keeps `%tcHH:MM:SS`. A date-formatted header cell labels its variable with the cell's serial number rather than the formatted date. With `allstring`, a date beyond year 9999 reads as the formatted date rather than as Stata's numeric text.

### Export

`cexport excel` creates a new workbook or, with `sheet(..., modify)` or `sheet(..., replace)`, adds, updates, or replaces one worksheet in an existing workbook. Everything else in the workbook is kept as it was, including other sheets, formulas, styles, and the date system.

- `modify` leaves cells outside the range you write untouched, and `keepcellfmt` keeps the styles of the cells it overwrites.
- If you overwrite part of a shared formula, the cells that remain keep their formula.
- After an XLSX edit, Excel recalculates formulas the next time the workbook is opened.
- XLS output uses the Excel 97–2003 (BIFF8) format. Older BIFF5 workbooks can't be edited.
- Charts and objects embedded in an edited worksheet are kept. A chart sheet can't be the target of `sheet()` (error 601).
- `sheet()`, `sheetmodify` and `sheetreplace` follow `export excel`: the separate options take precedence over a mode in `sheet()`, the first of `modify` and `replace` in `sheet()` applies, and without `sheet()` a mode applies to the first worksheet.

Some Excel output differs from `export excel`:

- Dates: `cexport` writes `%td`, `%-td` and other period dates as Excel dates, `%tb` dates through their business calendar, and `%tg` values as plain numbers. `export excel` writes a `%-td` value (and an unconverted `%tC` value in .xls) as the raw Stata number, `%tb` as empty text, and `%tg` as a datetime.
- Number formats: `cexport` writes ordinary numbers in Excel's General format. `export excel` assigns one of five built-in formats from the display format (`0.00` for `%9.2f`, `#,##0.00` for `%9.2fc`, `#,##0` for any `%gc`, `0.00E+00` for `%9.2e`, and `0` for every other numeric format, which shows `%9.0g` values such as 1234.57 as 1235). Values are identical.
- Long text: an Excel cell holds at most 32,767 characters. `cexport excel` refuses longer strL text with error 109; `export excel` truncates it to 32,766 characters.

## dBase

`cexport dbase` writes dBase IV by default and dBase III with `version(3)`. As with `export dbase`, the dataset and the exported variables are limited to 255 variables and the exported variables to an observation width (`c(width)`) of 4,000 bytes; the record written can be longer (4,001 bytes for 200 doubles). Numeric field widths match `export dbase`, with one exception. Stata cuts off the last digit of very negative integers: bytes from -100 to -127, ints from -10,000, and longs from -1,000,000,000. `cexport` widens those fields by one column so the values are written in full.

## SAS, SPSS, and SAS XPORT

- **SAS:** reads .sas7bdat files with `bcat()` catalogs and `encoding()`, and accepts a variable list and `if`/`in`.
- **SPSS:** reads .sav and .zsav files with `encoding()`, and accepts a variable list and `if`/`in`. User-defined missing values and ranges become Stata missing values.
- **XPORT v5:** supports `member()` for files with several datasets, `novallabels`, and `describe`. Value labels round-trip through a `formats.xpf` file.
- **XPORT v8:** like `import sasxport8`, a nonempty companion file isn't applied. The import stops with error 4 and leaves the companion's contents in memory.

Missing values and dates follow Stata's own handling for each format. For example, XPORT v5 imports keep SAS special missing values and SAS and XPORT v8 imports don't, just as with Stata's commands. SPSS and XPORT v8 exports write `%tC` values without leap seconds, and monthly, quarterly, weekly and half-yearly dates as the first day of the period with a matching date format, as `export spss` and `export sasxport8` do. The conversion uses one temporary variable for each such variable, so an export of a dataset at `c(maxvar)` variables fails with error 900 where Stata's commands succeed.

As with Stata's commands, SPSS numeric formats are the display width plus three, COMMA for comma formats and E for `%e`. An XPORT v5 member is named after the file name, upper case, up to its first period and eight bytes, as written; XPORT v8 members are named `DATASET`. XPORT v5 value-label names lose trailing digits (`lab1` is written as `LAB`).

These cases get different display formats than in Stata, because Stata's own choice is invalid:

- Some long-string SAS and SPSS variables come with formats like `%2104s`, which Stata can't assign. `cimport` stores them as strL with `%9s`.
- For numeric SPSS variables in Z format, Stata assigns the string format `%8s`. `cimport` uses `%10.0g`.
- For business-calendar (`%tb`) variables, `export spss` writes a string format on a numeric variable, so its file reads back as garbled strings, and `export sasxport8` writes a `$8` format. `cexport` writes the business-day numbers as plain numbers.

## Shapefiles

`cimport shp` produces the same layout as `import shp`: the binary `rec_header`, separator rows, `shape_order`, coordinates, sort order, and dataset characteristics. As in Stata, measure (M) values for shape types 28 and 31 are dropped on import, and exporting those shapes fails with error 3300.

## When an import or export fails

As in Stata, SAS, SPSS, XPORT, dBase, and shapefile imports clear the data in memory before reading the file, so a failed import leaves an empty dataset unless you specify `snapshot`. Excel imports check the file before replacing anything. See [Compatibility](COMPATIBILITY.md#when-a-command-fails) for `snapshot` and `fast`.

Exports write to temporary files next to their destinations and publish them after writing succeeds. As with Stata's commands, `replace` turns a destination that is a symbolic link into a regular file and leaves the link's target unchanged, while `sheet(..., modify)` writes through the link to its target. A write-protected destination is not replaced, with the error code Stata's command gives (603 for delimited and Excel, 602 for SPSS, XPORT v8, dBase and shapefiles, 608 for XPORT v5); unlike `export spss` and `export sasxport8`, `cexport` leaves that file in place. Multi-file exports, such as shapefiles and XPORT files with value labels, retain the originals and roll back ordinary publication failures. This is not a crash-atomic transaction across files. If rollback itself fails, the original backups are retained and their paths reported. On POSIX filesystems without exclusive rename or hard links, the no-overwrite fallback copies into an exclusively created destination; that name is visible during the copy and is removed if copying fails.

## Testing

Tests compare with Stata's `import` and `export` where the reference command is available, checking values, variable names, storage types, display formats, labels, characteristics, sort order, and the data left in memory after an error. SAS export uses reader-based validation because the validation installation lacks the native SAS writer. List the discovered cases with `python3 validation/run.py --list`; see the [testing guide](../validation/README.md).

## Third-party code

SAS and SPSS support builds on [ReadStat](https://github.com/WizardMac/ReadStat) 1.1.9, and XLS import uses [libxls](https://github.com/libxls/libxls) 1.6.3. Both are compiled into the plugin. See [THIRD_PARTY_NOTICES](licenses/THIRD_PARTY_NOTICES) for licenses.
