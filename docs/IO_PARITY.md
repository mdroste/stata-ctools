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

Both commands carry variable names, storage types, variable and value labels, display formats, and dataset characteristics between Stata and the file. Text longer than 2,045 bytes is stored as strL, and so are shapefile binary headers.

## Delimited text

`cimport delimited` supports:

- `delimiters()` with several characters, a multi-character string, or `collapse`
- explicit variable names, row and column ranges, `favorstrfixed`, and `maxquotedrows()`
- UTF-8, UTF-16, UTF-32, and named code pages
- automatic encoding detection, using the same ICU detectors as Stata
- `locale()` with all 1,016 locale names Stata recognizes, so numbers like `1.234,5` parse correctly

`cexport delimited` handles strL variables.

## Excel

### Import

`cimport excel` supports `sheet()`, `cellrange()`, `firstrow`, column lists (including renames like `identifier=A`), `allstring` with an optional numeric format, and `describe`, which reports on the workbook without touching the data in memory.

Dates and numbers are handled as in `import excel`:

- Both the 1900 and 1904 date systems work.
- Each column is a date or a datetime column. A single cell with a time of day makes the whole column datetime.
- When a column shares one built-in Excel number format, it gets the matching Stata display format. Custom formats and mixed columns get `%10.0g`.
- TRUE/FALSE cells import as 1/0, error cells as missing, and text formula results as strings.

`cimport excel` differs from `import excel` in a few ways, all deliberate:

- **Default sheet.** Without `sheet()`, it reads the first worksheet. Stata reads the first tab, so a workbook that starts with a chart, dialog, or macro sheet gives you an empty dataset.
- **Non-data sheets.** Chart, dialog, and macro sheets can't be imported (error 601), and `describe` lists them by name.
- **Content Stata loses.** It keeps emoji and other characters outside Unicode's Basic Multilingual Plane, cell text longer than 32,766 bytes, cells and rows that don't record their own position, and worksheets stored outside the workbook's usual folder.
- **Files Stata rejects.** It reads XLSX files that Stata rejects with error 603, such as files missing SpreadsheetML content types or with an empty phonetic run.

### Export

`cexport excel` creates a new workbook or, with `sheet(..., modify)` or `sheet(..., replace)`, adds, updates, or replaces one worksheet in an existing workbook. Everything else in the workbook is kept as it was, including other sheets, formulas, styles, and the date system.

- `modify` leaves cells outside the range you write untouched, and `keepcellfmt` keeps the styles of the cells it overwrites.
- If you overwrite part of a shared formula, the cells that remain keep their formula.
- After an XLSX edit, Excel recalculates formulas the next time the workbook is opened.
- XLS output uses the Excel 97–2003 (BIFF8) format. Older BIFF5 workbooks can't be edited.

## dBase

`cexport dbase` writes dBase IV by default and dBase III with `version(3)`. Numeric field widths match `export dbase`, with one exception. Stata cuts off the last digit of very negative integers: bytes from -100 to -127, ints from -10,000, and longs from -1,000,000,000. `cexport` widens those fields by one column so the values are written in full.

## SAS, SPSS, and SAS XPORT

- **SAS:** reads .sas7bdat files with `bcat()` catalogs and `encoding()`, and accepts a variable list and `if`/`in`.
- **SPSS:** reads .sav and .zsav files with `encoding()`, and accepts a variable list and `if`/`in`. User-defined missing values and ranges become Stata missing values.
- **XPORT v5:** supports `member()` for files with several datasets, `novallabels`, and `describe`. Value labels round-trip through a `formats.xpf` file.
- **XPORT v8:** like `import sasxport8`, a nonempty companion file isn't applied. The import stops with error 4 and leaves the companion's contents in memory.

Missing values and dates follow Stata's own handling for each format. For example, XPORT v5 imports keep SAS special missing values and SAS and XPORT v8 imports don't, just as with Stata's commands.

Two cases get different display formats than in Stata, because Stata's own choice is invalid:

- Some long-string SAS and SPSS variables come with formats like `%2104s`, which Stata can't assign. `cimport` stores them as strL with `%9s`.
- For numeric SPSS variables in Z format, Stata assigns the string format `%8s`. `cimport` uses `%10.0g`.

## Shapefiles

`cimport shp` produces the same layout as `import shp`: the binary `rec_header`, separator rows, `shape_order`, coordinates, sort order, and dataset characteristics. As in Stata, measure (M) values for shape types 28 and 31 are dropped on import, and exporting those shapes fails with error 3300.

## When an import or export fails

As in Stata, SAS, SPSS, XPORT, dBase, and shapefile imports clear the data in memory before reading the file, so a failed import leaves an empty dataset unless you specify `snapshot`. Excel imports check the file before replacing anything. See [Compatibility](COMPATIBILITY.md#when-a-command-fails) for `snapshot` and `fast`.

Exports write to a temporary file next to the destination and replace it only after the write succeeds, so a failed export never leaves a half-written file. Formats that write more than one file, such as shapefiles and XPORT files with value labels, replace each file in turn. If a later file fails, the earlier ones have already been replaced.

## Testing

Every format is tested against Stata's own `import` and `export`. The tests compare values, variable names, storage types, display formats, labels, characteristics, sort order, and the data left in memory after an error. The suites are listed in [`validation/suites.json`](../validation/suites.json).

## Third-party code

SAS and SPSS support builds on [ReadStat](https://github.com/WizardMac/ReadStat) 1.1.9, and XLS import uses [libxls](https://github.com/libxls/libxls) 1.6.3. Both are compiled into the plugin. See [THIRD_PARTY_NOTICES](../THIRD_PARTY_NOTICES) for licenses.
