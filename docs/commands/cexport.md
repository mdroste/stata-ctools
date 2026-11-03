# cexport

C writers for delimited text, Excel (.xls/.xlsx), dBase, SAS7BDAT, SAS XPORT5/8,
SPSS SAV, and ESRI shapefiles. See
[import and export formats](../IO_PARITY.md) for what's supported and how it differs from `export`.

## Other file formats

```stata
cexport excel using "output.xls", firstrow(variables) replace
cexport spss id amount using "output.sav", replace
cexport sasxport5 using "output.xpt", rename vallabfile(both) replace
cexport sasxport8 using "output.v8xpt", vallabfile replace
cexport dbase using "output.dbf", datafmt replace
cexport shp using "boundaries.shp", shx replace
cexport sas using "output.sas7bdat", replace
```

Statistical writers accept variable lists and if/in selection. SAS export is
tested by reading the generated files. The Stata installation used for testing
does not provide SAS export. Both XLS and XLSX export support text longer than 2,045 bytes, subject to
Excel's cell-length limit. Both Excel writers support
`sheet(..., modify/replace)`, adding worksheets, and `keepcellfmt`. Modification
retains untouched cells and workbook content; keepcellfmt preserves styles in
the written rectangle. When an XLSX edit overwrites part of a shared formula, the remaining cells
retain their formulas. Excel recalculates formulas when the workbook is next opened.
See [import and export formats](../IO_PARITY.md) for format-specific behavior and validation.

## Delimited text

## Overview

`cexport delimited` uses a C plugin with parallel data loading, long-string support and chunked formatting.

## Syntax

```stata
cexport delimited [varlist] using filename [if] [in] [, options]
```

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `delimiter(char)` | Field delimiter (default: comma). Use `tab` for tab-delimited |
| `replace` | Overwrite existing file |

### Formatting Options
| Option | Description |
|--------|-------------|
| `novarnames` | Do not write variable names as header row |
| `quote` | Quote all string fields |
| `noquoteif` | Never quote strings, even if they contain delimiters |
| `datafmt` | Apply numeric and date/time display formats |
| `nolabel` | Write numeric codes instead of value-label text |
| `datestring(fmt)` | Write date/time variables as strings using a common Stata format |

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information during export |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |

## Examples

```stata
* Export all variables to a CSV file
cexport delimited using output.csv, replace

* Export selected variables
cexport delimited id name value using output.csv, replace

* Export to tab-delimited file
cexport delimited using output.tsv, delimiter(tab) replace

* Export without header row
cexport delimited using output.csv, novarnames replace

* Export with verbose timing output
cexport delimited using output.csv, replace verbose

* Export subset of observations
cexport delimited using subset.csv if year > 2020, replace

* Export a specific range
cexport delimited using sample.csv in 1/1000, replace

* Quote all string fields
cexport delimited using output.csv, quote replace
```

## Performance

`cexport` achieves its speed advantage through:

### Implementation details

- **Vectored I/O**: Uses `pwritev()` on POSIX and batched async `WriteFile` on Windows to write multiple buffers in a single syscall, reducing kernel transitions
- **Async mmap flush**: Uses `MS_ASYNC` instead of `MS_SYNC` so the OS writes dirty pages in the background while formatting continues
- **File pre-allocation**: `posix_fallocate` (POSIX) or `SetFilePointerEx`+`SetEndOfFile` (Windows) reserves disk space upfront, avoiding expensive file-growth operations during writing
- **Zero-copy mmap access**: Data is formatted directly into memory-mapped output buffers, bypassing userspace buffering
- **Deferred disk write**: `nofsync` skips the final request to write buffered data to disk. The command may finish sooner, but a system crash or power loss can leave the file incomplete
- **Windows async batching**: On Windows, multiple output buffers can be written concurrently
- **Optional write-through**: `direct` opens the file with `FILE_FLAG_WRITE_THROUGH` on Windows; it does not bypass the page cache, because unbuffered writes would need sector-aligned rows. macOS and Linux builds ignore it
- **SIMD-accelerated data loading**: AVX2/SSE2/NEON bulk load from Stata with 16x loop unrolling in `ctools_data_io`
- **Parallel formatting**: Multiple OpenMP threads convert numeric values to text concurrently
- **Automatic quoting**: Quotes fields containing delimiters, quotation marks, or newlines

### Timing Breakdown

With `verbose`, you'll see:
- **Load time**: Time to read data from Stata
- **Write time**: Time to format and write the file
- **Throughput**: MB/s writing speed

## Stored Results

`cexport delimited` stores the following in `r()`:

### Scalars
| Result | Description |
|--------|-------------|
| `r(N)` | Number of observations exported |
| `r(k)` | Number of variables exported |
| `r(time)` | Elapsed time in seconds |

### Macros
| Result | Description |
|--------|-------------|
| `r(filename)` | Name of the output file |

## Quoting Behavior

| Mode | Behavior |
|------|----------|
| Default | Quote strings containing delimiter, quotes, or newlines |
| `quote` | Quote all string fields unconditionally |
| `noquoteif` | Never quote fields. Delimiters or newlines within a value may then be read as field or row boundaries |

Embedded quotes within strings are escaped by doubling them (`"` becomes `""`).

## Output Format

- Ordinary numeric missing values are written as empty fields; `.a`–`.z` retain their codes
- String missing values are written as empty strings
- Lines end with LF on every platform; the `crlf` option writes CRLF instead
- UTF-8 encoding is used for output

## Technical Notes

- Value labels can be decoded to their text representation
- Variables are exported in the order specified (or dataset order if `varlist` is omitted)
- The `replace` option is required to overwrite existing files

## Comparison with Native `export delimited`

| Feature | Stata `export delimited` | `cexport delimited` |
|---------|--------------------------|---------------------|
| Implementation | Stata | C with OpenMP |
| Parallelization | No | Yes |
| Chunked I/O | No | Yes |

## See Also

- [ctools Overview](../../README.md)
- [cimport](cimport.md) - Fast CSV import

## Excel and output fidelity

`cexport excel using "output file.xlsx", firstrow(variables) replace` writes an
XLSX workbook. `sheet("A sheet")` and filenames may contain spaces. Text resembling an option, such as `threads(2)`, is exported as text. `threads(#)` also limits Excel export workers.

Delimited and Excel exports check variable attributes, then write a temporary file
in the destination directory and move it to the requested filename after writing
succeeds. A failed export leaves an existing destination intact. Without `replace`,
the command returns an error if that filename exists, including if another process
creates it during export. The destination directory must allow temporary files.
As with `export`, `replace` turns a symbolic-link destination into a regular file
and leaves its target unchanged, sheet modification writes through the link, and a
write-protected destination is not replaced.

New Excel workbooks use the 1900 date system; edits preserve the existing workbook's date system. Daily dates retain their calendar day;
`%tc` retains milliseconds as fractional days, and `%tC` is converted to `%tc`
with `cofC()` first. Leap seconds cannot be represented in Excel. Monthly,
quarterly, weekly, half-yearly, and yearly dates use the first day of the period;
`%tb` dates are converted through their business calendar and `%tg` values are
written as numbers.
`datafmt` and `datestring()` instead export date/time values as text.

Cell formats follow `export excel`: `m/d/yy` for dates, `m/d/yy h:mm` for datetimes (`h:mm:ss AM/PM` when a `%tc` format begins with an hour, minute or am/pm code), and for numbers `0.00` (`%9.2f`), `#,##0.00` (`%9.2fc`), `#,##0` (any `%gc`), `0.00E+00` (`%9.2e`) or `0` (any other format, so `%9.0g` shows 1234.57 as 1235). See [import and export formats](../IO_PARITY.md) for the remaining differences.

`keepcellfmt` preserves the existing cell formats in the written rectangle when modifying a worksheet; variables with `%t` formats still receive their formats and other numbers receive none, as with `export excel`. Cells outside that rectangle and other worksheets retain their content and formatting. Use `sheet(..., modify)` for edits and `sheet(..., replace)` to replace the selected worksheet. Whole-workbook `replace` retains only style definitions with `keepcellfmt`, not the old cells, their style assignments, or other sheets.

## More details

Excel export writes data in the first row by default, matching Stata. Use `firstrow(variables)` for variable names or `firstrow(varlabels)` for labels. See [import and export formats](../IO_PARITY.md) for supported formats, differences from `export`, and testing.
