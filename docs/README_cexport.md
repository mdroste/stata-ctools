# cexport

C writers for delimited text, Excel (.xls/.xlsx), dBase, SAS7BDAT, SAS XPORT5/8,
SPSS SAV, and ESRI shapefiles. Full native parity is not yet established; see
[the compatibility audit](IO_PARITY.md).

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

Statistical writers accept variable lists and if/in selection. Native SAS export
is unavailable in the local installation, so that writer has reader-based
validation rather than a native writer comparison. XLS supports long text via
BIFF8 SST continuation. XLSX loads long text through a dedicated C path rather
than the general str2045 transfer buffer. Both Excel writers support
`sheet(..., modify/replace)`, adding worksheets, and `keepcellfmt`. Modification
retains untouched cells and workbook content; keepcellfmt preserves styles in
the written rectangle. XLSX edits move overwritten shared-formula anchors to
surviving cells, remove calculation chains, and request recalculation on opening.
See the audit for complex-workbook validation limits.

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

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information during export |

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

### Speedup Tricks

- **Vectored I/O**: Uses `pwritev()` on POSIX and batched async `WriteFile` on Windows to write multiple buffers in a single syscall, reducing kernel transitions
- **Async mmap flush**: Uses `MS_ASYNC` instead of `MS_SYNC` so the OS writes dirty pages in the background while formatting continues
- **File pre-allocation**: `posix_fallocate` (POSIX) or `SetFilePointerEx`+`SetEndOfFile` (Windows) reserves disk space upfront, avoiding expensive file-growth operations during writing
- **Zero-copy mmap access**: Data is formatted directly into memory-mapped output buffers, bypassing userspace buffering
- **Deferred fsync**: Optional `CEXPORT_IO_FLAG_NOFSYNC` skips the final `fsync`, trading durability for speed on non-critical exports
- **Windows async batching**: `WaitForMultipleObjects` allows concurrent overlapped writes, saturating disk bandwidth
- **Optional direct I/O**: `FILE_FLAG_NO_BUFFERING` on Windows bypasses the OS page cache for large exports, avoiding cache pollution
- **SIMD-accelerated data loading**: AVX2/SSE2/NEON bulk load from Stata with 16x loop unrolling in `ctools_data_io`
- **Parallel formatting**: Multiple OpenMP threads convert numeric values to text concurrently
- **Smart quoting**: Only quotes fields containing delimiters, quotes, or newlines—avoids unnecessary work

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
| `noquoteif` | Never quote (use with caution) |

Embedded quotes within strings are escaped by doubling them (`"` becomes `""`).

## Output Format

- Ordinary numeric missing values are written as empty fields; `.a`–`.z` retain their codes
- String missing values are written as empty strings
- Line endings use the system default (LF on Unix/Mac, CRLF on Windows)
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

- [ctools Overview](../README.md)
- [cimport](README_cimport.md) - Fast CSV import

## Excel and output fidelity

`cexport excel using "output file.xlsx", firstrow(variables) replace` writes an
XLSX workbook. `sheet("A sheet")` and filenames may contain spaces. Text fields
are transferred separately from plugin options, including text resembling
`threads(2)`. `threads(#)` also limits Excel export workers.

Delimited and Excel exports validate per-column metadata and write a sibling temporary
file before publishing it. A failed export leaves an existing destination intact;
without `replace`, publication refuses an existing file, including one created
while export is running. The destination directory must allow temporary files.

Excel export uses the 1900 date system. Daily dates retain their calendar day;
`%tc` retains milliseconds as fractional days, and `%tC` is converted to `%tc`
with `cofC()` first. Leap seconds cannot be represented in Excel. Monthly,
quarterly, weekly, half-yearly, and yearly dates use the first day of the period.
`datafmt` and `datestring()` instead export date/time values as text.

`keepcellfmt` only copies an existing style-definition table. It does not retain worksheet cell-style assignments or other sheets; existing date style indices may be incompatible with newly written date cells. It is not a formatted-template update facility.

## Parity audit

Excel export now writes data in the first row by default, matching Stata. Use `firstrow(variables)` for variable names or `firstrow(varlabels)` for labels. See [the full import/export audit](IO_PARITY.md) for exact regression coverage and remaining gaps.
