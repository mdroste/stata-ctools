# cimport

C readers for delimited text, Excel (.xls/.xlsx), dBase, SAS7BDAT and catalogs,
SAS XPORT5/8, SPSS SAV/ZSAV, and ESRI shapefiles. See
[import and export formats](../IO_PARITY.md) for what's supported and how it differs from `import`.

## Other file formats

```stata
cimport excel using "data.xls", firstrow clear
cimport excel identifier=A note=C using "data.xlsx", clear
cimport excel using "data.xlsx", describe
cimport sas using "data.sas7bdat", bcat("labels.sas7bcat") clear
cimport spss id amount if amount>0 using "data.sav", clear
cimport spss using "data.zsav", zsav clear
cimport sasxport5 using "data.xpt", member(TABLE) clear
cimport sasxport8 using "data.v8xpt", clear
cimport dbase using "data.dbf", clear
cimport shp using "boundaries.shp", clear
```

Excel supports sheet names, cell ranges, firstrow, allstring with an optional numeric format, column lists and
describe. SAS/SPSS support variable and observation selection; SAS supports
catalogs and explicit encoding. XPORT5 reads `formats.xpf` beside the transport
file unless `novallabels` is requested. Long delimited/Excel/statistical text and SHP
binary headers use strL. With `clear`, statistical, dBase and SHP readers clear
data before parsing, matching native failure behavior. Excel retains existing
data when its container cannot be opened. Other failures or Break restore the
data in memory from an in-memory copy (a frame, Stata 16 or later), not a
snapshot on disk; option `snapshot` uses `preserve` instead and option `fast`
clears the data first and undoes nothing. Excel `describe` is exclusive of
data-loading options.

## Delimited text

## Overview

`cimport delimited` uses C parsing and native-command differential tests to implement Stata's `import delimited` interface.

## Syntax

```stata
cimport delimited [using] filename [, options]
cimport delimited extvarlist using filename [, options]
```

As in `import delimited`, an `extvarlist` renames the first imported variables in order, exactly as typed, after the import; header detection is unchanged.

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `clear` | Clear data in memory before loading |
| `delimiters(chars)` | Delimiter set (default: automatic detection, as native: the most frequent of tab, comma, semicolon, colon and pipe outside quotes in the first 50 nonempty lines, header included); `asstring` treats a sequence literally and `collapse` merges consecutive separators |

### Variable Name Options
| Option | Description |
|--------|-------------|
| `varnames(rule)` | How to read variable names: automatic detection by default; `1` forces the first row, `nonames` disables headers |
| `case(option)` | Variable name case: `lower` (default), `preserve`, or `upper` |

### Parsing Options
| Option | Description |
|--------|-------------|
| `bindquotes(option)` | Quote handling: `loose` (default), `strict`, or `nobind` |
| `stripquotes(default\|yes\|no)` | Control quote removal; `no` retains literal quotes |
| `encoding(encoding)` | Automatic detection or explicit UTF-8/16/32 and system code pages |
| `rowrange([start][:end])` / `colrange([start][:end])` | Select rows and columns before inference; rows are file lines as in native (the header, blank lines and quoted line breaks count) |
| `maxquotedrows(#)` | Strict quoted-row limit; accepts `unlimited` |
| `parselocale(locale)` | C numeric profiles for installed Java locale names |
| `favorstrfixed` | Prefer fixed strings up to Stata's 2,045-byte limit |

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information and throughput (MB/s) |

## Examples

```stata
* Import a CSV file
cimport delimited using data.csv, clear

* Import a tab-delimited file
cimport delimited using data.tsv, clear delimiters(tab)

* Import with verbose output and lowercase variable names
cimport delimited using data.csv, clear case(lower) verbose

* Import only lines 1000-2000
cimport delimited using bigdata.csv, clear rowrange(1000:2000)

* Import from line 500 to end
cimport delimited using bigdata.csv, clear rowrange(500:)

* Import with timing output
cimport delimited using bigdata.csv, clear verbose

* First row is data, not variable names
cimport delimited using data.csv, clear varnames(nonames)
```

## Performance

`cimport` uses a three-phase approach for maximum efficiency:

### Phase 1: Scan
- Parse the file to determine column types and widths
- Automatic type inference (numeric vs. string)
- Detect maximum string lengths

### Phase 2: Create
- Create Stata variables with appropriate types
- Allocate memory efficiently based on scan results

### Phase 3: Load
- Load data into variables using parallel processing
- OpenMP-parallelized parsing

For delimited files of at least 256 MiB, `cimport` caches blocks of typed
column data. Integer columns use a block minimum plus one-, two-, or four-byte
differences; constant and all-missing blocks need no per-row numeric payload.
Decimals use float or double storage, and short strings are decoded once into
the cache. This avoids retaining a double, a field length and a row pointer
for every cell or row.

Actual block allocations share a budget of one quarter of installed RAM
(512 MiB if the OS memory query is unavailable). When the budget fills, each
worker keeps its cached prefix and reparses only its remaining rows during
loading. Type inference still checks every input value. If a later value
changes a column's type, affected cached columns are recovered from the
original CSV so numeric precision and original text are preserved.

When every row is cached and no column still needs the source, `cimport`
releases the mapped CSV and any converted input buffer before Stata allocates
its rollback copy and destination. Small files keep the existing row cache.
`rowrange()`, files adding columns after the header, and `strL` also use that
path. Long strings and late type changes can retain the source mapping.
The cache budget excludes worker scratch space, allocator overhead, input
buffers and Stata datasets; it is not a limit on total process memory or a
measurement of free RAM. `verbose` reports cached bytes and row coverage.

### Speedup Tricks

- **Memory-mapped I/O**: The entire CSV file is `mmap`'d for zero-copy access, avoiding `read()` syscall overhead and letting the OS handle page-level caching
- **OS-level prefetch hints**: `madvise(MADV_SEQUENTIAL | MADV_WILLNEED)` on POSIX tells the kernel to read ahead aggressively
- **SIMD CSV parsing**: SSE2 (x86) and NEON (ARM64) intrinsics scan 16 bytes at a time for delimiters and newlines, accelerating field boundary detection
- **8-way unrolled quote detection**: Inner loop processes 8 characters per iteration to validate quoted fields with minimal branch overhead
- **Parallel row processing**: Each OpenMP thread processes an independent chunk of rows; chunk boundaries are determined by a fast newline pre-scan
- **Compact block storage**: Large inputs use typed column payloads and a shared allocation budget; cached prefixes remain usable if the budget fills
- **Arena allocator for strings**: String values are bulk-allocated in a thread-safe arena with atomic CAS, avoiding per-string `malloc` calls
- **Persistent thread pool**: Worker threads are reused across scan and load phases

### Throughput

With `verbose` output, you'll see:
- File size
- Number of rows and columns
- Import throughput in MB/s
- Timing breakdown

## Stored Results

`cimport delimited` stores the following in `r()`:

### Scalars
| Result | Description |
|--------|-------------|
| `r(N)` | Number of observations imported |
| `r(k)` | Number of variables created |
| `r(time)` | Elapsed time in seconds |

### Macros
| Result | Description |
|--------|-------------|
| `r(filename)` | Name of the imported file |
| `r(encoding)` | Source encoding selected or supplied |
| `r(delimiters)` | Field delimiters used, with tabs reported as `\t` |

## Supported Delimiters

| Delimiter | Syntax |
|-----------|--------|
| Comma | `delimiters(",")` |
| Tab | `delimiters(tab)` or `delimiters("\t")` |
| Semicolon | `delimiters(";")` |
| Pipe | `delimiters("|")` |
| Custom | `delimiters("X")` for any character X |

## Variable Type Inference

`cimport` automatically determines variable types:
- If all values are numeric, creates a numeric variable
- If any value contains non-numeric characters, creates a string variable
- String storage follows decoded widths and the native strL selection rule

## Technical Notes

- Automatic encoding detection uses a C adaptation of ICU 67.1's statistical recognizers and its 8,000-byte sample. Explicit Unicode encodings include UTF-8, UTF-16LE/BE and UTF-32LE/BE; named legacy encodings use system C conversion. Charset alias coverage differs by platform.
- Fields exceeding 2,045 UTF-8 bytes use strL. Shorter fields follow native width/average selection unless favorstrfixed is specified.
- Empty cells are imported as missing (`.` for numeric, `""` for string)
- Quoted fields handle embedded delimiters and newlines correctly
- Variable names are sanitized to be valid Stata names

## Comparison with Native `import delimited`

| Feature | Stata `import delimited` | `cimport delimited` |
|---------|--------------------------|---------------------|
| Implementation | Stata | C with OpenMP |
| Parallelization | No | Yes |
| Memory Efficiency | Standard | Optimized |

## See Also

- [ctools Overview](../../README.md)
- [cexport](cexport.md) - Fast CSV export

## Validation and failure behavior

Excel import honors the 1900 and 1904 date systems. Daily dates become Stata
daily dates through the native millisecond rounding step; timestamps become
Stata milliseconds. Built-in Excel formats and standard/custom date and time
formats retain their native display metadata and allstring output. Serial 60 in the 1900 system
maps to 28 February 1900, following native import. Arbitrary Excel display
formats are not fully mapped to Stata formats.

## More details

See [import and export formats](../IO_PARITY.md) for supported formats, differences from `import`, and testing.
