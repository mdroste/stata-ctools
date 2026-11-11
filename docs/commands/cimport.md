# cimport

C readers for delimited text, Excel (.xls/.xlsx), Apache Parquet, dBase,
SAS7BDAT and catalogs, SAS XPORT5/8, SPSS SAV/ZSAV, and ESRI shapefiles; FRED
and Haver imports run Stata's own commands. See
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
cimport parquet using "data.parquet", describe
cimport parquet make price using "auto.parquet", rowrange(64:74) clear
```

## Parquet

`cimport parquet` has the syntax of StataNow's `import parquet`:

```stata
cimport parquet [using] filename [, clear rowrange([start][:end]) favormemory snapshot fast]
cimport parquet columnlist using filename [, options]
cimport parquet [using] filename, describe
```

It creates the variables `import parquet` creates, with the same names, types,
formats, labels and values, prints the same notes and errors, and stores
`r(N)` and `r(k)` (`r(N_rows)` and `r(n_columns)` with `describe`). It reads
every Parquet type, row groups, data pages v1 and v2, every encoding, and
Snappy, Gzip, Brotli, Zstandard, LZ4 and LZ4_RAW compression; nested columns
are skipped with a note, as in Stata. `favormemory` has no effect. Dictionary (categorical)
columns written by pandas or pyarrow become labeled numeric variables, as in
Stata. See [Parquet](../IO_PARITY.md#parquet) for the type mapping.

## FRED and Haver

`cimport fred`, `cimport haver` and `cimport haverdirect` pass their arguments
unchanged to `import fred`, `import haver` and `import haverdirect`, so they
behave exactly as Stata's commands; `snapshot` and `fast` do not apply.

Excel supports sheet names, cell ranges, firstrow, allstring with an optional numeric format, column lists and
describe. As with `import excel`, `all` abbreviates `allstring`, only the full `allstring(format)` carries a
format, and the workbook kind (.xlsx package or binary .xls) comes from the file's content, so `.xlsm`,
`.xltx` or any other name reads. SAS/SPSS support variable and observation selection, read as `import sas`
and `import spss` read them: exact, case-sensitive names (a hyphen separates names, so `a-b` is `a` and `b`),
names the file lacks ignored, and "no matching names found in data" (r(111)) without a match. The names and an
`in` range are applied while the file is read, so only the selection is held in memory; `if` is applied
afterwards to the selected variables. SAS supports
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

`cimport delimited` reads delimited files in C using Stata's `import delimited` syntax. Tests compare its results with `import delimited`.

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
| `encoding(encoding)` | Automatic detection, or a Java charset name or alias as `import delimited` takes them (`latin1`, `cyrillic`, `KOI8-U`, `x-MacCyrillic`, `IBM037`, `windows-31j`...) |
| `rowrange([start][:end])` / `colrange([start][:end])` | Select rows and columns before determining variable types. Row numbers count file lines, including headers, blank lines, and line breaks within quoted fields |
| `maxquotedrows(#)` | Strict quoted-row limit; accepts `unlimited` |
| `parselocale(locale)` | C numeric profiles from the pinned OpenJDK locale tables; no installed Java runtime is needed |
| `favorstrfixed` | Prefer fixed strings up to Stata's 2,045-byte limit |
| `stringcols(numlist \| _all)` | Force selected columns to string |
| `numericcols(numlist \| _all)` | Force selected columns to numeric |
| `asdouble`, `asfloat` | Select numeric storage precision |

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information and throughput (MB/s) |
| `threads(#)` | Maximum threads; zero uses the ctools thread setting |
| `snapshot`, `fast` | Use a whole-dataset recovery snapshot, or skip recovery; mutually exclusive |
| `tempdta` | Load strL-heavy data through a temporary `.dta` file (faster, but the file may be written to disk) |

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

Import proceeds in three stages: scan the file, create variables, and load values.

### Phase 1: Scan
- Parse the file to determine column types and widths
- Automatic type inference (numeric vs. string)
- Detect maximum string lengths

### Phase 2: Create
- Create Stata variables with appropriate types
- Allocate memory using the column types and widths found during the scan

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
`rowrange()` and files adding columns after the header also use that path.
Large files with `strL` columns keep the block cache: their values are read
from the cached field references, or by reparsing the blocks that hold no
text references. Long strings and late type changes can retain the source
mapping.
The cache budget excludes worker scratch space, allocator overhead, input
buffers and Stata datasets; it is not a limit on total process memory or a
measurement of free RAM. `verbose` reports cached bytes and row coverage.

### Implementation details

- **Memory-mapped I/O**: The entire CSV file is `mmap`'d for zero-copy access, avoiding `read()` syscall overhead and letting the OS handle page-level caching
- **OS-level prefetch hints**: `madvise(MADV_SEQUENTIAL | MADV_WILLNEED)` on POSIX requests that the operating system read upcoming file data in advance
- **SIMD CSV parsing**: SSE2 (x86) and NEON (ARM64) intrinsics scan 16 bytes at a time for delimiters and newlines, accelerating field boundary detection
- **8-way unrolled quote detection**: Checks eight characters per iteration for quotation marks
- **Parallel row processing**: Each OpenMP thread processes an independent chunk of rows; chunk boundaries are determined by a fast newline pre-scan
- **Compact block storage**: Large inputs are cached in blocks using storage appropriate to each column's type. When the cache reaches its memory limit, previously cached rows are retained and the remaining rows are read again during loading
- **Arena allocator for strings**: String values are bulk-allocated in a thread-safe arena with atomic CAS, avoiding per-string `malloc` calls
- **Persistent thread pool**: Worker threads are reused across scan and load phases
- **UTF-8 in place**: A file named or detected as UTF-8 that is well-formed UTF-8 is read from the mapping directly instead of being copied through a converter. One parallel pass over the input checks UTF-8 validity and looks for binary zeros and lone carriage returns, so input without them is not rewritten
- **Quoted text skipped in runs**: Inside quotes the parser jumps from one quote (or line end) to the next instead of reading byte by byte, and decoding copies the text between quotes in runs. With `bindquotes(strict)`, the `maxquotedrows()` check runs inside the parallel parse, on the rows that span lines
- **strL transfer**: The plugin API cannot write strL values (`SF_sstore()` on a strL variable stores nothing). By default Mata stores them in blocks of up to 65,536 rows with one `st_sstore()` per column, and nothing is written to disk: on Stata 16 and later the plugin keeps each block in memory and hands Mata the positions of its values through a temporary frame. When `cimport` reads the file's own bytes (no encoding, delimiter, line-ending or binary-zero rewrite), a value found in the file is sent only as its position, and Mata cuts it from the block's byte range of the CSV, undoing doubled quotes or CR LF in it with `subinstr()` when that reproduces the decoded value; the few values it cannot reproduce travel as text in global macros. Empty values of mostly empty columns are not sent at all. Stata 14 and 15, and rewritten input, use a short-lived transfer file instead, marked temporary on Windows so that the system keeps it in memory if it can. With option `tempdta`, when the strL text is at least as large as the other variables' data and the import may replace the dataset outright, the plugin instead writes the whole imported dataset, filled in parallel by the same store workers, as a temporary format-118 `.dta` that one `use` loads, much faster than Mata can store strLs. If that cannot be done (an uncached suffix too long to hold in memory, for example), it writes only the strL columns and the other variables are added and put in file order. `c(filename)`, `c(filedate)` and `c(changed)` end as they do without the option. The option is off by default because the temporary file, about the size of the imported data, may be written to disk
- **Batched variable creation**: Names, types, labels and formats are settled in Mata and all variables are created with one `st_addvar()` call

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

- Automatic encoding detection uses a C adaptation of ICU 67.1's statistical recognizers and its 8,000-byte sample.
- `encoding()` takes the names and aliases of the Java charsets `import delimited` uses (StataNow 19's Java 21 runtime), in any case. The 104 single-byte charsets (ISO-8859-x, windows-125x, KOI8-R/U, the Mac and IBM PC code pages, EBCDIC pages, TIS-620...) are decoded from Java's own tables on every platform, so each byte reads as it does natively and a byte the charset does not map is dropped. Multi-byte charsets (Shift_JIS, windows-31j, GBK, GB18030, Big5, EUC-JP/KR, x-windows-949, the UTF encodings...) are decoded by the system converter (iconv; Windows code pages on Windows). On Windows a few have no code page and give r(198): `Big5-HKSCS` and its variants, `x-EUC-TW`, `x-ISCII91`, `ISO-2022-JP-2` and the IBM EBCDIC CJK pages. Invalid multi-byte sequences become replacement characters on Windows, where native drops them.
- Fields exceeding 2,045 UTF-8 bytes use strL. For shorter text, storage type follows `import delimited`'s rules based on maximum and average text length, unless `favorstrfixed` is specified.
- Empty cells are imported as missing (`.` for numeric, `""` for string)
- Quoted fields handle embedded delimiters and newlines correctly
- Variable names are sanitized to be valid Stata names

## Comparison with Native `import delimited`

| Feature | Stata `import delimited` | `cimport delimited` |
|---------|--------------------------|---------------------|
| Implementation | Stata | C with OpenMP |
| Parallelization | No | Yes |

Memory use depends on the dataset, options, and number of threads.

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
