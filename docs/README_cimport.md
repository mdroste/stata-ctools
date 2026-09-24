# cimport delimited

High-performance C-accelerated CSV/delimited text import.

## Overview

`cimport delimited` is a high-performance replacement for Stata's `import delimited` command that uses multi-threaded parallel parsing.

## Syntax

```stata
cimport delimited [using] filename [, options]
```

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `clear` | Clear data in memory before loading |
| `delimiters(chars)` | Field delimiter (default: automatic detection). Use `tab` or `\t` for tab-delimited |

### Variable Name Options
| Option | Description |
|--------|-------------|
| `varnames(rule)` | How to read variable names: `1` (first row, default) or `nonames` |
| `case(option)` | Variable name case: `preserve` (default), `lower`, or `upper` |

### Parsing Options
| Option | Description |
|--------|-------------|
| `bindquotes(option)` | Quote handling: `loose` (default) or `strict` |
| `stripquotes` | Accepted for compatibility; no additional stripping is implemented |
| `encoding(encoding)` | Automatic detection or explicit supported encoding; UTF-32 is rejected |
| `rowrange([start][:end])` | Range of rows to import |

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

* Import only rows 1000-2000
cimport delimited using bigdata.csv, clear rowrange(1000:2000)

* Import from row 500 to end
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

### Speedup Tricks

- **Memory-mapped I/O**: The entire CSV file is `mmap`'d for zero-copy access, avoiding `read()` syscall overhead and letting the OS handle page-level caching
- **OS-level prefetch hints**: `madvise(MADV_SEQUENTIAL | MADV_WILLNEED)` on POSIX tells the kernel to read ahead aggressively
- **SIMD CSV parsing**: SSE2 (x86) and NEON (ARM64) intrinsics scan 16 bytes at a time for delimiters and newlines, accelerating field boundary detection
- **8-way unrolled quote detection**: Inner loop processes 8 characters per iteration to validate quoted fields with minimal branch overhead
- **Parallel row processing**: Each OpenMP thread processes an independent chunk of rows; chunk boundaries are determined by a fast newline pre-scan
- **Compact field references**: Each parsed field is stored as a 64-bit offset + 32-bit length (12 bytes total), minimizing metadata memory for files with millions of fields
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
- String variables are sized to fit the longest value

## Technical Notes

- Encoding is detected automatically. Explicit encodings include UTF-8, UTF-16LE/BE, ASCII, Latin-1/9, Windows-1252, and Mac Roman. UTF-32 is rejected; convert it to UTF-8 first.
- Fields exceeding 2045 UTF-8 bytes are rejected before clearing the current data.
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

- [ctools Overview](../README.md)
- [cexport](README_cexport.md) - Fast CSV export

## Validation and failure behavior

Excel import honors workbookPr date1904. Daily date cells become Stata daily dates, and datetime cells become Stata milliseconds. In the 1900 system, serial 60 (the nonexistent 29 February 1900) is missing; serials below and above it use their correct offsets. Inline-string headers and values are retained in both the serial and parallel parser paths. Imported display formats are not inferred from arbitrary Excel formatting.
