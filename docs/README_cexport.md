# cexport delimited

High-performance C-accelerated delimited text export.

## Overview

`cexport delimited` is a high-performance replacement for Stata's `export delimited` command that uses a C plugin with parallel data loading and chunked formatting. It can achieve up to 25x speedup over native Stata export.

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
| `datafmt` | Use display formats for numeric variables (not yet implemented) |

### Reporting Options
| Option | Description |
|--------|-------------|
| `verbose` | Display progress information during export |
| `timeit` | Display timing breakdown |

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
cexport delimited using output.csv, replace verbose timeit

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

With `timeit`, you'll see:
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

- Numeric missing values are written as empty fields
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
| Typical Speedup | 1x (baseline) | 10-25x |
| Chunked I/O | No | Yes |

## See Also

- [ctools Overview](../README.md)
- [cimport](README_cimport.md) - Fast CSV import
