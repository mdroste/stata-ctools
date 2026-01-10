# cimport delimited

High-performance C-accelerated CSV/delimited text import.

## Overview

`cimport delimited` is a high-performance replacement for Stata's `import delimited` command that uses multi-threaded parallel parsing. It can achieve up to 50x speedup on large files through efficient memory management and parallel data loading.

## Syntax

```stata
cimport delimited using filename [, options]
```

## Options

### Main Options
| Option | Description |
|--------|-------------|
| `clear` | Clear data in memory before loading |
| `delimiters(chars)` | Field delimiter (default: comma). Use `tab` or `\t` for tab-delimited |

### Variable Name Options
| Option | Description |
|--------|-------------|
| `varnames(rule)` | How to read variable names: `1` (first row, default) or `nonames` |
| `case(option)` | Variable name case: `preserve` (default), `lower`, or `upper` |

### Parsing Options
| Option | Description |
|--------|-------------|
| `bindquotes(option)` | Quote handling: `strict` (default) or `loose` |
| `stripquotes` | Remove surrounding quotes from string values |
| `encoding(encoding)` | File encoding (currently only UTF-8 supported) |
| `rowrange([start][:end])` | Range of rows to import |

### Performance Options
| Option | Description |
|--------|-------------|
| `fast` | Use fast DTA mode (writes temp DTA file, then loads) |

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

* Fast import for large files
cimport delimited using bigdata.csv, clear fast verbose

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

### Fast Mode

The `fast` option provides additional speedups for very large files by:
1. Writing a temporary DTA file using C
2. Loading the DTA file with Stata's native `use` command
3. Cleaning up the temporary file

This bypasses some Stata overhead and can be significantly faster for files with millions of rows.

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
| Comma | `delimiters(",")` or default |
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

- UTF-8 encoding is assumed (other encodings not yet supported)
- Empty cells are imported as missing (`.` for numeric, `""` for string)
- Quoted fields handle embedded delimiters and newlines correctly
- Variable names are sanitized to be valid Stata names

## Comparison with Native `import delimited`

| Feature | Stata `import delimited` | `cimport delimited` |
|---------|--------------------------|---------------------|
| Implementation | Stata | C with OpenMP |
| Parallelization | No | Yes |
| Typical Speedup | 1x (baseline) | 10-50x |
| Memory Efficiency | Standard | Optimized |
| Fast Mode | No | Yes |

## See Also

- [ctools Overview](../README.md)
- [cexport](README_cexport.md) - Fast CSV export
