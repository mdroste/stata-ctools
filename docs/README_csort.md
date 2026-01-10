# csort

High-performance C-accelerated sorting for Stata datasets.

## Overview

`csort` is a drop-in replacement for Stata's `sort` command that provides 2-4x speedup over native Stata sorting. It supports multiple sorting algorithms optimized for different use cases:

- **LSD Radix Sort** (default): Best for fixed-width keys and general-purpose sorting
- **MSD Radix Sort**: Best for variable-length strings with common prefixes
- **Timsort**: Best for partially sorted data (panel data, time series)

## Syntax

```stata
csort varlist [if] [in] [, options]
```

## Options

| Option | Description |
|--------|-------------|
| `algorithm(string)` | Sort algorithm: `lsd` (default), `msd`, or `timsort` |
| `verbose` | Display progress information including variable names and sort key indices |
| `timeit` | Display timing breakdown (load, sort, store, total) |

## Examples

```stata
* Basic usage - sort by single variable (uses LSD radix sort by default)
csort id

* Sort by multiple variables
csort state year firm_id

* Sort with timing output
csort myvar, verbose timeit

* Sort with if/in conditions
csort price if foreign == 1

* Sort a subset of observations
csort value in 1/10000

* Use MSD radix sort (better for variable-length strings)
csort company_name, algorithm(msd)

* Use Timsort (better for partially sorted data)
csort id time, algorithm(timsort)

* Sort panel data by ID then time (Timsort excels here)
csort panelid year, algorithm(timsort) verbose
```

## Algorithm Selection Guide

### LSD Radix Sort (Default)
**Best for:**
- Integer or fixed-width numeric keys
- Data with uniform key distribution
- General-purpose sorting when you're unsure

**How it works:**
- Processes bytes from least significant to most significant
- 8-bit radix (256 buckets per pass)
- 8 passes for 64-bit doubles
- Stable sort: maintains relative order of equal elements

### MSD Radix Sort
**Best for:**
- Variable-length strings (names, addresses, identifiers)
- Strings with common prefixes (state names, company names)
- Data where elements can be distinguished by their first few characters

**How it works:**
- Processes characters from most significant (leftmost) first
- Short-circuits when elements become distinguishable
- Uses insertion sort for small buckets
- O(N*k) where k is the average distinguishing prefix length

**Example use cases:**
```stata
* Sorting names - MSD stops early when prefixes differ
csort last_name first_name, algorithm(msd)

* Sorting categorical strings with common prefixes
csort state_name, algorithm(msd)
```

### Timsort
**Best for:**
- Partially sorted data (common in panel/time series data)
- Data with natural "runs" of sorted elements
- When data might already be mostly sorted

**How it works:**
- Finds natural ascending/descending runs in the data
- Extends short runs using binary insertion sort
- Merges runs using an optimized merge strategy
- O(N) best case for already-sorted data, O(N log N) worst case

**Example use cases:**
```stata
* Panel data already sorted by ID, needs time sort
sort panelid  // Already sorted by ID
csort panelid year, algorithm(timsort)  // Timsort exploits existing order

* Re-sorting data that's mostly in order
csort date_var, algorithm(timsort)
```

## Performance

All three algorithms are optimized for performance:

### LSD Radix Sort Performance
- **Complexity**: O(N*k) where k is key width (8 for doubles)
- **Parallelization**: OpenMP parallel histogram and scatter
- **Best case**: 2-4x faster than Stata's sort

### MSD Radix Sort Performance
- **Complexity**: O(N*k) where k is average distinguishing prefix length
- **Parallelization**: OpenMP task parallelism for independent buckets
- **Best case**: 2-5x faster for variable-length strings

### Timsort Performance
- **Complexity**: O(N) best case, O(N log N) worst case
- **Parallelization**: Parallel permutation application
- **Best case**: Near-instant for already-sorted data

### Common Optimizations
1. **Cache-Efficient Memory Access**: 8-way loop unrolling for data loading/storing
2. **Stable Sort**: All three algorithms maintain relative order of equal elements
3. **Parallel Permutation**: Variable data reordering done in parallel

### Timing Breakdown

When using the `timeit` option, you'll see:
- **Load**: Time to load data from Stata into C memory
- **Sort**: Time for the actual radix sort operation
- **Store**: Time to write sorted data back to Stata
- **Total**: Overall elapsed time

## Supported Variable Types

- Numeric variables (all Stata numeric types)
- String variables
- Multiple sort keys (sorts by first variable, then second, etc.)

## Stored Results

`csort` stores timing information in the following scalars (when `timeit` is specified):

| Scalar | Description |
|--------|-------------|
| `_csort_time_load` | Data loading time (seconds) |
| `_csort_time_sort` | Sort time (seconds) |
| `_csort_time_store` | Data storing time (seconds) |
| `_csort_time_total` | Total elapsed time (seconds) |

## Technical Notes

- The sort is stable, meaning observations with equal sort keys maintain their original relative order
- After sorting, Stata's internal `sortedby` characteristic is updated appropriately
- Missing values are sorted to the end (consistent with Stata's convention)

## Comparison with Native `sort`

| Feature | Stata `sort` | `csort` (LSD) | `csort` (MSD) | `csort` (Timsort) |
|---------|--------------|---------------|---------------|-------------------|
| Algorithm | Modified merge sort | LSD radix | MSD radix | Hybrid merge/insert |
| Complexity | O(N log N) | O(N*k) | O(N*k) | O(N) to O(N log N) |
| Parallelization | No | Yes (OpenMP) | Yes (OpenMP) | Partial |
| Best For | General | Fixed-width keys | Variable strings | Partially sorted |
| Stable Sort | Yes (with `stable`) | Yes | Yes | Yes |
| Typical Speedup | 1x | 2-4x | 2-5x (strings) | 2-10x (sorted) |

## See Also

- [ctools Overview](../README.md)
- [cmerge](README_cmerge.md) - Fast dataset merging (uses csort internally)
