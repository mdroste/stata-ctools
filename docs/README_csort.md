# csort

High-performance C-accelerated sorting for Stata datasets.

## Overview

`csort` is a drop-in replacement for Stata's `sort` command that uses a parallel Least Significant Digit (LSD) radix sort algorithm implemented in C. It provides 2-4x speedup over native Stata sorting, especially on large datasets.

## Syntax

```stata
csort varlist [if] [in] [, options]
```

## Options

| Option | Description |
|--------|-------------|
| `verbose` | Display progress information including variable names and sort key indices |
| `timeit` | Display timing breakdown (load, sort, store, total) |

## Examples

```stata
* Basic usage - sort by single variable
csort id

* Sort by multiple variables
csort state year firm_id

* Sort with timing output
csort myvar, verbose timeit

* Sort with if/in conditions
csort price if foreign == 1

* Sort a subset of observations
csort value in 1/10000
```

## Performance

`csort` achieves its speed advantage through:

1. **Parallel LSD Radix Sort**: O(N*k) complexity where k is the key width, compared to O(N log N) for comparison-based sorts
2. **OpenMP Parallelization**: Utilizes multiple CPU cores for parallel bucket operations
3. **Cache-Efficient Memory Access**: 8-way loop unrolling for optimized data loading and storing
4. **Stable Sort**: Maintains relative order of equal elements (same as Stata's default)

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

| Feature | Stata `sort` | `csort` |
|---------|--------------|---------|
| Algorithm | Modified merge sort | Parallel LSD radix sort |
| Parallelization | No | Yes (OpenMP) |
| Typical Speedup | 1x (baseline) | 2-4x |
| Stable Sort | Yes (with `stable`) | Yes (always) |

## See Also

- [ctools Overview](../README.md)
- [cmerge](README_cmerge.md) - Fast dataset merging (uses csort internally)
