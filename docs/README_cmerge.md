# cmerge

High-performance C-accelerated merge for Stata datasets.

## Overview

`cmerge` is a drop-in replacement for Stata's `merge` command that performs all merge operations entirely in C. It uses parallel data loading, radix sort, and optimized merge algorithms to achieve 2-5x speedup over native Stata merging.

## Syntax

```stata
cmerge merge_type varlist using filename [, options]
```

Where `merge_type` is one of:
- `1:1` - One-to-one merge
- `m:1` - Many-to-one merge
- `1:m` - One-to-many merge
- `m:m` - Many-to-many merge (not recommended)

## Options

| Option | Description |
|--------|-------------|
| `keep(results)` | Which observations to keep: `match`, `master`, `using`, `1`, `2`, `3` |
| `assert(results)` | Verify merge results match expected values |
| `generate(varname)` | Name for merge indicator variable (default: `_merge`) |
| `nogenerate` | Do not create merge indicator variable |
| `keepusing(varlist)` | Variables to keep from using dataset |
| `sorted` | Assert data is already sorted (skip internal sorting) |
| `force` | Allow string/numeric type mismatches in key variables |
| `noreport` | Suppress merge result table |
| `verbose` | Display timing and progress information |
| `nolabel` | Do not copy value labels from using data |
| `timeit` | Display total elapsed time |

## Examples

```stata
* One-to-one merge
cmerge 1:1 id using data2.dta

* Many-to-one merge keeping only matched observations
cmerge m:1 state year using statedata.dta, keep(match)

* Keep only specific variables from using dataset
cmerge 1:1 id using fulldata.dta, keepusing(var1 var2)

* Merge with verbose output and no merge variable
cmerge m:1 id using lookup.dta, nogenerate verbose

* Assert all observations should match
cmerge 1:1 id using data2.dta, assert(match)

* Force merge when key variable types differ
cmerge m:1 code using lookup.dta, force
```

## Performance

`cmerge` achieves its speed advantage through a three-phase approach:

### Phase 1: Load Using Dataset
- Loads the using dataset into C memory
- Uses parallel data loading with 8-way loop unrolling

### Phase 2: Execute Merge
- Performs high-performance radix sort on key variables (if not pre-sorted)
- Uses single-pass sorted merge join
- Parallel output construction

### Phase 3: Apply Options
- Applies `keep()` and `assert()` filters
- Renames placeholder variables to match using dataset

### Optimization Features

| Feature | Description |
|---------|-------------|
| Hash-based lookup | For m:1/1:1 merges when `keep()` excludes using-only rows |
| Pre-sorted detection | Automatically detects if data is already sorted |
| Parallel radix sort | OpenMP-parallelized LSD radix sort |

## Merge Indicator Values

The `_merge` variable (or variable specified by `generate()`) indicates:

| Value | Label | Meaning |
|-------|-------|---------|
| 1 | master only | Observation from master, no match in using |
| 2 | using only | Observation from using, no match in master |
| 3 | matched | Observation matched in both datasets |

## Stored Results

`cmerge` stores the following in `r()`:

### Scalars
| Result | Description |
|--------|-------------|
| `r(N)` | Number of observations after merge |
| `r(N_1)` | Number of observations only in master |
| `r(N_2)` | Number of observations only in using |
| `r(N_3)` | Number of matched observations |
| `r(time)` | Elapsed time in seconds |

### Macros
| Result | Description |
|--------|-------------|
| `r(using)` | Name of using file |
| `r(keyvars)` | Key variable names |

## Technical Notes

- Key variable types must match between master and using (unless `force` is specified)
- The using file must be a Stata `.dta` file (extension can be omitted)
- Variable name conflicts: master values are kept by default
- Missing values in key variables are treated as valid match keys

## Comparison with Native `merge`

| Feature | Stata `merge` | `cmerge` |
|---------|---------------|----------|
| Implementation | Stata + Mata | Pure C |
| Parallelization | No | Yes (OpenMP) |
| Sort Algorithm | Stata default | Radix sort |
| Typical Speedup | 1x (baseline) | 2-5x |
| Memory Efficiency | Standard | Optimized |

## See Also

- [ctools Overview](../README.md)
- [csort](README_csort.md) - High-performance sorting
