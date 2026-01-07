# ctools Architecture

This document describes the internal architecture of the ctools Stata plugin suite.

## Overview

ctools is a high-performance C plugin for Stata that provides accelerated versions of common data operations. The architecture is designed around three key principles:

1. **Modularity** - Core routines are reusable across commands
2. **Parallelism** - Operations use multithreading where beneficial
3. **Memory efficiency** - Column-major storage for cache-friendly access

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              Stata Session                                   │
│                                                                             │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────────┐           │
│  │ csort   │ │ cmerge  │ │cimport  │ │cexport  │ │ creghdfe    │  .ado     │
│  │  .ado   │ │  .ado   │ │  .ado   │ │  .ado   │ │   .ado      │  files    │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └──────┬──────┘           │
│       │           │           │           │             │                   │
│       └───────────┴───────────┼───────────┴─────────────┘                   │
│                               │                                             │
│                    plugin call ctools_plugin                                │
│                               │                                             │
└───────────────────────────────┼─────────────────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────────────────┐
│                          ctools.plugin (C)                                    │
│                                                                               │
│  ┌────────────────────────────────────────────────────────────────────────┐  │
│  │                      ctools_plugin.c (Dispatcher)                       │  │
│  │                                                                         │  │
│  │   stata_call() → parses command → dispatches to *_main()               │  │
│  └────────────────────────────────────────────────────────────────────────┘  │
│                                    │                                          │
│       ┌────────────────────────────┼────────────────────────────┐            │
│       │              │             │              │              │            │
│       ▼              ▼             ▼              ▼              ▼            │
│  ┌─────────┐   ┌─────────┐   ┌─────────┐   ┌─────────┐   ┌───────────┐      │
│  │ csort   │   │ cmerge  │   │cimport  │   │cexport  │   │ creghdfe  │      │
│  │ _main() │   │ _main() │   │ _main() │   │ _main() │   │  _main()  │      │
│  └────┬────┘   └────┬────┘   └────┬────┘   └────┬────┘   └─────┬─────┘      │
│       │             │             │             │               │            │
│       └─────────────┼─────────────┼─────────────┼───────────────┘            │
│                     │             │             │                             │
│                     ▼             ▼             ▼                             │
│  ┌────────────────────────────────────────────────────────────────────────┐  │
│  │                     Core Utilities (ctools_*)                          │  │
│  │                                                                         │  │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌────────────┐ │  │
│  │  │ data_load.c  │  │ data_store.c │  │sort_radix.c  │  │ threads.c  │ │  │
│  │  │              │  │              │  │              │  │            │ │  │
│  │  │ Stata → C    │  │ C → Stata    │  │ LSD radix    │  │Thread pool │ │  │
│  │  │ parallel I/O │  │ parallel I/O │  │ parallel     │  │ utilities  │ │  │
│  │  └──────────────┘  └──────────────┘  └──────────────┘  └────────────┘ │  │
│  │                                                                         │  │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐                  │  │
│  │  │  types.c/h   │  │  error.c/h   │  │  timer.c/h   │                  │  │
│  │  │              │  │              │  │              │                  │  │
│  │  │ stata_data   │  │ ctools_msg() │  │ High-res     │                  │  │
│  │  │ stata_var    │  │ ctools_err() │  │ timing       │                  │  │
│  │  └──────────────┘  └──────────────┘  └──────────────┘                  │  │
│  └────────────────────────────────────────────────────────────────────────┘  │
│                                    │                                          │
│                                    ▼                                          │
│  ┌────────────────────────────────────────────────────────────────────────┐  │
│  │               Stata Plugin Interface (SPI) - stplugin.c/h              │  │
│  │                                                                         │  │
│  │   SF_vdata() SF_sdata() SF_vstore() SF_sstore() SF_in1() SF_in2() ... │  │
│  └────────────────────────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────────────────────┘
```

## Memory Layout

### Column-Major Storage

Data is stored in column-major format for cache efficiency during variable operations:

```
Stata Dataset:          C Memory (stata_data):
┌────┬────┬────┐       vars[0].data.dbl = [1.0, 2.0, 3.0]  (contiguous)
│ v1 │ v2 │ v3 │       vars[1].data.dbl = [4.0, 5.0, 6.0]  (contiguous)
├────┼────┼────┤       vars[2].data.dbl = [7.0, 8.0, 9.0]  (contiguous)
│1.0 │4.0 │7.0 │
│2.0 │5.0 │8.0 │       sort_order = [0, 1, 2]  (identity permutation)
│3.0 │6.0 │9.0 │
└────┴────┴────┘
```

### String Variables

String variables use a pointer array with heap-allocated strings:

```
vars[i].data.str:
┌─────────┬─────────┬─────────┐
│ ptr[0]  │ ptr[1]  │ ptr[2]  │
└────┬────┴────┬────┴────┬────┘
     │         │         │
     ▼         ▼         ▼
 "Alice"   "Bob"    "Carol"   (heap-allocated, freed by stata_data_free)
```

## Data Flow Patterns

### Pattern 1: csort (Load → Sort → Store)

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│    Stata     │     │   C Memory   │     │    Stata     │
│   (unsorted) │────►│    (sort)    │────►│   (sorted)   │
└──────────────┘     └──────────────┘     └──────────────┘
                           │
         ctools_data_load()│ctools_sort_radix_lsd()
                           │ctools_data_store()
```

1. `ctools_data_load()` - Parallel read via SF_vdata/SF_sdata
2. `ctools_sort_radix_lsd()` - LSD radix sort with permutation
3. `apply_permutation_to_all_vars()` - Physical reordering
4. `ctools_data_store()` - Parallel write via SF_vstore/SF_sstore

### Pattern 2: cimport (File → Stata)

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│   CSV File   │────►│ Parse & Type │────►│    Stata     │
│              │     │   Detection  │     │  (new data)  │
└──────────────┘     └──────────────┘     └──────────────┘
      │                    │                    │
      │   mmap + parallel  │   SF_vstore()      │
      │   line scanning    │   SF_sstore()      │
```

Modes:
- **Standard**: Scan → Create vars → Load via SPI
- **Fast (DTA)**: Write DTA directly → `use` in Stata

### Pattern 3: cmerge (Two Datasets → One)

```
┌──────────────┐
│    Master    │──┐
│   (in mem)   │  │     ┌──────────────┐     ┌──────────────┐
└──────────────┘  ├────►│   C Merge    │────►│   Merged     │
┌──────────────┐  │     │   Engine     │     │   Dataset    │
│    Using     │──┘     └──────────────┘     └──────────────┘
│   (.dta)     │
└──────────────┘
```

Phase 1: Load using dataset into C memory
Phase 2: Execute merge (sort-merge or hash-merge depending on options)
Phase 3: Write merged result back to Stata

### Pattern 4: creghdfe (HDFE Regression)

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│  y, X, FE    │────►│  HDFE Init   │────►│  Partial Out │
│   vars       │     │  (factors)   │     │  (CG solver) │
└──────────────┘     └──────────────┘     └──────────────┘
                                                 │
                           ┌─────────────────────┘
                           │
                           ▼
                    ┌──────────────┐     ┌──────────────┐
                    │   OLS on     │────►│   VCE +      │
                    │  residuals   │     │   Results    │
                    └──────────────┘     └──────────────┘
```

1. Initialize FE factors, detect singletons
2. Partial out FE via conjugate gradient
3. OLS on demeaned data
4. Compute variance-covariance matrix (robust/cluster)

## Parallel Execution Model

### Thread Pool Architecture

```
ctools_thread_pool:
┌────────────────────────────────────────────────────────────┐
│                                                            │
│   pthread_t threads[N]     void *args      ctools_thread_func
│   ┌───┬───┬───┬───┐       ┌─────────┐     ┌─────────────┐
│   │ 0 │ 1 │ 2 │...│       │args[0]  │────►│ func(arg)   │
│   └───┴───┴───┴───┘       │args[1]  │     │             │
│                           │args[2]  │     │ return NULL │
│                           │...      │     │ on success  │
│                           └─────────┘     └─────────────┘
│                                                            │
└────────────────────────────────────────────────────────────┘
```

### Parallel Data I/O

When `nvars >= 2`, data load/store uses one thread per variable:

```
Thread 0: Load/Store var[0] ─────────────────────►
Thread 1: Load/Store var[1] ─────────────────────►
Thread 2: Load/Store var[2] ─────────────────────►
...
Thread N: Load/Store var[N] ─────────────────────►
                                                  │
                                          pthread_join()
```

### Parallel Radix Sort

For large datasets (`nobs >= MIN_OBS_PER_THREAD * 2`):

```
Phase 1: Parallel Histogram
┌─────────────────────────────────────────────────────┐
│ T0: count[0..N/4]    → local_counts[0][256]        │
│ T1: count[N/4..N/2]  → local_counts[1][256]        │
│ T2: count[N/2..3N/4] → local_counts[2][256]        │
│ T3: count[3N/4..N]   → local_counts[3][256]        │
└─────────────────────────────────────────────────────┘
                         │
                         ▼ combine histograms
                         │
Phase 2: Compute Offsets
┌─────────────────────────────────────────────────────┐
│ global_counts = Σ local_counts                      │
│ global_offsets = prefix_sum(global_counts)          │
│ thread_offsets[t][b] = per-thread starting position │
└─────────────────────────────────────────────────────┘
                         │
                         ▼
Phase 3: Parallel Scatter
┌─────────────────────────────────────────────────────┐
│ T0: scatter[0..N/4]    using thread_offsets[0]     │
│ T1: scatter[N/4..N/2]  using thread_offsets[1]     │
│ T2: scatter[N/2..3N/4] using thread_offsets[2]     │
│ T3: scatter[3N/4..N]   using thread_offsets[3]     │
└─────────────────────────────────────────────────────┘
```

## Radix Sort Details

### Numeric Sorting

IEEE 754 doubles are converted to uint64 for correct unsigned comparison:

```
Original double bits:     Transformed for sorting:
─────────────────────     ───────────────────────
Positive: 0xxx...xxx  →   1xxx...xxx  (flip sign bit)
Negative: 1xxx...xxx  →   ~(1xxx...xxx)  (flip all bits)
Missing:  (any NaN)   →   UINT64_MAX  (sort to end)
```

### String Sorting

LSD radix sort from rightmost character:

```
Pass 0: Sort by char[max_len-1]
Pass 1: Sort by char[max_len-2]
...
Pass N: Sort by char[0]

Shorter strings padded with 0x00 (sort before longer strings with same prefix)
```

### Multi-Key Sorting

Keys processed from least to most significant (LSD principle):

```
sort by (name, age):
  1. Sort by age (least significant)
  2. Sort by name (most significant)

Result: Ties in name are broken by age (stable sort)
```

## Error Handling

### Return Codes (stata_retcode)

```c
STATA_OK = 0               // Success
STATA_ERR_MEMORY = 1       // malloc/calloc failed
STATA_ERR_INVALID_INPUT = 2 // NULL pointer or bad parameter
STATA_ERR_STATA_READ = 3   // SF_vdata/SF_sdata failed
STATA_ERR_STATA_WRITE = 4  // SF_vstore/SF_sstore failed
STATA_ERR_UNSUPPORTED_TYPE = 5 // Unknown variable type
```

### Stata Error Codes

Common codes returned to Stata:
- `0` - Success
- `198` - Syntax error / invalid arguments
- `601` - File not found
- `920` - Memory allocation failure
- `2000` - No observations

### Error Macros

```c
CTOOLS_CHECK_ALLOC(ptr, module, retval)
CTOOLS_CHECK_ALLOC_CLEANUP(ptr, module, cleanup, retval)
```

## creghdfe Architecture

### Module Structure

```
creghdfe/
├── creghdfe_types.h    - HDFE_State, FE_Factor, UnionFind
├── creghdfe_impl.c     - Entry point, dispatches subcommands
├── creghdfe_hdfe.c     - Factor creation, singleton detection
├── creghdfe_solver.c   - Conjugate gradient solver
├── creghdfe_main.c     - Full regression orchestration
├── creghdfe_ols.c      - OLS computation on demeaned data
├── creghdfe_vce.c      - Variance estimation (robust, cluster)
└── creghdfe_utils.c    - Helper functions
```

### HDFE Algorithm

1. **Initialize Factors**: Create level mappings for each FE variable
2. **Detect Singletons**: Iteratively drop obs appearing only once in any FE
3. **Partial Out FE**: For each variable, demean using CG solver:
   ```
   r = y - P*y  (where P projects onto FE space)
   Iteratively: y ← y - mean(y|FE_g) for each g
   ```
4. **OLS**: Standard OLS on demeaned y, X
5. **VCE**: Compute standard errors (unadjusted, robust, or clustered)

### Conjugate Gradient Solver

Used to solve the normal equations for FE projection:

```
Goal: Find x such that (I - P)y = x
Method: Iterative CG with tolerance check

For each iteration:
  1. Project residual onto FE space
  2. Update search direction
  3. Line search for step size
  4. Update solution and residual
  5. Check convergence: ||r|| < tolerance
```

## Performance Considerations

### When to Use Parallel Sort

- `nobs >= MIN_OBS_PER_THREAD * 2` (default: 200,000)
- Below threshold: sequential sort has better cache locality

### When to Use Parallel I/O

- `nvars >= 2` (one thread per variable)
- Below threshold: sequential I/O avoids thread overhead

### Memory Alignment

- All large allocations use 64-byte alignment (cache line)
- Prevents false sharing in parallel code
- Enables SIMD operations where applicable

### Optimizations Summary

| Optimization | Location | Benefit |
|--------------|----------|---------|
| 8x loop unrolling | data_load, data_store | ~20% faster I/O |
| Pointer swapping | radix sort | Avoids memcpy per pass |
| Early exit | radix sort | Skips uniform distributions |
| Pre-cached strlen | string sort | Avoids repeated strlen calls |
| Aligned allocation | everywhere | Better cache performance |
| SD_FASTMODE | compile flag | Skips SPI bounds checks |
