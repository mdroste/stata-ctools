# ctools Developer Guide

This guide documents the internal C infrastructure for developers contributing to ctools.

## Table of Contents

- [Architecture Overview](#architecture-overview)
- [Core Data Structures](#core-data-structures)
- [Data I/O](#data-io)
- [Sorting Algorithms](#sorting-algorithms)
- [Thread Pool](#thread-pool)
- [Timing and Profiling](#timing-and-profiling)
- [Argument Parsing](#argument-parsing)
- [Error Handling](#error-handling)
- [Memory Management](#memory-management)
- [Adding a New Command](#adding-a-new-command)
- [Performance Guidelines](#performance-guidelines)

---

## Architecture Overview

### Plugin Lifecycle

Every ctools command follows this pattern:

```
Stata → C Plugin → Stata
  1. Load data from Stata into C memory
  2. Process (sort, merge, regress, etc.)
  3. Store results back to Stata
```

### Memory Model

- **Column-major storage**: Each variable is a contiguous array
- **Numeric variables**: `double[]` (8 bytes per observation)
- **String variables**: `char*[]` (pointer array, heap-allocated strings)
- **Aligned allocations**: 64-byte cache line alignment for SIMD/prefetch efficiency

### Key Files

| File | Purpose |
|------|---------|
| `src/ctools_plugin.c` | Main dispatcher - routes commands to implementations |
| `src/ctools_data_io.c` | Parallel data load/store between Stata and C |
| `src/ctools_types.h` | Core data structures and function declarations |
| `src/ctools_config.h` | Performance tuning constants and macros |
| `src/stplugin.c/h` | Stata Plugin Interface (DO NOT MODIFY) |

---

## Core Data Structures

### stata_variable

Represents a single Stata variable in C memory:

```c
typedef struct {
    stata_vartype type;      // STATA_TYPE_DOUBLE or STATA_TYPE_STRING
    size_t nobs;             // Number of observations
    union {
        double *dbl;         // Numeric: contiguous double[nobs]
        char **str;          // String: char*[nobs], heap-allocated
    } data;
    size_t str_maxlen;       // Max string length (string vars only)
    void *_arena;            // Internal: string arena for bulk free
} stata_variable;
```

### stata_data

Complete dataset in C memory:

```c
typedef struct {
    size_t nobs;             // Number of observations
    size_t nvars;            // Number of variables
    stata_variable *vars;    // Array of variables [nvars]
    size_t *sort_order;      // Permutation array [nobs]
} stata_data;
```

### Return Codes

```c
typedef enum {
    STATA_OK = 0,               // Success
    STATA_ERR_MEMORY = 1,       // Memory allocation failed
    STATA_ERR_INVALID_INPUT = 2,// NULL pointer or invalid parameter
    STATA_ERR_STATA_READ = 3,   // Stata plugin read error
    STATA_ERR_STATA_WRITE = 4,  // Stata plugin write error
    STATA_ERR_UNSUPPORTED_TYPE = 5
} stata_retcode;
```

---

## Data I/O

Include: `#include "ctools_types.h"`

### Loading Data from Stata

```c
stata_data data;
stata_data_init(&data);

// Load all variables
stata_retcode rc = ctools_data_load(&data, nvars);

// Or load specific variables (1-based indices)
int var_indices[] = {1, 3, 5};
rc = ctools_data_load_selective(&data, var_indices, 3, 0, 0);
```

### Storing Data to Stata

```c
// Write all variables back (obs1 = first observation, 1-based)
stata_retcode rc = ctools_data_store(&data, obs1);

// Write with permutation (for merge operations)
rc = ctools_stream_var_permuted(var_idx, source_rows, output_nobs, obs1);
```

### Cleanup

```c
stata_data_free(&data);  // Safe to call multiple times
```

---

## Sorting Algorithms

Include: `#include "ctools_types.h"`

### Available Algorithms

| Algorithm | Best For | Function |
|-----------|----------|----------|
| LSD Radix | Fixed-width keys (default) | `ctools_sort_radix_lsd()` |
| MSD Radix | Variable-length strings | `ctools_sort_radix_msd()` |
| Timsort | Partially sorted data | `ctools_sort_timsort()` |
| Sample Sort | Large datasets, many cores | `ctools_sort_sample()` |
| Counting Sort | Integers with small range | `ctools_sort_counting()` |
| Merge Sort | Stable, predictable O(n log n) | `ctools_sort_merge()` |
| IPS4o | Memory-efficient parallel | `ctools_sort_ips4o()` |

### Usage

```c
// Sort by variables 1 and 2 (1-based indices)
int sort_vars[] = {1, 2};
stata_retcode rc = ctools_sort_radix_lsd(&data, sort_vars, 2);

// Or use the unified interface
rc = ctools_sort(&data, sort_vars, 2, SORT_ALG_LSD);

// Get permutation without applying it
size_t *perm = malloc(data.nobs * sizeof(size_t));
rc = ctools_sort_radix_lsd_with_perm(&data, sort_vars, 2, perm);
```

### Choosing an Algorithm

```c
// Check if counting sort is suitable
if (ctools_counting_sort_suitable(&data, var_idx)) {
    ctools_sort_counting(&data, sort_vars, nsort);
} else {
    ctools_sort_radix_lsd(&data, sort_vars, nsort);
}
```

---

## Thread Pool

Include: `#include "ctools_threads.h"`

### Basic Usage

```c
// Define thread arguments
typedef struct {
    int thread_id;
    double *input;
    double *output;
    size_t start, end;
} my_thread_args;

// Thread function
void *my_worker(void *arg) {
    my_thread_args *a = (my_thread_args *)arg;
    for (size_t i = a->start; i < a->end; i++) {
        a->output[i] = process(a->input[i]);
    }
    return NULL;  // NULL = success
}

// Launch threads
my_thread_args args[8];
// ... fill in args ...

int result = ctools_parallel_run(my_worker, args, 8, sizeof(my_thread_args));
if (result != 0) { /* handle error */ }
```

### Thread Pool API

```c
ctools_thread_pool pool;

// Initialize
ctools_pool_init(&pool, num_threads, args, sizeof(args[0]));

// Launch and wait
ctools_pool_launch(&pool, my_worker);
int result = ctools_pool_join(&pool);

// Or combined
int result = ctools_pool_run(&pool, my_worker);

// Cleanup
ctools_pool_free(&pool);
```

### Convenience Macro

```c
int result;
CTOOLS_PARALLEL_FOR(my_worker, args, 8, sizeof(args[0]), result);
```

---

## Threading Model (OpenMP vs Persistent Pool)

ctools uses a mixed threading model. OpenMP is preferred for tight, regular loops or divide-and-conquer tasks. The persistent pool is preferred for batches of heterogeneous work items (e.g., per-variable I/O) where thread creation overhead would dominate.

### OpenMP usage (parallel loops/tasks)

- Sorting algorithms: `ctools_sort_merge.c`, `ctools_sort_counting.c`, `ctools_sort_ips4o.c`, `ctools_sort_sample.c`, `ctools_sort_timsort.c`, `ctools_sort_radix_msd.c` (recursive tasks)
- Merge: `cmerge/cmerge_impl.c` (parallel loop for merge steps)
- Regression/HDFE: `creghdfe/*.c`, `civreghdfe/*.c` (parallel loops and SIMD for linear algebra/partialling)
- Quantile regression: `cqreg/*.c` (parallel loops, tasking, and SIMD)

### Global persistent pool usage (ctools_persistent_pool)

- Stata I/O: `ctools_data_io.c` (load/store variables in parallel)
- Merge helpers: `cmerge/cmerge_impl.c` (histogram/scatter for radix-based merge paths)
- Apply-permutation stages: `ctools_sort_radix_lsd.c`, `ctools_sort_radix_msd.c`, `ctools_sort_timsort.c`

### When to pick which

- **OpenMP**: tight, uniform loops where per-iteration work is similar and data access is contiguous.
- **Persistent pool**: many small/medium independent tasks, per-variable work, or repeated batches across a command where thread reuse saves overhead.
- Avoid nested parallelism (e.g., OpenMP inside code already using the persistent pool) to reduce oversubscription.

### Potential swap opportunities

- `ctools_sort_radix_msd.c` top-level histogram/scatter currently create pthreads per call; this could move to the global persistent pool to reduce thread creation overhead when sorting is called repeatedly in a session.
- If OpenMP regions are inside phases that already use the persistent pool (e.g., sort + apply-permutation + I/O), consider serializing one layer or lowering OpenMP threads to avoid oversubscribing cores.

### OpenMP tuning tips

- `OMP_NUM_THREADS`: set to physical cores for compute-bound code; keep below `CTOOLS_IO_MAX_THREADS + OMP_NUM_THREADS` to avoid oversubscription.
- `OMP_PROC_BIND=spread` and `OMP_PLACES=cores`: keep threads pinned and reduce cache thrash.
- `OMP_DYNAMIC=FALSE`: avoid runtime thread count changes that interfere with sizing heuristics.
- `OMP_WAIT_POLICY=PASSIVE` (lower CPU use when waiting) or `ACTIVE` (lower latency for short tasks).
- `GOMP_CPU_AFFINITY` (GNU OpenMP) or `KMP_AFFINITY` (Intel/OpenMP) for explicit pinning when needed.

---

## Timing and Profiling

Include: `#include "ctools_timer.h"`

### Initialize Once

```c
ctools_timer_init();  // Call at plugin start
```

### Timing a Code Block

```c
double elapsed;
CTOOLS_TIME_BLOCK(elapsed) {
    // code to time
}
printf("Elapsed: %.3f ms\n", elapsed);
```

### Phase Timing

```c
CTOOLS_TIMER_START(load);
// ... loading code ...
CTOOLS_TIMER_END(load);
printf("Load took %.3f sec\n", _timer_load_elapsed);
```

### Store to Variable

```c
double t_load, t_sort, t_store;

CTOOLS_TIMER_BEGIN(load);
ctools_data_load(&data, nvars);
CTOOLS_TIMER_STORE(load, t_load);

CTOOLS_TIMER_BEGIN(sort);
ctools_sort_radix_lsd(&data, sort_vars, nsort);
CTOOLS_TIMER_STORE(sort, t_sort);
```

### Raw Functions

```c
double start = ctools_timer_seconds();
// ... work ...
double elapsed = ctools_timer_seconds() - start;
```

---

## Argument Parsing

Include: `#include "ctools_args.h"`

Arguments are passed from Stata as a space-separated string.

### Check for Flags

```c
if (ctools_args_has_flag(args, "verbose")) {
    verbose = 1;
}

if (ctools_args_has_flag(args, "stable")) {
    use_stable_sort = 1;
}
```

### Get Option Values

```c
// Get integer option with default
int nthreads = ctools_args_get_int(args, "threads=", 4);

// Find option string
const char *alg = ctools_args_find_option(args, "alg=");
if (alg && strncmp(alg, "timsort", 7) == 0) {
    algorithm = SORT_ALG_TIMSORT;
}
```

### Parse Integer Arrays

```c
// Parse "1 2 3 alg=lsd" -> vars = {1, 2, 3}
int *vars;
size_t nvars;
if (ctools_args_parse_int_array(args, "alg=", &vars, &nvars) == 0) {
    // use vars[0..nvars-1]
    free(vars);
}
```

---

## Error Handling

Include: `#include "ctools_error.h"`

### Display Messages

```c
// Informational (goes to SF_display)
ctools_msg("csort", "Loaded %zu observations", nobs);

// Error (goes to SF_error)
ctools_error("csort", "Failed to open file: %s", filename);

// Verbose (only if flag set)
ctools_verbose("csort", verbose, "Sorted in %.1f ms", elapsed);
```

### Allocation Errors

```c
ctools_error_alloc("csort");
// Output: "csort: memory allocation failed"

ctools_error_alloc_ctx("csort", "thread buffers");
// Output: "csort: failed to allocate thread buffers"
```

### Check Macros

```c
// Check and return on failure
double *buf = malloc(n * sizeof(double));
CTOOLS_CHECK_ALLOC(buf, "csort", STATA_ERR_MEMORY);

// Check with cleanup
double *buf = malloc(n * sizeof(double));
CTOOLS_CHECK_ALLOC_CLEANUP(buf, "csort", { free(other); }, STATA_ERR_MEMORY);
```

---

## Memory Management

Include: `#include "ctools_config.h"`

### Aligned Allocation

Always use aligned allocation for data arrays (enables SIMD, better cache behavior):

```c
// Allocate (64-byte aligned)
double *data = ctools_aligned_alloc(CACHE_LINE_SIZE, n * sizeof(double));
if (!data) { /* handle error */ }

// Free (required - regular free() crashes on Windows)
ctools_aligned_free(data);  // Safe with NULL
```

### Prefetching

```c
// Prefetch for reading
CTOOLS_PREFETCH(&data[i + PREFETCH_DISTANCE]);

// Prefetch for writing
CTOOLS_PREFETCH_W(&output[i + PREFETCH_DISTANCE]);
```

### Branch Hints

```c
if (CTOOLS_LIKELY(ptr != NULL)) {
    // common path
}

if (CTOOLS_UNLIKELY(error)) {
    // rare path
}
```

### Configuration Constants

| Constant | Default | Purpose |
|----------|---------|---------|
| `MIN_PARALLEL_BYTES` | 1 MB | Minimum data size for parallel I/O |
| `MIN_OBS_PER_THREAD` | 100,000 | Minimum obs per thread for parallel sort |
| `NUM_THREADS` | 8 | Maximum threads for parallel operations |
| `CACHE_LINE_SIZE` | 64 | Alignment for allocations |
| `RADIX_BITS` | 8 | Radix sort bucket size (256 buckets) |
| `PREFETCH_DISTANCE` | 16 | Elements ahead to prefetch |

---

## Adding a New Command

### 1. Create Implementation Files

```
src/newcmd/
  newcmd_impl.c
  newcmd_impl.h
```

### 2. Implement Main Entry Point

```c
// newcmd_impl.h
#ifndef NEWCMD_IMPL_H
#define NEWCMD_IMPL_H
#include "stplugin.h"
ST_retcode newcmd_main(const char *args);
#endif

// newcmd_impl.c
#include "newcmd_impl.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
#include "ctools_error.h"
#include "ctools_args.h"

ST_retcode newcmd_main(const char *args)
{
    ctools_timer_init();
    CTOOLS_TIMER_START(total);

    int verbose = ctools_args_has_flag(args, "verbose");

    // Parse arguments
    int *var_indices;
    size_t nvars;
    if (ctools_args_parse_int_array(args, NULL, &var_indices, &nvars) != 0) {
        ctools_error("newcmd", "Failed to parse arguments");
        return 198;
    }

    // Load data
    stata_data data;
    stata_data_init(&data);

    stata_retcode rc = ctools_data_load_selective(&data, var_indices, nvars, 0, 0);
    free(var_indices);
    if (rc != STATA_OK) {
        ctools_error("newcmd", "Failed to load data");
        return 920;
    }

    // Process...

    // Store results
    rc = ctools_data_store(&data, SF_in1());
    stata_data_free(&data);

    if (verbose) {
        CTOOLS_TIMER_END(total);
        ctools_verbose("newcmd", verbose, "Total time: %.3f sec",
                       _timer_total_elapsed);
    }

    return rc == STATA_OK ? 0 : 920;
}
```

### 3. Add Dispatch Case

In `src/ctools_plugin.c`:

```c
#include "newcmd/newcmd_impl.h"

// In stata_call():
if (strncmp(cmd, "newcmd", 6) == 0) {
    return newcmd_main(args);
}
```

### 4. Update Makefile

```makefile
NEWCMD_SRCS = src/newcmd/newcmd_impl.c
NEWCMD_HEADERS = src/newcmd/newcmd_impl.h

SRCS += $(NEWCMD_SRCS)
HEADERS += $(NEWCMD_HEADERS)
```

### 5. Create Stata Files

```
build/newcmd.ado    - Stata wrapper
build/newcmd.sthlp  - Help file
```

---

## Performance Guidelines

### Parallelization

1. **Check data size before parallelizing**:
   ```c
   size_t total_bytes = nobs * nvars * sizeof(double);
   if (total_bytes < MIN_PARALLEL_BYTES) {
       // Use sequential path
   }
   ```

2. **Use OpenMP for simple loops**:
   ```c
   #pragma omp parallel for num_threads(NUM_THREADS)
   for (size_t i = 0; i < n; i++) {
       output[i] = process(input[i]);
   }
   ```

3. **Use thread pool for complex work**:
   ```c
   ctools_parallel_run(worker_func, args, num_threads, sizeof(args[0]));
   ```

### Memory Access

1. **Process data in cache-friendly order** (sequential access):
   ```c
   // Good: sequential access
   for (size_t i = 0; i < n; i++) {
       sum += data[i];
   }

   // Bad: strided access
   for (size_t i = 0; i < n; i++) {
       sum += data[i * stride];
   }
   ```

2. **Prefetch for indirect access**:
   ```c
   for (size_t i = 0; i < n; i++) {
       CTOOLS_PREFETCH(&data[order[i + PREFETCH_DISTANCE]]);
       output[i] = data[order[i]];
   }
   ```

3. **Use aligned allocations**:
   ```c
   double *buf = ctools_aligned_alloc(CACHE_LINE_SIZE, n * sizeof(double));
   ```

### Avoid Common Pitfalls

1. **Don't mix `malloc`/`ctools_aligned_free`** - causes crashes on Windows
2. **Check allocations** - use `CTOOLS_CHECK_ALLOC` macro
3. **Free in reverse order** - prevents dangling pointer issues
4. **Use `verbose` option** - always support timing output for debugging

---

## Quick Reference

### Common Includes

```c
#include "stplugin.h"        // Stata plugin interface
#include "ctools_types.h"    // Data structures, I/O, sorting
#include "ctools_config.h"   // Aligned alloc, prefetch, constants
#include "ctools_threads.h"  // Thread pool
#include "ctools_timer.h"    // Timing utilities
#include "ctools_args.h"     // Argument parsing
#include "ctools_error.h"    // Error reporting
```

### Typical Command Structure

```c
ST_retcode mycmd_main(const char *args)
{
    ctools_timer_init();
    int verbose = ctools_args_has_flag(args, "verbose");

    // 1. Parse arguments
    // 2. Load data: ctools_data_load()
    // 3. Process data
    // 4. Store results: ctools_data_store()
    // 5. Cleanup: stata_data_free()
    // 6. Report timing if verbose

    return 0;
}
```
