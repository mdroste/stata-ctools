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
- [Arena Allocators](#arena-allocators)
- [Hash Tables](#hash-tables)
- [OLS and Regression Infrastructure](#ols-and-regression-infrastructure)
- [Adding a New Command](#adding-a-new-command)
- [Performance Guidelines](#performance-guidelines)

---

## Architecture Overview

### Plugin Lifecycle

Every ctools command follows this pattern:

```
Stata -> C Plugin -> Stata
  1. Load data from Stata into C memory
  2. Process (sort, merge, regress, etc.)
  3. Store results back to Stata
```

On each invocation the dispatcher in `ctools_plugin.c`:
1. Parses the command name and optional `threads()` setting
2. Calls `ctools_cleanup_stale_state()` to free leftover state from previously interrupted commands (while preserving state for multi-phase operations like cmerge)
3. Dispatches to the appropriate `*_main()` function
4. Destroys the global thread pool on exit

### Memory Model

- **Column-major storage**: Each variable is a contiguous array
- **Numeric variables**: `double[]` (8 bytes per observation)
- **String variables**: `char*[]` (pointer array, heap-allocated strings)
- **Aligned allocations**: 64-byte cache line alignment for SIMD/prefetch efficiency

### Key Files

| File | Purpose |
|------|---------|
| `src/ctools_plugin.c` | Main dispatcher - routes commands to implementations |
| `src/ctools_types.h/c` | Core data structures, data I/O, sorting declarations |
| `src/ctools_config.h` | Aligned alloc, prefetch, thread management, constants |
| `src/ctools_runtime.h/c` | Error reporting, timing, lifecycle cleanup |
| `src/ctools_parse.h/c` | Argument parsing utilities |
| `src/ctools_threads.h/c` | Persistent thread pool |
| `src/ctools_arena.h/c` | Arena allocators (growing + string) |
| `src/ctools_hash.h/c` | Hash tables and label utilities |
| `src/ctools_ols.h/c` | Cholesky, precision arithmetic, dot products, collinearity |
| `src/ctools_matrix.h/c` | Matrix multiplication and sandwich VCE |
| `src/ctools_hdfe_utils.h/c` | Singleton detection, FE remapping, union-find, DOF |
| `src/ctools_spi.h` | Stata Plugin Interface convenience wrappers |
| `src/ctools_simd.h` | SIMD intrinsic utilities |
| `src/ctools_select.h` | Selection algorithms (nth element) |
| `src/ctools_sort_pairs.h` | Pair sorting utilities |
| `src/ctools_unroll.h` | Loop unrolling macros |
| `src/ctools_eisel_lemire.h` | Fast float parsing (Eisel-Lemire algorithm) |
| `src/stplugin.c/h` | Stata Plugin Interface (DO NOT MODIFY) |

### Dispatched Commands

The following commands are registered in `ctools_plugin.c`:

| Dispatch name | Handler | Source directory |
|---------------|---------|------------------|
| `csort` | `csort_main` | `src/csort/` |
| `creghdfe` | `creghdfe_main` | `src/creghdfe/` |
| `cimport` | `cimport_main` | `src/cimport/` |
| `cexport` | `cexport_main` | `src/cexport/` |
| `cexport_xlsx` | `cexport_xlsx_main` | `src/cexport/` |
| `cmerge` | `cmerge_main` | `src/cmerge/` |
| `cqreg` | `cqreg_main` | `src/cqreg/` |
| `cbinscatter` | `cbinscatter_main` | `src/cbinscatter/` |
| `civreghdfe` | `civreghdfe_main` | `src/civreghdfe/` |
| `cencode` | `cencode_main` | `src/cencode/` |
| `cwinsor` | `cwinsor_main` | `src/cwinsor/` |
| `cdestring` | `cdestring_main` | `src/cdestring/` |
| `cdecode` | `cdecode_main` | `src/cdecode/` |
| `cdecode_scan` | `cdecode_scan_main` | `src/cdecode/` |
| `csample` | `csample_main` | `src/csample/` |
| `cbsample` | `cbsample_main` | `src/cbsample/` |
| `crangestat` | `crangestat_main` | `src/crangestat/` |
| `cpsmatch` | `cpsmatch_main` | `src/cpsmatch/` |

---

## Core Data Structures

Include: `#include "ctools_types.h"`

### perm_idx_t

Permutation index type. Uses `uint32_t` instead of `size_t` for 50% memory savings on 64-bit systems. Stata's SPI limits observations to < 2^31.

```c
typedef uint32_t perm_idx_t;
#define PERM_IDX_MAX UINT32_MAX
```

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
    size_t nvars;            // Number of variables (columns)
    stata_variable *vars;    // Array of all variables [nvars]
    perm_idx_t *sort_order;  // Permutation array (0-based) [nobs] - uint32_t
} stata_data;
```

### ctools_filtered_data

Dataset loaded with if/in filtering applied at load time:

```c
typedef struct {
    stata_data data;         // Standard data (N = N_filtered)
    perm_idx_t *obs_map;     // obs_map[i] = 1-based Stata obs for filtered index i
    size_t n_range;          // Original range size (before filtering)
    int was_filtered;        // 1 if any filtering occurred, 0 if identity
} ctools_filtered_data;
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

The primary loading function is `ctools_data_load()`, which handles selective variable loading and if/in filtering in a single call:

```c
ctools_filtered_data fd;
ctools_filtered_data_init(&fd);

// Load specific variables (1-based indices) with if/in filtering
int var_indices[] = {1, 3, 5};
stata_retcode rc = ctools_data_load(&fd, var_indices, 3, 0, 0, CTOOLS_LOAD_CHECK_IF);

// Or load ALL variables (pass NULL for var_indices)
rc = ctools_data_load(&fd, NULL, 0, 0, 0, CTOOLS_LOAD_CHECK_IF);

// Skip if/in filtering (load everything in range)
rc = ctools_data_load(&fd, var_indices, 3, 0, 0, CTOOLS_LOAD_SKIP_IF);

// Access data: fd.data.vars[k].data.dbl[i]
// Observation count: fd.data.nobs (= N_filtered)
```

Parameters:
- `var_indices`: Array of 1-based Stata variable indices, or `NULL` to load all
- `nvars`: Number of variables (ignored if `var_indices` is NULL)
- `obs_start`, `obs_end`: 1-based range (0 = use `SF_in1()`/`SF_in2()`)
- `flags`: `CTOOLS_LOAD_CHECK_IF` (default) or `CTOOLS_LOAD_SKIP_IF`

Loading performs two passes:
1. Count observations passing `SF_ifobs()` and build `obs_map`
2. Load only filtered observations (memory = O(N_filtered x K))

### Storing Data to Stata

```c
// Write all variables back (obs1 = first observation, 1-based)
stata_retcode rc = ctools_data_store(&fd.data, obs1);

// Write specific variables to specific Stata variable indices
rc = ctools_data_store_selective(&fd.data, var_indices, nvars, obs1);

// Write filtered results back using obs_map
rc = ctools_store_filtered(values, N, var_idx, fd.obs_map);

// Write with permutation mapping (for merge operations)
rc = ctools_stream_var_permuted(var_idx, source_rows, output_nobs, obs1);
```

### Cleanup

```c
ctools_filtered_data_free(&fd);  // Safe to call multiple times

// Or for raw stata_data:
stata_data_free(&data);
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
// Sort by variables 1 and 2 (1-based indices), applies permutation immediately
int sort_vars[] = {1, 2};
stata_retcode rc = ctools_sort_radix_lsd(&data, sort_vars, 2);

// Sort + get permutation mapping (sorted_idx -> original_idx)
size_t *perm = malloc(data.nobs * sizeof(size_t));
rc = ctools_sort_radix_lsd_with_perm(&data, sort_vars, 2, perm);

// Or use the unified dispatcher (order-only: computes sort_order, does NOT apply)
rc = ctools_sort_dispatch(&data, sort_vars, 2, SORT_ALG_LSD);
// ... then apply separately:
rc = ctools_apply_permutation(&data);
```

### Algorithm Enum

```c
typedef enum {
    SORT_ALG_LSD = 0,      // LSD radix sort
    SORT_ALG_MSD = 1,      // MSD radix sort
    SORT_ALG_TIMSORT = 2,  // Timsort
    SORT_ALG_SAMPLE = 3,   // Sample sort
    SORT_ALG_COUNTING = 4, // Counting sort
    SORT_ALG_MERGE = 5,    // Parallel merge sort
    SORT_ALG_IPS4O = 6,    // IPS4o
    SORT_ALG_AUTO = 7      // Auto-select based on data
} sort_algorithm_t;
```

### Order-Only Variants

Every sort algorithm has an `_order_only` variant that computes `data->sort_order` without applying the permutation. Use `ctools_apply_permutation()` to apply afterward:

```c
rc = ctools_sort_radix_lsd_order_only(&data, sort_vars, 2);
// data->sort_order is set, data is unchanged
rc = ctools_apply_permutation(&data);
// data is now reordered
```

### Choosing an Algorithm

```c
// Check if counting sort is suitable
if (ctools_counting_sort_suitable(&data, var_idx)) {
    ctools_sort_counting(&data, sort_vars, nsort);
} else {
    ctools_sort_radix_lsd(&data, sort_vars, nsort);
}

// Or let the dispatcher handle it (SORT_ALG_COUNTING falls back to LSD if unsuitable)
ctools_sort_dispatch(&data, sort_vars, nsort, SORT_ALG_COUNTING);
```

---

## Thread Pool

Include: `#include "ctools_threads.h"`

ctools uses a **persistent thread pool** that keeps workers alive between tasks, eliminating thread creation overhead for repeated parallel operations.

### Global Pool (Most Common)

```c
// Get the singleton pool (lazily initialized)
ctools_persistent_pool *pool = ctools_get_global_pool();
if (!pool) { /* handle error */ }

// Submit a single work item
ctools_persistent_pool_submit(pool, my_worker, &my_arg);

// Submit a batch of work items (like parallel-for)
my_args args[8];
// ... fill in args ...
ctools_persistent_pool_submit_batch(pool, my_worker, args, 8, sizeof(args[0]));

// Wait for all submitted work to complete
int result = ctools_persistent_pool_wait(pool);
if (result != 0) { /* a worker failed */ }
```

### Manual Pool Management

```c
ctools_persistent_pool pool;

// Initialize with specific worker count
ctools_persistent_pool_init(&pool, num_threads);

// Submit and wait (same API as global pool)
ctools_persistent_pool_submit_batch(&pool, my_worker, args, 8, sizeof(args[0]));
int result = ctools_persistent_pool_wait(&pool);

// Destroy when done
ctools_persistent_pool_destroy(&pool);
```

### Thread Function Signature

```c
typedef void *(*ctools_thread_func)(void *arg);
// Return NULL for success, non-NULL for failure
```

### Thread Count Management

```c
int n = ctools_get_max_threads();      // Current thread limit
ctools_set_max_threads(4);             // Override
ctools_reset_max_threads();            // Reset to default (omp_get_max_threads)
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

### OpenMP tuning tips

- `OMP_NUM_THREADS`: set to physical cores for compute-bound code; keep below `CTOOLS_IO_MAX_THREADS + OMP_NUM_THREADS` to avoid oversubscription.
- `OMP_PROC_BIND=spread` and `OMP_PLACES=cores`: keep threads pinned and reduce cache thrash.
- `OMP_DYNAMIC=FALSE`: avoid runtime thread count changes that interfere with sizing heuristics.
- `OMP_WAIT_POLICY=PASSIVE` (lower CPU use when waiting) or `ACTIVE` (lower latency for short tasks).
- `GOMP_CPU_AFFINITY` (GNU OpenMP) or `KMP_AFFINITY` (Intel/OpenMP) for explicit pinning when needed.

---

## Timing and Profiling

Include: `#include "ctools_runtime.h"`

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
ctools_data_load(&fd, NULL, 0, 0, 0, 0);
CTOOLS_TIMER_STORE(load, t_load);

CTOOLS_TIMER_BEGIN(sort);
ctools_sort_radix_lsd(&fd.data, sort_vars, nsort);
CTOOLS_TIMER_STORE(sort, t_sort);
```

### Raw Functions

```c
double start = ctools_timer_seconds();
// ... work ...
double elapsed = ctools_timer_seconds() - start;

double ms = ctools_timer_ms();  // Millisecond variant
```

---

## Argument Parsing

Include: `#include "ctools_parse.h"`

Arguments are passed from Stata as a space-separated string.

### Check for Flags

```c
if (ctools_parse_bool_option(args, "verbose")) {
    verbose = 1;
}

if (ctools_parse_bool_option(args, "stable")) {
    use_stable_sort = 1;
}
```

### Get Option Values

```c
// Get integer option ("key=value" format)
int nthreads = 4;  // default
ctools_parse_int_option(args, "threads", &nthreads);

// Get double option with default
double pctl = ctools_parse_double_option(args, "p", 1.0);

// Get size_t option with default
size_t chunk = ctools_parse_size_option(args, "chunk", 10000);

// Get string option
char name[64];
if (ctools_parse_string_option(args, "name", name, sizeof(name))) {
    // name = value from "name=value"
}
```

### Parse Integer Arrays

```c
// Parse leading integers from cursor position
const char *cursor = args;  // e.g., "1 2 3 verbose"
int arr[3];
if (ctools_parse_int_array(arr, 3, &cursor) == 0) {
    // arr = {1, 2, 3}, cursor points past them
}

// Parse one value at a time
size_t val;
ctools_parse_next_size(&cursor, &val);

int ival;
ctools_parse_next_int(&cursor, &ival);
```

---

## Error Handling

Include: `#include "ctools_runtime.h"`

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

### Safe Allocation (Overflow-Checked)

```c
// Two-factor: ptr = malloc(a * b) with overflow check
double *buf = ctools_safe_malloc2(n, sizeof(double));

// Three-factor: ptr = malloc(a * b * c) with overflow check
double *mat = ctools_safe_malloc3(N, K, sizeof(double));

// Zero-initialized variants
double *zbuf = ctools_safe_calloc2(n, sizeof(double));

// Aligned variant
double *abuf = ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, n, sizeof(double));
```

All return `NULL` on overflow or allocation failure.

### Prefetching

```c
// Prefetch for reading
CTOOLS_PREFETCH(&data[i + 16]);

// Prefetch for writing
CTOOLS_PREFETCH_W(&output[i + 16]);
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
| `MIN_OBS_PER_THREAD` | 100,000 | Minimum obs per thread for parallel sort |
| `CACHE_LINE_SIZE` | 64 | Alignment for allocations |
| `RADIX_BITS` | 8 | Radix sort bucket size (256 buckets) |
| `CTOOLS_IO_BUFFER_SIZE` | 64 KB | File I/O buffer size |
| `CTOOLS_IMPORT_CHUNK_SIZE` | 8 MB | Bytes per CSV parsing chunk |
| `CTOOLS_EXPORT_CHUNK_SIZE` | 10,000 | Rows per CSV formatting chunk |
| `CTOOLS_ARENA_BLOCK_SIZE` | 1 MB | Arena allocator block size |
| `CTOOLS_MAX_VARNAME_LEN` | 32 | Maximum Stata variable name length |
| `CTOOLS_MAX_STRING_LEN` | 2045 | Maximum Stata string length |
| `CTOOLS_MAX_COLUMNS` | 32,767 | Maximum Stata columns |

### Thread Count

Thread limits are determined at runtime via `omp_get_max_threads()` and can be overridden per-command with `threads()`:

```c
int n = ctools_get_max_threads();
ctools_set_max_threads(4);      // Override
ctools_reset_max_threads();     // Back to default
```

---

## Arena Allocators

Include: `#include "ctools_arena.h"`

Two arena types for different allocation patterns.

### Growing Arena (unknown total size)

A chain of fixed-size blocks. When the current block is exhausted, a new one is allocated. All blocks freed together.

```c
ctools_arena arena;
ctools_arena_init(&arena, CTOOLS_ARENA_DEFAULT_BLOCK_SIZE);  // 1MB blocks

// Allocate (8-byte aligned)
void *ptr = ctools_arena_alloc(&arena, 256);

// Aligned allocation (e.g., cache-line)
void *aligned = ctools_arena_alloc_aligned(&arena, 1024, 64);

// String duplication
char *s = ctools_arena_strdup(&arena, "hello");

// Reset (keep blocks, mark as empty - for reuse)
ctools_arena_reset(&arena);

// Free everything
ctools_arena_free(&arena);
```

### String Arena (known/estimated total size)

A single contiguous block optimized for string pooling. Three fallback modes when full:

```c
// Create with estimated capacity and fallback mode
ctools_string_arena *sa = ctools_string_arena_create(
    1024 * 1024,                           // 1MB capacity
    CTOOLS_STRING_ARENA_STRDUP_FALLBACK    // fall back to strdup if full
);

// Duplicate strings into the arena
char *s1 = ctools_string_arena_strdup(sa, "value1");
char *s2 = ctools_string_arena_strdup(sa, "value2");

// Check ownership (useful for selective freeing)
if (ctools_string_arena_owns(sa, s1)) { /* arena-owned */ }

// Free (strings from STRDUP_FALLBACK must be freed separately)
ctools_string_arena_free(sa);
```

Fallback modes:
- `CTOOLS_STRING_ARENA_NO_FALLBACK` - return NULL when full
- `CTOOLS_STRING_ARENA_STRDUP_FALLBACK` - fall back to `strdup()` (caller must track and free)
- `CTOOLS_STRING_ARENA_STATIC_FALLBACK` - return a static empty string (never NULL)

---

## Hash Tables

Include: `#include "ctools_hash.h"`

Open-addressing hash tables with linear probing, automatic resizing at 75% load factor. Used primarily by `cencode` and `cdecode`.

### String -> Integer (for cencode)

```c
ctools_str_hash_table ht;
ctools_str_hash_init(&ht, CTOOLS_HASH_INIT_SIZE);

// Insert with auto-assigned value (1-based)
uint32_t hash = ctools_str_hash_compute("label_text");
int code = ctools_str_hash_insert(&ht, "label_text", hash);  // returns assigned code

// Insert with specific value
int code2 = ctools_str_hash_insert_value(&ht, "other", 42);

// Lookup (returns 0 if not found)
int val = ctools_str_hash_lookup(&ht, "label_text");

ctools_str_hash_free(&ht);
```

### Integer -> String (for cdecode)

```c
ctools_int_hash_table ht;
ctools_int_hash_init(&ht, CTOOLS_HASH_INIT_SIZE);

// Insert
ctools_int_hash_insert(&ht, 1, "Male");
ctools_int_hash_insert(&ht, 2, "Female");

// Lookup (returns NULL if not found)
const char *label = ctools_int_hash_lookup(&ht, 1);

ctools_int_hash_free(&ht);
```

### Label Utilities

Shared label escape/unescape/parse/serialize functions for `cdecode` and `cencode`:

```c
// Escape/unescape labels for serialization
char dst[2048];
ctools_label_unescape("escaped\\|text", dst, sizeof(dst));
size_t len = ctools_label_escape("raw|text", dst, sizeof(dst));

// Parse label files directly
ctools_int_hash_table ht;
ctools_int_hash_init(&ht, 1024);
int max_len;
ctools_label_parse_stata_file_int("labels.do", &ht, &max_len);

// Write labels as Stata .do file
ctools_label_write_stata_file(strings, codes, n_labels, "myvar", "output.do");
```

---

## OLS and Regression Infrastructure

### Linear Algebra (`ctools_ols.h`)

```c
// Cholesky decomposition: A = L * L' (in-place, returns L in lower triangle)
ST_int ctools_cholesky(ST_double *A, ST_int n);

// Matrix inversion via Cholesky
ST_int ctools_invert_from_cholesky(const ST_double *L, ST_int n, ST_double *inv);

// Solve Ax = b
ST_int ctools_solve_cholesky(const ST_double *A, const ST_double *b, ST_int n, ST_double *x);

// Quad-precision dot products (match Stata's quadcross)
ST_double kahan_dot(const ST_double *x, const ST_double *y, ST_int N);
ST_double fast_dot(const ST_double *x, const ST_double *y, ST_int N);

// Quad-precision sum of squares
ST_double dd_sum_sq(const ST_double *x, ST_int N);
ST_double dd_sum_sq_weighted(const ST_double *x, const ST_double *w, ST_int N);

// Compute X'X and X'y (Kahan-compensated)
void compute_xtx_xty(const ST_double *data, ST_int N, ST_int K,
                      ST_double *xtx, ST_double *xty);
void compute_xtx_xty_weighted(const ST_double *data, const ST_double *weights,
                               ST_int weight_type, ST_int N, ST_int K,
                               ST_double *xtx, ST_double *xty);

// Collinearity detection via modified Cholesky
ST_int detect_collinearity(const ST_double *xx, ST_int K,
                           ST_int *is_collinear, ST_int verbose);
```

### Matrix Operations (`ctools_matrix.h`)

```c
// C = A' * B  (A: N x K1, B: N x K2, C: K1 x K2)
void ctools_matmul_atb(const ST_double *A, const ST_double *B,
                       ST_int N, ST_int K1, ST_int K2, ST_double *C);

// C = A * B  (A: K1 x K2, B: K2 x K3, C: K1 x K3)
void ctools_matmul_ab(const ST_double *A, const ST_double *B,
                      ST_int K1, ST_int K2, ST_int K3, ST_double *C);

// C = A' * diag(w) * B  (weighted)
void ctools_matmul_atdb(const ST_double *A, const ST_double *B, const ST_double *w,
                        ST_int N, ST_int K1, ST_int K2, ST_double *C);
```

### Sandwich VCE (`ctools_matrix.h`)

```c
// Data bundle for VCE computation
typedef struct {
    const ST_double *X_eff;   // Effective regressors (X for OLS, P_Z*X for IV)
    const ST_double *D;       // Bread: (X'X)^-1 or (XkX)^-1
    const ST_double *resid;   // Residuals
    const ST_double *weights; // NULL if unweighted
    ST_int weight_type;       // 0=none, 1=aweight, 2=fweight, 3=pweight
    ST_int N, K;
    ST_int normalize_weights; // 1=normalize aw/pw (OLS), 0=raw (IV)
} ctools_vce_data;

void ctools_vce_robust(const ctools_vce_data *d, ST_double dof_adj, ST_double *V);
void ctools_vce_cluster(const ctools_vce_data *d, const ST_int *cluster_ids,
                        ST_int num_clusters, ST_double dof_adj, ST_double *V);
```

### HDFE Utilities (`ctools_hdfe_utils.h`)

```c
// FE factor and global HDFE state
typedef struct { ST_int num_levels, max_level; ST_int *levels; ... } FE_Factor;
typedef struct { ST_int G, N, K; FE_Factor *factors; ... } HDFE_State;

// Singleton detection (iterative, handles cascading)
ST_int ctools_remove_singletons(ST_int **fe_levels, ST_int G, ST_int N,
                                 ST_int *mask, ST_int max_iter, ST_int verbose);

// Cluster/FE remapping to contiguous indices
ST_int ctools_remap_cluster_ids(ST_int *cluster_ids, ST_int N, ST_int *num_clusters);

// Array compaction (remove flagged observations)
ST_int ctools_compact_array_double(const ST_double *src, ST_double *dest,
                                    const ST_int *mask, ST_int N_src, ST_int N_dest);
ST_int ctools_compact_matrix_double(const ST_double *src, ST_double *dest,
                                     const ST_int *mask, ST_int N_src, ST_int N_dest, ST_int K);

// Connected components (for mobility groups / DOF)
ST_int ctools_count_connected_components(const ST_int *fe1_levels, const ST_int *fe2_levels,
                                          ST_int N, ST_int num_levels1, ST_int num_levels2);

// FE nesting check
ST_int ctools_fe_nested_in_cluster(const ST_int *fe_levels, ST_int num_fe_levels,
                                    const ST_int *cluster_ids, ST_int N);

// DOF calculation
ST_int ctools_compute_hdfe_dof(const FE_Factor *factors, ST_int G, ST_int N,
                                ST_int *df_a, ST_int *mobility_groups);

// Cache-friendly FE projection via sorted permutation
ST_int ctools_build_sorted_permutation(FE_Factor *f, ST_int N);

// Allocate CG solver buffers
ST_int ctools_hdfe_alloc_buffers(HDFE_State *state, ST_int alloc_proj, ST_int max_columns);

// Cleanup
void ctools_hdfe_state_cleanup(HDFE_State *state);
```

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
#include "ctools_runtime.h"
#include "ctools_parse.h"

ST_retcode newcmd_main(const char *args)
{
    ctools_timer_init();
    CTOOLS_TIMER_START(total);

    int verbose = ctools_parse_bool_option(args, "verbose");

    // Parse arguments
    const char *cursor = args;
    int var_indices[10];
    size_t nvars = 3;
    if (ctools_parse_int_array(var_indices, nvars, &cursor) != 0) {
        ctools_error("newcmd", "Failed to parse variable indices");
        return 198;
    }

    // Load data with if/in filtering
    ctools_filtered_data fd;
    ctools_filtered_data_init(&fd);

    stata_retcode rc = ctools_data_load(&fd, var_indices, nvars, 0, 0,
                                         CTOOLS_LOAD_CHECK_IF);
    if (rc != STATA_OK) {
        ctools_error("newcmd", "Failed to load data");
        return 920;
    }

    // Process... (fd.data.vars[k].data.dbl[i])

    // Store results
    rc = ctools_data_store(&fd.data, SF_in1());
    ctools_filtered_data_free(&fd);

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
else if (strcmp(cmd_name, "newcmd") == 0) {
    rc = newcmd_main(cmd_args);
}
```

Also add a `CTOOLS_CMD_NEWCMD` entry to the `ctools_command_t` enum in `ctools_runtime.h` and handle it in `get_command_type()`.

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
   if (nobs < MIN_OBS_PER_THREAD) {
       // Use sequential path
   }
   ```

2. **Use OpenMP for simple loops**:
   ```c
   int nt = ctools_get_max_threads();
   #pragma omp parallel for num_threads(nt)
   for (size_t i = 0; i < n; i++) {
       output[i] = process(input[i]);
   }
   ```

3. **Use persistent pool for variable-level work**:
   ```c
   ctools_persistent_pool *pool = ctools_get_global_pool();
   ctools_persistent_pool_submit_batch(pool, worker_func, args, count, sizeof(args[0]));
   ctools_persistent_pool_wait(pool);
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
       CTOOLS_PREFETCH(&data[order[i + 16]]);
       output[i] = data[order[i]];
   }
   ```

3. **Use aligned allocations**:
   ```c
   double *buf = ctools_aligned_alloc(CACHE_LINE_SIZE, n * sizeof(double));
   ```

4. **Use overflow-safe allocations for user-controlled sizes**:
   ```c
   double *mat = ctools_safe_malloc3(N, K, sizeof(double));
   CTOOLS_CHECK_ALLOC(mat, "mymod", STATA_ERR_MEMORY);
   ```

### Avoid Common Pitfalls

1. **Don't mix `malloc`/`ctools_aligned_free`** - causes crashes on Windows
2. **Check allocations** - use `CTOOLS_CHECK_ALLOC` macro
3. **Free in reverse order** - prevents dangling pointer issues
4. **Use `verbose` option** - always support timing output for debugging
5. **Use `ctools_safe_malloc2/3`** for sizes derived from user data (prevents integer overflow)

---

## Quick Reference

### Common Includes

```c
#include "stplugin.h"        // Stata plugin interface
#include "ctools_types.h"    // Data structures, I/O, sorting
#include "ctools_config.h"   // Aligned alloc, prefetch, constants, thread count
#include "ctools_runtime.h"  // Error reporting, timing, lifecycle cleanup
#include "ctools_parse.h"    // Argument parsing
#include "ctools_threads.h"  // Persistent thread pool
#include "ctools_arena.h"    // Arena allocators
#include "ctools_hash.h"     // Hash tables and label utilities
#include "ctools_ols.h"      // Cholesky, dot products, collinearity
#include "ctools_matrix.h"   // Matrix multiplication, sandwich VCE
#include "ctools_hdfe_utils.h" // FE factors, singletons, DOF
```

### Typical Command Structure

```c
ST_retcode mycmd_main(const char *args)
{
    ctools_timer_init();
    int verbose = ctools_parse_bool_option(args, "verbose");

    // 1. Parse arguments (ctools_parse.h)
    // 2. Load data: ctools_data_load()
    // 3. Process data
    // 4. Store results: ctools_data_store() or ctools_store_filtered()
    // 5. Cleanup: ctools_filtered_data_free()
    // 6. Report timing if verbose

    return 0;
}
```
