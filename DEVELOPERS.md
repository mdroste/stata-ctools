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

Most commands follow this pattern. `cdecode` delegates literal value-label decoding to Stata and commits staged outputs only after every variable succeeds.

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
- **String variables**: `char*[]` pointing into owned flat buffers or arenas, with individually allocated fallback strings
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
| `src/ctools_order.h/c` | Stable lexicographic row-index ordering with distinct extended missing keys |
| `src/ctools_sort_pairs.h` | Pair sorting utilities |
| `src/ctools_unroll.h` | Loop unrolling macros |
| `src/ctools_eisel_lemire.h` | Fast float parsing (Eisel-Lemire algorithm) |
| `src/stplugin.c/h` | Stata Plugin Interface (DO NOT MODIFY) |

### Dispatched Commands

The following commands are registered in `ctools_plugin.c`:

| Dispatch name | Handler | Source directory |
|---------------|---------|------------------|
| `csort` | `csort_main` | `src/csort/` |
| `cipolate` | `cipolate_main` | `src/cipolate/` |
| `csplit` | `csplit_main` | `src/csplit/` |
| `crangejoin` | `crangejoin_main` | `src/crangejoin/` |
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
| `cdecode`, `cdecode_scan` | Retired entry points; request matching ado update | `src/cdecode/` |
| `csample` | `csample_main` | `src/csample/` |
| `cbsample` | `cbsample_main` | `src/cbsample/` |
| `crangestat` | `crangestat_main` | `src/crangestat/` |
| `cpsmatch` | `cpsmatch_main` | `src/cpsmatch/` |
| `cpplmhdfe` | `cpplmhdfe_main` | `src/cpplmhdfe/` |

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
        char **str;          // String: char*[nobs], column-owned storage
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
    perm_idx_t *sort_order;  // Identity/order [nobs], or NULL with NO_SORT_ORDER
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

// Access data: fd.data.vars[k].data.dbl[i]
// Observation count: fd.data.nobs (= N_filtered)
```

Parameters:

- `var_indices`: Array of 1-based plugin-visible Stata variable indices, or `NULL` to load all
- `nvars`: Number of variables (ignored if `var_indices` is NULL)
- `obs_start`, `obs_end`: 1-based range (0 = use `SF_in1()`/`SF_in2()`)
- `flags`: Bitwise OR of the options below

| Flag | Value | Behavior |
|------|-------|----------|
| `CTOOLS_LOAD_CHECK_IF` | `0x00` | Evaluate `SF_ifobs()` in the resolved observation range. |
| `CTOOLS_LOAD_SKIP_IF` | `0x01` | Load every observation in the range; use only when the caller knows filtering is unnecessary. |
| `CTOOLS_LOAD_NO_SORT_ORDER` | `0x02` | Omit the identity permutation and leave `fd.data.sort_order == NULL`. |

Loading first resolves the selection into `obs_map`, then loads only selected
observations. The map always contains the original 1-based Stata observation
indices. Memory for the loaded values is O(N_selected x K).

The default load creates an identity `sort_order`. Callers that only inspect or
transform columns can add `CTOOLS_LOAD_NO_SORT_ORDER`; filtering, column order,
and `obs_map` are unchanged. This saves `sizeof(perm_idx_t) * N_selected` bytes
and their initialization. Timing benefits depend on the data shape.
The flag is used by the regression families, binscatter, matching, encode,
destring, split, and the bulk CSV/XLSX export paths. Sort, sampling, and other
order-dependent paths retain the default.

A caller using this flag must not pass the result to sorting or permutation
routines that require `sort_order`. Ordinary and filtered stores remain valid;
`ctools_data_store_sorted()` rejects a nonempty dataset with no order. If later
processing needs to sort, load with the default flag instead.

`ctools_data_load_ex()` optionally accepts string-width hints. Its array is
indexed by **plugin-visible variable index minus one** and must have
`SF_nvars()` entries, even when the loaded varlist selects or reorders a subset.
Passing `NULL` reads the caller's `_ctools_strw` metadata when available.
Hints select an allocation strategy; they do not relax read bounds. Known
`str2045` columns use packed arenas with a reservation capped at 2 GiB per
column, then retry the original 64-byte-per-row estimate if allocation fails.
Actual strings are stored densely. Virtual/committed memory accounting depends
on the platform; short values are not expanded into 2046-byte slots.

`strL` reads stay on the calling thread. Fixed strings and numeric columns may
use the transport schedulers. Shared numeric writes retain the checked SPI
callback; switching to the unchecked store has produced cross-column corruption
in real Stata. See [the transport measurements](docs/PERFORMANCE_TRANSPORT.md).

Small transfers avoid worker-pool setup using a width-weighted work estimate.
All-numeric loads use 512-row tiles when there are at least 128 columns and
50,000 selected observations. See [the adaptive transport benchmarks](docs/PERFORMANCE_TRANSPORT_ADAPTIVE.md)
for the measured thresholds, allocation changes, and performance tradeoffs.

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

The filtered result owns its column buffers, string arenas/fallback strings,
observation map, and optional sort order. Free the previous result before
loading into the same object again. Do not free individual loaded string
pointers: they may point inside an arena. A caller that takes ownership of a
buffer must clear the corresponding pointer in the result before cleanup.
String-column ownership includes both the pointer array and its arena metadata.

```c
ctools_filtered_data_free(&fd);  // Safe to call multiple times

// Or for raw stata_data:
stata_data_free(&data);
```

### Transport validation

The native transport suites use independent SPI mocks and build both serial
and OpenMP variants:

- `validation/test_transport_native.py`: selection, values, width bounds, and ownership.
- `validation/test_transport_scheduling.py`: variable/row scheduling, permutations, destination maps, and SPI failures.
- `validation/test_transport_adaptive_native.py`: wide numeric gate boundaries, allocation retries, packed-string ownership, optional sort order, and reordered fixed-string/strL/numeric loads with errors.
- `validation/test_transport_store_native.py`: duplicate and invalid destination maps, sparse/OOM fallbacks, checked writes and cancellation, and the sorted-gather size boundary.

Run a suite directly with Python, or run the complete native inventory with
`python3 validation/suite_registry.py --native`. `CTOOLS_TEST_SOURCE` selects an
isolated source directory; `CTOOLS_SANITIZERS=address,undefined` enables ASan plus
UBSan where supported. Linux CI includes the suites and explicitly checks the
new adaptive-load and sorted-store ownership paths with both sanitizers and
leak detection. Mock tests complement the real-Stata command and transport regressions;
they do not establish SPI thread safety or end-to-end performance.

---

## Sorting Algorithms

Include: `#include "ctools_types.h"`

### Available Algorithms

The public `csort` default is `algorithm(auto)`. It selects an engine based on the data; an explicit low-level LSD entry point does not define the public default.

Every engine computes `data->sort_order` without moving the loaded columns (an "order-only" sort). Most callers should go through `ctools_sort_dispatch()`.

| Algorithm | Best For | Order-only function |
|-----------|----------|---------------------|
| LSD Radix | Fixed-width keys | `ctools_sort_radix_lsd_order_only()` |
| MSD Radix | Variable-length strings | `ctools_sort_radix_msd_order_only()` |
| Timsort | Partially sorted data | `ctools_sort_timsort_order_only()` |
| Sample Sort | Large datasets, many cores | `ctools_sort_sample_order_only()` |
| Counting Sort | Integers with small range | `ctools_sort_counting_order_only()` |
| Merge Sort | Stable, predictable O(n log n) | `ctools_sort_merge_order_only()` |
| IPS4o | Memory-efficient parallel | `ctools_sort_ips4o_order_only()` |

`ctools_sort_ips4o_with_perm()` also reorders `data->vars` and returns the permutation; cmerge uses it.

### Usage

```c
// Sort by variables 1 and 2 (1-based indices); computes sort_order only
int sort_vars[] = {1, 2};
stata_retcode rc = ctools_sort_dispatch(&data, sort_vars, 2, SORT_ALG_AUTO);

// Then either write straight to Stata in sorted order...
rc = ctools_data_store_sorted(&data, 1);
// ...or physically reorder the loaded columns first
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

### Choosing an Algorithm

`SORT_ALG_AUTO` uses MSD when any key is a string and counting sort otherwise, falling back to LSD when the data are unsuitable for counting sort. An explicit `SORT_ALG_COUNTING` also falls back to LSD:

```c
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
ctools_sort_dispatch(&fd.data, sort_vars, nsort, SORT_ALG_AUTO);
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

An open-addressing string -> integer hash table with linear probing and automatic resizing at 75% load factor, used by cencode. The public `cdecode` command uses native Stata decoding.

### String -> Integer (for cencode)

```c
ctools_str_hash_table ht;
ctools_str_hash_init(&ht, 1024);

// Insert with auto-assigned value (1-based)
uint32_t hash = ctools_str_hash_compute("label_text");
int code = ctools_str_hash_insert(&ht, "label_text", hash);  // returns assigned code

// Insert with a specific value (0 on success, -1 on allocation failure)
ctools_str_hash_insert_value(&ht, "other", 42);

// Lookup (returns 1 and sets val if found, 0 otherwise)
int val;
if (ctools_str_hash_lookup(&ht, "label_text", &val)) { /* ... */ }

ctools_str_hash_free(&ht);
```

### Value Label Output

```c
// Write a .do file that rebuilds the labels with Mata st_vlmodify()
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

// Singleton detection: peels to the fixed point (cascading chains included);
// with fweights a level is a singleton only when its total fweight is 1
ST_int ctools_remove_singletons(ST_int **fe_levels, ST_int G, ST_int N,
                                 ST_int *mask, const ST_double *fweights, ST_int verbose);

// Cluster/FE remapping to contiguous indices
ST_int ctools_remap_cluster_ids(ST_int *cluster_ids, ST_int N, ST_int *num_clusters);

// Connected components (for mobility groups / DOF)
ST_int ctools_count_connected_components(const ST_int *fe1_levels, const ST_int *fe2_levels,
                                          ST_int N, ST_int num_levels1, ST_int num_levels2);

// FE nesting check
ST_int ctools_fe_nested_in_cluster(const ST_int *fe_levels, ST_int num_fe_levels,
                                    const ST_int *cluster_ids, ST_int N);

// DOF calculation
ST_int ctools_compute_hdfe_dof(const FE_Factor *factors, ST_int G, ST_int N,
                                ST_int *df_a, ST_int *mobility_groups);

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

### 4. Source discovery

The Makefile discovers every C/header file recursively under `src/`, tracks
headers as dependencies, and adds source directories to the include path.
Do not add per-command `*_SRCS`/`*_HEADERS` lists. Third-party libdeflate sources
are compiled separately with their own warning flags.

### 5. Create Stata Files

```
build/newcmd.ado    - Stata wrapper
build/newcmd.sthlp  - Help file
```

Add both files and every prediction/helper ado to `build/ctools.pkg`. Register
the command in the help and compatibility inventories, add its component to
`validation/validate_all.do` and `validation/run_stata_audit.py`, and put the
completion marker at the actual end of the component. Define cleanup for every
owned allocation and any persistent state. `validation/test_release_gate.py`
checks agreement between the public help, source, package, and test inventories.

---

## Performance Guidelines

The [README command benchmarks](docs/BENCHMARKS_README_COMMANDS.md) compare
`cipolate`, `csplit`, and `crangejoin` against native/SSC commands at 1M, 5M, and
20M rows. Use `validation/prepare_readme_benchmarks.py` with a frozen build,
then `validation/summarize_readme_benchmarks.py` to validate complete logs and
retain every trial. The harness uses a monotonic clock, alternates execution
order, and checks every output cell outside the timed region.

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

## Correctness and distribution

Use `scripts/fetch_validation_data.py` once to cache the checksum-pinned official
Stata fixtures, then `validation/run_stata_audit.py` for the complete offline
suite. The driver uses the machine's `stata` shell alias and exits cleanly.
Every component must reach its completion marker. Missing references, unexpected
skips, and failed assertions block publication. ATT-SE comparisons with psmatch2
are separately counted as documented method differences, never as passes.

Numeric variable indices are positions in the plugin varlist. Use shared checked
I/O or `SD_SAFEMODE` for direct SPI calls; raw `SD_FASTMODE` stores can use a
different index space and overwrite inputs. Never change `src/stplugin.c/h`.

`make package` stages a host-platform manifest. The complete checker verifies
exactly its declared platform scope and all helper files. CI release packages
must declare all four platforms and require successful licensed Stata validation.

## Shared failure and wrapper contracts

`_ctools_load.ado` owns platform selection and plugin identity checks. Stata
scopes plugin registration to the calling ado, so each wrapper registers the
helper-selected binary in its own scope. Missing platform files may fall back to
`ctools.plugin`; an incompatible file returns an error. `_ctools_newvars.ado`
validates complete output-name lists before allocation, and `_ctools_weight.ado`
stages expression weights while leaving command-specific sample rules in callers.
`cqreg` uses native `qreg_p`; the unused custom prediction ado is no longer shipped.

`partial_out_columns()` returns status, convergence, and iteration count.
Callers must propagate a failed projection and must not post successful estimates.
The shared robust/cluster covariance APIs also return status: allocation or
cluster-sort failure is 920, invalid covariance configuration is 498. Output
matrices remain unpublished after an error. Legitimate zero residual variance
still produces a successful zero covariance matrix.

OpenMP sort capacity is distinct from pthread/hardware capacity. Logical work
partitions use work-sharing loops so reduced runtime teams still cover every row.
Without OpenMP, sorting uses one usable OpenMP worker. Arena reset retains and
reuses its block chain; aligned and ordinary allocation reject size overflow.

Platform builds compile dependency sources in a fresh private directory, stop
on the first compiler error, link only that directory's explicit object list,
and publish the plugin only after a successful link. Failed builds preserve the
previous plugin and remove temporary objects.

The September 24 regressions are `validation/validate_sep24.do` and
`validation/test_sep24_native.py`. The native runner accepts
`CTOOLS_TEST_OPENMP_PREFIX` for constrained-team tests and
`CTOOLS_SANITIZERS=address,undefined` for ASan plus UBSan. The Stata cases are part
of the full offline release gate. See `SEP24_ASTRA_FIXES.md` for repair evidence.

## Interpolation, splitting, and interval joins

`cipolate` performs one stable index sort and scans groups. `csplit` has a
scan/size phase followed by a write to temporary output variables. `crangejoin`
caches using columns, prepares matches and checked output counts, then writes
the expanded dataset. Both multi-phase commands register cache cleanup with
`ctools_runtime`; their ado wrappers clear caches after errors.

The three command suites are part of `validate_all.do`. The interval-join
reference tests additionally require SSC `rangejoin` 1.1.3 (and `rangestat`) on
the test runner's adopath. Production `crangejoin` has neither dependency.
`python3 validation/test_newcommands_native.py` checks the kernels under UBSan
and fails each allocation in turn; set `CTOOLS_SANITIZERS=address,undefined` to
include ASan. Linux CI runs both sanitizers. `validation/benchmark_newcommands.do`
provides reproducible end-to-end timings and reference-result checks.

The shared stable order helper selects radix passes for large numeric inputs,
parallel merge tiles for string/mixed inputs, and a no-sort path for ordered data.
It preserves extended missing values and stable ties. Run
`python3 validation/test_order_parallel.py` with `CTOOLS_LIBOMP_PREFIX` pointing
to the macOS OpenMP prefix (Linux uses `-fopenmp`).
`python3 validation/test_split_tokens_native.py` checks in-place token storage
with the existing SPI/allocator mock; it also accepts `CTOOLS_SANITIZERS`.
`validation/validate_command_optimizations.do` adds reference checks for the
optimized paths; run it through an error-capturing Stata driver.

For repeated timing and phase profiles, use
`validation/prepare_command_performance.py` and
`validation/summarize_command_performance.py`. The benchmark uses separate
processes for frozen baseline and candidate plugins and a monotonic timer.
See [the measurement report](docs/PERFORMANCE_NEWCOMMANDS.md) for workload
shapes, raw observations, build identities, results, and reproduction commands.
