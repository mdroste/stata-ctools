/*
    ctools_data_io.c
    Stata <-> C Data Transfer Module

    Handles bidirectional data transfer between Stata's data space and C memory
    using the Stata Plugin Interface (SPI). This module is designed to be
    reusable for any operation requiring bulk data transfer, not just sorting.

    Architecture:
    - Column-major storage: Each variable is stored in a contiguous array
    - Parallel I/O: One thread per variable when nvars >= 2
    - Uses SF_vdata()/SF_sdata() for reading, SF_vstore()/SF_sstore() for writing
    - SPI convention: variable index first, then observation index (both 1-based)

    Thread Safety:
    - Each thread operates on a different Stata variable (column)
    - No shared state during I/O phase
    - Assumes Stata's SPI is thread-safe for column-level parallelism
    - Note: Custom streaming code (e.g., in cmerge) should use sequential I/O
      to avoid potential race conditions with complex access patterns

    Performance Optimizations:
    - Parallel variable I/O overlaps operations across columns
    - 16x loop unrolling reduces loop overhead and improves instruction pipelining
    - Batched SPI calls via SF_VDATA_BATCH16/SF_VSTORE_BATCH16 macros
    - SD_FASTMODE (compile flag) disables SPI bounds checking
    - Cache-line aligned allocations for optimal memory access
    - Software prefetching for improved memory latency hiding
    - SIMD operations for bulk memory copies and initialization
    - Non-temporal stores for large write operations
    - String arena allocator for reduced malloc overhead
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_threads.h"
#include "ctools_config.h"
#include "ctools_spi.h"

/* SIMD and prefetch intrinsics - platform-specific */
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #define CTOOLS_X86 1
    #include <immintrin.h>
    #include <xmmintrin.h>
    #include <emmintrin.h>
    #ifdef __AVX2__
        #define CTOOLS_AVX2 1
    #endif
#elif defined(__aarch64__) || defined(_M_ARM64)
    #define CTOOLS_ARM64 1
    #include <arm_neon.h>
    /* ARM has hardware prefetch, but we can hint */
    #define _mm_prefetch(addr, hint) __builtin_prefetch(addr, 0, 3)
    #define _MM_HINT_T0 3
#else
    /* Fallback: no-op prefetch */
    #define _mm_prefetch(addr, hint) ((void)0)
    #define _MM_HINT_T0 0
#endif

/* Maximum string buffer size for Stata string variables */
#define STATA_STR_MAXLEN 2045

/* Threshold for non-temporal stores (bypass cache for very large writes) */
#define NONTEMPORAL_THRESHOLD 1000000

/* ===========================================================================
   Memory Allocation Utilities
   =========================================================================== */

/*
    Allocate memory aligned to cache line boundary.
    Returns NULL on failure.
*/
static inline void *aligned_alloc_cacheline(size_t size)
{
    void *ptr = NULL;
#if defined(_WIN32)
    ptr = _aligned_malloc(size, CACHE_LINE_SIZE);
#else
    if (posix_memalign(&ptr, CACHE_LINE_SIZE, size) != 0) {
        return NULL;
    }
#endif
    return ptr;
}

/*
    Free cache-line aligned memory.
*/
static inline void aligned_free(void *ptr)
{
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/* ===========================================================================
   String Arena Allocator
   =========================================================================== */

/*
    Arena allocator for string variables.
    Reduces malloc overhead by allocating strings from a contiguous block.
*/
typedef struct {
    char *base;         /* Base pointer to arena memory */
    size_t capacity;    /* Total arena capacity in bytes */
    size_t used;        /* Currently used bytes */
    int has_fallback;   /* Non-zero if any strings were allocated via strdup fallback */
} string_arena;

/*
    Create a string arena with given capacity.
    Returns NULL on allocation failure.
*/
static string_arena *arena_create(size_t capacity)
{
    string_arena *arena = (string_arena *)malloc(sizeof(string_arena));
    if (arena == NULL) return NULL;

    arena->base = (char *)malloc(capacity);
    if (arena->base == NULL) {
        free(arena);
        return NULL;
    }

    arena->capacity = capacity;
    arena->used = 0;
    arena->has_fallback = 0;
    return arena;
}

/*
    Allocate a string copy from the arena.
    Falls back to strdup if arena is full (and sets has_fallback flag).
    Returns NULL on complete failure.
*/
static char *arena_strdup(string_arena *arena, const char *s)
{
    size_t len = strlen(s) + 1;

    /* Check if arena has space - with overflow protection */
    if (arena != NULL && arena->used <= arena->capacity - len && len <= arena->capacity) {
        char *ptr = arena->base + arena->used;
        memcpy(ptr, s, len);
        arena->used += len;
        return ptr;
    }

    /* Fallback to regular strdup - mark that we have fallback allocations */
    if (arena != NULL) {
        arena->has_fallback = 1;
    }
    return strdup(s);
}

/*
    Free arena memory.
    Note: Individual strings allocated via fallback strdup must be freed separately.
*/
static void arena_free(string_arena *arena)
{
    if (arena != NULL) {
        free(arena->base);
        free(arena);
    }
}

/* ===========================================================================
   SIMD-Accelerated Utilities
   =========================================================================== */

/*
    Copy 8 doubles from local buffer to destination array.
    Uses SIMD when available for improved throughput.
*/
static inline void copy_doubles_8(double * restrict dst, const double * restrict src)
{
#if defined(CTOOLS_AVX2) || defined(CTOOLS_X86)
    #ifdef __AVX__
    _mm256_storeu_pd(&dst[0], _mm256_loadu_pd(&src[0]));
    _mm256_storeu_pd(&dst[4], _mm256_loadu_pd(&src[4]));
    #else
    _mm_storeu_pd(&dst[0], _mm_loadu_pd(&src[0]));
    _mm_storeu_pd(&dst[2], _mm_loadu_pd(&src[2]));
    _mm_storeu_pd(&dst[4], _mm_loadu_pd(&src[4]));
    _mm_storeu_pd(&dst[6], _mm_loadu_pd(&src[6]));
    #endif
#elif defined(CTOOLS_ARM64)
    vst1q_f64(&dst[0], vld1q_f64(&src[0]));
    vst1q_f64(&dst[2], vld1q_f64(&src[2]));
    vst1q_f64(&dst[4], vld1q_f64(&src[4]));
    vst1q_f64(&dst[6], vld1q_f64(&src[6]));
#else
    dst[0] = src[0]; dst[1] = src[1]; dst[2] = src[2]; dst[3] = src[3];
    dst[4] = src[4]; dst[5] = src[5]; dst[6] = src[6]; dst[7] = src[7];
#endif
}

/*
    Copy 16 doubles using SIMD (aggressive unrolling).
*/
static inline void copy_doubles_16(double * restrict dst, const double * restrict src)
{
#if defined(CTOOLS_AVX2) || defined(CTOOLS_X86)
    #ifdef __AVX__
    _mm256_storeu_pd(&dst[0], _mm256_loadu_pd(&src[0]));
    _mm256_storeu_pd(&dst[4], _mm256_loadu_pd(&src[4]));
    _mm256_storeu_pd(&dst[8], _mm256_loadu_pd(&src[8]));
    _mm256_storeu_pd(&dst[12], _mm256_loadu_pd(&src[12]));
    #else
    _mm_storeu_pd(&dst[0], _mm_loadu_pd(&src[0]));
    _mm_storeu_pd(&dst[2], _mm_loadu_pd(&src[2]));
    _mm_storeu_pd(&dst[4], _mm_loadu_pd(&src[4]));
    _mm_storeu_pd(&dst[6], _mm_loadu_pd(&src[6]));
    _mm_storeu_pd(&dst[8], _mm_loadu_pd(&src[8]));
    _mm_storeu_pd(&dst[10], _mm_loadu_pd(&src[10]));
    _mm_storeu_pd(&dst[12], _mm_loadu_pd(&src[12]));
    _mm_storeu_pd(&dst[14], _mm_loadu_pd(&src[14]));
    #endif
#elif defined(CTOOLS_ARM64)
    vst1q_f64(&dst[0], vld1q_f64(&src[0]));
    vst1q_f64(&dst[2], vld1q_f64(&src[2]));
    vst1q_f64(&dst[4], vld1q_f64(&src[4]));
    vst1q_f64(&dst[6], vld1q_f64(&src[6]));
    vst1q_f64(&dst[8], vld1q_f64(&src[8]));
    vst1q_f64(&dst[10], vld1q_f64(&src[10]));
    vst1q_f64(&dst[12], vld1q_f64(&src[12]));
    vst1q_f64(&dst[14], vld1q_f64(&src[14]));
#else
    copy_doubles_8(dst, src);
    copy_doubles_8(dst + 8, src + 8);
#endif
}

/*
    Non-temporal store for 8 doubles (bypasses cache).
    Use for large writes that won't be read again soon.
*/
static inline void stream_doubles_8(double * restrict dst, const double * restrict src)
{
#if defined(CTOOLS_X86) && defined(__SSE2__)
    _mm_stream_pd(&dst[0], _mm_loadu_pd(&src[0]));
    _mm_stream_pd(&dst[2], _mm_loadu_pd(&src[2]));
    _mm_stream_pd(&dst[4], _mm_loadu_pd(&src[4]));
    _mm_stream_pd(&dst[6], _mm_loadu_pd(&src[6]));
#elif defined(CTOOLS_ARM64)
    vst1q_f64(&dst[0], vld1q_f64(&src[0]));
    vst1q_f64(&dst[2], vld1q_f64(&src[2]));
    vst1q_f64(&dst[4], vld1q_f64(&src[4]));
    vst1q_f64(&dst[6], vld1q_f64(&src[6]));
#else
    copy_doubles_8(dst, src);
#endif
}

/*
    Non-temporal store for 16 doubles (bypasses cache).
    Note: Currently unused but kept for future optimization opportunities.
*/
static inline __attribute__((unused)) void stream_doubles_16(double * restrict dst, const double * restrict src)
{
#if defined(CTOOLS_X86) && defined(__SSE2__)
    _mm_stream_pd(&dst[0], _mm_loadu_pd(&src[0]));
    _mm_stream_pd(&dst[2], _mm_loadu_pd(&src[2]));
    _mm_stream_pd(&dst[4], _mm_loadu_pd(&src[4]));
    _mm_stream_pd(&dst[6], _mm_loadu_pd(&src[6]));
    _mm_stream_pd(&dst[8], _mm_loadu_pd(&src[8]));
    _mm_stream_pd(&dst[10], _mm_loadu_pd(&src[10]));
    _mm_stream_pd(&dst[12], _mm_loadu_pd(&src[12]));
    _mm_stream_pd(&dst[14], _mm_loadu_pd(&src[14]));
#elif defined(CTOOLS_ARM64)
    copy_doubles_16(dst, src);
#else
    copy_doubles_16(dst, src);
#endif
}

/* ===========================================================================
   Data Loading (Stata -> C)
   =========================================================================== */

/*
    Thread argument structure for multi-variable I/O.
    Each worker thread processes a range of variables [var_start, var_end).
*/
typedef struct {
    stata_variable *vars;       /* Array of variable structures */
    int *var_indices;           /* Array of 1-based Stata variable indices */
    size_t var_start;           /* First variable index (0-based in arrays) */
    size_t var_end;             /* End variable index (exclusive) */
    size_t obs1;                /* First observation (1-based) */
    size_t nobs;                /* Number of observations */
    int *is_string;             /* Array indicating if each var is string */
    int success;                /* 1 on success, 0 on failure */
} ctools_multi_var_io_args;

/*
    Load a single variable from Stata into C memory (internal helper).
    Used by both single-var and multi-var loading.
*/
static int load_single_variable(stata_variable *var, int var_idx, size_t obs1,
                                 size_t nobs, int is_string)
{
    size_t i;
    char strbuf[STATA_STR_MAXLEN + 1];

    var->nobs = nobs;

    if (is_string) {
        /* String variable - use arena allocator for fast bulk free */
        var->type = STATA_TYPE_STRING;
        var->str_maxlen = STATA_STR_MAXLEN;
        var->_arena = NULL;

        /* Allocate pointer array (cache-line aligned, zero-initialized for safe cleanup) */
        size_t str_array_size;
        if (ctools_safe_mul_size(nobs, sizeof(char *), &str_array_size) != 0) {
            return -1;  /* Overflow */
        }
        var->data.str = (char **)aligned_alloc_cacheline(str_array_size);
        if (var->data.str == NULL) {
            return -1;
        }
        memset(var->data.str, 0, str_array_size);

        char **str_ptr = var->data.str;

        /* Create arena for all strings - estimate avg 64 bytes per string */
        /* Overflow check: nobs * 64 */
        size_t arena_capacity = 0;
        if (nobs <= SIZE_MAX / 64) {
            arena_capacity = nobs * 64;
        } else {
            arena_capacity = SIZE_MAX;  /* Will likely fail, triggering strdup fallback */
        }
        string_arena *arena = arena_create(arena_capacity);
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings with prefetching */
        for (i = 0; i < nobs; i++) {
            if (i + PREFETCH_DISTANCE < nobs) {
                _mm_prefetch((const char *)&str_ptr[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }

            SF_sdata((ST_int)var_idx, (ST_int)(i + obs1), strbuf);
            str_ptr[i] = arena_strdup(arena, strbuf);
            if (str_ptr[i] == NULL) {
                return -1;
            }
        }
    } else {
        /* Numeric variable - use cache-line aligned allocation */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;
        size_t dbl_array_size;
        if (ctools_safe_mul_size(nobs, sizeof(double), &dbl_array_size) != 0) {
            return -1;  /* Overflow */
        }
        var->data.dbl = (double *)aligned_alloc_cacheline(dbl_array_size);

        if (var->data.dbl == NULL) {
            return -1;
        }

        double * restrict dbl_ptr = var->data.dbl;

        /* Double buffers for software pipelining - use pointer swap */
        double buf_a[16], buf_b[16];
        double *v0 = buf_a, *v1 = buf_b;

        /* Calculate loop bounds */
        size_t i_end_16 = nobs - (nobs % 16);
        size_t prefetch_end = (nobs > 32) ? (i_end_16 - 32) : 0;

        /* Phase 1: Main loop with prefetching (bulk of iterations) */
        i = 0;
        if (prefetch_end > 0) {
            /* Prime the pipeline */
            SF_VDATA_BATCH16((ST_int)var_idx, obs1, v0);

            for (; i < prefetch_end; i += 16) {
                /* Prefetch 2 cache lines ahead */
                _mm_prefetch((const char *)&dbl_ptr[i + 32], _MM_HINT_T0);
                _mm_prefetch((const char *)&dbl_ptr[i + 40], _MM_HINT_T0);

                /* Load next batch into v1 */
                SF_VDATA_BATCH16((ST_int)var_idx, i + 16 + obs1, v1);

                /* Store previous batch */
                copy_doubles_16(&dbl_ptr[i], v0);

                /* Swap buffer pointers instead of memcpy */
                double *tmp = v0;
                v0 = v1;
                v1 = tmp;
            }
        }

        /* Phase 2: Remaining 16-element blocks without prefetch */
        for (; i < i_end_16; i += 16) {
            SF_VDATA_BATCH16((ST_int)var_idx, i + obs1, v0);
            copy_doubles_16(&dbl_ptr[i], v0);
        }

        /* Phase 3: Handle remaining elements (< 16) */
        for (; i < nobs; i++) {
            SF_vdata((ST_int)var_idx, (ST_int)(i + obs1), &dbl_ptr[i]);
        }
    }

    return 0;
}

/*
    Thread function: Load multiple variables from Stata into C memory.
    Each thread processes variables in range [var_start, var_end).
*/
static void *load_multi_variable_thread(void *arg)
{
    ctools_multi_var_io_args *args = (ctools_multi_var_io_args *)arg;
    args->success = 0;

    for (size_t v = args->var_start; v < args->var_end; v++) {
        if (load_single_variable(&args->vars[v], args->var_indices[v],
                                  args->obs1, args->nobs, args->is_string[v]) != 0) {
            return (void *)1;
        }
    }

    args->success = 1;
    return NULL;
}

/*
    Thread function: Load a single variable from Stata into C memory.
    (Legacy interface for backward compatibility)

    Allocates memory and reads all observations for one variable using
    SF_vdata (numeric) or SF_sdata (string). Uses optimized loops with:
    - 8x unrolling for reduced loop overhead
    - Batched SPI calls via macros
    - Separated load/store phases for better pipelining
    - Software prefetching for write destination
    - Cache-line aligned allocations
    - String arena for reduced malloc overhead

    @param arg  Pointer to ctools_var_io_args with input/output parameters
    @return     NULL on success, non-NULL on failure (for ctools_threads)
*/
static void *load_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    args->success = 0;

    if (load_single_variable(args->var, args->var_idx, args->obs1,
                              args->nobs, args->is_string) != 0) {
        return (void *)1;
    }

    args->success = 1;
    return NULL;
}

/*
    Load all variables from Stata's data space into C memory.

    This is the core data transfer function - it performs pure data loading
    without any sort-specific setup. Use this for general-purpose data transfer
    between Stata and C.

    @param data   [out] stata_data structure to populate (caller frees)
    @param nvars  [in]  Total number of variables in dataset

    @return STATA_OK on success, or:
            STATA_ERR_INVALID_INPUT if data is NULL or nvars is 0
            STATA_ERR_MEMORY on allocation failure

    Memory Layout After Return:
    - data->vars[j].data.dbl: contiguous array of nobs doubles (numeric)
    - data->vars[j].data.str: array of nobs char* pointers (string)
    - data->sort_order: identity permutation [0, 1, 2, ..., nobs-1]
*/
stata_retcode ctools_data_load(stata_data *data, size_t nvars)
{
    size_t nobs;
    size_t obs1;
    ctools_var_io_args *thread_args = NULL;
    size_t j;
    int use_parallel;

    if (data == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Cache observation bounds */
    obs1 = (size_t)SF_in1();
    nobs = (size_t)(SF_in2() - SF_in1() + 1);

    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = (size_t)nobs;
    data->nvars = nvars;

    /* Allocate array of variables (cache-line aligned for better access) */
    data->vars = (stata_variable *)aligned_alloc_cacheline(nvars * sizeof(stata_variable));
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }
    memset(data->vars, 0, nvars * sizeof(stata_variable));

    /* Allocate sort order array (cache-line aligned, uses perm_idx_t for 50% memory savings) */
    data->sort_order = (perm_idx_t *)aligned_alloc_cacheline(nobs * sizeof(perm_idx_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation */
    for (size_t i = 0; i < nobs; i++) {
        data->sort_order[i] = (perm_idx_t)i;
    }

    /* Allocate thread arguments (used for both parallel and sequential) */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = (int)(j + 1);
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = SF_var_is_string((ST_int)(j + 1));
        thread_args[j].success = 0;
    }

    /* Parallel loading with capped thread count.
     * Cap threads to min(nvars, CTOOLS_IO_MAX_THREADS) to avoid oversubscription.
     * Each thread processes a chunk of variables for better efficiency. */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Determine optimal thread count */
        size_t num_threads = nvars;
        if (num_threads > CTOOLS_IO_MAX_THREADS) {
            num_threads = CTOOLS_IO_MAX_THREADS;
        }

        /* If few variables, use legacy one-thread-per-variable approach */
        if (nvars <= CTOOLS_IO_MAX_THREADS) {
            /* Legacy path: one thread per variable */
            ctools_thread_pool pool;
            if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
                free(thread_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }

            int pool_result = ctools_pool_run(&pool, load_variable_thread);
            ctools_pool_free(&pool);
            free(thread_args);

            if (pool_result != 0) {
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }
        } else {
            /* Multi-variable per thread: distribute variables among capped threads */
            free(thread_args);  /* Don't need single-var args */

            /* Build arrays for multi-var threading */
            int *var_indices = (int *)malloc(nvars * sizeof(int));
            int *is_string_arr = (int *)malloc(nvars * sizeof(int));
            ctools_multi_var_io_args *multi_args = (ctools_multi_var_io_args *)malloc(
                num_threads * sizeof(ctools_multi_var_io_args));

            if (var_indices == NULL || is_string_arr == NULL || multi_args == NULL) {
                free(var_indices);
                free(is_string_arr);
                free(multi_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }

            /* Initialize variable info */
            for (j = 0; j < nvars; j++) {
                var_indices[j] = (int)(j + 1);
                is_string_arr[j] = SF_var_is_string((ST_int)(j + 1));
            }

            /* Distribute variables among threads */
            size_t vars_per_thread = nvars / num_threads;
            size_t extra_vars = nvars % num_threads;

            size_t var_offset = 0;
            for (size_t t = 0; t < num_threads; t++) {
                size_t this_thread_vars = vars_per_thread + (t < extra_vars ? 1 : 0);
                multi_args[t].vars = data->vars;
                multi_args[t].var_indices = var_indices;
                multi_args[t].var_start = var_offset;
                multi_args[t].var_end = var_offset + this_thread_vars;
                multi_args[t].obs1 = obs1;
                multi_args[t].nobs = nobs;
                multi_args[t].is_string = is_string_arr;
                multi_args[t].success = 0;
                var_offset += this_thread_vars;
            }

            /* Run multi-variable threads */
            ctools_thread_pool pool;
            if (ctools_pool_init(&pool, num_threads, multi_args,
                                  sizeof(ctools_multi_var_io_args)) != 0) {
                free(var_indices);
                free(is_string_arr);
                free(multi_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }

            int pool_result = ctools_pool_run(&pool, load_multi_variable_thread);
            ctools_pool_free(&pool);
            free(var_indices);
            free(is_string_arr);
            free(multi_args);

            if (pool_result != 0) {
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }
        }
    } else {
        /* Sequential loading for single variable */
        for (j = 0; j < nvars; j++) {
            if (load_variable_thread(&thread_args[j]) != NULL) {
                free(thread_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }
        }
        free(thread_args);
    }

    /* Memory barrier to ensure all thread-side effects are visible before returning */
    ctools_memory_barrier();

    return STATA_OK;
}

/* ===========================================================================
   Data Storing (C -> Stata)
   =========================================================================== */

/*
    Store a single variable from C memory to Stata (internal helper).
    Used by both single-var and multi-var storing.
*/
static void store_single_variable(stata_variable *var, int var_idx, size_t obs1, size_t nobs)
{
    size_t i;
    ST_int stata_var_idx = (ST_int)var_idx;

    if (var->type == STATA_TYPE_DOUBLE) {
        /* Numeric variable - aggressive optimization */
        const double * restrict dbl_data = var->data.dbl;

        /* Double buffers for software pipelining - use pointer swap */
        double buf_a[16], buf_b[16];
        double *v0 = buf_a, *v1 = buf_b;

        /* Calculate loop bounds */
        size_t i_end_16 = nobs - (nobs % 16);
        size_t prefetch_end = (nobs > 32) ? (i_end_16 - 32) : 0;

        /* Phase 1: Main loop with prefetching and pipelining */
        i = 0;
        if (prefetch_end > 0) {
            /* Prime the pipeline */
            copy_doubles_16(v0, &dbl_data[0]);

            for (; i < prefetch_end; i += 16) {
                /* Prefetch source data ahead */
                _mm_prefetch((const char *)&dbl_data[i + 32], _MM_HINT_T0);
                _mm_prefetch((const char *)&dbl_data[i + 40], _MM_HINT_T0);

                /* Load next batch while writing current */
                copy_doubles_16(v1, &dbl_data[i + 16]);

                /* Write current batch to Stata */
                SF_VSTORE_BATCH16(stata_var_idx, i + obs1, v0);

                /* Swap buffer pointers instead of memcpy */
                double *tmp = v0;
                v0 = v1;
                v1 = tmp;
            }
        }

        /* Phase 2: Remaining 16-element blocks without prefetch */
        for (; i < i_end_16; i += 16) {
            copy_doubles_16(v0, &dbl_data[i]);
            SF_VSTORE_BATCH16(stata_var_idx, i + obs1, v0);
        }

        /* Phase 3: Handle remaining elements */
        for (; i < nobs; i++) {
            SF_vstore(stata_var_idx, (ST_int)(i + obs1), dbl_data[i]);
        }

    } else {
        /* String variable - optimized with prefetching */
        char * const * restrict str_data = var->data.str;

        /* Calculate prefetch region */
        size_t prefetch_end = (nobs > PREFETCH_DISTANCE) ? (nobs - PREFETCH_DISTANCE) : 0;

        /* Phase 1: Main loop with prefetching */
        for (i = 0; i < prefetch_end; i++) {
            /* Prefetch pointer array and string content ahead */
            _mm_prefetch((const char *)&str_data[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            _mm_prefetch(str_data[i + PREFETCH_DISTANCE], _MM_HINT_T0);

            SF_sstore(stata_var_idx, (ST_int)(i + obs1), str_data[i]);
        }

        /* Phase 2: Tail without prefetch */
        for (; i < nobs; i++) {
            SF_sstore(stata_var_idx, (ST_int)(i + obs1), str_data[i]);
        }
    }
}

/*
    Thread function: Store multiple variables from C memory to Stata.
    Each thread processes variables in range [var_start, var_end).
*/
static void *store_multi_variable_thread(void *arg)
{
    ctools_multi_var_io_args *args = (ctools_multi_var_io_args *)arg;
    args->success = 0;

    for (size_t v = args->var_start; v < args->var_end; v++) {
        store_single_variable(&args->vars[v], args->var_indices[v],
                              args->obs1, args->nobs);
    }

    args->success = 1;
    return NULL;
}

/*
    Thread function: Store a single variable from C memory to Stata.
    (Legacy interface for backward compatibility)

    Writes all observations for one variable using SF_vstore (numeric) or
    SF_sstore (string). Aggressive optimizations include:
    - 16x unrolled loop with batched SPI writes
    - Software pipelining with double buffering
    - Branch-free prefetch regions
    - restrict keyword for aliasing optimization

    @param arg  Pointer to ctools_var_io_args with input/output parameters
    @return     NULL on success (for ctools_threads compatibility)
*/
static void *store_variable_thread(void *arg)
{
    ctools_var_io_args *args = (ctools_var_io_args *)arg;
    store_single_variable(args->var, args->var_idx, args->obs1, args->nobs);
    args->success = 1;
    return NULL;
}

/*
    Write all variables from C memory back to Stata's data space.

    This is the main entry point for the C-to-Stata data transfer module.
    It writes all variables back to Stata in their current order. The data
    must already be in final sorted order (permutation applied by the sort
    module).

    @param data  [in] stata_data structure with data to write; data arrays
                      must be in final order (sequential read, no permutation)
    @param obs1  [in] First observation in Stata (1-based, from SF_in1)

    @return STATA_OK on success, or:
            STATA_ERR_INVALID_INPUT if data is NULL
            STATA_ERR_MEMORY on allocation failure (parallel mode only)

    Parallelization:
    - Uses parallel writes when nvars >= 2 (one thread per variable)
    - Falls back to sequential writes for single-variable datasets
*/
stata_retcode ctools_data_store(stata_data *data, size_t obs1)
{
    size_t j;
    size_t nobs, nvars;
    ctools_var_io_args *thread_args = NULL;
    int use_parallel;

    if (data == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Cache frequently accessed values */
    nobs = data->nobs;
    nvars = data->nvars;

    /* Allocate thread arguments (used for both parallel and sequential) */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = (int)(j + 1);
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = 0;  /* Not used by store */
        thread_args[j].success = 0;
    }

    /* Parallel storing with capped thread count.
     * Cap threads to min(nvars, CTOOLS_IO_MAX_THREADS) to avoid oversubscription. */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Determine optimal thread count */
        size_t num_threads = nvars;
        if (num_threads > CTOOLS_IO_MAX_THREADS) {
            num_threads = CTOOLS_IO_MAX_THREADS;
        }

        /* If few variables, use legacy one-thread-per-variable approach */
        if (nvars <= CTOOLS_IO_MAX_THREADS) {
            /* Legacy path: one thread per variable */
            ctools_thread_pool pool;
            if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
                free(thread_args);
                return STATA_ERR_MEMORY;
            }

            ctools_pool_run(&pool, store_variable_thread);
            ctools_pool_free(&pool);
            free(thread_args);
        } else {
            /* Multi-variable per thread: distribute variables among capped threads */
            free(thread_args);  /* Don't need single-var args */

            /* Build arrays for multi-var threading */
            int *var_indices = (int *)malloc(nvars * sizeof(int));
            int *is_string_arr = (int *)malloc(nvars * sizeof(int));
            ctools_multi_var_io_args *multi_args = (ctools_multi_var_io_args *)malloc(
                num_threads * sizeof(ctools_multi_var_io_args));

            if (var_indices == NULL || is_string_arr == NULL || multi_args == NULL) {
                free(var_indices);
                free(is_string_arr);
                free(multi_args);
                return STATA_ERR_MEMORY;
            }

            /* Initialize variable info */
            for (j = 0; j < nvars; j++) {
                var_indices[j] = (int)(j + 1);
                is_string_arr[j] = 0;  /* Not used by store */
            }

            /* Distribute variables among threads */
            size_t vars_per_thread = nvars / num_threads;
            size_t extra_vars = nvars % num_threads;

            size_t var_offset = 0;
            for (size_t t = 0; t < num_threads; t++) {
                size_t this_thread_vars = vars_per_thread + (t < extra_vars ? 1 : 0);
                multi_args[t].vars = data->vars;
                multi_args[t].var_indices = var_indices;
                multi_args[t].var_start = var_offset;
                multi_args[t].var_end = var_offset + this_thread_vars;
                multi_args[t].obs1 = obs1;
                multi_args[t].nobs = nobs;
                multi_args[t].is_string = is_string_arr;
                multi_args[t].success = 0;
                var_offset += this_thread_vars;
            }

            /* Run multi-variable threads */
            ctools_thread_pool pool;
            if (ctools_pool_init(&pool, num_threads, multi_args,
                                  sizeof(ctools_multi_var_io_args)) != 0) {
                free(var_indices);
                free(is_string_arr);
                free(multi_args);
                return STATA_ERR_MEMORY;
            }

            ctools_pool_run(&pool, store_multi_variable_thread);
            ctools_pool_free(&pool);
            free(var_indices);
            free(is_string_arr);
            free(multi_args);
        }
    } else {
        /* Sequential storing for single variable */
        for (j = 0; j < nvars; j++) {
            store_variable_thread(&thread_args[j]);
        }
        free(thread_args);
    }

    /* Memory barrier to ensure all thread-side effects are visible before returning */
    ctools_memory_barrier();

    return STATA_OK;
}

/* ===========================================================================
   Selective Data Loading (Stata -> C) - Load only specified variables
   =========================================================================== */

/*
    Load only specified variables from Stata into C memory.

    This function loads a subset of variables, identified by their 1-based
    Stata indices. This is more memory-efficient when only specific columns
    are needed (e.g., key variables for merge operations).

    @param data        [out] stata_data structure to populate
    @param var_indices [in]  Array of 1-based Stata variable indices
    @param nvars       [in]  Number of variables to load
    @param obs_start   [in]  First observation (1-based), 0 = use SF_in1()
    @param obs_end     [in]  Last observation (1-based), 0 = use SF_in2()

    @return STATA_OK on success, error code otherwise
*/
stata_retcode ctools_data_load_selective(stata_data *data, int *var_indices,
                                          size_t nvars, size_t obs_start, size_t obs_end)
{
    size_t nobs;
    size_t obs1;
    ctools_var_io_args *thread_args = NULL;
    size_t j;
    int use_parallel;

    if (data == NULL || var_indices == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Determine observation range */
    if (obs_start == 0 || obs_end == 0) {
        obs1 = (size_t)SF_in1();
        nobs = (size_t)(SF_in2() - SF_in1() + 1);
    } else {
        obs1 = obs_start;
        nobs = obs_end - obs_start + 1;
    }

    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = nobs;
    data->nvars = nvars;

    /* Allocate array of variables (cache-line aligned) */
    size_t vars_size;
    if (ctools_safe_mul_size(nvars, sizeof(stata_variable), &vars_size) != 0) {
        return STATA_ERR_MEMORY;  /* Overflow */
    }
    data->vars = (stata_variable *)aligned_alloc_cacheline(vars_size);
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }
    memset(data->vars, 0, vars_size);

    /* Allocate sort order array (cache-line aligned, uses perm_idx_t for 50% memory savings) */
    size_t sort_order_size;
    if (ctools_safe_mul_size(nobs, sizeof(perm_idx_t), &sort_order_size) != 0) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;  /* Overflow */
    }
    data->sort_order = (perm_idx_t *)aligned_alloc_cacheline(sort_order_size);
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation */
    for (size_t i = 0; i < nobs; i++) {
        data->sort_order[i] = (perm_idx_t)i;
    }

    /* Allocate thread arguments */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments with SPECIFIED variable indices */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = var_indices[j];  /* Use specified index (1-based) */
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = SF_var_is_string((ST_int)var_indices[j]);
        thread_args[j].success = 0;
    }

    /* Parallel loading when multiple variables - one thread per variable.
     * Each thread operates on a different Stata variable (column).
     * Note: This assumes Stata's SPI is thread-safe for column-level access. */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Parallel loading using thread pool */
        ctools_thread_pool pool;
        if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
            free(thread_args);
            stata_data_free(data);
            return STATA_ERR_MEMORY;
        }

        int pool_result = ctools_pool_run(&pool, load_variable_thread);
        ctools_pool_free(&pool);
        free(thread_args);

        if (pool_result != 0) {
            stata_data_free(data);
            return STATA_ERR_MEMORY;
        }
    } else {
        /* Sequential loading for single variable */
        for (j = 0; j < nvars; j++) {
            if (load_variable_thread(&thread_args[j]) != NULL) {
                free(thread_args);
                stata_data_free(data);
                return STATA_ERR_MEMORY;
            }
        }
        free(thread_args);
    }

    /* Memory barrier to ensure all thread-side effects are visible before returning */
    ctools_memory_barrier();

    return STATA_OK;
}

/* ===========================================================================
   Selective Data Storing (C -> Stata) - Store to specified variables
   =========================================================================== */

/*
    Store variables from C memory to specified Stata variable indices.

    This function writes data to specific Stata variables, identified by their
    1-based indices. This is needed when the C data structure doesn't map 1:1
    to Stata's variable order (e.g., after merge operations).

    @param data        [in]  stata_data structure containing data to write
    @param var_indices [in]  Array of 1-based Stata variable indices
                             var_indices[j] is the Stata variable for data->vars[j]
    @param nvars       [in]  Number of variables to store
    @param obs1        [in]  First observation in Stata (1-based)

    @return STATA_OK on success, error code otherwise
*/
stata_retcode ctools_data_store_selective(stata_data *data, int *var_indices,
                                           size_t nvars, size_t obs1)
{
    size_t j;
    size_t nobs;
    ctools_var_io_args *thread_args = NULL;
    int use_parallel;

    if (data == NULL || var_indices == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    nobs = data->nobs;

    /* Allocate thread arguments */
    thread_args = (ctools_var_io_args *)malloc(nvars * sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments with SPECIFIED variable indices */
    for (j = 0; j < nvars; j++) {
        thread_args[j].var = &data->vars[j];
        thread_args[j].var_idx = var_indices[j];  /* Use specified index (1-based) */
        thread_args[j].obs1 = obs1;
        thread_args[j].nobs = nobs;
        thread_args[j].is_string = 0;  /* Not used by store */
        thread_args[j].success = 0;
    }

    /* Parallel storing when multiple variables */
    use_parallel = (nvars >= 2);

    if (use_parallel) {
        /* Parallel storing using thread pool */
        ctools_thread_pool pool;
        if (ctools_pool_init(&pool, nvars, thread_args, sizeof(ctools_var_io_args)) != 0) {
            free(thread_args);
            return STATA_ERR_MEMORY;
        }

        ctools_pool_run(&pool, store_variable_thread);
        ctools_pool_free(&pool);
    } else {
        /* Sequential storing for single variable */
        for (j = 0; j < nvars; j++) {
            store_variable_thread(&thread_args[j]);
        }
    }

    free(thread_args);

    /* Memory barrier to ensure all thread-side effects are visible */
    ctools_memory_barrier();

    return STATA_OK;
}

/* ===========================================================================
   Streaming Permuted Variable Write (C -> Stata with reordering)
   =========================================================================== */

/*
    Write a single variable to Stata with row permutation.

    For each output row i, reads from Stata row source_rows[i] and writes
    to Stata row i. This enables in-place reordering without loading the
    entire variable into C memory first.

    Optimizations:
    - Cache-line aligned buffer allocation
    - Software prefetching for source_rows lookups
    - String arena for reduced malloc overhead
    - Batched SPI writes with SIMD buffer reads
    - Non-temporal stores for large datasets

    @param var_idx      [in] 1-based Stata variable index
    @param source_rows  [in] Array mapping output row -> source row (0-based)
                             -1 indicates missing value should be written
    @param output_nobs  [in] Number of output rows
    @param obs1         [in] First output observation in Stata (1-based)

    @return STATA_OK on success, error code otherwise
*/
stata_retcode ctools_stream_var_permuted(int var_idx, int64_t *source_rows,
                                          size_t output_nobs, size_t obs1)
{
    size_t i;
    ST_int stata_var = (ST_int)var_idx;
    int is_string = SF_var_is_string(stata_var);

    if (is_string) {
        /* String variable: allocate buffer, gather, scatter */
        char **buf = (char **)aligned_alloc_cacheline(output_nobs * sizeof(char *));
        if (!buf) return STATA_ERR_MEMORY;
        memset(buf, 0, output_nobs * sizeof(char *));  /* Zero-init for safe cleanup */

        char strbuf[2048];

        /* Create arena for string storage (estimate 64 bytes avg per string) */
        /* Overflow check: output_nobs * 64 */
        size_t arena_capacity = 0;
        if (output_nobs <= SIZE_MAX / 64) {
            arena_capacity = output_nobs * 64;
        } else {
            arena_capacity = SIZE_MAX;  /* Will likely fail, triggering strdup fallback */
        }
        string_arena *arena = arena_create(arena_capacity);
        /* Arena failure is ok - we fall back to strdup */

        /* GATHER: Read from source positions with prefetching */
        for (i = 0; i < output_nobs; i++) {
            /* Prefetch source_rows array ahead */
            if (i + PREFETCH_DISTANCE < output_nobs) {
                _mm_prefetch((const char *)&source_rows[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }

            if (source_rows[i] >= 0) {
                SF_sdata(stata_var, (ST_int)(source_rows[i] + obs1), strbuf);
                buf[i] = arena_strdup(arena, strbuf);
            } else {
                buf[i] = arena_strdup(arena, "");  /* Missing -> empty string */
            }
            if (!buf[i]) {
                /* Cleanup on allocation failure - free any fallback strings first */
                if (arena != NULL && arena->has_fallback) {
                    for (size_t j = 0; j < i; j++) {
                        if (buf[j] != NULL &&
                            (buf[j] < arena->base || buf[j] >= arena->base + arena->capacity)) {
                            free(buf[j]);
                        }
                    }
                } else if (arena == NULL) {
                    /* No arena - free all strdup'd strings */
                    for (size_t j = 0; j < i; j++) {
                        free(buf[j]);
                    }
                }
                arena_free(arena);
                aligned_free(buf);
                return STATA_ERR_MEMORY;
            }
        }

        /* SCATTER: Write to sequential output positions with prefetching */
        for (i = 0; i < output_nobs; i++) {
            /* Prefetch buffer pointers and string content ahead */
            if (i + PREFETCH_DISTANCE < output_nobs) {
                _mm_prefetch((const char *)&buf[i + PREFETCH_DISTANCE], _MM_HINT_T0);
                if (buf[i + PREFETCH_DISTANCE / 2] != NULL) {
                    _mm_prefetch(buf[i + PREFETCH_DISTANCE / 2], _MM_HINT_T0);
                }
            }
            SF_sstore(stata_var, (ST_int)(i + obs1), buf[i]);
        }

        /* Free fallback strings (those outside arena range), then free arena */
        if (arena != NULL && arena->has_fallback) {
            for (i = 0; i < output_nobs; i++) {
                if (buf[i] != NULL &&
                    (buf[i] < arena->base || buf[i] >= arena->base + arena->capacity)) {
                    free(buf[i]);
                }
            }
        } else if (arena == NULL) {
            /* No arena - free all strdup'd strings */
            for (i = 0; i < output_nobs; i++) {
                if (buf[i] != NULL) {
                    free(buf[i]);
                }
            }
        }
        arena_free(arena);
        aligned_free(buf);

    } else {
        /* Numeric variable: allocate aligned buffer, gather, scatter */
        double *buf = (double *)aligned_alloc_cacheline(output_nobs * sizeof(double));
        if (!buf) return STATA_ERR_MEMORY;

        /* GATHER: Read from source positions with prefetching */
        for (i = 0; i < output_nobs; i++) {
            /* Prefetch source_rows array ahead for next lookup */
            if (i + PREFETCH_DISTANCE < output_nobs) {
                _mm_prefetch((const char *)&source_rows[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }
            /* Prefetch destination buffer ahead for writes */
            if (i + PREFETCH_DISTANCE < output_nobs) {
                _mm_prefetch((const char *)&buf[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }

            if (source_rows[i] >= 0) {
                SF_vdata(stata_var, (ST_int)(source_rows[i] + obs1), &buf[i]);
            } else {
                buf[i] = SV_missval;  /* Missing value */
            }
        }

        /* SCATTER: Write to sequential output positions */
        size_t i_end = output_nobs - (output_nobs % 8);

        /* Use non-temporal stores for very large datasets */
        if (output_nobs >= NONTEMPORAL_THRESHOLD) {
            /* Local buffer for batched reads */
            double v[8];

            for (i = 0; i < i_end; i += 8) {
                /* Stream read from buffer (bypasses cache on way back) */
                stream_doubles_8(v, &buf[i]);
                /* Batch write to Stata */
                SF_VSTORE_BATCH8(stata_var, i + obs1, v);
            }

#if defined(CTOOLS_X86) && defined(__SSE2__)
            _mm_sfence();  /* Ensure all streaming stores complete */
#endif
        } else {
            /* Standard path with prefetching */
            double v[8];

            for (i = 0; i < i_end; i += 8) {
                /* Prefetch source data ahead */
                if (i + PREFETCH_DISTANCE < output_nobs) {
                    _mm_prefetch((const char *)&buf[i + PREFETCH_DISTANCE], _MM_HINT_T0);
                }

                /* Batch read using SIMD */
                copy_doubles_8(v, &buf[i]);
                /* Batch write to Stata */
                SF_VSTORE_BATCH8(stata_var, i + obs1, v);
            }
        }

        /* Handle remaining elements */
        for (; i < output_nobs; i++) {
            SF_vstore(stata_var, (ST_int)(i + obs1), buf[i]);
        }

        aligned_free(buf);
    }

    return STATA_OK;
}
