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
#include "ctools_arena.h"
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
    /* ARM64: Use STNP (Store Pair Non-temporal) for cache bypass */
    __asm__ __volatile__(
        "stnp %d[v0], %d[v1], [%[dst]]\n\t"
        "stnp %d[v2], %d[v3], [%[dst], #16]\n\t"
        "stnp %d[v4], %d[v5], [%[dst], #32]\n\t"
        "stnp %d[v6], %d[v7], [%[dst], #48]"
        :
        : [v0] "w" (src[0]), [v1] "w" (src[1]),
          [v2] "w" (src[2]), [v3] "w" (src[3]),
          [v4] "w" (src[4]), [v5] "w" (src[5]),
          [v6] "w" (src[6]), [v7] "w" (src[7]),
          [dst] "r" (dst)
        : "memory"
    );
#else
    copy_doubles_8(dst, src);
#endif
}

/*
    Non-temporal store for 16 doubles (bypasses cache).
    Use for large writes that won't be read again soon.
*/
__attribute__((unused))
static inline void stream_doubles_16(double * restrict dst, const double * restrict src)
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
    /* ARM64: Use STNP (Store Pair Non-temporal) for cache bypass */
    __asm__ __volatile__(
        "stnp %d[v0], %d[v1], [%[dst]]\n\t"
        "stnp %d[v2], %d[v3], [%[dst], #16]\n\t"
        "stnp %d[v4], %d[v5], [%[dst], #32]\n\t"
        "stnp %d[v6], %d[v7], [%[dst], #48]\n\t"
        "stnp %d[v8], %d[v9], [%[dst], #64]\n\t"
        "stnp %d[v10], %d[v11], [%[dst], #80]\n\t"
        "stnp %d[v12], %d[v13], [%[dst], #96]\n\t"
        "stnp %d[v14], %d[v15], [%[dst], #112]"
        :
        : [v0] "w" (src[0]), [v1] "w" (src[1]),
          [v2] "w" (src[2]), [v3] "w" (src[3]),
          [v4] "w" (src[4]), [v5] "w" (src[5]),
          [v6] "w" (src[6]), [v7] "w" (src[7]),
          [v8] "w" (src[8]), [v9] "w" (src[9]),
          [v10] "w" (src[10]), [v11] "w" (src[11]),
          [v12] "w" (src[12]), [v13] "w" (src[13]),
          [v14] "w" (src[14]), [v15] "w" (src[15]),
          [dst] "r" (dst)
        : "memory"
    );
#else
    copy_doubles_16(dst, src);
#endif
}

/* ===========================================================================
   Data Loading (Stata -> C)
   =========================================================================== */

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

    /* Edge case: if no observations, allocate minimal buffer to avoid NULL pointer issues.
     * This prevents crashes when downstream code assumes data.dbl/data.str is non-NULL
     * after successful load. Allocating 1 element is minimal overhead for safety. */
    if (nobs == 0) {
        if (is_string) {
            var->type = STATA_TYPE_STRING;
            var->str_maxlen = 0;
            var->_arena = NULL;
            /* Allocate minimal buffer to avoid NULL - calloc ensures it's zeroed */
            var->data.str = (char **)calloc(1, sizeof(char *));
            if (!var->data.str) return -1;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            /* Allocate minimal buffer to avoid NULL */
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return -1;
            var->data.dbl[0] = SV_missval;  /* Initialize to missing */
        }
        return 0;  /* Success with empty data */
    }

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
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            return -1;
        }
        memset(var->data.str, 0, str_array_size);

        char **str_ptr = var->data.str;

        /* Create arena for all strings - estimate avg 64 bytes per string */
        /* Overflow check: nobs * 64 - skip arena if overflow, rely on strdup fallback */
        ctools_string_arena *arena = NULL;
        if (nobs <= SIZE_MAX / 64) {
            size_t arena_capacity = nobs * 64;
            arena = ctools_string_arena_create(arena_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings with prefetching */
        for (i = 0; i < nobs; i++) {
            if (i + PREFETCH_DISTANCE < nobs) {
                _mm_prefetch((const char *)&str_ptr[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }

            SF_sdata((ST_int)var_idx, (ST_int)(i + obs1), strbuf);
            str_ptr[i] = ctools_string_arena_strdup(arena, strbuf);
            if (str_ptr[i] == NULL) {
                /* Cleanup on allocation failure:
                 * - Free fallback strings (not owned by arena)
                 * - Free the string pointer array
                 * - Free the arena itself */
                if (arena != NULL) {
                    for (size_t j = 0; j < i; j++) {
                        if (str_ptr[j] != NULL && !ctools_string_arena_owns(arena, str_ptr[j])) {
                            free(str_ptr[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    var->_arena = NULL;
                } else {
                    /* No arena - all strings were strdup'd */
                    for (size_t j = 0; j < i; j++) {
                        free(str_ptr[j]);
                    }
                }
                ctools_aligned_free(var->data.str);
                var->data.str = NULL;
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
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);

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

/* ===========================================================================
   Internal Helper Functions
   =========================================================================== */

/*
    Data structure initialization helper.
    Allocates vars array and sort_order with proper alignment.
    Returns STATA_OK on success, error code on failure.
*/
static stata_retcode init_data_structure(stata_data *data, size_t nvars, size_t nobs)
{
    /* Initialize the data structure */
    stata_data_init(data);
    data->nobs = nobs;
    data->nvars = nvars;

    /* Allocate array of variables (cache-line aligned, overflow-safe) */
    data->vars = (stata_variable *)ctools_safe_cacheline_alloc2(nvars, sizeof(stata_variable));
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }
    memset(data->vars, 0, nvars * sizeof(stata_variable));

    /* Allocate sort order array (cache-line aligned, overflow-safe) */
    data->sort_order = (perm_idx_t *)ctools_safe_cacheline_alloc2(nobs, sizeof(perm_idx_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation (parallel for large datasets) */
    #pragma omp parallel for schedule(static) if(nobs >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < nobs; i++) {
        data->sort_order[i] = (perm_idx_t)i;
    }

    return STATA_OK;
}

/*
    I/O mode enum for thread args initialization.
*/
typedef enum { IO_MODE_LOAD, IO_MODE_STORE } io_mode_t;

/*
    Thread args batch initialization helper.
    Initializes an array of ctools_var_io_args for load or store operations.

    @param args         [out] Pre-allocated array of nvars args
    @param data         [in]  stata_data structure
    @param var_indices  [in]  Array of 1-based variable indices, or NULL for sequential
    @param nvars        [in]  Number of variables
    @param obs1         [in]  First observation (1-based)
    @param nobs         [in]  Number of observations
    @param mode         [in]  IO_MODE_LOAD or IO_MODE_STORE
*/
static void init_io_thread_args(ctools_var_io_args *args, stata_data *data,
                                int *var_indices, size_t nvars,
                                size_t obs1, size_t nobs, io_mode_t mode)
{
    for (size_t j = 0; j < nvars; j++) {
        args[j].var = &data->vars[j];
        args[j].var_idx = var_indices ? var_indices[j] : (int)(j + 1);
        args[j].obs1 = obs1;
        args[j].nobs = nobs;
        args[j].is_string = (mode == IO_MODE_LOAD)
            ? SF_var_is_string((ST_int)args[j].var_idx)
            : 0;  /* Not used by store */
        args[j].success = 0;
    }
}

/*
    Thread function type for I/O operations.
*/
typedef void *(*io_thread_func)(void *);

/*
    Thread pool orchestration helper.
    Executes I/O operations either in parallel (using thread pool) or sequentially.

    @param args         [in]  Array of thread arguments
    @param nvars        [in]  Number of variables (and args)
    @param func         [in]  Thread function to execute
    @param check_error  [in]  Non-zero to check for errors in load mode

    @return 0 on success, non-zero on failure
*/
static int execute_io_parallel(ctools_var_io_args *args, size_t nvars,
                               io_thread_func func, int check_error)
{
    int use_parallel = (nvars >= 2);

    if (use_parallel) {
        ctools_persistent_pool *pool = ctools_get_global_pool();

        if (pool != NULL) {
            /* Submit batch to thread pool */
            if (ctools_persistent_pool_submit_batch(pool, func, args, nvars,
                                                     sizeof(ctools_var_io_args)) != 0) {
                return -1;
            }

            int pool_result = ctools_persistent_pool_wait(pool);
            if (check_error && pool_result != 0) {
                return -1;
            }
        } else {
            /* Fallback to sequential */
            for (size_t j = 0; j < nvars; j++) {
                void *result = func(&args[j]);
                if (check_error && result != NULL) {
                    return -1;
                }
            }
        }
    } else {
        /* Sequential execution */
        for (size_t j = 0; j < nvars; j++) {
            void *result = func(&args[j]);
            if (check_error && result != NULL) {
                return -1;
            }
        }
    }

    return 0;
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
    ctools_var_io_args *thread_args = NULL;

    if (data == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    size_t nobs = data->nobs;
    size_t nvars = data->nvars;

    /* Allocate thread arguments */
    thread_args = (ctools_var_io_args *)ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments */
    init_io_thread_args(thread_args, data, NULL, nvars, obs1, nobs, IO_MODE_STORE);

    /* Execute parallel/sequential store (no error checking for store) */
    execute_io_parallel(thread_args, nvars, store_variable_thread, 0);

    free(thread_args);

    /* Memory barrier to ensure all thread-side effects are visible */
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
    ctools_var_io_args *thread_args = NULL;

    if (data == NULL || var_indices == NULL || nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    size_t nobs = data->nobs;

    /* Allocate thread arguments */
    thread_args = (ctools_var_io_args *)ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
    if (thread_args == NULL) {
        return STATA_ERR_MEMORY;
    }

    /* Initialize thread arguments with specified variable indices */
    init_io_thread_args(thread_args, data, var_indices, nvars, obs1, nobs, IO_MODE_STORE);

    /* Execute parallel/sequential store (no error checking for store) */
    execute_io_parallel(thread_args, nvars, store_variable_thread, 0);

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
    ST_int nobs_max = SF_nobs();

    if (is_string) {
        /* String variable: allocate buffer, gather, scatter (overflow-safe) */
        char **buf = (char **)ctools_safe_cacheline_alloc2(output_nobs, sizeof(char *));
        if (!buf) return STATA_ERR_MEMORY;
        memset(buf, 0, output_nobs * sizeof(char *));  /* Zero-init for safe cleanup */

        char strbuf[2048];

        /* Create arena for string storage (estimate 64 bytes avg per string) */
        /* Overflow check: output_nobs * 64 - skip arena if overflow, rely on strdup fallback */
        ctools_string_arena *arena = NULL;
        if (output_nobs <= SIZE_MAX / 64) {
            size_t arena_capacity = output_nobs * 64;
            arena = ctools_string_arena_create(arena_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        /* Arena failure is ok - we fall back to strdup */

        /* GATHER: Read from source positions with prefetching */
        for (i = 0; i < output_nobs; i++) {
            /* Prefetch source_rows array ahead */
            if (i + PREFETCH_DISTANCE < output_nobs) {
                _mm_prefetch((const char *)&source_rows[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }

            if (source_rows[i] >= 0) {
                /* Overflow check: ensure addition doesn't overflow int64_t */
                if (source_rows[i] > INT64_MAX - (int64_t)obs1) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld + obs1=%zu would overflow\n",
                             i, (long long)source_rows[i], obs1);
                    SF_error(errbuf);
                    /* Cleanup */
                    for (size_t j = 0; j < i; j++) {
                        if (buf[j] != NULL && !ctools_string_arena_owns(arena, buf[j])) {
                            free(buf[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                /* Bounds check: source observation must be within valid Stata range */
                int64_t src_obs = source_rows[i] + (int64_t)obs1;
                if (src_obs < 1 || src_obs > nobs_max) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld -> obs %lld out of bounds [1,%d]\n",
                             i, (long long)source_rows[i], (long long)src_obs, (int)nobs_max);
                    SF_error(errbuf);
                    /* Cleanup */
                    for (size_t j = 0; j < i; j++) {
                        if (buf[j] != NULL && !ctools_string_arena_owns(arena, buf[j])) {
                            free(buf[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                SF_sdata(stata_var, (ST_int)src_obs, strbuf);
                buf[i] = ctools_string_arena_strdup(arena, strbuf);
            } else {
                buf[i] = ctools_string_arena_strdup(arena, "");  /* Missing -> empty string */
            }
            if (!buf[i]) {
                /* Cleanup on allocation failure - free any fallback strings first */
                for (size_t j = 0; j < i; j++) {
                    if (buf[j] != NULL && !ctools_string_arena_owns(arena, buf[j])) {
                        free(buf[j]);
                    }
                }
                ctools_string_arena_free(arena);
                ctools_aligned_free(buf);
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

        /* Free fallback strings (those not owned by arena), then free arena */
        for (i = 0; i < output_nobs; i++) {
            if (buf[i] != NULL && !ctools_string_arena_owns(arena, buf[i])) {
                free(buf[i]);
            }
        }
        ctools_string_arena_free(arena);
        ctools_aligned_free(buf);

    } else {
        /* Numeric variable: allocate aligned buffer, gather, scatter (overflow-safe) */
        double *buf = (double *)ctools_safe_cacheline_alloc2(output_nobs, sizeof(double));
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
                /* Overflow check: ensure addition doesn't overflow int64_t */
                if (source_rows[i] > INT64_MAX - (int64_t)obs1) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld + obs1=%zu would overflow\n",
                             i, (long long)source_rows[i], obs1);
                    SF_error(errbuf);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                /* Bounds check: source observation must be within valid Stata range */
                int64_t src_obs = source_rows[i] + (int64_t)obs1;
                if (src_obs < 1 || src_obs > nobs_max) {
                    char errbuf[128];
                    snprintf(errbuf, sizeof(errbuf),
                             "ctools: source_rows[%zu]=%lld -> obs %lld out of bounds [1,%d]\n",
                             i, (long long)source_rows[i], (long long)src_obs, (int)nobs_max);
                    SF_error(errbuf);
                    ctools_aligned_free(buf);
                    return STATA_ERR_INVALID_INPUT;
                }
                SF_vdata(stata_var, (ST_int)src_obs, &buf[i]);
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

        ctools_aligned_free(buf);
    }

    return STATA_OK;
}

/* ===========================================================================
   Filtered Data Loading (Stata -> C with if/in filtering)
   =========================================================================== */

/*
    Initialize filtered data structure to safe empty state.
*/
void ctools_filtered_data_init(ctools_filtered_data *fd)
{
    if (fd == NULL) return;
    stata_data_init(&fd->data);
    fd->obs_map = NULL;
    fd->n_range = 0;
    fd->was_filtered = 0;
}

/*
    Free all memory associated with filtered data.
*/
void ctools_filtered_data_free(ctools_filtered_data *fd)
{
    if (fd == NULL) return;
    stata_data_free(&fd->data);
    if (fd->obs_map != NULL) {
        ctools_aligned_free(fd->obs_map);
        fd->obs_map = NULL;
    }
    fd->n_range = 0;
    fd->was_filtered = 0;
}

/*
    Quick identity detection: probe first and last observations.
    Returns 1 if likely no filtering needed, 0 if filtering detected.
*/
static int probe_for_identity(ST_int obs1, ST_int obs2, size_t n_range)
{
    /* Probe size: 16 at each end, or full range if small */
    size_t probe_size = 16;
    if (n_range <= 32) {
        /* Small range: check all observations */
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (!SF_ifobs(obs)) {
                return 0;  /* Filtering detected */
            }
        }
        return 1;  /* Identity confirmed */
    }

    /* Probe first 16 observations */
    for (size_t i = 0; i < probe_size; i++) {
        if (!SF_ifobs(obs1 + (ST_int)i)) {
            return 0;  /* Filtering detected */
        }
    }

    /* Probe last 16 observations */
    for (size_t i = 0; i < probe_size; i++) {
        if (!SF_ifobs(obs2 - (ST_int)i)) {
            return 0;  /* Filtering detected */
        }
    }

    return 1;  /* Likely identity (will verify during full pass if needed) */
}

/*
    Load a single variable from Stata for filtered observations only.
    Uses obs_map to read only the observations that passed filtering.
*/
static int load_filtered_variable(stata_variable *var, int var_idx,
                                   perm_idx_t *obs_map, size_t n_filtered,
                                   int is_string)
{
    size_t i;
    char strbuf[STATA_STR_MAXLEN + 1];

    var->nobs = n_filtered;

    /* Edge case: empty filtered set */
    if (n_filtered == 0) {
        if (is_string) {
            var->type = STATA_TYPE_STRING;
            var->str_maxlen = 0;
            var->_arena = NULL;
            var->data.str = (char **)calloc(1, sizeof(char *));
            if (!var->data.str) return -1;
        } else {
            var->type = STATA_TYPE_DOUBLE;
            var->_arena = NULL;
            var->data.dbl = (double *)ctools_cacheline_alloc(sizeof(double));
            if (!var->data.dbl) return -1;
            var->data.dbl[0] = SV_missval;
        }
        return 0;
    }

    if (is_string) {
        /* String variable - use arena allocator */
        var->type = STATA_TYPE_STRING;
        var->str_maxlen = STATA_STR_MAXLEN;
        var->_arena = NULL;

        /* Allocate pointer array */
        size_t str_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(char *), &str_array_size) != 0) {
            return -1;
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            return -1;
        }
        memset(var->data.str, 0, str_array_size);

        char **str_ptr = var->data.str;

        /* Create arena for string storage */
        ctools_string_arena *arena = NULL;
        if (n_filtered <= SIZE_MAX / 64) {
            size_t arena_capacity = n_filtered * 64;
            arena = ctools_string_arena_create(arena_capacity, CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        /* Load strings using obs_map */
        for (i = 0; i < n_filtered; i++) {
            SF_sdata((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
            str_ptr[i] = ctools_string_arena_strdup(arena, strbuf);
            if (str_ptr[i] == NULL) {
                /* Cleanup on failure */
                if (arena != NULL) {
                    for (size_t j = 0; j < i; j++) {
                        if (str_ptr[j] != NULL && !ctools_string_arena_owns(arena, str_ptr[j])) {
                            free(str_ptr[j]);
                        }
                    }
                    ctools_string_arena_free(arena);
                    var->_arena = NULL;
                } else {
                    for (size_t j = 0; j < i; j++) {
                        free(str_ptr[j]);
                    }
                }
                ctools_aligned_free(var->data.str);
                var->data.str = NULL;
                return -1;
            }
        }
    } else {
        /* Numeric variable */
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;

        size_t dbl_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(double), &dbl_array_size) != 0) {
            return -1;
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);
        if (var->data.dbl == NULL) {
            return -1;
        }

        double * restrict dbl_ptr = var->data.dbl;

        /* Load using obs_map - obs_map contains 1-based Stata obs numbers */
        for (i = 0; i < n_filtered; i++) {
            SF_vdata((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i]);
        }
    }

    return 0;
}

/*
    Thread function for filtered variable loading.
*/
typedef struct {
    stata_variable *var;
    int var_idx;
    perm_idx_t *obs_map;
    size_t n_filtered;
    int is_string;
    int success;
} filtered_var_io_args;

static void *load_filtered_variable_thread(void *arg)
{
    filtered_var_io_args *args = (filtered_var_io_args *)arg;
    args->success = 0;

    if (load_filtered_variable(args->var, args->var_idx, args->obs_map,
                                args->n_filtered, args->is_string) != 0) {
        return (void *)1;
    }

    args->success = 1;
    return NULL;
}

/*
    Load variables with if/in filtering applied at load time.
*/
stata_retcode ctools_data_load(ctools_filtered_data *result,
                                         int *var_indices, size_t nvars,
                                         size_t obs_start, size_t obs_end,
                                         int flags)
{
    ST_int obs1, obs2;
    size_t n_range, n_filtered;
    perm_idx_t *obs_map = NULL;
    int is_identity = 0;
    int skip_if_check = (flags & CTOOLS_LOAD_SKIP_IF) != 0;
    int *auto_indices = NULL;  /* For "load all" mode when var_indices is NULL */

    if (result == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Initialize result structure */
    ctools_filtered_data_init(result);

    /*
        "Load all variables" mode: if var_indices is NULL, build indices for all variables.
        This replaces the old ctools_data_load() function.
    */
    if (var_indices == NULL) {
        nvars = (size_t)SF_nvars();
        if (nvars == 0) {
            return STATA_OK;  /* No variables to load */
        }
        auto_indices = (int *)malloc(nvars * sizeof(int));
        if (auto_indices == NULL) {
            return STATA_ERR_MEMORY;
        }
        for (size_t i = 0; i < nvars; i++) {
            auto_indices[i] = (int)(i + 1);  /* 1-based Stata indices */
        }
        var_indices = auto_indices;
    } else if (nvars == 0) {
        return STATA_ERR_INVALID_INPUT;
    }

    /* Resolve observation range */
    if (obs_start == 0 || obs_end == 0) {
        obs1 = SF_in1();
        obs2 = SF_in2();
        if (obs1 < 1 || obs2 < obs1) {
            result->n_range = 0;
            result->was_filtered = 0;
            if (auto_indices) free(auto_indices);
            return STATA_OK;  /* Empty but valid */
        }
    } else {
        obs1 = (ST_int)obs_start;
        obs2 = (ST_int)obs_end;
    }

    n_range = (size_t)(obs2 - obs1 + 1);
    result->n_range = n_range;

    /*
        Fast path: skip SF_ifobs checks entirely when CTOOLS_LOAD_SKIP_IF is set.
        This saves N SPI calls when there's no `if` condition.
    */
    if (skip_if_check) {
        is_identity = 1;
        n_filtered = n_range;
        result->was_filtered = 0;

        /* Build identity obs_map */
        obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
        if (obs_map == NULL) {
            if (auto_indices) free(auto_indices);
            return STATA_ERR_MEMORY;
        }
        for (size_t i = 0; i < n_filtered; i++) {
            obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
        }
        result->obs_map = obs_map;
    }
    /* Quick identity detection for small ranges */
    else if (n_range <= 32) {
        is_identity = probe_for_identity(obs1, obs2, n_range);
        if (is_identity) {
            /* Small range confirmed as identity - use direct load */
            result->was_filtered = 0;
            n_filtered = n_range;

            /* Allocate obs_map as identity for write-back consistency */
            obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
            if (obs_map == NULL) {
                if (auto_indices) free(auto_indices);
                return STATA_ERR_MEMORY;
            }
            for (size_t i = 0; i < n_filtered; i++) {
                obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
            }
            result->obs_map = obs_map;
        }
    }

    /*
        Standard path: check SF_ifobs for each observation.
        Two-pass approach: count filtered, then build obs_map.
    */
    if (!skip_if_check && !is_identity) {
        /* Pass 1: Count filtered observations */
        n_filtered = 0;
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (SF_ifobs(obs)) {
                n_filtered++;
            }
        }

        if (n_filtered == 0) {
            /* No observations pass filter */
            result->was_filtered = 1;
            result->obs_map = NULL;
            /* Initialize empty data structure */
            stata_retcode rc = init_data_structure(&result->data, nvars, 0);
            if (auto_indices) free(auto_indices);
            return rc;
        }

        /* Check if identity (all passed) */
        is_identity = (n_filtered == n_range);
        result->was_filtered = !is_identity;

        /* Allocate obs_map */
        obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
        if (obs_map == NULL) {
            if (auto_indices) free(auto_indices);
            return STATA_ERR_MEMORY;
        }

        /* Pass 2: build obs_map */
        size_t idx = 0;
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (SF_ifobs(obs)) {
                obs_map[idx++] = (perm_idx_t)obs;
            }
        }

        result->obs_map = obs_map;
    }

    /* Load variables - use different paths for identity vs filtered cases */
    stata_retcode rc;

    if (is_identity && n_filtered > 0) {
        /*
            Identity case: all observations pass, use standard contiguous load.
            This is faster because it can use sequential memory access patterns.
        */
        rc = init_data_structure(&result->data, nvars, n_filtered);
        if (rc != STATA_OK) {
            ctools_aligned_free(obs_map);
            result->obs_map = NULL;
            if (auto_indices) free(auto_indices);
            return rc;
        }

        /* Allocate and initialize thread arguments for contiguous load */
        ctools_var_io_args *thread_args = (ctools_var_io_args *)
            ctools_safe_malloc2(nvars, sizeof(ctools_var_io_args));
        if (thread_args == NULL) {
            ctools_filtered_data_free(result);
            if (auto_indices) free(auto_indices);
            return STATA_ERR_MEMORY;
        }

        init_io_thread_args(thread_args, &result->data, var_indices, nvars, obs1, n_filtered, IO_MODE_LOAD);

        if (execute_io_parallel(thread_args, nvars, load_variable_thread, 1) != 0) {
            free(thread_args);
            ctools_filtered_data_free(result);
            if (auto_indices) free(auto_indices);
            return STATA_ERR_MEMORY;
        }
        free(thread_args);
    } else if (n_filtered > 0) {
        /*
            Filtered case: only some observations pass, use obs_map for gather.
            This uses random access but only loads the filtered observations.
        */
        rc = init_data_structure(&result->data, nvars, n_filtered);
        if (rc != STATA_OK) {
            ctools_aligned_free(obs_map);
            result->obs_map = NULL;
            if (auto_indices) free(auto_indices);
            return rc;
        }

        filtered_var_io_args *thread_args = (filtered_var_io_args *)
            ctools_safe_malloc2(nvars, sizeof(filtered_var_io_args));
        if (thread_args == NULL) {
            ctools_filtered_data_free(result);
            if (auto_indices) free(auto_indices);
            return STATA_ERR_MEMORY;
        }

        /* Initialize thread arguments */
        for (size_t j = 0; j < nvars; j++) {
            thread_args[j].var = &result->data.vars[j];
            thread_args[j].var_idx = var_indices[j];
            thread_args[j].obs_map = obs_map;
            thread_args[j].n_filtered = n_filtered;
            thread_args[j].is_string = SF_var_is_string((ST_int)var_indices[j]);
            thread_args[j].success = 0;
        }

        /* Execute load (parallel for multiple vars) */
        int use_parallel = (nvars >= 2);

        if (use_parallel) {
            ctools_persistent_pool *pool = ctools_get_global_pool();

            if (pool != NULL) {
                if (ctools_persistent_pool_submit_batch(pool, load_filtered_variable_thread,
                                                         thread_args, nvars,
                                                         sizeof(filtered_var_io_args)) != 0) {
                    free(thread_args);
                    ctools_filtered_data_free(result);
                    if (auto_indices) free(auto_indices);
                    return STATA_ERR_MEMORY;
                }

                int pool_result = ctools_persistent_pool_wait(pool);
                if (pool_result != 0) {
                    free(thread_args);
                    ctools_filtered_data_free(result);
                    if (auto_indices) free(auto_indices);
                    return STATA_ERR_MEMORY;
                }
            } else {
                /* Fallback to sequential */
                for (size_t j = 0; j < nvars; j++) {
                    void *res = load_filtered_variable_thread(&thread_args[j]);
                    if (res != NULL) {
                        free(thread_args);
                        ctools_filtered_data_free(result);
                        if (auto_indices) free(auto_indices);
                        return STATA_ERR_MEMORY;
                    }
                }
            }
        } else {
            /* Sequential execution */
            for (size_t j = 0; j < nvars; j++) {
                void *res = load_filtered_variable_thread(&thread_args[j]);
                if (res != NULL) {
                    free(thread_args);
                    ctools_filtered_data_free(result);
                    if (auto_indices) free(auto_indices);
                    return STATA_ERR_MEMORY;
                }
            }
        }

        free(thread_args);
    }

    /* Memory barrier */
    ctools_memory_barrier();

    if (auto_indices) free(auto_indices);
    return STATA_OK;
}

/* ===========================================================================
   Row-Parallel Single-Variable Loader (experimental)
   =========================================================================== */

/*
    Load a single variable using row-parallel OpenMP threads.
    For single-variable loads (the common cdestring case), the standard
    ctools_data_load() is entirely sequential since it parallelizes across
    columns. This function splits observations across threads instead.

    Strings use strdup (no arena) since arenas aren't thread-safe.
    The _arena field is set to NULL so stata_data_free uses the per-string
    free path.
*/
stata_retcode ctools_data_load_single_var_rowpar(
    ctools_filtered_data *result,
    int var_idx,
    size_t obs_start,
    size_t obs_end,
    int flags)
{
    ST_int obs1, obs2;
    size_t n_range, n_filtered;
    perm_idx_t *obs_map = NULL;
    int is_identity = 0;
    int skip_if_check = (flags & CTOOLS_LOAD_SKIP_IF) != 0;
    int is_string;

    if (result == NULL) {
        return STATA_ERR_INVALID_INPUT;
    }

    ctools_filtered_data_init(result);

    /* Resolve observation range */
    obs1 = (obs_start > 0) ? (ST_int)obs_start : SF_in1();
    obs2 = (obs_end > 0)   ? (ST_int)obs_end   : SF_in2();
    if (obs1 < 1) obs1 = 1;
    if (obs2 < obs1) {
        /* Empty range */
        return init_data_structure(&result->data, 1, 0);
    }

    n_range = (size_t)(obs2 - obs1 + 1);
    result->n_range = n_range;

    is_string = SF_var_is_string((ST_int)var_idx);

    /* --- Filter pass (sequential — requires SPI) --- */
    if (skip_if_check) {
        is_identity = 1;
        n_filtered = n_range;
    } else {
        /* Pass 1: count */
        n_filtered = 0;
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (SF_ifobs(obs)) n_filtered++;
        }
        is_identity = (n_filtered == n_range);
    }

    if (n_filtered == 0) {
        return init_data_structure(&result->data, 1, 0);
    }

    /* Build obs_map */
    if (is_identity) {
        obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
        if (obs_map == NULL) return STATA_ERR_MEMORY;
        for (size_t i = 0; i < n_filtered; i++) {
            obs_map[i] = (perm_idx_t)(obs1 + (ST_int)i);
        }
    } else {
        obs_map = (perm_idx_t *)ctools_safe_cacheline_alloc2(n_filtered, sizeof(perm_idx_t));
        if (obs_map == NULL) return STATA_ERR_MEMORY;
        size_t idx = 0;
        for (ST_int obs = obs1; obs <= obs2; obs++) {
            if (SF_ifobs(obs)) obs_map[idx++] = (perm_idx_t)obs;
        }
    }

    result->obs_map = obs_map;
    result->was_filtered = !is_identity;

    /* --- Allocate data structure for 1 variable --- */
    stata_retcode rc = init_data_structure(&result->data, 1, n_filtered);
    if (rc != STATA_OK) {
        ctools_aligned_free(obs_map);
        result->obs_map = NULL;
        return rc;
    }

    stata_variable *var = &result->data.vars[0];
    var->nobs = n_filtered;

    /* --- Row-parallel load --- */
    if (is_string) {
        var->type = STATA_TYPE_STRING;
        var->str_maxlen = STATA_STR_MAXLEN;
        var->_arena = NULL;

        size_t str_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(char *), &str_array_size) != 0) {
            ctools_filtered_data_free(result);
            return STATA_ERR_MEMORY;
        }
        var->data.str = (char **)ctools_cacheline_alloc(str_array_size);
        if (var->data.str == NULL) {
            ctools_filtered_data_free(result);
            return STATA_ERR_MEMORY;
        }
        memset(var->data.str, 0, str_array_size);

        /* Create arena and partition into per-thread slabs so each thread
           does plain pointer bumping with no atomic contention. */
        ctools_string_arena *arena = NULL;
        if (n_filtered <= SIZE_MAX / 64) {
            size_t arena_capacity = n_filtered * 64;
            arena = ctools_string_arena_create(arena_capacity,
                                               CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
        }
        if (arena != NULL) {
            var->_arena = arena;
        }

        volatile int load_error = 0;

        #pragma omp parallel if(n_filtered >= MIN_OBS_PER_THREAD * 2)
        {
            /* Partition arena into per-thread slabs */
            char *slab_ptr = NULL;
            char *slab_end = NULL;
            if (arena != NULL) {
                int nthreads = 1;
                int tid = 0;
                #ifdef _OPENMP
                nthreads = omp_get_num_threads();
                tid = omp_get_thread_num();
                #endif
                size_t slab_size = arena->capacity / (size_t)nthreads;
                slab_ptr = arena->base + (size_t)tid * slab_size;
                slab_end = (tid == nthreads - 1)
                    ? arena->base + arena->capacity
                    : slab_ptr + slab_size;
            }

            #pragma omp for schedule(static)
            for (size_t i = 0; i < n_filtered; i++) {
                if (load_error) continue;
                char strbuf[STATA_STR_MAXLEN + 1];
                SF_sdata((ST_int)var_idx, (ST_int)obs_map[i], strbuf);
                size_t len = strlen(strbuf) + 1;

                char *s;
                if (slab_ptr != NULL && slab_ptr + len <= slab_end) {
                    memcpy(slab_ptr, strbuf, len);
                    s = slab_ptr;
                    slab_ptr += len;
                } else {
                    /* Slab full or no arena — fall back to strdup */
                    s = strdup(strbuf);
                    if (arena != NULL) arena->has_fallback = 1;
                    if (s == NULL) {
                        load_error = 1;
                        continue;
                    }
                }
                var->data.str[i] = s;
            }

            /* Update arena->used so arena_owns() covers all slabs */
            #pragma omp single
            {
                if (arena != NULL) arena->used = arena->capacity;
            }
        }

        if (load_error) {
            ctools_filtered_data_free(result);
            return STATA_ERR_MEMORY;
        }
    } else {
        var->type = STATA_TYPE_DOUBLE;
        var->_arena = NULL;

        size_t dbl_array_size;
        if (ctools_safe_mul_size(n_filtered, sizeof(double), &dbl_array_size) != 0) {
            ctools_filtered_data_free(result);
            return STATA_ERR_MEMORY;
        }
        var->data.dbl = (double *)ctools_cacheline_alloc(dbl_array_size);
        if (var->data.dbl == NULL) {
            ctools_filtered_data_free(result);
            return STATA_ERR_MEMORY;
        }

        double * restrict dbl_ptr = var->data.dbl;

        #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
        for (size_t i = 0; i < n_filtered; i++) {
            SF_vdata((ST_int)var_idx, (ST_int)obs_map[i], &dbl_ptr[i]);
        }
    }

    ctools_memory_barrier();
    return STATA_OK;
}

/* ===========================================================================
   Filtered Data Storing (C -> Stata with obs_map)
   =========================================================================== */

/*
    Write filtered numeric values back to Stata using obs_map.
*/
stata_retcode ctools_store_filtered(double *values, size_t n_filtered,
                                     int var_idx, perm_idx_t *obs_map)
{
    if (values == NULL || obs_map == NULL || n_filtered == 0) {
        return (values == NULL && n_filtered == 0) ? STATA_OK : STATA_ERR_INVALID_INPUT;
    }

    ST_int stata_var = (ST_int)var_idx;
    ST_int nobs_max = SF_nobs();

    /* Write values using obs_map for indexing */
    for (size_t i = 0; i < n_filtered; i++) {
        /* Bounds check: obs_map values must be valid 1-based Stata observation indices */
        perm_idx_t obs = obs_map[i];
        if (obs < 1 || obs > (perm_idx_t)nobs_max) {
            char buf[128];
            snprintf(buf, sizeof(buf),
                     "ctools: obs_map[%zu]=%u out of bounds [1,%d]\n",
                     i, (unsigned)obs, (int)nobs_max);
            SF_error(buf);
            return STATA_ERR_INVALID_INPUT;
        }
        SF_vstore(stata_var, (ST_int)obs, values[i]);
    }

    return STATA_OK;
}

/*
    Row-parallel store for a single numeric variable.
    Counterpart to ctools_data_load_single_var_rowpar().
    Skips per-element bounds checks (obs_map is trusted from our own loader).
*/
stata_retcode ctools_store_filtered_rowpar(double *values, size_t n_filtered,
                                            int var_idx, perm_idx_t *obs_map)
{
    if (values == NULL || obs_map == NULL || n_filtered == 0) {
        return (values == NULL && n_filtered == 0) ? STATA_OK : STATA_ERR_INVALID_INPUT;
    }

    ST_int stata_var = (ST_int)var_idx;

    #pragma omp parallel for schedule(static) if(n_filtered >= MIN_OBS_PER_THREAD * 2)
    for (size_t i = 0; i < n_filtered; i++) {
        SF_vstore(stata_var, (ST_int)obs_map[i], values[i]);
    }

    return STATA_OK;
}

/*
    Write filtered string values back to Stata using obs_map.
*/
stata_retcode ctools_store_filtered_str(char **strings, size_t n_filtered,
                                         int var_idx, perm_idx_t *obs_map)
{
    if (strings == NULL || obs_map == NULL || n_filtered == 0) {
        return (strings == NULL && n_filtered == 0) ? STATA_OK : STATA_ERR_INVALID_INPUT;
    }

    ST_int stata_var = (ST_int)var_idx;
    ST_int nobs_max = SF_nobs();

    /* Write strings using obs_map for indexing */
    for (size_t i = 0; i < n_filtered; i++) {
        /* Bounds check: obs_map values must be valid 1-based Stata observation indices */
        perm_idx_t obs = obs_map[i];
        if (obs < 1 || obs > (perm_idx_t)nobs_max) {
            char buf[128];
            snprintf(buf, sizeof(buf),
                     "ctools: obs_map[%zu]=%u out of bounds [1,%d]\n",
                     i, (unsigned)obs, (int)nobs_max);
            SF_error(buf);
            return STATA_ERR_INVALID_INPUT;
        }
        SF_sstore(stata_var, (ST_int)obs, strings[i] ? strings[i] : "");
    }

    return STATA_OK;
}
