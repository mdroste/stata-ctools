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

    Performance Optimizations:
    - 8x loop unrolling reduces loop overhead and improves instruction pipelining
    - Parallel variable I/O overlaps operations across columns
    - SD_FASTMODE (compile flag) disables SPI bounds checking
    - Cache-line aligned allocations for optimal memory access
    - Software prefetching for improved memory latency hiding
    - SIMD operations for bulk memory initialization
    - Non-temporal stores for large write operations
    - String arena allocator for reduced malloc overhead

    Thread Safety:
    - Each thread operates on a different Stata variable (column)
    - No shared state during I/O phase
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
    return arena;
}

/*
    Allocate a string copy from the arena.
    Falls back to strdup if arena is full.
    Returns NULL on complete failure.
*/
static char *arena_strdup(string_arena *arena, const char *s)
{
    size_t len = strlen(s) + 1;

    /* Check if arena has space */
    if (arena != NULL && arena->used + len <= arena->capacity) {
        char *ptr = arena->base + arena->used;
        memcpy(ptr, s, len);
        arena->used += len;
        return ptr;
    }

    /* Fallback to regular strdup */
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
    Initialize an array with identity permutation [0, 1, 2, ..., n-1]
    using SIMD acceleration where available.
*/
static void init_identity_permutation(size_t *arr, size_t n)
{
    size_t i = 0;

#if defined(CTOOLS_AVX2)
    /* AVX2: Initialize 4 size_t values at a time (256 bits = 4 x 64-bit) */
    __m256i increment = _mm256_set1_epi64x(4);
    __m256i indices = _mm256_set_epi64x(3, 2, 1, 0);

    size_t vec_end = n - (n % 4);
    for (; i < vec_end; i += 4) {
        _mm256_storeu_si256((__m256i *)&arr[i], indices);
        indices = _mm256_add_epi64(indices, increment);
    }
#elif defined(CTOOLS_ARM64)
    /* NEON: Initialize 2 size_t values at a time (128 bits = 2 x 64-bit) */
    uint64x2_t increment = vdupq_n_u64(2);
    uint64x2_t indices = {0, 1};

    size_t vec_end = n - (n % 2);
    for (; i < vec_end; i += 2) {
        vst1q_u64((uint64_t *)&arr[i], indices);
        indices = vaddq_u64(indices, increment);
    }
#else
    /* Scalar fallback with 8x unrolling */
    size_t vec_end = n - (n % 8);
    for (; i < vec_end; i += 8) {
        arr[i]     = i;
        arr[i + 1] = i + 1;
        arr[i + 2] = i + 2;
        arr[i + 3] = i + 3;
        arr[i + 4] = i + 4;
        arr[i + 5] = i + 5;
        arr[i + 6] = i + 6;
        arr[i + 7] = i + 7;
    }
#endif

    /* Scalar cleanup for remaining elements */
    for (; i < n; i++) {
        arr[i] = i;
    }
}

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
    Thread function: Load a single variable from Stata into C memory.

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
    size_t i;
    size_t nobs = args->nobs;
    size_t obs1 = args->obs1;
    ST_int var_idx = (ST_int)args->var_idx;
    char strbuf[STATA_STR_MAXLEN + 1];

    args->var->nobs = nobs;
    args->success = 0;  /* Assume failure until complete */

    if (args->is_string) {
        /* String variable - use arena allocator for fast bulk free */
        args->var->type = STATA_TYPE_STRING;
        args->var->str_maxlen = STATA_STR_MAXLEN;
        args->var->_arena = NULL;

        /* Allocate pointer array (cache-line aligned) */
        args->var->data.str = (char **)aligned_alloc_cacheline(nobs * sizeof(char *));
        if (args->var->data.str == NULL) {
            return (void *)1;  /* Signal failure */
        }

        char **str_ptr = args->var->data.str;

        /* Create arena for all strings - estimate avg 64 bytes per string.
           Arena allows O(1) bulk free instead of O(n) individual frees. */
        size_t arena_capacity = nobs * 64;
        string_arena *arena = arena_create(arena_capacity);
        if (arena != NULL) {
            args->var->_arena = arena;  /* Store for later bulk free */
        }

        /* Load strings with prefetching */
        for (i = 0; i < nobs; i++) {
            /* Prefetch destination pointer array ahead */
            if (i + PREFETCH_DISTANCE < nobs) {
                _mm_prefetch((const char *)&str_ptr[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            }

            SF_sdata(var_idx, (ST_int)(i + obs1), strbuf);
            str_ptr[i] = arena_strdup(arena, strbuf);
            if (str_ptr[i] == NULL) {
                /* Arena full or allocation failed - this shouldn't happen with
                   proper arena sizing, but we handle it gracefully */
                return (void *)1;  /* Signal failure */
            }
        }

    } else {
        /* Numeric variable - use cache-line aligned allocation */
        args->var->type = STATA_TYPE_DOUBLE;
        args->var->_arena = NULL;  /* No arena for numeric data */
        args->var->data.dbl = (double *)aligned_alloc_cacheline(nobs * sizeof(double));

        if (args->var->data.dbl == NULL) {
            return (void *)1;  /* Signal failure */
        }

        double * restrict dbl_ptr = args->var->data.dbl;

        /*
         * Aggressive optimization strategy:
         * - 16x unrolling for better instruction-level parallelism
         * - Double buffering for software pipelining (overlap load/store)
         * - Branch-free prefetch by splitting into prefetch and non-prefetch regions
         * - Use restrict to help compiler optimize
         */

        /* Double buffers for software pipelining */
        double v0[16], v1[16];

        /* Calculate loop bounds - separate prefetch region from tail */
        size_t i_end_16 = nobs - (nobs % 16);
        size_t prefetch_end = (nobs > 32) ? (i_end_16 - 32) : 0;

        /* Phase 1: Main loop with prefetching (bulk of iterations) */
        i = 0;
        if (prefetch_end > 0) {
            /* Prime the pipeline: load first batch */
            SF_VDATA_BATCH16(var_idx, obs1, v0);

            for (; i < prefetch_end; i += 16) {
                /* Prefetch 2 cache lines ahead (32 doubles = 256 bytes) */
                _mm_prefetch((const char *)&dbl_ptr[i + 32], _MM_HINT_T0);
                _mm_prefetch((const char *)&dbl_ptr[i + 40], _MM_HINT_T0);

                /* Load next batch into v1 while we still have v0 */
                SF_VDATA_BATCH16(var_idx, i + 16 + obs1, v1);

                /* Store previous batch */
                copy_doubles_16(&dbl_ptr[i], v0);

                /* Swap buffers */
                double *tmp = (double *)memcpy(v0, v1, sizeof(v0));
                (void)tmp;
            }
        }

        /* Phase 2: Remaining 16-element blocks without prefetch */
        for (; i < i_end_16; i += 16) {
            SF_VDATA_BATCH16(var_idx, i + obs1, v0);
            copy_doubles_16(&dbl_ptr[i], v0);
        }

        /* Phase 3: Handle remaining elements (< 16) */
        for (; i < nobs; i++) {
            SF_vdata(var_idx, (ST_int)(i + obs1), &dbl_ptr[i]);
        }
    }

    args->success = 1;
    return NULL;  /* Success */
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

    /* Allocate sort order array (cache-line aligned) */
    data->sort_order = (size_t *)aligned_alloc_cacheline(nobs * sizeof(size_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation using SIMD */
    init_identity_permutation(data->sort_order, nobs);

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

    /* Decide whether to use parallel loading */
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
        /* Sequential loading - reuse thread function directly */
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
    Thread function: Store a single variable from C memory to Stata.

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
    size_t i;
    size_t nobs = args->nobs;
    size_t obs1 = args->obs1;
    ST_int var_idx = (ST_int)args->var_idx;

    if (args->var->type == STATA_TYPE_DOUBLE) {
        /* Numeric variable - aggressive optimization */
        const double * restrict dbl_data = args->var->data.dbl;

        /* Double buffers for software pipelining */
        double v0[16], v1[16];

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
                SF_VSTORE_BATCH16(var_idx, i + obs1, v0);

                /* Swap buffers */
                double *tmp = (double *)memcpy(v0, v1, sizeof(v0));
                (void)tmp;
            }
        }

        /* Phase 2: Remaining 16-element blocks without prefetch */
        for (; i < i_end_16; i += 16) {
            copy_doubles_16(v0, &dbl_data[i]);
            SF_VSTORE_BATCH16(var_idx, i + obs1, v0);
        }

        /* Phase 3: Handle remaining elements */
        for (; i < nobs; i++) {
            SF_vstore(var_idx, (ST_int)(i + obs1), dbl_data[i]);
        }

    } else {
        /* String variable - optimized with prefetching */
        char * const * restrict str_data = args->var->data.str;

        /* Calculate prefetch region */
        size_t prefetch_end = (nobs > PREFETCH_DISTANCE) ? (nobs - PREFETCH_DISTANCE) : 0;

        /* Phase 1: Main loop with prefetching */
        for (i = 0; i < prefetch_end; i++) {
            /* Prefetch pointer array and string content ahead */
            _mm_prefetch((const char *)&str_data[i + PREFETCH_DISTANCE], _MM_HINT_T0);
            _mm_prefetch(str_data[i + PREFETCH_DISTANCE], _MM_HINT_T0);

            SF_sstore(var_idx, (ST_int)(i + obs1), str_data[i]);
        }

        /* Phase 2: Tail without prefetch */
        for (; i < nobs; i++) {
            SF_sstore(var_idx, (ST_int)(i + obs1), str_data[i]);
        }
    }

    args->success = 1;
    return NULL;  /* Success */
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

    /* Decide whether to use parallel storing */
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
        /* Sequential storing - reuse thread function directly */
        for (j = 0; j < nvars; j++) {
            store_variable_thread(&thread_args[j]);
        }
    }

    free(thread_args);

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
    data->vars = (stata_variable *)aligned_alloc_cacheline(nvars * sizeof(stata_variable));
    if (data->vars == NULL) {
        return STATA_ERR_MEMORY;
    }
    memset(data->vars, 0, nvars * sizeof(stata_variable));

    /* Allocate sort order array (cache-line aligned) */
    data->sort_order = (size_t *)aligned_alloc_cacheline(nobs * sizeof(size_t));
    if (data->sort_order == NULL) {
        stata_data_free(data);
        return STATA_ERR_MEMORY;
    }

    /* Initialize sort order to identity permutation using SIMD */
    init_identity_permutation(data->sort_order, nobs);

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

    /* Decide whether to use parallel loading */
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
        /* Sequential loading */
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

        char strbuf[2048];

        /* Create arena for string storage (estimate 64 bytes avg per string) */
        string_arena *arena = arena_create(output_nobs * 64);
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
                /* Cleanup on allocation failure */
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

        /* Note: We don't free individual strings - they're in the arena
           (or if fallback strdup was used, we accept the leak for simplicity).
           In practice, for very large strings that overflow the arena,
           this could be improved with tracking. */
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
