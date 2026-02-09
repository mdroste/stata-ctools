/*
    ctools_types.c
    Implementation of common utility functions for Stata-C data structures
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_threads.h"
#include "ctools_arena.h"
#include "ctools_eisel_lemire.h"

/* SIMD intrinsics for gather operations */
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #define CTOOLS_TYPES_X86 1
    #include <immintrin.h>
    #ifdef __AVX2__
        #define CTOOLS_TYPES_AVX2 1
    #endif
#elif defined(__aarch64__) || defined(_M_ARM64)
    #define CTOOLS_TYPES_ARM64 1
    #include <arm_neon.h>
#endif

/*
    NOTE: Use ctools_aligned_free() from ctools_config.h for freeing aligned memory.
    This ensures consistent behavior across all platforms including Windows.
*/
#define aligned_free_internal ctools_aligned_free

void stata_data_init(stata_data *data)
{
    if (data == NULL) return;

    data->nobs = 0;
    data->nvars = 0;
    data->vars = NULL;
    data->sort_order = NULL;
}

void stata_data_free(stata_data *data)
{
    size_t i;

    if (data == NULL) return;

    /* Sanity check: if nvars > 0 but vars is NULL, we have corruption */
    if (data->nvars > 0 && data->vars == NULL) {
        /* This should never happen - indicates memory corruption */
        /* Reset to safe state and return without crashing */
        data->nvars = 0;
        data->nobs = 0;
        data->sort_order = NULL;
        return;
    }

    /* Free each variable's data (allocated with aligned_alloc_cacheline) */
    if (data->vars != NULL) {
        for (i = 0; i < data->nvars; i++) {
            if (data->vars[i].type == STATA_TYPE_DOUBLE) {
                /* Numeric data uses aligned allocation */
                aligned_free_internal(data->vars[i].data.dbl);
            } else if (data->vars[i].type == STATA_TYPE_STRING) {
                if (data->vars[i].data.str != NULL) {
                    ctools_string_arena *arena = (ctools_string_arena *)data->vars[i]._arena;

                    if (arena != NULL && !arena->has_fallback) {
                        /* Fast path: all strings are in the arena, O(1) bulk free */
                        ctools_string_arena_free(arena);
                        data->vars[i]._arena = NULL;
                    } else if (arena != NULL && arena->has_fallback) {
                        /* Mixed path: some strings in arena, some from strdup.
                           Free fallback strings individually (those outside arena range). */
                        for (size_t j = 0; j < data->vars[i].nobs; j++) {
                            char *str = data->vars[i].data.str[j];
                            if (str != NULL && !ctools_string_arena_owns(arena, str)) {
                                free(str);
                            }
                        }
                        ctools_string_arena_free(arena);
                        data->vars[i]._arena = NULL;
                    } else {
                        /* No arena: all strings are from strdup (or NULL).
                           Free each non-NULL pointer individually. */
                        for (size_t j = 0; j < data->vars[i].nobs; j++) {
                            if (data->vars[i].data.str[j] != NULL) {
                                free(data->vars[i].data.str[j]);
                            }
                        }
                    }
                    /* Free the pointer array (aligned allocation) */
                    aligned_free_internal(data->vars[i].data.str);
                }
            }
        }
        /* Free vars array (aligned) */
        aligned_free_internal(data->vars);
    }

    /* Free sort order array (aligned) */
    aligned_free_internal(data->sort_order);

    /* Reset the structure */
    stata_data_init(data);
}

/* ===========================================================================
   Type Conversion Utilities Implementation
   =========================================================================== */

/*
    Power of 10 lookup table for fast float parsing/formatting.
*/
const double ctools_pow10_table[23] = {
    1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,
    1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
    1e20, 1e21, 1e22
};

/*
    Two-digit lookup table for fast digit pair output.
    "00", "01", "02", ... "99"
    Exported for use by cexport and other modules.
*/
const char CTOOLS_DIGIT_PAIRS[200] = {
    '0','0', '0','1', '0','2', '0','3', '0','4', '0','5', '0','6', '0','7', '0','8', '0','9',
    '1','0', '1','1', '1','2', '1','3', '1','4', '1','5', '1','6', '1','7', '1','8', '1','9',
    '2','0', '2','1', '2','2', '2','3', '2','4', '2','5', '2','6', '2','7', '2','8', '2','9',
    '3','0', '3','1', '3','2', '3','3', '3','4', '3','5', '3','6', '3','7', '3','8', '3','9',
    '4','0', '4','1', '4','2', '4','3', '4','4', '4','5', '4','6', '4','7', '4','8', '4','9',
    '5','0', '5','1', '5','2', '5','3', '5','4', '5','5', '5','6', '5','7', '5','8', '5','9',
    '6','0', '6','1', '6','2', '6','3', '6','4', '6','5', '6','6', '6','7', '6','8', '6','9',
    '7','0', '7','1', '7','2', '7','3', '7','4', '7','5', '7','6', '7','7', '7','8', '7','9',
    '8','0', '8','1', '8','2', '8','3', '8','4', '8','5', '8','6', '8','7', '8','8', '8','9',
    '9','0', '9','1', '9','2', '9','3', '9','4', '9','5', '9','6', '9','7', '9','8', '9','9'
};

int ctools_uint64_to_str(uint64_t val, char *buf)
{
    if (val == 0) {
        buf[0] = '0';
        buf[1] = '\0';
        return 1;
    }

    /* Determine number of digits using binary search */
    int len;
    if (val < 10000000000ULL) {
        if (val < 100000) {
            if (val < 100) { len = (val < 10) ? 1 : 2; }
            else { len = (val < 1000) ? 3 : (val < 10000) ? 4 : 5; }
        } else {
            if (val < 10000000) { len = (val < 1000000) ? 6 : 7; }
            else { len = (val < 100000000) ? 8 : (val < 1000000000ULL) ? 9 : 10; }
        }
    } else {
        if (val < 1000000000000000ULL) {
            if (val < 1000000000000ULL) { len = (val < 100000000000ULL) ? 11 : 12; }
            else { len = (val < 10000000000000ULL) ? 13 : (val < 100000000000000ULL) ? 14 : 15; }
        } else {
            if (val < 100000000000000000ULL) { len = (val < 10000000000000000ULL) ? 16 : 17; }
            else { len = (val < 1000000000000000000ULL) ? 18 : (val < 10000000000000000000ULL) ? 19 : 20; }
        }
    }

    buf[len] = '\0';

    /* Write digits from right to left directly into output buffer */
    char *p = buf + len;
    while (val >= 100) {
        int idx = (int)(val % 100) * 2;
        val /= 100;
        *(--p) = CTOOLS_DIGIT_PAIRS[idx + 1];
        *(--p) = CTOOLS_DIGIT_PAIRS[idx];
    }
    if (val >= 10) {
        int idx = (int)val * 2;
        *(--p) = CTOOLS_DIGIT_PAIRS[idx + 1];
        *(--p) = CTOOLS_DIGIT_PAIRS[idx];
    } else {
        *(--p) = '0' + (int)val;
    }

    return len;
}

int ctools_int64_to_str(int64_t val, char *buf)
{
    if (val < 0) {
        buf[0] = '-';
        return 1 + ctools_uint64_to_str((uint64_t)(-val), buf + 1);
    }
    return ctools_uint64_to_str((uint64_t)val, buf);
}

int ctools_uint64_to_str_fast(uint64_t val, char *buf)
{
    /* Now identical to ctools_uint64_to_str which already uses digit pairs */
    return ctools_uint64_to_str(val, buf);
}

/* ===========================================================================
   Permutation Application - Optimized for Large Datasets
   =========================================================================== */

#ifdef _OPENMP
#include <omp.h>
#endif

/*
    Prefetch distances for permutation operations.
    Uses adaptive distances based on detected cache hierarchy.
    - NEAR: prefetch into L1 cache
    - FAR: prefetch into L2/L3 cache
*/
static int _perm_prefetch_near = 0;
static int _perm_prefetch_far = 0;

static inline void init_perm_prefetch_distances(void)
{
    if (_perm_prefetch_near == 0) {
        ctools_prefetch_distances dist = ctools_get_prefetch_distances();
        _perm_prefetch_near = dist.near_dist;
        _perm_prefetch_far = dist.far_dist;
    }
}

#define PERM_PREFETCH_NEAR  (_perm_prefetch_near ? _perm_prefetch_near : 8)
#define PERM_PREFETCH_FAR   (_perm_prefetch_far ? _perm_prefetch_far : 32)

/*
    Apply permutation using GATHER pattern with multi-level prefetching.
    new_data[i] = old_data[perm[i]] -- random read from old_data, sequential write

    GATHER is preferred over SCATTER for large datasets because:
    1. Sequential writes can use streaming stores (bypass cache)
    2. Random reads with prefetching are cheaper than random writes
    3. No overhead to compute inverse permutation

    NOTE: These functions do NOT use internal parallelism since they're called
    from an outer parallel region that parallelizes across variables.

    SIMD Optimization (AVX2):
    Uses _mm256_i32gather_pd to gather 4 doubles per instruction when available.
    Falls back to scalar loop for ARM64 and non-AVX2 x86.
*/

#ifdef CTOOLS_TYPES_AVX2
/*
    AVX2 SIMD gather: processes 4 doubles per iteration using hardware gather.
    The _mm256_i32gather_pd intrinsic gathers 4 doubles using 32-bit indices.
*/
static inline void apply_perm_gather_double_simd(double * CTOOLS_RESTRICT new_data,
                                                  const double * CTOOLS_RESTRICT old_data,
                                                  const perm_idx_t * CTOOLS_RESTRICT perm,
                                                  size_t nobs)
{
    size_t i;
    size_t nobs_vec4 = nobs & ~(size_t)3;  /* Round down to multiple of 4 */

    /* Main SIMD loop: gather 4 doubles at a time */
    for (i = 0; i < nobs_vec4; i += 4) {
        /* Prefetch ahead for next iterations */
        if (i + PERM_PREFETCH_FAR < nobs) {
            CTOOLS_PREFETCH(&perm[i + PERM_PREFETCH_FAR]);
        }

        /* Load 4 indices into 128-bit register */
        __m128i indices = _mm_loadu_si128((const __m128i *)&perm[i]);

        /* Gather 4 doubles using hardware gather instruction
           Scale = 8 because we're gathering doubles (8 bytes each) */
        __m256d gathered = _mm256_i32gather_pd(old_data, indices, 8);

        /* Store 4 doubles sequentially */
        _mm256_storeu_pd(&new_data[i], gathered);
    }

    /* Scalar tail: handle remaining elements */
    for (; i < nobs; i++) {
        new_data[i] = old_data[perm[i]];
    }
}
#endif /* CTOOLS_TYPES_AVX2 */

/*
    Scalar gather with multi-level prefetching.
    Used on ARM64 and non-AVX2 x86 platforms.
*/
#ifndef CTOOLS_TYPES_AVX2
static inline void apply_perm_gather_double_scalar(double * CTOOLS_RESTRICT new_data,
                                                    const double * CTOOLS_RESTRICT old_data,
                                                    const perm_idx_t * CTOOLS_RESTRICT perm,
                                                    size_t nobs)
{
    size_t i;
    /* Process with multi-level prefetching */
    for (i = 0; i < nobs; i++) {
        /* Far prefetch - into L2/L3 */
        if (i + PERM_PREFETCH_FAR < nobs) {
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_FAR]]);
        }
        /* Near prefetch - into L1 */
        if (i + PERM_PREFETCH_NEAR < nobs) {
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_NEAR]]);
        }
        new_data[i] = old_data[perm[i]];
    }
}
#endif /* !CTOOLS_TYPES_AVX2 */

/* ============================================================================
   Non-Temporal Store Variants for Large Datasets

   For very large datasets (> 1M elements), use non-temporal (streaming) stores
   to bypass the cache. This avoids cache pollution when writing data that
   won't be read again soon. The data goes directly to main memory.
   ============================================================================ */

/* Threshold for non-temporal stores: 1M doubles = 8MB */
#define PERM_NONTEMPORAL_THRESHOLD 1000000

#ifdef CTOOLS_TYPES_X86
/*
    x86 streaming gather: uses non-temporal stores to bypass cache.
    For very large datasets where new_data won't fit in cache.
*/
__attribute__((unused))
static void apply_perm_gather_double_streaming(double * CTOOLS_RESTRICT new_data,
                                                const double * CTOOLS_RESTRICT old_data,
                                                const perm_idx_t * CTOOLS_RESTRICT perm,
                                                size_t nobs)
{
    size_t i;
    size_t nobs_vec2 = nobs & ~(size_t)1;  /* Round down to multiple of 2 */

    /* Main loop: gather and stream-store 2 doubles at a time */
    for (i = 0; i < nobs_vec2; i += 2) {
        /* Prefetch source data */
        if (i + PERM_PREFETCH_FAR < nobs) {
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_FAR]]);
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_FAR + 1]]);
        }

        /* Gather two values */
        double v0 = old_data[perm[i]];
        double v1 = old_data[perm[i + 1]];

        /* Non-temporal store (bypasses cache) */
        _mm_stream_pd(&new_data[i], _mm_set_pd(v1, v0));
    }

    /* Handle odd element */
    if (i < nobs) {
        new_data[i] = old_data[perm[i]];
    }

    /* Memory fence to ensure all streaming stores complete before continuing */
    _mm_sfence();
}
#endif /* CTOOLS_TYPES_X86 */

#ifdef CTOOLS_TYPES_ARM64
/*
    ARM64 streaming gather: uses STNP (Store Pair Non-temporal) to bypass cache.
    STNP is a hint that the data won't be accessed again soon, allowing the
    memory system to optimize for write-through behavior.
*/
__attribute__((unused))
static void apply_perm_gather_double_streaming(double * CTOOLS_RESTRICT new_data,
                                                const double * CTOOLS_RESTRICT old_data,
                                                const perm_idx_t * CTOOLS_RESTRICT perm,
                                                size_t nobs)
{
    size_t i;
    size_t nobs_vec2 = nobs & ~(size_t)1;  /* Round down to multiple of 2 */

    /* Main loop: gather and stream-store 2 doubles at a time using STNP */
    for (i = 0; i < nobs_vec2; i += 2) {
        /* Prefetch source data */
        if (i + PERM_PREFETCH_FAR < nobs) {
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_FAR]]);
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_FAR + 1]]);
        }

        /* Gather two values */
        double v0 = old_data[perm[i]];
        double v1 = old_data[perm[i + 1]];

        /* Non-temporal store pair using STNP instruction */
        __asm__ __volatile__(
            "stnp %d[val0], %d[val1], [%[dst]]"
            :
            : [val0] "w" (v0), [val1] "w" (v1), [dst] "r" (&new_data[i])
            : "memory"
        );
    }

    /* Handle odd element */
    if (i < nobs) {
        new_data[i] = old_data[perm[i]];
    }

    /* Data memory barrier to ensure stores are visible */
    __asm__ __volatile__("dmb ishst" ::: "memory");
}
#endif /* CTOOLS_TYPES_ARM64 */

/*
    Main dispatch function: selects SIMD or scalar path based on platform.
*/
static inline void apply_perm_gather_double(double * CTOOLS_RESTRICT new_data,
                                             const double * CTOOLS_RESTRICT old_data,
                                             const perm_idx_t * CTOOLS_RESTRICT perm,
                                             size_t nobs)
{
#ifdef CTOOLS_TYPES_AVX2
    apply_perm_gather_double_simd(new_data, old_data, perm, nobs);
#else
    apply_perm_gather_double_scalar(new_data, old_data, perm, nobs);
#endif
}

static inline void apply_perm_gather_ptr(char ** CTOOLS_RESTRICT new_data,
                                          char ** CTOOLS_RESTRICT old_data,
                                          const perm_idx_t * CTOOLS_RESTRICT perm,
                                          size_t nobs)
{
    size_t i;
    for (i = 0; i < nobs; i++) {
        if (i + PERM_PREFETCH_FAR < nobs) {
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_FAR]]);
        }
        if (i + PERM_PREFETCH_NEAR < nobs) {
            CTOOLS_PREFETCH(&old_data[perm[i + PERM_PREFETCH_NEAR]]);
        }
        new_data[i] = old_data[perm[i]];
    }
}

/*
    Apply permutation (sort_order) to all variables in the dataset.
    After calling, data is physically reordered and sort_order is reset to identity.

    Uses GATHER pattern (random reads, sequential writes) which is optimal because:
    - Sequential writes are cache-friendly and can be optimized by the CPU
    - Random reads benefit from prefetching
    - No memory overhead for inverse permutation

    MEMORY OPTIMIZATION: Uses a shared buffer pool instead of allocating a new
    buffer for each variable in parallel. This reduces peak memory from
    2x data size to: data size + (num_threads × variable_size).

    @param data  [in/out] stata_data with computed sort_order
    @return      STATA_OK on success, STATA_ERR_MEMORY on allocation failure
*/
stata_retcode ctools_apply_permutation(stata_data *data)
{
    size_t nvars = data->nvars;
    size_t nobs = data->nobs;
    perm_idx_t *perm = data->sort_order;
    size_t j;

    if (nvars == 0 || nobs == 0) return STATA_OK;

    /* Initialize adaptive prefetch distances (once per process) */
    init_perm_prefetch_distances();

    #ifdef _OPENMP
    /*
        BUFFER POOL APPROACH: Pre-allocate one buffer per thread.
        Each thread reuses its buffer across multiple variables.
        This limits peak memory to: original_data + (nthreads × buffer_size)
        instead of: original_data × 2
    */
    int nthreads = omp_get_max_threads();
    if (nthreads > (int)nvars) nthreads = (int)nvars;
    if (nthreads < 1) nthreads = 1;

    /* Allocate buffer pool - one buffer per thread */
    double **dbl_buffers = (double **)malloc(nthreads * sizeof(double *));
    char ***ptr_buffers = (char ***)malloc(nthreads * sizeof(char **));
    if (!dbl_buffers || !ptr_buffers) {
        free(dbl_buffers);
        free(ptr_buffers);
        return STATA_ERR_MEMORY;
    }

    int alloc_success = 1;
    for (int t = 0; t < nthreads; t++) {
        dbl_buffers[t] = (double *)ctools_safe_aligned_alloc2(64, nobs, sizeof(double));
        ptr_buffers[t] = (char **)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(char *));
        if (!dbl_buffers[t] || !ptr_buffers[t]) {
            alloc_success = 0;
            break;
        }
    }

    if (!alloc_success) {
        for (int t = 0; t < nthreads; t++) {
            ctools_aligned_free(dbl_buffers[t]);
            ctools_aligned_free(ptr_buffers[t]);
        }
        free(dbl_buffers);
        free(ptr_buffers);
        return STATA_ERR_MEMORY;
    }

    /* Apply permutation using thread-local buffers */
    int success = 1;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
    for (j = 0; j < nvars; j++) {
        if (!success) continue;

        int tid = omp_get_thread_num();
        stata_variable *var = &data->vars[j];

        if (var->type == STATA_TYPE_DOUBLE) {
            double *old_data = var->data.dbl;
            double *new_data = dbl_buffers[tid];

            apply_perm_gather_double(new_data, old_data, perm, nobs);

            /* Swap: old buffer goes to pool, new data becomes variable's data */
            dbl_buffers[tid] = old_data;
            var->data.dbl = new_data;
        } else {
            char **old_data = var->data.str;
            char **new_data = ptr_buffers[tid];

            apply_perm_gather_ptr(new_data, old_data, perm, nobs);

            /* Swap: old buffer goes to pool, new data becomes variable's data */
            ptr_buffers[tid] = old_data;
            var->data.str = new_data;
        }
    }

    /* Free the buffer pool (now contains the old data arrays) */
    for (int t = 0; t < nthreads; t++) {
        ctools_aligned_free(dbl_buffers[t]);
        ctools_aligned_free(ptr_buffers[t]);
    }
    free(dbl_buffers);
    free(ptr_buffers);

    if (!success) return STATA_ERR_MEMORY;

    #else
    /* Non-OpenMP: process sequentially with single reused buffer */
    double *dbl_buffer = (double *)ctools_safe_aligned_alloc2(64, nobs, sizeof(double));
    char **ptr_buffer = (char **)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(char *));

    if (!dbl_buffer || !ptr_buffer) {
        ctools_aligned_free(dbl_buffer);
        ctools_aligned_free(ptr_buffer);
        return STATA_ERR_MEMORY;
    }

    for (j = 0; j < nvars; j++) {
        stata_variable *var = &data->vars[j];

        if (var->type == STATA_TYPE_DOUBLE) {
            double *old_data = var->data.dbl;
            apply_perm_gather_double(dbl_buffer, old_data, perm, nobs);
            /* Swap buffers */
            var->data.dbl = dbl_buffer;
            dbl_buffer = old_data;
        } else {
            char **old_data = var->data.str;
            apply_perm_gather_ptr(ptr_buffer, old_data, perm, nobs);
            /* Swap buffers */
            var->data.str = ptr_buffer;
            ptr_buffer = old_data;
        }
    }

    ctools_aligned_free(dbl_buffer);
    ctools_aligned_free(ptr_buffer);
    #endif

    /* Reset sort_order to identity */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (j = 0; j < nobs; j++) {
        perm[j] = (perm_idx_t)j;
    }

    return STATA_OK;
}

/* ===========================================================================
   Unified Sort Dispatcher
   =========================================================================== */

/*
    Auto-select the best sort algorithm based on data characteristics.

    Selection logic:
    1. If any sort key is a string variable → MSD radix (optimized for strings)
    2. For numeric-only data → Counting sort with LSD radix fallback
       - Counting sort: O(n+k), optimal for small-range integers (dates, codes)
       - LSD radix fallback: O(n*d), still faster than comparison sorts

    This provides massive speedups for common Stata workloads:
    - Panel data with date/time keys: counting sort is often 100x+ faster
    - ID variables with moderate range: counting sort or LSD radix
    - String keys (names, codes): MSD radix handles variable-length efficiently
*/
static sort_algorithm_t auto_select_algorithm(stata_data *data, int *sort_vars, size_t nsort)
{
    /* Check if any sort variable is a string */
    for (size_t i = 0; i < nsort; i++) {
        int var_idx = sort_vars[i] - 1;  /* Convert to 0-based */
        if (var_idx >= 0 && (size_t)var_idx < data->nvars) {
            if (data->vars[var_idx].type == STATA_TYPE_STRING) {
                /* String variable detected - use MSD radix */
                return SORT_ALG_MSD;
            }
        }
    }

    /* All numeric - use counting sort (with automatic LSD radix fallback) */
    return SORT_ALG_COUNTING;
}

stata_retcode ctools_sort_dispatch(stata_data *data, int *sort_vars, size_t nsort,
                                    sort_algorithm_t algorithm)
{
    stata_retcode rc;

    /* Handle auto-selection */
    if (algorithm == SORT_ALG_AUTO) {
        algorithm = auto_select_algorithm(data, sort_vars, nsort);
    }

    switch (algorithm) {
        case SORT_ALG_MSD:
            rc = ctools_sort_radix_msd_order_only(data, sort_vars, nsort);
            break;
        case SORT_ALG_TIMSORT:
            rc = ctools_sort_timsort_order_only(data, sort_vars, nsort);
            break;
        case SORT_ALG_SAMPLE:
            rc = ctools_sort_sample_order_only(data, sort_vars, nsort);
            break;
        case SORT_ALG_COUNTING:
            rc = ctools_sort_counting_order_only(data, sort_vars, nsort);
            /* If counting sort not suitable, fall back to LSD radix */
            if (rc == STATA_ERR_UNSUPPORTED_TYPE) {
                rc = ctools_sort_radix_lsd_order_only(data, sort_vars, nsort);
            }
            break;
        case SORT_ALG_MERGE:
            rc = ctools_sort_merge_order_only(data, sort_vars, nsort);
            break;
        case SORT_ALG_LSD:
            rc = ctools_sort_radix_lsd_order_only(data, sort_vars, nsort);
            break;
        case SORT_ALG_IPS4O:
        default:
            rc = ctools_sort_ips4o_order_only(data, sort_vars, nsort);
            break;
    }

    return rc;
}

/* ===========================================================================
   Type Conversion Utilities Implementation
   =========================================================================== */

bool ctools_parse_double_fast(const char *str, int len, double *result, double missval)
{
    const char *p = str;
    const char *end = str + len;

    /* Skip leading whitespace */
    while (p < end && (*p == ' ' || *p == '\t')) p++;
    /* Skip trailing whitespace */
    while (end > p && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\r' || end[-1] == '\n')) end--;

    len = (int)(end - p);

    /* Empty string -> missing */
    if (len == 0) {
        *result = missval;
        return true;
    }

    /* Single "." -> Stata missing */
    if (len == 1 && *p == '.') {
        *result = missval;
        return true;
    }

    /* "NA" -> missing */
    if (len == 2 && (p[0] == 'N' || p[0] == 'n') && (p[1] == 'A' || p[1] == 'a')) {
        *result = missval;
        return true;
    }

    /* "NaN" -> missing */
    if (len == 3 && (p[0] == 'N' || p[0] == 'n') && (p[1] == 'a' || p[1] == 'A') && (p[2] == 'N' || p[2] == 'n')) {
        *result = missval;
        return true;
    }

    /* Parse sign */
    bool negative = false;
    if (*p == '-') {
        negative = true;
        p++;
    } else if (*p == '+') {
        p++;
    }

    if (p >= end) return false;

    /* Parse mantissa digits into uint64 */
    uint64_t mantissa = 0;
    int decimal_pos = -1;
    int digit_count = 0;
    int total_digits = 0;
    bool has_digits = false;
    int exponent = 0;

    while (p < end) {
        char c = *p;
        if (c >= '0' && c <= '9') {
            has_digits = true;
            total_digits++;
            if (digit_count < 19) {
                mantissa = mantissa * 10 + (uint64_t)(c - '0');
                digit_count++;
            }
            p++;
        } else if (c == '.' && decimal_pos < 0) {
            decimal_pos = total_digits;
            p++;
        } else if (c == 'e' || c == 'E') {
            p++;
            if (p >= end) return false;
            bool exp_negative = false;
            if (*p == '-') {
                exp_negative = true;
                p++;
            } else if (*p == '+') {
                p++;
            }

            while (p < end && *p >= '0' && *p <= '9') {
                if (exponent < 10000) /* prevent overflow in exponent accumulator */
                    exponent = exponent * 10 + (*p - '0');
                p++;
            }

            if (exp_negative) exponent = -exponent;
            break;
        } else {
            break;
        }
    }

    if (!has_digits) return false;
    if (p != end) return false;

    /* Handle zero mantissa */
    if (mantissa == 0) {
        *result = negative ? -0.0 : 0.0;
        return true;
    }

    /* Compute the effective decimal exponent.
     * The mantissa holds the first digit_count digits as an integer.
     * With a decimal point at position decimal_pos (count of integer digits),
     * value = mantissa * 10^(exponent + decimal_pos - digit_count).
     * Without a decimal point, unconsumed trailing digits need positive shift. */
    int final_exp;
    if (decimal_pos >= 0) {
        final_exp = exponent + decimal_pos - digit_count;
    } else {
        final_exp = exponent + (total_digits - digit_count);
    }

    /* Eisel-Lemire fast path: handles ~99% of inputs exactly */
    if (digit_count <= 19 && eisel_lemire_compute(mantissa, final_exp, negative, result)) {
        return true;
    }

    /* Fallback: use strtod for ambiguous / very long inputs */
    {
        /* p was advanced past the parsed portion; use the trimmed range */
        const char *start = str;
        while (start < end && (*start == ' ' || *start == '\t')) start++;

        /* strtod needs a null-terminated string. If the input isn't
         * null-terminated at 'end', copy to a small buffer. */
        int trimmed_len = (int)(end - start);
        char stack_buf[128];
        char *buf;
        bool allocated = false;

        if (trimmed_len < (int)sizeof(stack_buf)) {
            buf = stack_buf;
        } else {
            buf = (char *)malloc(trimmed_len + 1);
            if (!buf) return false;
            allocated = true;
        }
        memcpy(buf, start, trimmed_len);
        buf[trimmed_len] = '\0';

        char *endptr;
        *result = strtod(buf, &endptr);
        bool ok = (endptr == buf + trimmed_len);

        if (allocated) free(buf);
        return ok;
    }
}

/*
    Parse double with custom decimal and group separators.
    Handles European formats like "1.234,56" (decimal=',', group='.').

    @param str       Input string
    @param len       Length of input
    @param result    Output value
    @param missval   Value to use for missing
    @param dec_sep   Decimal separator character (default '.')
    @param grp_sep   Group separator character (default '\0' = none)
    @return true on success, false on parse error
*/
bool ctools_parse_double_with_separators(const char *str, int len, double *result, double missval,
                                          char dec_sep, char grp_sep)
{
    /* If using default separators, use the fast path */
    if (dec_sep == '.' && grp_sep == '\0') {
        return ctools_parse_double_fast(str, len, result, missval);
    }

    const char *p = str;
    const char *end = str + len;

    /* Skip leading whitespace */
    while (p < end && (*p == ' ' || *p == '\t')) p++;
    /* Skip trailing whitespace */
    while (end > p && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\r' || end[-1] == '\n')) end--;

    len = (int)(end - p);

    /* Empty string -> missing */
    if (len == 0) {
        *result = missval;
        return true;
    }

    /* Single decimal separator -> Stata missing (like single ".") */
    if (len == 1 && *p == dec_sep) {
        *result = missval;
        return true;
    }

    /* "NA" -> missing */
    if (len == 2 && (p[0] == 'N' || p[0] == 'n') && (p[1] == 'A' || p[1] == 'a')) {
        *result = missval;
        return true;
    }

    /* "NaN" -> missing */
    if (len == 3 && (p[0] == 'N' || p[0] == 'n') && (p[1] == 'a' || p[1] == 'A') && (p[2] == 'N' || p[2] == 'n')) {
        *result = missval;
        return true;
    }

    /* Parse sign */
    bool negative = false;
    if (*p == '-') {
        negative = true;
        p++;
    } else if (*p == '+') {
        p++;
    }

    if (p >= end) return false;

    /* Parse mantissa, handling group separators */
    uint64_t mantissa = 0;
    int decimal_pos = -1;
    int digit_count = 0;
    int total_digits = 0;
    bool has_digits = false;

    while (p < end) {
        char c = *p;
        if (c >= '0' && c <= '9') {
            has_digits = true;
            total_digits++;
            if (digit_count < 18) {
                mantissa = mantissa * 10 + (c - '0');
                digit_count++;
            }
            p++;
        } else if (c == dec_sep && decimal_pos < 0) {
            /* Decimal separator */
            decimal_pos = total_digits;
            p++;
        } else if (grp_sep != '\0' && c == grp_sep) {
            /* Group separator - just skip it */
            p++;
        } else if (c == 'e' || c == 'E') {
            /* Scientific notation */
            p++;
            if (p >= end) return false;

            bool exp_negative = false;
            if (*p == '-') {
                exp_negative = true;
                p++;
            } else if (*p == '+') {
                p++;
            }

            int exponent = 0;
            while (p < end && *p >= '0' && *p <= '9') {
                exponent = exponent * 10 + (*p - '0');
                p++;
            }

            if (exp_negative) exponent = -exponent;

            double value = (double)mantissa;
            int decimal_shift = (decimal_pos >= 0) ? (total_digits - decimal_pos) : 0;
            int extra_digits = (total_digits > 18) ? (total_digits - 18) : 0;
            if (decimal_pos >= 0 && decimal_pos < total_digits) {
                extra_digits -= (total_digits - decimal_pos);
                if (extra_digits < 0) extra_digits = 0;
            }

            int final_exp = exponent - decimal_shift + extra_digits;

            if (final_exp > 0 && final_exp <= 22) {
                value *= ctools_pow10_table[final_exp];
            } else if (final_exp < 0 && final_exp >= -22) {
                value /= ctools_pow10_table[-final_exp];
            } else if (final_exp != 0) {
                value *= pow(10.0, final_exp);
            }

            *result = negative ? -value : value;
            return (p == end);
        } else {
            break;
        }
    }

    if (!has_digits) return false;
    if (p != end) return false;

    double value = (double)mantissa;

    /* Apply decimal shift */
    if (decimal_pos >= 0) {
        int decimal_digits = total_digits - decimal_pos;
        if (total_digits > 18 && decimal_pos < 18) {
            decimal_digits -= (total_digits - 18);
            if (decimal_digits < 0) decimal_digits = 0;
        }
        if (decimal_digits > 0 && decimal_digits <= 22) {
            value /= ctools_pow10_table[decimal_digits];
        } else if (decimal_digits > 22) {
            value /= pow(10.0, decimal_digits);
        }
    } else if (total_digits > 18) {
        int extra = total_digits - 18;
        if (extra <= 22) {
            value *= ctools_pow10_table[extra];
        } else {
            value *= pow(10.0, extra);
        }
    }

    *result = negative ? -value : value;
    return true;
}
