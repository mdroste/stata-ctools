/*
    csort_stream.c
    Memory-efficient streaming permutation for large dataset sorting

    This module provides an alternative sorting approach that minimizes memory
    usage by only loading key variables into C, then applying the permutation
    to non-key variables one at a time.

    Optimizations:
    - Inverse permutation for sequential Stata reads (much faster than random)
    - String arena allocation to avoid per-row malloc/free
    - Batched SPI calls for sequential access patterns
    - Identity permutation detection for early exit

    Memory usage:
    - Standard mode: O(nobs * nvars * 8) bytes for all data
    - Memeff mode: O(nobs * nkeys * 8) + O(nobs * 8) bytes
      (keys + one variable buffer at a time)
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* SIMD support for scatter operations and non-temporal stores */
#if defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h>
    /* AVX-512 scatter: compile-time or runtime detection (see below) */
    #ifdef __AVX512F__
        #define HAVE_AVX512_SCATTER 1
    #endif
    /* Runtime AVX-512 detection on GCC/Clang without compile-time -mavx512f */
    #if !defined(HAVE_AVX512_SCATTER) && (defined(__GNUC__) || defined(__clang__))
        #define HAVE_AVX512_RUNTIME 1
    #endif
    /* Non-temporal stores via SSE2 (always available on x86-64) */
    #define HAVE_NT_STORE 1
#elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
    /* ARM64 has STNP for non-temporal pair stores */
    #define HAVE_NT_STORE 1
#endif

/* Non-temporal prefetch hint macro for reads that won't be reused */
#if defined(__GNUC__) || defined(__clang__)
    #define CTOOLS_PREFETCH_NTA(addr) __builtin_prefetch((addr), 0, 0)  /* Non-temporal access */
#else
    #define CTOOLS_PREFETCH_NTA(addr) ((void)0)
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_spi.h"
#include "ctools_arena.h"
#include "csort_stream.h"

/* ============================================================================
   Configuration
   ============================================================================ */

/* Maximum memory to use for block buffers (default 1GB) */
#define STREAM_MAX_BUFFER_MEMORY (1024ULL * 1024 * 1024)

/* String arena size per variable (1MB should handle most string columns) */
#define STRING_ARENA_SIZE (1024 * 1024)

/* ============================================================================
   SIMD Scatter for Phase 1 (AVX-512 only)

   Scatters 8 doubles to random locations via permutation indices.
   Falls back to scalar for non-AVX512 systems.
   ============================================================================ */

#ifdef HAVE_AVX512_SCATTER
/*
    AVX-512 scatter: write 8 doubles to scattered locations in one instruction.
    indices: 8 destination offsets (32-bit, will be scaled by 8 for double*)
    values:  8 double values to scatter
    base:    base pointer of destination array
*/
static inline void scatter_8_doubles_avx512(
    double *base,
    const perm_idx_t *indices,
    const double *values)
{
    /* Load 8 indices into 256-bit register (perm_idx_t is uint32_t) */
    __m256i idx = _mm256_loadu_si256((const __m256i *)indices);

    /* Load 8 values into 512-bit register */
    __m512d vals = _mm512_loadu_pd(values);

    /* Scatter: base[indices[i]] = values[i] for i=0..7 */
    /* Scale of 8 because each double is 8 bytes */
    _mm512_i32scatter_pd(base, idx, vals, 8);
}
#endif

#ifdef HAVE_AVX512_RUNTIME
/*
    Runtime-dispatched AVX-512 scatter: compiled with target attribute so it
    can coexist in a binary built with -mavx2. Dispatched via function pointer
    after checking __builtin_cpu_supports("avx512f") at init time.
*/
__attribute__((target("avx512f")))
static void scatter_8_avx512_rt(
    double *base,
    const perm_idx_t *indices,
    const double *values)
{
    __m256i idx = _mm256_loadu_si256((const __m256i *)indices);
    __m512d vals = _mm512_loadu_pd(values);
    _mm512_i32scatter_pd(base, idx, vals, 8);
}

static void scatter_8_scalar(
    double *base,
    const perm_idx_t *indices,
    const double *values)
{
    for (int i = 0; i < 8; i++) base[indices[i]] = values[i];
}

typedef void (*scatter_fn_t)(double *, const perm_idx_t *, const double *);
static scatter_fn_t scatter_fn = NULL;

static void init_scatter_fn(void)
{
    if (scatter_fn) return;
    if (__builtin_cpu_supports("avx512f"))
        scatter_fn = scatter_8_avx512_rt;
    else
        scatter_fn = scatter_8_scalar;
}
#endif

/* ============================================================================
   Non-Temporal Stores for Phase 2

   Stream 8 doubles sequentially using non-temporal stores.
   Bypasses cache since data won't be reused.
   ============================================================================ */

#ifdef HAVE_NT_STORE
#if defined(__x86_64__) || defined(_M_X64)
/*
    x86-64: Use _mm_stream_pd to write 2 doubles at a time (SSE2).
    Writes 8 doubles total using 4 streaming store instructions.
    Destination must be 16-byte aligned for best performance.
*/
static inline void stream_8_doubles_x86(double *dst, const double *src)
{
    _mm_stream_pd(dst + 0, _mm_loadu_pd(src + 0));
    _mm_stream_pd(dst + 2, _mm_loadu_pd(src + 2));
    _mm_stream_pd(dst + 4, _mm_loadu_pd(src + 4));
    _mm_stream_pd(dst + 6, _mm_loadu_pd(src + 6));
}

#elif defined(__aarch64__) || defined(_M_ARM64)
/*
    ARM64: Use STNP (Store Pair Non-temporal) instruction.
    Writes 8 doubles using 4 STNP instructions.
*/
static inline void stream_8_doubles_arm(double *dst, const double *src)
{
    /* STNP stores a pair of values with non-temporal hint */
    __asm__ __volatile__(
        "ldp q0, q1, [%[src]]\n\t"
        "ldp q2, q3, [%[src], #32]\n\t"
        "stnp q0, q1, [%[dst]]\n\t"
        "stnp q2, q3, [%[dst], #32]\n\t"
        :
        : [dst] "r" (dst), [src] "r" (src)
        : "v0", "v1", "v2", "v3", "memory"
    );
}
#endif

/* Unified interface for non-temporal 8-double store */
__attribute__((unused))
static inline void stream_store_8_doubles(double *dst, const double *src)
{
#if defined(__x86_64__) || defined(_M_X64)
    stream_8_doubles_x86(dst, src);
#elif defined(__aarch64__) || defined(_M_ARM64)
    stream_8_doubles_arm(dst, src);
#else
    /* Fallback: regular copy */
    for (int i = 0; i < 8; i++) dst[i] = src[i];
#endif
}
#endif /* HAVE_NT_STORE */

/* ============================================================================
   Inverse Permutation Builder + Identity Detection

   Given perm[] where sorted[i] = original[perm[i]], compute inv[] where
   inv[perm[i]] = i, so we can read Stata sequentially: read original[j],
   store at buffer[inv[j]].

   Returns 1 if permutation is identity (already sorted), 0 otherwise.

   Optimizations:
   - Parallel construction with OpenMP
   - Vectorized identity check
   ============================================================================ */

static int build_inverse_permutation(
    const perm_idx_t *perm,
    perm_idx_t *inv,
    size_t nobs)
{
    int is_identity = 1;

    #ifdef _OPENMP
    /* Parallel inverse permutation build - each thread handles a chunk */
    #pragma omp parallel
    {
        int local_is_identity = 1;

        #pragma omp for schedule(static) nowait
        for (size_t i = 0; i < nobs; i++) {
            inv[perm[i]] = (perm_idx_t)i;
            if (perm[i] != i) {
                local_is_identity = 0;
            }
        }

        /* Reduce identity flag */
        if (!local_is_identity) {
            #pragma omp atomic write
            is_identity = 0;
        }
    }
    #else
    /* Sequential fallback */
    for (size_t i = 0; i < nobs; i++) {
        inv[perm[i]] = (perm_idx_t)i;
        if (perm[i] != i) {
            is_identity = 0;
        }
    }
    #endif

    return is_identity;
}

/* ============================================================================
   Optimized Streaming Permutation

   Key optimizations:
   1. Parallel variable processing - each thread handles different variables
   2. Multi-variable batching - process multiple variables per block to amortize
      permutation lookups
   3. Prefetching for permutation array and destination buffers
   4. Cache-blocked processing for better memory locality

   Memory access pattern:
   - Sequential reads from Stata (cache-friendly on Stata side)
   - Scattered writes to local buffer via inverse permutation
   - Sequential writes back to Stata
   ============================================================================ */

/* Number of variables to process together in multi-var mode */
#define STREAM_MULTI_VAR_BATCH 4

/*
    Adaptive prefetch distance for streaming operations.
    Uses cache-aware distances based on detected hardware.
*/
static int _stream_prefetch_dist = 0;

__attribute__((unused))
static inline int get_stream_prefetch_dist(void)
{
    if (_stream_prefetch_dist == 0) {
        ctools_prefetch_distances dist = ctools_get_prefetch_distances();
        _stream_prefetch_dist = dist.stream;
    }
    return _stream_prefetch_dist;
}

#define STREAM_PREFETCH_DIST (get_stream_prefetch_dist())

/*
    Process a string variable with arena allocation.
    Uses ctools_string_arena with NO_FALLBACK mode - returns error if arena exhausted.
*/
static int stream_permute_string_var(
    ST_int stata_var,
    const perm_idx_t *inv_perm,
    size_t nobs,
    const perm_idx_t *obs_map)
{
    size_t i;

    /* String variable with arena allocation */
    size_t arena_size = nobs * 256;  /* Estimate avg 256 bytes per string */
    if (arena_size < STRING_ARENA_SIZE) arena_size = STRING_ARENA_SIZE;

    /* Use NO_FALLBACK mode - if arena is exhausted, we fall back to strdup manually */
    ctools_string_arena *arena = ctools_string_arena_create(arena_size,
                                    CTOOLS_STRING_ARENA_STRDUP_FALLBACK);
    char **str_ptrs = (char **)calloc(nobs, sizeof(char *));

    if (!arena || !str_ptrs) {
        ctools_string_arena_free(arena);
        free(str_ptrs);
        return -1;
    }

    char strbuf[2048];

    /* Sequential read from Stata using obs_map, scatter to buffer via inverse perm */
    for (i = 0; i < nobs; i++) {
        SF_sdata(stata_var, (ST_int)obs_map[i], strbuf);
        str_ptrs[inv_perm[i]] = ctools_string_arena_strdup(arena, strbuf);

        if (!str_ptrs[inv_perm[i]]) {
            /* Cleanup on allocation failure */
            for (size_t k = 0; k < nobs; k++) {
                if (str_ptrs[k] && !ctools_string_arena_owns(arena, str_ptrs[k])) {
                    free(str_ptrs[k]);
                }
            }
            ctools_string_arena_free(arena);
            free(str_ptrs);
            return -1;
        }
    }

    /* Sequential write back to Stata using obs_map */
    for (i = 0; i < nobs; i++) {
        SF_sstore(stata_var, (ST_int)obs_map[i], str_ptrs[i] ? str_ptrs[i] : "");
    }

    /* Free fallback allocations (those not owned by arena) */
    if (arena->has_fallback) {
        for (i = 0; i < nobs; i++) {
            if (str_ptrs[i] && !ctools_string_arena_owns(arena, str_ptrs[i])) {
                free(str_ptrs[i]);
            }
        }
    }
    ctools_string_arena_free(arena);
    free(str_ptrs);

    return 0;
}

stata_retcode csort_stream_apply_permutation(
    const perm_idx_t *perm,
    size_t nobs,
    int *nonkey_var_indices,
    size_t nvars_nonkey,
    const perm_idx_t *obs_map,
    size_t block_size,
    int vars_per_batch,
    csort_stream_timings *timings)
{
    size_t v;
    perm_idx_t *inv_perm = NULL;
    int is_identity;
    double t_phase;
    (void)block_size;

    /* Validate and constrain vars_per_batch */
    if (vars_per_batch < 1) vars_per_batch = 1;
    if (vars_per_batch > 16) vars_per_batch = 16;

    /* Initialize timing accumulators */
    double t_build_inv = 0.0;
    double t_scatter = 0.0;
    double t_writeback = 0.0;
    double t_string = 0.0;

    if (nvars_nonkey == 0 || nobs == 0) {
        return STATA_OK;
    }

#ifdef HAVE_AVX512_RUNTIME
    init_scatter_fn();
#endif

    /* Build inverse permutation once (shared across all variables) */
    t_phase = ctools_timer_seconds();

    inv_perm = (perm_idx_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE,
                                                         nobs, sizeof(perm_idx_t));
    if (!inv_perm) {
        return STATA_ERR_MEMORY;
    }

    is_identity = build_inverse_permutation(perm, inv_perm, nobs);

    t_build_inv = ctools_timer_seconds() - t_phase;

    /* Early exit if data is already sorted */
    if (is_identity) {
        ctools_aligned_free(inv_perm);
        if (timings) {
            timings->stream_build_inv_time = t_build_inv;
            timings->stream_scatter_time = 0.0;
            timings->stream_writeback_time = 0.0;
            timings->stream_string_time = 0.0;
        }
        return STATA_OK;
    }

    /* Separate string and numeric variables */
    int *numeric_vars = (int *)malloc(nvars_nonkey * sizeof(int));
    int *string_vars = (int *)malloc(nvars_nonkey * sizeof(int));
    size_t n_numeric = 0, n_string = 0;

    if (!numeric_vars || !string_vars) {
        free(numeric_vars);
        free(string_vars);
        ctools_aligned_free(inv_perm);
        return STATA_ERR_MEMORY;
    }

    for (v = 0; v < nvars_nonkey; v++) {
        if (SF_var_is_string((ST_int)nonkey_var_indices[v])) {
            string_vars[n_string++] = nonkey_var_indices[v];
        } else {
            numeric_vars[n_numeric++] = nonkey_var_indices[v];
        }
    }

    /* Store variable counts in timings */
    if (timings) {
        timings->n_numeric_vars = n_numeric;
        timings->n_string_vars = n_string;
    }

    /* Process numeric variables in batches with shared buffer pool
       vars_per_batch controls how many variables are loaded at a time */
    #ifdef _OPENMP
    int success = 1;
    double **batch_bufs = NULL;
    int nthreads = 0;
    int actual_batch_size = vars_per_batch;

    if (n_numeric > 0) {
        /* Limit batch size to number of numeric variables */
        if (actual_batch_size > (int)n_numeric) {
            actual_batch_size = (int)n_numeric;
        }

        nthreads = omp_get_max_threads();
        if (nthreads > actual_batch_size) nthreads = actual_batch_size;
        if (nthreads < 1) nthreads = 1;

        /* Pre-allocate buffers for the batch (one per variable in batch) */
        batch_bufs = (double **)malloc(actual_batch_size * sizeof(double *));
        if (!batch_bufs) {
            free(numeric_vars);
            free(string_vars);
            ctools_aligned_free(inv_perm);
            return STATA_ERR_MEMORY;
        }
        for (int b = 0; b < actual_batch_size; b++) {
            batch_bufs[b] = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(double));
            if (!batch_bufs[b]) {
                for (int k = 0; k < b; k++) ctools_aligned_free(batch_bufs[k]);
                free(batch_bufs);
                free(numeric_vars);
                free(string_vars);
                ctools_aligned_free(inv_perm);
                return STATA_ERR_MEMORY;
            }
        }
    }

    /* Track scatter and writeback times */
    if (n_numeric > 0) {
        /* Process variables in batches of actual_batch_size */
        size_t batch_start = 0;
        while (batch_start < n_numeric) {
            /* Determine how many variables in this batch */
            int this_batch = actual_batch_size;
            if (batch_start + this_batch > n_numeric) {
                this_batch = (int)(n_numeric - batch_start);
            }

            /* Determine threads for this batch */
            int batch_threads = nthreads;
            if (batch_threads > this_batch) batch_threads = this_batch;

            double batch_scatter_time = 0.0;
            double batch_writeback_time = 0.0;

            /* Phase 1: Read and scatter for all variables in batch (parallel) */
            double t_batch_scatter = ctools_timer_seconds();

            #pragma omp parallel for schedule(static) num_threads(batch_threads)
            for (int b = 0; b < this_batch; b++) {
                size_t var_idx = batch_start + b;
                double *num_buf = batch_bufs[b];
                ST_int stata_var = (ST_int)numeric_vars[var_idx];
                size_t i;

                size_t nobs_batch16 = (nobs / 16) * 16;
                double batch_vals[16] __attribute__((aligned(64)));

                for (i = 0; i < nobs_batch16; i += 16) {
                    /* Prefetch future permutation indices for next iteration */
                    if (i + STREAM_PREFETCH_DIST < nobs) {
                        CTOOLS_PREFETCH(&inv_perm[i + STREAM_PREFETCH_DIST]);
                    }

                    /* Read data from Stata */
                    SF_VDATA_BATCH16(stata_var, (ST_int)obs_map[i], batch_vals);

#ifdef HAVE_AVX512_SCATTER
                    /* AVX-512: Scatter 8 doubles at a time using SIMD scatter instruction.
                       Two scatter calls for 16 values. */
                    scatter_8_doubles_avx512(num_buf, &inv_perm[i + 0], &batch_vals[0]);
                    scatter_8_doubles_avx512(num_buf, &inv_perm[i + 8], &batch_vals[8]);
#elif defined(HAVE_AVX512_RUNTIME)
                    /* Runtime-dispatched AVX-512 or scalar fallback */
                    scatter_fn(num_buf, &inv_perm[i + 0], &batch_vals[0]);
                    scatter_fn(num_buf, &inv_perm[i + 8], &batch_vals[8]);
#else
                    /* Scalar fallback: Prefetch destinations then scatter */
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 0]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 1]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 2]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 3]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 4]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 5]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 6]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 7]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 8]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 9]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 10]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 11]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 12]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 13]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 14]]);
                    CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 15]]);

                    num_buf[inv_perm[i + 0]]  = batch_vals[0];
                    num_buf[inv_perm[i + 1]]  = batch_vals[1];
                    num_buf[inv_perm[i + 2]]  = batch_vals[2];
                    num_buf[inv_perm[i + 3]]  = batch_vals[3];
                    num_buf[inv_perm[i + 4]]  = batch_vals[4];
                    num_buf[inv_perm[i + 5]]  = batch_vals[5];
                    num_buf[inv_perm[i + 6]]  = batch_vals[6];
                    num_buf[inv_perm[i + 7]]  = batch_vals[7];
                    num_buf[inv_perm[i + 8]]  = batch_vals[8];
                    num_buf[inv_perm[i + 9]]  = batch_vals[9];
                    num_buf[inv_perm[i + 10]] = batch_vals[10];
                    num_buf[inv_perm[i + 11]] = batch_vals[11];
                    num_buf[inv_perm[i + 12]] = batch_vals[12];
                    num_buf[inv_perm[i + 13]] = batch_vals[13];
                    num_buf[inv_perm[i + 14]] = batch_vals[14];
                    num_buf[inv_perm[i + 15]] = batch_vals[15];
#endif
                }
                for (i = nobs_batch16; i < nobs; i++) {
                    double val;
                    SF_vdata(stata_var, (ST_int)obs_map[i], &val);
                    num_buf[inv_perm[i]] = val;
                }
            }

            batch_scatter_time = ctools_timer_seconds() - t_batch_scatter;
            t_scatter += batch_scatter_time;

            /* Phase 2: Write back all variables in batch (parallel) */
            double t_batch_writeback = ctools_timer_seconds();

            #pragma omp parallel for schedule(static) num_threads(batch_threads)
            for (int b = 0; b < this_batch; b++) {
                size_t var_idx = batch_start + b;
                double *num_buf = batch_bufs[b];
                ST_int stata_var = (ST_int)numeric_vars[var_idx];
                size_t i;

                size_t nobs_batch8 = (nobs / 8) * 8;
                for (i = 0; i < nobs_batch8; i += 8) {
                    /* Prefetch next chunk with non-temporal hint (won't be reused) */
                    CTOOLS_PREFETCH_NTA(&num_buf[i + 64]);

                    /* Store 8 values to Stata */
                    SF_vstore(stata_var, (ST_int)obs_map[i + 0], num_buf[i + 0]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 1], num_buf[i + 1]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 2], num_buf[i + 2]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 3], num_buf[i + 3]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 4], num_buf[i + 4]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 5], num_buf[i + 5]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 6], num_buf[i + 6]);
                    SF_vstore(stata_var, (ST_int)obs_map[i + 7], num_buf[i + 7]);
                }
                for (i = nobs_batch8; i < nobs; i++) {
                    SF_vstore(stata_var, (ST_int)obs_map[i], num_buf[i]);
                }
            }

            batch_writeback_time = ctools_timer_seconds() - t_batch_writeback;
            t_writeback += batch_writeback_time;

            batch_start += this_batch;
        }

        /* Free batch buffers */
        for (int b = 0; b < actual_batch_size; b++) {
            ctools_aligned_free(batch_bufs[b]);
        }
        free(batch_bufs);
    }

    #else
    /* Non-OpenMP fallback: process sequentially with single reused buffer */
    (void)vars_per_batch;  /* Unused in non-OpenMP build */
    double *num_buf = NULL;
    if (n_numeric > 0) {
        num_buf = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(double));
        if (!num_buf) {
            free(numeric_vars);
            free(string_vars);
            ctools_aligned_free(inv_perm);
            return STATA_ERR_MEMORY;
        }
    }
    if (num_buf) {
        for (v = 0; v < n_numeric; v++) {
            ST_int stata_var = (ST_int)numeric_vars[v];
            size_t i;
            size_t nobs_batch16 = (nobs / 16) * 16;
            double batch_vals[16] __attribute__((aligned(64)));
            double t_var_scatter, t_var_writeback;

            /* Phase 1: Sequential read, scatter via inverse perm */
            t_var_scatter = ctools_timer_seconds();

            for (i = 0; i < nobs_batch16; i += 16) {
                /* Read data from Stata */
                SF_VDATA_BATCH16(stata_var, (ST_int)obs_map[i], batch_vals);

#ifdef HAVE_AVX512_SCATTER
                /* AVX-512: Scatter 8 doubles at a time */
                scatter_8_doubles_avx512(num_buf, &inv_perm[i + 0], &batch_vals[0]);
                scatter_8_doubles_avx512(num_buf, &inv_perm[i + 8], &batch_vals[8]);
#elif defined(HAVE_AVX512_RUNTIME)
                /* Runtime-dispatched AVX-512 or scalar fallback */
                scatter_fn(num_buf, &inv_perm[i + 0], &batch_vals[0]);
                scatter_fn(num_buf, &inv_perm[i + 8], &batch_vals[8]);
#else
                /* Scalar fallback: Prefetch destinations then scatter */
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 0]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 1]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 2]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 3]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 4]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 5]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 6]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 7]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 8]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 9]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 10]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 11]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 12]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 13]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 14]]);
                CTOOLS_PREFETCH_W(&num_buf[inv_perm[i + 15]]);

                num_buf[inv_perm[i + 0]]  = batch_vals[0];
                num_buf[inv_perm[i + 1]]  = batch_vals[1];
                num_buf[inv_perm[i + 2]]  = batch_vals[2];
                num_buf[inv_perm[i + 3]]  = batch_vals[3];
                num_buf[inv_perm[i + 4]]  = batch_vals[4];
                num_buf[inv_perm[i + 5]]  = batch_vals[5];
                num_buf[inv_perm[i + 6]]  = batch_vals[6];
                num_buf[inv_perm[i + 7]]  = batch_vals[7];
                num_buf[inv_perm[i + 8]]  = batch_vals[8];
                num_buf[inv_perm[i + 9]]  = batch_vals[9];
                num_buf[inv_perm[i + 10]] = batch_vals[10];
                num_buf[inv_perm[i + 11]] = batch_vals[11];
                num_buf[inv_perm[i + 12]] = batch_vals[12];
                num_buf[inv_perm[i + 13]] = batch_vals[13];
                num_buf[inv_perm[i + 14]] = batch_vals[14];
                num_buf[inv_perm[i + 15]] = batch_vals[15];
#endif
            }
            for (i = nobs_batch16; i < nobs; i++) {
                double val;
                SF_vdata(stata_var, (ST_int)obs_map[i], &val);
                num_buf[inv_perm[i]] = val;
            }

            t_var_scatter = ctools_timer_seconds() - t_var_scatter;
            t_scatter += t_var_scatter;

            /* Phase 2: Sequential write back to Stata
               Use non-temporal prefetch hints since buffer won't be reused. */
            t_var_writeback = ctools_timer_seconds();

            size_t nobs_batch8 = (nobs / 8) * 8;
            for (i = 0; i < nobs_batch8; i += 8) {
                /* Prefetch next chunk with non-temporal hint (won't be reused) */
                CTOOLS_PREFETCH_NTA(&num_buf[i + 64]);

                /* Store 8 values to Stata */
                SF_vstore(stata_var, (ST_int)obs_map[i + 0], num_buf[i + 0]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 1], num_buf[i + 1]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 2], num_buf[i + 2]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 3], num_buf[i + 3]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 4], num_buf[i + 4]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 5], num_buf[i + 5]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 6], num_buf[i + 6]);
                SF_vstore(stata_var, (ST_int)obs_map[i + 7], num_buf[i + 7]);
            }
            for (i = nobs_batch8; i < nobs; i++) {
                SF_vstore(stata_var, (ST_int)obs_map[i], num_buf[i]);
            }

            t_var_writeback = ctools_timer_seconds() - t_var_writeback;
            t_writeback += t_var_writeback;
        }
        ctools_aligned_free(num_buf);
    }
    #endif

    /* Process string variables (parallel - each var has independent arena/buffers) */
    t_phase = ctools_timer_seconds();

    #ifdef _OPENMP
    {
        int string_threads = omp_get_max_threads();
        if (string_threads > (int)n_string) string_threads = (int)n_string;
        if (string_threads < 1) string_threads = 1;

        #pragma omp parallel for schedule(static) num_threads(string_threads)
        for (int sv = 0; sv < (int)n_string; sv++) {
            if (stream_permute_string_var((ST_int)string_vars[sv], inv_perm, nobs, obs_map) != 0) {
                #pragma omp atomic write
                success = 0;
            }
        }
    }
    #else
    for (v = 0; v < n_string; v++) {
        if (stream_permute_string_var((ST_int)string_vars[v], inv_perm, nobs, obs_map) != 0) {
            /* error */
        }
    }
    #endif

    t_string = ctools_timer_seconds() - t_phase;

    /* Store detailed timings */
    if (timings) {
        timings->stream_build_inv_time = t_build_inv;
        timings->stream_scatter_time = t_scatter;
        timings->stream_writeback_time = t_writeback;
        timings->stream_string_time = t_string;
        timings->n_numeric_vars = n_numeric;
        timings->n_string_vars = n_string;
    }

    free(numeric_vars);
    free(string_vars);
    ctools_aligned_free(inv_perm);

    #ifdef _OPENMP
    if (!success) {
        return STATA_ERR_MEMORY;
    }
    #endif

    return STATA_OK;
}

/* ============================================================================
   Block Size Selection (kept for API compatibility)
   ============================================================================ */

size_t csort_stream_choose_block_size(
    size_t nobs,
    size_t nvars_nonkey,
    size_t available_memory)
{
    size_t bytes_per_row;
    size_t max_block;

    if (available_memory == 0) {
        available_memory = STREAM_MAX_BUFFER_MEMORY;
    }

    bytes_per_row = nvars_nonkey * sizeof(double) + 64;
    max_block = available_memory / bytes_per_row;

    if (max_block < CSORT_STREAM_MIN_BLOCK_SIZE) {
        max_block = CSORT_STREAM_MIN_BLOCK_SIZE;
    }
    if (max_block > CSORT_STREAM_DEFAULT_BLOCK_SIZE * 4) {
        max_block = CSORT_STREAM_DEFAULT_BLOCK_SIZE * 4;
    }
    if (max_block > nobs) {
        max_block = nobs;
    }

    return max_block;
}

int csort_stream_recommended(size_t nobs, size_t nvars, size_t nkeys)
{
    size_t total_bytes;
    size_t key_bytes;

    total_bytes = nobs * nvars * sizeof(double);
    total_bytes += nobs * sizeof(perm_idx_t);
    total_bytes += nobs * nvars * sizeof(double);

    key_bytes = nobs * nkeys * sizeof(double);
    key_bytes += nobs * sizeof(perm_idx_t);

    if (total_bytes > CSORT_STREAM_THRESHOLD_BYTES &&
        key_bytes < total_bytes / 4) {
        return 1;
    }

    return 0;
}

/* ============================================================================
   Main Entry Point
   ============================================================================ */

stata_retcode csort_stream_sort(
    int *key_var_indices,
    size_t nkeys,
    int *all_var_indices,
    size_t nvars,
    sort_algorithm_t algorithm,
    size_t block_size,
    int vars_per_batch,
    csort_stream_timings *timings)
{
    ctools_filtered_data key_filtered;
    perm_idx_t *obs_map = NULL;
    stata_retcode rc;
    double t_start, t_phase;
    size_t nobs = 0;
    size_t i, j;
    perm_idx_t *saved_perm = NULL;
    int *nonkey_var_indices = NULL;
    int *local_sort_vars = NULL;

    csort_stream_timings local_timings = {0};
    t_start = ctools_timer_seconds();

    /* ================================================================
       Phase 1: Load only key variables using filtered loading
       ================================================================ */
    t_phase = ctools_timer_seconds();

    ctools_filtered_data_init(&key_filtered);
    rc = ctools_data_load(&key_filtered, key_var_indices, nkeys, 0, 0, 0);

    local_timings.load_keys_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
    }

    /* Get filtered observation count and obs_map */
    nobs = key_filtered.data.nobs;
    obs_map = key_filtered.obs_map;

    if (nobs == 0) {
        /* No observations to sort */
        ctools_filtered_data_free(&key_filtered);
        if (timings) *timings = local_timings;
        return STATA_OK;
    }

    /* ================================================================
       Phase 2: Sort key variables (compute permutation only)
       ================================================================ */
    t_phase = ctools_timer_seconds();

    local_sort_vars = (int *)malloc(nkeys * sizeof(int));
    if (!local_sort_vars) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }
    for (i = 0; i < nkeys; i++) {
        local_sort_vars[i] = (int)(i + 1);
    }

    /* Call unified sort dispatcher (computes sort_order only) */
    rc = ctools_sort_dispatch(&key_filtered.data, local_sort_vars, nkeys, algorithm);

    local_timings.sort_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
    }

    /* ================================================================
       Phase 2.5: Save permutation before applying (it gets reset)
       ================================================================ */
    saved_perm = (perm_idx_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE,
                                                           nobs, sizeof(perm_idx_t));
    if (!saved_perm) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }
    memcpy(saved_perm, key_filtered.data.sort_order, nobs * sizeof(perm_idx_t));

    /* ================================================================
       Phase 3: Apply permutation to key variables (in C memory)
       ================================================================ */
    t_phase = ctools_timer_seconds();

    rc = ctools_apply_permutation(&key_filtered.data);

    local_timings.permute_keys_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
    }

    /* ================================================================
       Phase 4: Write sorted key variables back to Stata using obs_map
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Store each key variable using obs_map */
    for (size_t k = 0; k < nkeys && rc == STATA_OK; k++) {
        rc = ctools_store_filtered(key_filtered.data.vars[k].data.dbl,
                                   nobs, key_var_indices[k], obs_map);
    }

    local_timings.store_keys_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
    }

    /* DEBUG: Verify first key is sorted in Stata after store */
    #ifdef CSORT_DEBUG_VERIFY
    {
        double prev_val, curr_val;
        int sorted_ok = 1;
        SF_vdata(key_var_indices[0], (ST_int)obs_map[0], &prev_val);
        for (size_t check_i = 1; check_i < nobs && sorted_ok; check_i++) {
            SF_vdata(key_var_indices[0], (ST_int)obs_map[check_i], &curr_val);
            if (curr_val < prev_val) {
                char msg[256];
                snprintf(msg, sizeof(msg),
                    "csort DEBUG: Key not sorted at row %zu: prev=%.2f, curr=%.2f\n",
                    check_i, prev_val, curr_val);
                SF_display(msg);
                sorted_ok = 0;
            }
            prev_val = curr_val;
        }
        if (sorted_ok) {
            SF_display("csort DEBUG: Keys verified sorted after store\n");
        }
    }
    #endif

    /* Take ownership of obs_map before freeing key_filtered, since
     * ctools_filtered_data_free would free it but we still need it for
     * stream-permuting non-key variables. */
    key_filtered.obs_map = NULL;  /* Prevent free */
    ctools_filtered_data_free(&key_filtered);
    stata_data_init(&key_filtered.data);

    /* ================================================================
       Phase 5: Stream-permute non-key variables
       ================================================================ */
    t_phase = ctools_timer_seconds();

    size_t nvars_nonkey = nvars - nkeys;

    if (nvars_nonkey > 0) {
        nonkey_var_indices = (int *)malloc(nvars_nonkey * sizeof(int));
        if (!nonkey_var_indices) {
            rc = STATA_ERR_MEMORY;
            goto cleanup;
        }

        size_t nk = 0;
        for (i = 0; i < nvars; i++) {
            int is_key = 0;
            for (j = 0; j < nkeys; j++) {
                if (all_var_indices[i] == key_var_indices[j]) {
                    is_key = 1;
                    break;
                }
            }
            if (!is_key) {
                nonkey_var_indices[nk++] = all_var_indices[i];
            }
        }

        /* Sanity check: nk should equal nvars_nonkey */
        if (nk != nvars_nonkey) {
            rc = STATA_ERR_INVALID_INPUT;
            goto cleanup;
        }

        rc = csort_stream_apply_permutation(saved_perm, nobs, nonkey_var_indices,
                                             nvars_nonkey, obs_map, block_size,
                                             vars_per_batch, &local_timings);
    }

    local_timings.stream_nonkeys_time = ctools_timer_seconds() - t_phase;

cleanup:
    free(local_sort_vars);
    free(nonkey_var_indices);
    ctools_aligned_free(saved_perm);
    if (obs_map) ctools_aligned_free(obs_map);  /* We took ownership from key_filtered */
    ctools_filtered_data_free(&key_filtered);

    local_timings.total_time = ctools_timer_seconds() - t_start;
    local_timings.block_size = block_size ? block_size : CSORT_STREAM_DEFAULT_BLOCK_SIZE;
    local_timings.num_blocks = (nobs + local_timings.block_size - 1) / local_timings.block_size;

    if (timings) {
        *timings = local_timings;
    }

    /* Save timing scalars for Stata */
    SF_scal_save("_csort_stream_time_load_keys", local_timings.load_keys_time);
    SF_scal_save("_csort_stream_time_sort", local_timings.sort_time);
    SF_scal_save("_csort_stream_time_permute_keys", local_timings.permute_keys_time);
    SF_scal_save("_csort_stream_time_store_keys", local_timings.store_keys_time);
    SF_scal_save("_csort_stream_time_stream_nonkeys", local_timings.stream_nonkeys_time);
    /* Detailed breakdown of stream_nonkeys phase */
    SF_scal_save("_csort_stream_time_build_inv", local_timings.stream_build_inv_time);
    SF_scal_save("_csort_stream_time_scatter", local_timings.stream_scatter_time);
    SF_scal_save("_csort_stream_time_writeback", local_timings.stream_writeback_time);
    SF_scal_save("_csort_stream_time_strings", local_timings.stream_string_time);
    SF_scal_save("_csort_stream_n_numeric", (double)local_timings.n_numeric_vars);
    SF_scal_save("_csort_stream_n_string", (double)local_timings.n_string_vars);
    SF_scal_save("_csort_stream_time_total", local_timings.total_time);

    return rc;
}
