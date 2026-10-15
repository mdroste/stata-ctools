/*
    csort_stream.c
    Memory-efficient streaming permutation for large dataset sorting

    This module provides an alternative sorting approach that minimizes memory
    usage by only loading key variables into C, then applying the permutation
    to non-key variables one at a time.

    Optimizations:
    - Inverse permutation for sequential Stata reads (much faster than random)
    - Flat fixed-width string buffers, reused across variables
    - Batched SPI calls for sequential access patterns
    - Identity permutation detection for early exit

    Memory usage:
    - Standard mode: O(nobs * nvars * 8) bytes for all data
    - Memeff mode: O(nobs * nkeys * 8) + O(nobs * 8) bytes
      (keys + one variable buffer at a time)

    Every buffer the non-key pass needs is allocated before the sorted keys
    are written, so an allocation failure leaves the dataset untouched.
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
#elif defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
#endif

/* Non-temporal prefetch hint macro for reads that won't be reused */
#if defined(__GNUC__) || defined(__clang__)
    #define CTOOLS_PREFETCH_NTA(addr) __builtin_prefetch((addr), 0, 0)  /* Non-temporal access */
#else
    #define CTOOLS_PREFETCH_NTA(addr) ((void)0)
#endif

/* Direct SPI reads and stores resolve positions in the plugin varlist, even
 * after estimation has introduced internal Stata variables. */
#define SD_SAFEMODE
#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_spi.h"
#include "csort_stream.h"

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

/* A string variable and the slot width its flat buffer needs. */
typedef struct {
    int var;      /* 1-based Stata variable index */
    int width;    /* str# width, from the ado's hint or measured */
} stream_string_var;

/*
    Everything the non-key pass needs. csort_stream_sort allocates all of it
    (stream_plan_inverse, stream_plan_buffers) before it writes anything to
    Stata, so once the first sorted key is stored only SPI calls remain.

    The buffer pool serves both phases, which run one after the other:
    numeric batch slot b uses bufs[b] as nobs doubles, and string worker t
    uses bufs[t] as nobs slots of strings[t].width + 1 bytes. Strings are
    sorted widest first and worker t handles variables t, t + team, ..., so
    its buffer fits every variable it meets.
*/
typedef struct {
    perm_idx_t *inv_perm;          /* inv_perm[j] = sorted position of row j */
    int is_identity;               /* already sorted: nothing to move */
    double build_inv_time;
    int *numeric_vars;
    size_t n_numeric;
    stream_string_var *strings;    /* widest first */
    size_t n_string;
    char **bufs;
    int nbufs;
    int batch_size;                /* numeric variables per batch */
    int string_threads;            /* string workers */
} stream_plan;

static void stream_plan_free(stream_plan *plan)
{
    for (int b = 0; plan->bufs && b < plan->nbufs; b++) {
        ctools_aligned_free(plan->bufs[b]);
    }
    free(plan->bufs);
    free(plan->numeric_vars);
    free(plan->strings);
    ctools_aligned_free(plan->inv_perm);
    memset(plan, 0, sizeof(*plan));
}

/* Build the inverse permutation once (shared across all variables). */
static stata_retcode stream_plan_inverse(stream_plan *plan, const perm_idx_t *perm, size_t nobs)
{
    double t_phase = ctools_timer_seconds();

    plan->inv_perm = (perm_idx_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE,
                                                               nobs, sizeof(perm_idx_t));
    if (!plan->inv_perm) {
        return STATA_ERR_MEMORY;
    }
    plan->is_identity = build_inverse_permutation(perm, plan->inv_perm, nobs);
    plan->build_inv_time = ctools_timer_seconds() - t_phase;
    return STATA_OK;
}

/*
    Width of the longest value of a string variable that came without a
    usable width hint (strL was rejected before any of this). Reads only.
*/
static int stream_measure_string_width(ST_int stata_var, size_t nobs, const perm_idx_t *obs_map)
{
    char strbuf[2048];
    size_t width = 1;

    for (size_t i = 0; i < nobs; i++) {
        strbuf[0] = '\0';
        SF_sdata(stata_var, (ST_int)obs_map[i], strbuf);
        size_t len = strlen(strbuf);
        if (len > width) width = len;
    }
    return (int)width;
}

/* Widest first; equal widths in variable order. */
static int stream_compare_string_vars(const void *a, const void *b)
{
    const stream_string_var *x = (const stream_string_var *)a;
    const stream_string_var *y = (const stream_string_var *)b;
    if (x->width != y->width) return x->width > y->width ? -1 : 1;
    return (x->var > y->var) - (x->var < y->var);
}

/*
    Split the non-key variables into numeric and string lists and allocate the
    buffer pool. If not every buffer fits, batches hold fewer numeric variables
    and fewer string workers run; only a pool without a single buffer fails.
*/
static stata_retcode stream_plan_buffers(
    stream_plan *plan,
    size_t nobs,
    const int *nonkey_var_indices,
    size_t nvars_nonkey,
    const perm_idx_t *obs_map,
    int vars_per_batch,
    const int *nonkey_str_widths)
{
    size_t v;
    int b;

    plan->numeric_vars = (int *)ctools_safe_malloc2(nvars_nonkey, sizeof(int));
    plan->strings = (stream_string_var *)ctools_safe_malloc2(nvars_nonkey,
                                                             sizeof(stream_string_var));
    if (!plan->numeric_vars || !plan->strings) {
        return STATA_ERR_MEMORY;
    }

    for (v = 0; v < nvars_nonkey; v++) {
        ST_int stata_var = (ST_int)nonkey_var_indices[v];
        if (SF_var_is_string(stata_var)) {
            int w = nonkey_str_widths ? nonkey_str_widths[v] : 0;
            if (w < 1 || w > CTOOLS_MAX_STRING_LEN) {
                w = stream_measure_string_width(stata_var, nobs, obs_map);
            }
            plan->strings[plan->n_string].var = (int)stata_var;
            plan->strings[plan->n_string].width = w;
            plan->n_string++;
        } else {
            plan->numeric_vars[plan->n_numeric++] = (int)stata_var;
        }
    }
    qsort(plan->strings, plan->n_string, sizeof(stream_string_var), stream_compare_string_vars);

    /* vars_per_batch controls how many numeric variables are loaded at a time */
    int batch_size = vars_per_batch < 1 ? 1 : (vars_per_batch > 16 ? 16 : vars_per_batch);
    if ((size_t)batch_size > plan->n_numeric) batch_size = (int)plan->n_numeric;
    int string_threads = 0;
    if (plan->n_string > 0) {
        string_threads = ctools_get_max_threads();
        if ((size_t)string_threads > plan->n_string) string_threads = (int)plan->n_string;
        if (string_threads < 1) string_threads = 1;
    }
#ifndef _OPENMP
    /* Non-OpenMP build: one buffer, reused for every variable in turn */
    if (batch_size > 1) batch_size = 1;
    if (string_threads > 1) string_threads = 1;
#endif

    int nbufs = batch_size > string_threads ? batch_size : string_threads;
    if (nbufs == 0) {
        return STATA_OK;
    }
    plan->bufs = (char **)calloc((size_t)nbufs, sizeof(char *));
    if (!plan->bufs) {
        return STATA_ERR_MEMORY;
    }
    for (b = 0; b < nbufs; b++) {
        size_t bytes = 0, string_bytes = 0;
        if (b < batch_size && ctools_safe_mul_size(nobs, sizeof(double), &bytes) != 0) break;
        if (b < string_threads &&
            ctools_safe_mul_size(nobs, (size_t)plan->strings[b].width + 1, &string_bytes) != 0) break;
        if (string_bytes > bytes) bytes = string_bytes;
        plan->bufs[b] = (char *)ctools_aligned_alloc(CACHE_LINE_SIZE, bytes);
        if (!plan->bufs[b]) break;
    }
    plan->nbufs = b;
    if (b == 0) {
        return STATA_ERR_MEMORY;
    }
    plan->batch_size = batch_size < b ? batch_size : b;
    plan->string_threads = string_threads < b ? string_threads : b;
    return STATA_OK;
}

/*
    Permute a fixed-width string variable through a flat buffer of nobs slots
    of max_width + 1 bytes, scattering and gathering slots directly (no arena,
    strlen or pointer indirection). Stata's SPI reports only whether a
    variable is a string, so max_width is the ado's width hint or measured.
*/
static void stream_permute_string_var_flat(
    ST_int stata_var,
    const perm_idx_t *inv_perm,
    size_t nobs,
    const perm_idx_t *obs_map,
    int max_width,
    char *buf)
{
    size_t i;
    size_t stride = (size_t)max_width + 1;  /* +1 for null terminator */

    /* Read from Stata sequentially, scatter to permuted position in flat buffer */
    for (i = 0; i < nobs; i++) {
        char *slot = buf + inv_perm[i] * stride;
        slot[0] = '\0';  /* a failed read leaves "" */
        SF_sdata(stata_var, (ST_int)obs_map[i], slot);
    }

    /* Write back sequentially from flat buffer (cache-friendly) */
    for (i = 0; i < nobs; i++) {
        SF_sstore(stata_var, (ST_int)obs_map[i], buf + i * stride);
    }
}

/*
    Move every non-key variable into sorted order with the planned buffers.
    Nothing here allocates, so nothing here can fail.
*/
static void stream_plan_run(
    const stream_plan *plan,
    size_t nobs,
    const perm_idx_t *obs_map,
    csort_stream_timings *timings)
{
    const perm_idx_t *inv_perm = plan->inv_perm;
    const int *numeric_vars = plan->numeric_vars;
    const size_t n_numeric = plan->n_numeric;
    const size_t n_string = plan->n_string;
    double t_phase;

    /* Initialize timing accumulators */
    double t_scatter = 0.0;
    double t_writeback = 0.0;
    double t_string = 0.0;

    /* Early exit if data is already sorted */
    if (plan->is_identity) {
        if (timings) {
            timings->stream_build_inv_time = plan->build_inv_time;
            timings->stream_scatter_time = 0.0;
            timings->stream_writeback_time = 0.0;
            timings->stream_string_time = 0.0;
        }
        return;
    }

#ifdef HAVE_AVX512_RUNTIME
    init_scatter_fn();
#endif

    /* Process numeric variables in batches; batch slot b uses pool buffer b.
       plan->batch_size controls how many variables are loaded at a time */
    #ifdef _OPENMP
    int nthreads = 0;
    const int actual_batch_size = plan->batch_size;

    if (n_numeric > 0) {
        nthreads = ctools_get_max_threads();
        if (nthreads > actual_batch_size) nthreads = actual_batch_size;
        if (nthreads < 1) nthreads = 1;
    }

    const int stream_prefetch = ctools_get_prefetch_distances().stream;

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
                double *num_buf = (double *)(void *)plan->bufs[b];
                ST_int stata_var = (ST_int)numeric_vars[var_idx];
                size_t i;

                size_t nobs_batch16 = (nobs / 16) * 16;
                double batch_vals[16] __attribute__((aligned(64)));

                for (i = 0; i < nobs_batch16; i += 16) {
                    /* Prefetch future permutation indices for next iteration */
                    if (i + stream_prefetch < nobs) {
                        CTOOLS_PREFETCH(&inv_perm[i + stream_prefetch]);
                    }

                    /* Read data from Stata */
                    for (int bb = 0; bb < 16; bb++)
                        SF_vdata(stata_var, (ST_int)(obs_map[i] + bb), &batch_vals[bb]);

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
                double *num_buf = (double *)(void *)plan->bufs[b];
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
    }

    #else
    /* Non-OpenMP fallback: process sequentially with the single pool buffer */
    if (n_numeric > 0) {
        double *num_buf = (double *)(void *)plan->bufs[0];
        for (size_t v = 0; v < n_numeric; v++) {
            ST_int stata_var = (ST_int)numeric_vars[v];
            size_t i;
            size_t nobs_batch16 = (nobs / 16) * 16;
            double batch_vals[16] __attribute__((aligned(64)));
            double t_var_scatter, t_var_writeback;

            /* Phase 1: Sequential read, scatter via inverse perm */
            t_var_scatter = ctools_timer_seconds();

            for (i = 0; i < nobs_batch16; i += 16) {
                /* Read data from Stata */
                for (int bb = 0; bb < 16; bb++)
                    SF_vdata(stata_var, (ST_int)(obs_map[i] + bb), &batch_vals[bb]);

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
    }
    #endif

    /* Process string variables (parallel - each worker reuses its own pool buffer) */
    t_phase = ctools_timer_seconds();

    #ifdef _OPENMP
    if (n_string > 0) {
        #pragma omp parallel num_threads(plan->string_threads)
        {
            size_t team = (size_t)omp_get_num_threads();
            size_t tid = (size_t)omp_get_thread_num();
            /* Widest first: variable tid is the widest this worker meets */
            for (size_t sv = tid; sv < n_string; sv += team) {
                stream_permute_string_var_flat((ST_int)plan->strings[sv].var, inv_perm, nobs,
                                               obs_map, plan->strings[sv].width, plan->bufs[tid]);
            }
        }
    }
    #else
    for (size_t sv = 0; sv < n_string; sv++) {
        stream_permute_string_var_flat((ST_int)plan->strings[sv].var, inv_perm, nobs, obs_map,
                                       plan->strings[sv].width, plan->bufs[0]);
    }
    #endif

    t_string = ctools_timer_seconds() - t_phase;

    /* Store detailed timings */
    if (timings) {
        timings->stream_build_inv_time = plan->build_inv_time;
        timings->stream_scatter_time = t_scatter;
        timings->stream_writeback_time = t_writeback;
        timings->stream_string_time = t_string;
        timings->n_numeric_vars = n_numeric;
        timings->n_string_vars = n_string;
    }
}

stata_retcode csort_stream_apply_permutation(
    const perm_idx_t *perm,
    size_t nobs,
    int *nonkey_var_indices,
    size_t nvars_nonkey,
    const perm_idx_t *obs_map,
    size_t block_size,
    int vars_per_batch,
    const int *nonkey_str_widths,
    csort_stream_timings *timings)
{
    stream_plan plan;
    stata_retcode rc;
    (void)block_size;

    if (nvars_nonkey == 0 || nobs == 0) {
        return STATA_OK;
    }

    /* Allocate everything before the first write, then permute */
    memset(&plan, 0, sizeof(plan));
    rc = stream_plan_inverse(&plan, perm, nobs);
    if (rc == STATA_OK && !plan.is_identity) {
        rc = stream_plan_buffers(&plan, nobs, nonkey_var_indices, nvars_nonkey, obs_map,
                                 vars_per_batch, nonkey_str_widths);
    }
    if (rc == STATA_OK) {
        stream_plan_run(&plan, nobs, obs_map, timings);
    }
    stream_plan_free(&plan);
    return rc;
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
    const int *str_widths,
    csort_stream_timings *timings)
{
    for (ST_int v = 1; v <= SF_nvars(); v++) {
        if (SF_var_is_string(v) && SF_var_is_strl(v)) {
            SF_error("csort: strL I/O is unsupported; recast to str2045 when lossless\n");
            return STATA_ERR_UNSUPPORTED_TYPE;
        }
    }

    ctools_filtered_data key_filtered;
    perm_idx_t *obs_map = NULL;
    stata_retcode rc;
    double t_start, t_phase;
    double t_plan = 0.0;
    size_t nobs = 0;
    size_t nvars_nonkey = 0;
    size_t i, j;
    int *nonkey_var_indices = NULL;
    int *nonkey_str_widths = NULL;
    int *local_sort_vars = NULL;
    stream_plan plan;

    csort_stream_timings local_timings = {0};
    t_start = ctools_timer_seconds();
    memset(&plan, 0, sizeof(plan));

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

    /* Non-key variables: the listed variables that are not keys. They are
       counted, not taken as nvars - nkeys, which a key listed twice (or more
       keys than variables) would make too large or wrap. */
    nonkey_var_indices = (int *)ctools_safe_malloc2(nvars ? nvars : 1, sizeof(int));
    nonkey_str_widths = (int *)ctools_safe_calloc2(nvars ? nvars : 1, sizeof(int));
    if (!nonkey_var_indices || !nonkey_str_widths) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }
    for (i = 0; i < nvars; i++) {
        int is_key = 0;
        for (j = 0; j < nkeys; j++) {
            if (all_var_indices[i] == key_var_indices[j]) {
                is_key = 1;
                break;
            }
        }
        if (!is_key) {
            nonkey_var_indices[nvars_nonkey] = all_var_indices[i];
            /* Map global var index to str_widths: all_var_indices[i] is
               1-based, str_widths is 0-based indexed by variable position */
            if (str_widths) {
                nonkey_str_widths[nvars_nonkey] = str_widths[all_var_indices[i] - 1];
            }
            nvars_nonkey++;
        }
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
       Phase 2.5: Invert the permutation before applying it (Phase 3
       resets sort_order)
       ================================================================ */
    if (nvars_nonkey > 0) {
        t_phase = ctools_timer_seconds();
        rc = stream_plan_inverse(&plan, key_filtered.data.sort_order, nobs);
        t_plan += ctools_timer_seconds() - t_phase;
        if (rc != STATA_OK) {
            goto cleanup;
        }
    }

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
       Phase 3.5: Allocate every buffer Phase 5 needs. Nothing has been
       written to Stata yet, so a failure up to here leaves the dataset
       untouched; from Phase 4 on only SPI calls remain.
       ================================================================ */
    if (nvars_nonkey > 0 && !plan.is_identity) {
        t_phase = ctools_timer_seconds();
        rc = stream_plan_buffers(&plan, nobs, nonkey_var_indices, nvars_nonkey, obs_map,
                                 vars_per_batch, nonkey_str_widths);
        t_plan += ctools_timer_seconds() - t_phase;
        if (rc != STATA_OK) {
            goto cleanup;
        }
    }

    /* ================================================================
       Phase 4: Write sorted key variables back to Stata using obs_map
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Store each key variable using obs_map */
    for (size_t k = 0; k < nkeys && rc == STATA_OK; k++) {
        if (key_filtered.data.vars[k].type == STATA_TYPE_STRING) {
            rc = ctools_store_filtered_str(key_filtered.data.vars[k].data.str,
                                           nobs, key_var_indices[k], obs_map);
        } else {
            rc = ctools_store_filtered(key_filtered.data.vars[k].data.dbl,
                                       nobs, key_var_indices[k], obs_map);
        }
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

    if (nvars_nonkey > 0) {
        stream_plan_run(&plan, nobs, obs_map, &local_timings);
    }

    /* Includes planning (inverse permutation and buffers) done before Phase 4 */
    local_timings.stream_nonkeys_time = t_plan + (ctools_timer_seconds() - t_phase);

cleanup:
    stream_plan_free(&plan);
    free(local_sort_vars);
    free(nonkey_var_indices);
    free(nonkey_str_widths);
    if (obs_map) ctools_aligned_free(obs_map);  /* We took ownership from key_filtered */
    key_filtered.obs_map = NULL;  /* Prevent double-free by ctools_filtered_data_free */
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
