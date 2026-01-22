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

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_timer.h"
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
    size_t obs1)
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

    /* Sequential read from Stata, scatter to buffer via inverse perm */
    for (i = 0; i < nobs; i++) {
        SF_sdata(stata_var, (ST_int)(i + obs1), strbuf);
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

    /* Sequential write back to Stata */
    for (i = 0; i < nobs; i++) {
        SF_sstore(stata_var, (ST_int)(i + obs1), str_ptrs[i] ? str_ptrs[i] : "");
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
    size_t obs1,
    size_t block_size)
{
    size_t v;
    perm_idx_t *inv_perm = NULL;
    int is_identity;
    (void)block_size;

    if (nvars_nonkey == 0 || nobs == 0) {
        return STATA_OK;
    }

    /* Build inverse permutation once (shared across all variables) */
    inv_perm = (perm_idx_t *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE,
                                                         nobs, sizeof(perm_idx_t));
    if (!inv_perm) {
        return STATA_ERR_MEMORY;
    }

    is_identity = build_inverse_permutation(perm, inv_perm, nobs);

    /* Early exit if data is already sorted */
    if (is_identity) {
        ctools_aligned_free(inv_perm);
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

    /* Process numeric variables in parallel with shared buffer pool
       Each thread gets its own buffer to avoid allocation overhead */
    #ifdef _OPENMP
    int success = 1;
    double **thread_bufs = NULL;
    int nthreads = 0;

    if (n_numeric > 0) {
        nthreads = omp_get_max_threads();
        if (nthreads > (int)n_numeric) nthreads = (int)n_numeric;
        if (nthreads < 1) nthreads = 1;

        /* Pre-allocate buffers for each thread */
        thread_bufs = (double **)malloc(nthreads * sizeof(double *));
        if (!thread_bufs) {
            free(numeric_vars);
            free(string_vars);
            ctools_aligned_free(inv_perm);
            return STATA_ERR_MEMORY;
        }
        for (int t = 0; t < nthreads; t++) {
            thread_bufs[t] = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(double));
            if (!thread_bufs[t]) {
                for (int k = 0; k < t; k++) ctools_aligned_free(thread_bufs[k]);
                free(thread_bufs);
                free(numeric_vars);
                free(string_vars);
                ctools_aligned_free(inv_perm);
                return STATA_ERR_MEMORY;
            }
        }
    }

    if (n_numeric > 0) {
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (v = 0; v < n_numeric; v++) {
            int tid = omp_get_thread_num();
            double *num_buf = thread_bufs[tid];
            ST_int stata_var = (ST_int)numeric_vars[v];
            size_t i;

            /* Phase 1: Sequential read, scatter via inverse perm */
            size_t nobs_batch16 = (nobs / 16) * 16;
            double batch[16];

            for (i = 0; i < nobs_batch16; i += 16) {
                /* Prefetch future permutation indices for next iteration */
                if (i + STREAM_PREFETCH_DIST < nobs) {
                    CTOOLS_PREFETCH(&inv_perm[i + STREAM_PREFETCH_DIST]);
                }

                /* Read data from Stata */
                SF_VDATA_BATCH16(stata_var, (ST_int)(i + obs1), batch);

                /* Prefetch ALL 16 scatter destinations with WRITE intent before scattering.
                   This tells the CPU we intend to write to these cache lines, enabling
                   better write combining and avoiding read-for-ownership stalls. */
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

                /* Scatter writes to buffer */
                num_buf[inv_perm[i + 0]]  = batch[0];
                num_buf[inv_perm[i + 1]]  = batch[1];
                num_buf[inv_perm[i + 2]]  = batch[2];
                num_buf[inv_perm[i + 3]]  = batch[3];
                num_buf[inv_perm[i + 4]]  = batch[4];
                num_buf[inv_perm[i + 5]]  = batch[5];
                num_buf[inv_perm[i + 6]]  = batch[6];
                num_buf[inv_perm[i + 7]]  = batch[7];
                num_buf[inv_perm[i + 8]]  = batch[8];
                num_buf[inv_perm[i + 9]]  = batch[9];
                num_buf[inv_perm[i + 10]] = batch[10];
                num_buf[inv_perm[i + 11]] = batch[11];
                num_buf[inv_perm[i + 12]] = batch[12];
                num_buf[inv_perm[i + 13]] = batch[13];
                num_buf[inv_perm[i + 14]] = batch[14];
                num_buf[inv_perm[i + 15]] = batch[15];
            }
            for (i = nobs_batch16; i < nobs; i++) {
                double val;
                SF_vdata(stata_var, (ST_int)(i + obs1), &val);
                num_buf[inv_perm[i]] = val;
            }

            /* Phase 2: Sequential write back */
            for (i = 0; i < nobs_batch16; i += 16) {
                SF_VSTORE_BATCH16(stata_var, (ST_int)(i + obs1), &num_buf[i]);
            }
            for (i = nobs_batch16; i < nobs; i++) {
                SF_vstore(stata_var, (ST_int)(i + obs1), num_buf[i]);
            }
        }

        /* Free thread buffers */
        for (int t = 0; t < nthreads; t++) {
            ctools_aligned_free(thread_bufs[t]);
        }
        free(thread_bufs);
    }

    #else
    /* Non-OpenMP fallback: process sequentially with single reused buffer */
    double *num_buf = NULL;
    if (n_numeric > 0) {
        num_buf = (double *)ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(double));
    }
    if (num_buf) {
        for (v = 0; v < n_numeric; v++) {
            ST_int stata_var = (ST_int)numeric_vars[v];
            size_t i;
            size_t nobs_batch16 = (nobs / 16) * 16;
            double batch[16];

            /* Phase 1: Sequential read, scatter via inverse perm */
            for (i = 0; i < nobs_batch16; i += 16) {
                /* Read data from Stata */
                SF_VDATA_BATCH16(stata_var, (ST_int)(i + obs1), batch);

                /* Prefetch ALL 16 scatter destinations with WRITE intent */
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

                /* Scatter writes to buffer */
                num_buf[inv_perm[i + 0]]  = batch[0];
                num_buf[inv_perm[i + 1]]  = batch[1];
                num_buf[inv_perm[i + 2]]  = batch[2];
                num_buf[inv_perm[i + 3]]  = batch[3];
                num_buf[inv_perm[i + 4]]  = batch[4];
                num_buf[inv_perm[i + 5]]  = batch[5];
                num_buf[inv_perm[i + 6]]  = batch[6];
                num_buf[inv_perm[i + 7]]  = batch[7];
                num_buf[inv_perm[i + 8]]  = batch[8];
                num_buf[inv_perm[i + 9]]  = batch[9];
                num_buf[inv_perm[i + 10]] = batch[10];
                num_buf[inv_perm[i + 11]] = batch[11];
                num_buf[inv_perm[i + 12]] = batch[12];
                num_buf[inv_perm[i + 13]] = batch[13];
                num_buf[inv_perm[i + 14]] = batch[14];
                num_buf[inv_perm[i + 15]] = batch[15];
            }
            for (i = nobs_batch16; i < nobs; i++) {
                double val;
                SF_vdata(stata_var, (ST_int)(i + obs1), &val);
                num_buf[inv_perm[i]] = val;
            }

            /* Phase 2: Sequential write back */
            for (i = 0; i < nobs_batch16; i += 16) {
                SF_VSTORE_BATCH16(stata_var, (ST_int)(i + obs1), &num_buf[i]);
            }
            for (i = nobs_batch16; i < nobs; i++) {
                SF_vstore(stata_var, (ST_int)(i + obs1), num_buf[i]);
            }
        }
        ctools_aligned_free(num_buf);
    }
    #endif

    /* Process string variables (sequentially - strings are harder to parallelize) */
    for (v = 0; v < n_string; v++) {
        if (stream_permute_string_var((ST_int)string_vars[v], inv_perm, nobs, obs1) != 0) {
            #ifdef _OPENMP
            success = 0;
            #endif
        }
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
    csort_stream_timings *timings)
{
    stata_data key_data;
    stata_retcode rc;
    double t_start, t_phase;
    size_t obs1, nobs;
    size_t i, j;
    perm_idx_t *saved_perm = NULL;
    int *nonkey_var_indices = NULL;
    int *local_sort_vars = NULL;

    csort_stream_timings local_timings = {0};
    t_start = ctools_timer_seconds();

    /* Validate observation bounds */
    {
        ST_int in1 = SF_in1();
        ST_int in2 = SF_in2();
        if (in1 < 1 || in2 < in1) {
            /* Invalid range - treat as zero rows (nothing to sort) */
            if (timings) *timings = local_timings;
            return STATA_OK;
        }
        obs1 = (size_t)in1;
        nobs = (size_t)(in2 - in1 + 1);
    }

    /* ================================================================
       Phase 1: Load only key variables
       ================================================================ */
    t_phase = ctools_timer_seconds();

    stata_data_init(&key_data);
    rc = ctools_data_load_selective(&key_data, key_var_indices, nkeys, 0, 0);

    local_timings.load_keys_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
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

    switch (algorithm) {
        case SORT_ALG_MSD:
            rc = ctools_sort_radix_msd_order_only(&key_data, local_sort_vars, nkeys);
            break;
        case SORT_ALG_TIMSORT:
            rc = ctools_sort_timsort_order_only(&key_data, local_sort_vars, nkeys);
            break;
        case SORT_ALG_SAMPLE:
            rc = ctools_sort_sample_order_only(&key_data, local_sort_vars, nkeys);
            break;
        case SORT_ALG_COUNTING:
            rc = ctools_sort_counting_order_only(&key_data, local_sort_vars, nkeys);
            if (rc == STATA_ERR_UNSUPPORTED_TYPE) {
                rc = ctools_sort_radix_lsd_order_only(&key_data, local_sort_vars, nkeys);
            }
            break;
        case SORT_ALG_MERGE:
            rc = ctools_sort_merge_order_only(&key_data, local_sort_vars, nkeys);
            break;
        case SORT_ALG_LSD:
            rc = ctools_sort_radix_lsd_order_only(&key_data, local_sort_vars, nkeys);
            break;
        case SORT_ALG_IPS4O:
        default:
            rc = ctools_sort_ips4o_order_only(&key_data, local_sort_vars, nkeys);
            break;
    }

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
    memcpy(saved_perm, key_data.sort_order, nobs * sizeof(perm_idx_t));

    /* ================================================================
       Phase 3: Apply permutation to key variables (in C memory)
       ================================================================ */
    t_phase = ctools_timer_seconds();

    rc = ctools_apply_permutation(&key_data);

    local_timings.permute_keys_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
    }

    /* ================================================================
       Phase 4: Write sorted key variables back to Stata
       ================================================================ */
    t_phase = ctools_timer_seconds();

    rc = ctools_data_store_selective(&key_data, key_var_indices, nkeys, obs1);

    local_timings.store_keys_time = ctools_timer_seconds() - t_phase;

    if (rc != STATA_OK) {
        goto cleanup;
    }

    /* DEBUG: Verify first key is sorted in Stata after store */
    #ifdef CSORT_DEBUG_VERIFY
    {
        double prev_val, curr_val;
        int sorted_ok = 1;
        SF_vdata(key_var_indices[0], (ST_int)obs1, &prev_val);
        for (size_t check_i = 1; check_i < nobs && sorted_ok; check_i++) {
            SF_vdata(key_var_indices[0], (ST_int)(check_i + obs1), &curr_val);
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

    stata_data_free(&key_data);
    stata_data_init(&key_data);

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
                                             nvars_nonkey, obs1, block_size);
    }

    local_timings.stream_nonkeys_time = ctools_timer_seconds() - t_phase;

cleanup:
    free(local_sort_vars);
    free(nonkey_var_indices);
    ctools_aligned_free(saved_perm);
    stata_data_free(&key_data);

    local_timings.total_time = ctools_timer_seconds() - t_start;
    local_timings.block_size = block_size ? block_size : CSORT_STREAM_DEFAULT_BLOCK_SIZE;
    local_timings.num_blocks = (nobs + local_timings.block_size - 1) / local_timings.block_size;

    if (timings) {
        *timings = local_timings;
    }

    SF_scal_save("_csort_stream_time_load_keys", local_timings.load_keys_time);
    SF_scal_save("_csort_stream_time_sort", local_timings.sort_time);
    SF_scal_save("_csort_stream_time_permute_keys", local_timings.permute_keys_time);
    SF_scal_save("_csort_stream_time_store_keys", local_timings.store_keys_time);
    SF_scal_save("_csort_stream_time_stream_nonkeys", local_timings.stream_nonkeys_time);
    SF_scal_save("_csort_stream_time_total", local_timings.total_time);

    return rc;
}
