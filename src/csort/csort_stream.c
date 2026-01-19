/*
    csort_stream.c
    Memory-efficient streaming permutation for large dataset sorting

    This module provides an alternative sorting approach that minimizes memory
    usage by only loading key variables into C, then applying the permutation
    to non-key variables one at a time.

    Memory usage:
    - Standard mode: O(nobs * nvars * 8) bytes for all data
    - Memeff mode: O(nobs * nkeys * 8) + O(nobs * 8) bytes
      (keys + one variable buffer at a time)

    Best when:
    - Wide datasets with many columns
    - Limited memory situations
    - Large datasets that don't fit in RAM with standard approach
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
#include "csort_stream.h"

/* ============================================================================
   Configuration and Tuning Constants
   ============================================================================ */

/* Maximum memory to use for block buffers (default 1GB) */
#define STREAM_MAX_BUFFER_MEMORY (1024ULL * 1024 * 1024)

/* ============================================================================
   Main Streaming Permutation - Variable-at-a-time approach

   We process one variable at a time completely:
   1. Allocate O(nobs) buffer for one variable
   2. Read all values in permuted order: buffer[i] = Stata[perm[i]]
   3. Write all values back: Stata[i] = buffer[i]
   4. Free buffer, move to next variable

   This uses O(nobs) memory (for one variable) instead of O(nobs * nvars),
   which is a significant savings for wide datasets.

   Variables are processed in parallel using OpenMP when available.
   ============================================================================ */

stata_retcode csort_stream_apply_permutation(
    const perm_idx_t *perm,
    size_t nobs,
    int *nonkey_var_indices,
    size_t nvars_nonkey,
    size_t obs1,
    size_t block_size)
{
    size_t v;
    (void)block_size;  /* Not used in variable-at-a-time approach */

    if (nvars_nonkey == 0 || nobs == 0) {
        return STATA_OK;
    }

    /* Process each non-key variable completely before moving to next */
    #ifdef _OPENMP
    int success = 1;
    #pragma omp parallel for schedule(dynamic, 1) if(nvars_nonkey >= 2)
    #endif
    for (v = 0; v < nvars_nonkey; v++) {
        #ifdef _OPENMP
        if (!success) continue;
        #endif

        ST_int stata_var = (ST_int)nonkey_var_indices[v];
        int is_string = SF_var_is_string(stata_var);

        if (is_string) {
            /* String variable: allocate char* buffer */
            char **str_buf = (char **)calloc(nobs, sizeof(char *));
            if (!str_buf) {
                #ifdef _OPENMP
                #pragma omp atomic write
                success = 0;
                continue;
                #else
                return STATA_ERR_MEMORY;
                #endif
            }

            char strbuf[2048];
            /* Read all values in permuted order */
            for (size_t i = 0; i < nobs; i++) {
                SF_sdata(stata_var, (ST_int)(perm[i] + obs1), strbuf);
                str_buf[i] = strdup(strbuf);
                if (!str_buf[i]) {
                    /* Cleanup and fail */
                    for (size_t k = 0; k < i; k++) free(str_buf[k]);
                    free(str_buf);
                    #ifdef _OPENMP
                    #pragma omp atomic write
                    success = 0;
                    continue;
                    #else
                    return STATA_ERR_MEMORY;
                    #endif
                }
            }

            /* Write all values back in order */
            for (size_t i = 0; i < nobs; i++) {
                SF_sstore(stata_var, (ST_int)(i + obs1), str_buf[i] ? str_buf[i] : "");
            }

            /* Free string buffer */
            for (size_t i = 0; i < nobs; i++) free(str_buf[i]);
            free(str_buf);

        } else {
            /* Numeric variable: allocate double buffer */
            double *num_buf = (double *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                                               nobs * sizeof(double));
            if (!num_buf) {
                #ifdef _OPENMP
                #pragma omp atomic write
                success = 0;
                continue;
                #else
                return STATA_ERR_MEMORY;
                #endif
            }

            /* Read all values in permuted order */
            for (size_t i = 0; i < nobs; i++) {
                SF_vdata(stata_var, (ST_int)(perm[i] + obs1), &num_buf[i]);
            }

            /* Write all values back in order */
            for (size_t i = 0; i < nobs; i++) {
                SF_vstore(stata_var, (ST_int)(i + obs1), num_buf[i]);
            }

            /* Free buffer */
            ctools_aligned_free(num_buf);
        }
    }

    #ifdef _OPENMP
    if (!success) {
        return STATA_ERR_MEMORY;
    }
    #endif

    return STATA_OK;
}

/* ============================================================================
   Block Size Selection (kept for API compatibility, not currently used)
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

    /* Estimate bytes per row: doubles + overhead */
    bytes_per_row = nvars_nonkey * sizeof(double) + 64;  /* 64 bytes overhead per row */

    /* Maximum block size that fits in available memory */
    max_block = available_memory / bytes_per_row;

    /* Clamp to reasonable range */
    if (max_block < CSORT_STREAM_MIN_BLOCK_SIZE) {
        max_block = CSORT_STREAM_MIN_BLOCK_SIZE;
    }
    if (max_block > CSORT_STREAM_DEFAULT_BLOCK_SIZE * 4) {
        max_block = CSORT_STREAM_DEFAULT_BLOCK_SIZE * 4;  /* Cap at 256K */
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

    /* Estimate total memory for standard approach */
    total_bytes = nobs * nvars * sizeof(double);  /* Data */
    total_bytes += nobs * sizeof(perm_idx_t);     /* Permutation */
    total_bytes += nobs * nvars * sizeof(double); /* Permutation temp */

    /* Estimate memory for streaming approach */
    key_bytes = nobs * nkeys * sizeof(double);
    key_bytes += nobs * sizeof(perm_idx_t);

    /* Recommend streaming if:
       1. Total data exceeds threshold AND
       2. Streaming saves significant memory (keys are small fraction) */
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

    /* Initialize timings */
    csort_stream_timings local_timings = {0};
    t_start = ctools_timer_seconds();

    /* Get observation range */
    obs1 = SF_in1();
    nobs = SF_in2() - obs1 + 1;

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

    /* Create sort_var indices relative to key_data (1-based within keys) */
    local_sort_vars = (int *)malloc(nkeys * sizeof(int));
    if (!local_sort_vars) {
        rc = STATA_ERR_MEMORY;
        goto cleanup;
    }
    for (i = 0; i < nkeys; i++) {
        local_sort_vars[i] = (int)(i + 1);  /* Keys are vars 1..nkeys in key_data */
    }

    /* Sort using the selected algorithm (order_only to just compute permutation) */
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
    saved_perm = (perm_idx_t *)ctools_aligned_alloc(CACHE_LINE_SIZE,
                                                      nobs * sizeof(perm_idx_t));
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

    /* Free key data now - we don't need it anymore */
    stata_data_free(&key_data);
    stata_data_init(&key_data);  /* Reset to safe state */

    /* ================================================================
       Phase 5: Stream-permute non-key variables
       ================================================================ */
    t_phase = ctools_timer_seconds();

    /* Build list of non-key variable indices */
    size_t nvars_nonkey = nvars - nkeys;

    if (nvars_nonkey > 0) {
        nonkey_var_indices = (int *)malloc(nvars_nonkey * sizeof(int));
        if (!nonkey_var_indices) {
            rc = STATA_ERR_MEMORY;
            goto cleanup;
        }

        /* Find non-key variables: those in all_var_indices but not in key_var_indices */
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

        /* Apply permutation to non-key variables using streaming */
        rc = csort_stream_apply_permutation(saved_perm, nobs, nonkey_var_indices,
                                             nvars_nonkey, obs1, block_size);
    }

    local_timings.stream_nonkeys_time = ctools_timer_seconds() - t_phase;

cleanup:
    /* Clean up all allocated resources */
    free(local_sort_vars);
    free(nonkey_var_indices);
    ctools_aligned_free(saved_perm);
    stata_data_free(&key_data);

    /* Finalize timings */
    local_timings.total_time = ctools_timer_seconds() - t_start;
    local_timings.block_size = block_size ? block_size : CSORT_STREAM_DEFAULT_BLOCK_SIZE;
    local_timings.num_blocks = (nobs + local_timings.block_size - 1) / local_timings.block_size;

    if (timings) {
        *timings = local_timings;
    }

    /* Store timing results in Stata scalars */
    SF_scal_save("_csort_stream_time_load_keys", local_timings.load_keys_time);
    SF_scal_save("_csort_stream_time_sort", local_timings.sort_time);
    SF_scal_save("_csort_stream_time_permute_keys", local_timings.permute_keys_time);
    SF_scal_save("_csort_stream_time_store_keys", local_timings.store_keys_time);
    SF_scal_save("_csort_stream_time_stream_nonkeys", local_timings.stream_nonkeys_time);
    SF_scal_save("_csort_stream_time_total", local_timings.total_time);

    return rc;
}
