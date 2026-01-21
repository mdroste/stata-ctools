/*
    cmerge_impl.c
    Optimized C-Accelerated Merge for Stata

    High-performance merge with MINIMAL data transfer:
    - Only loads key variables + keepusing variables into C
    - Master non-key variables are STREAMED (read-permute-write) without loading
    - Memory usage: O(nobs * (nkeys + nkeepusing)) instead of O(nobs * all_vars)

    Architecture (Two-Phase Plugin Call):
    Phase 1 - "load_using":
        - Load ONLY key + keepusing variables from using dataset
        - Sort using data on key variables
        - Store in static cache for phase 2

    Phase 2 - "execute":
        - Load ONLY master key variables + _orig_row index
        - Perform sorted merge join to produce output row mapping
        - Write _merge and keepusing variables directly
        - STREAM master non-key variables (parallel gather-scatter)

    Merge Types: 1:1, m:1, 1:m, m:m - all supported with streaming
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <pthread.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_timer.h"
#include "ctools_threads.h"
#include "ctools_spi.h"
#include "ctools_config.h"
#include "cmerge_impl.h"

/* Helper modules */
#include "cmerge_memory.h"
#include "cmerge_radix_sort.h"
#include "cmerge_keys.h"
#include "cmerge_group_search.h"
#include "cmerge_join.h"
#include "cmerge_io.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CMERGE_MAX_KEYVARS 32
#define CMERGE_MAX_VARS 1024

#define cmerge_get_time_ms() ctools_timer_ms()

/* ============================================================================
 * Global using data cache (used by cmerge_io.c)
 * ============================================================================ */

cmerge_using_cache_t g_using_cache = {0};

/* Spinlock for protecting cache access during concurrent operations.
 * Prevents race conditions if Stata interrupts/restarts during execution. */
static volatile int g_cache_lock = 0;
static volatile int g_cache_in_use = 0;  /* Set while execute is running */

/* Platform-specific spinlock acquire/release */
static void cache_lock_acquire(void)
{
#if defined(_WIN32)
    while (InterlockedCompareExchange((volatile long *)&g_cache_lock, 1, 0) != 0) {
        /* Spin - brief contention only */
    }
#elif defined(__GNUC__) || defined(__clang__)
    while (__sync_lock_test_and_set(&g_cache_lock, 1)) {
        /* Spin */
    }
#else
    /* Fallback: no locking */
#endif
}

static void cache_lock_release(void)
{
#if defined(_WIN32)
    InterlockedExchange((volatile long *)&g_cache_lock, 0);
#elif defined(__GNUC__) || defined(__clang__)
    __sync_lock_release(&g_cache_lock);
#else
    g_cache_lock = 0;
#endif
}

/* Helper to safely clear the cache in-use flag on error paths */
static void cache_clear_in_use(void)
{
    cache_lock_acquire();
    g_cache_in_use = 0;
    cache_lock_release();
}

/* ============================================================================
 * Phase 1: Load using dataset (keys + keepusing only)
 * ============================================================================ */

static ST_retcode cmerge_load_using(const char *args)
{
    ST_retcode rc;

    /* Parse: "load_using <nkeys> <key_idx1>...<key_idxN> n_keepusing <count>
              keepusing_indices <idx1>...<idxM> [verbose]" */
    char args_copy[4096];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    int nkeys = 0;
    int key_indices[CMERGE_MAX_KEYVARS];
    int n_keepusing = 0;
    int keepusing_indices[CMERGE_MAX_VARS];
    int verbose = 0;
    int sorted = 0;
    int merge_by_n = 0;

    char *token = strtok(args_copy, " ");
    int arg_idx = 0;
    int in_keepusing_indices = 0;
    int keepusing_idx = 0;

    while (token != NULL) {
        if (arg_idx == 0) {
            nkeys = atoi(token);
            if (nkeys > CMERGE_MAX_KEYVARS) nkeys = CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx <= nkeys) {
            key_indices[arg_idx - 1] = atoi(token);
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            if (token) n_keepusing = atoi(token);
        }
        else if (strcmp(token, "keepusing_indices") == 0) {
            in_keepusing_indices = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_indices && keepusing_idx < n_keepusing) {
            keepusing_indices[keepusing_idx++] = atoi(token);
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "merge_by_n") == 0) {
            merge_by_n = 1;
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    (void)verbose;  /* Timing output moved to ado file */

    /* Acquire lock before modifying global cache */
    cache_lock_acquire();

    /* Check if cache is currently in use by execute */
    if (g_cache_in_use) {
        cache_lock_release();
        SF_error("cmerge: Cannot load using data while merge is in progress\n");
        return 459;
    }

    /* Free any previous cache */
    if (g_using_cache.loaded) {
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }

    cache_lock_release();

    double t_phase1_start = cmerge_get_time_ms();

    double t_load_keys = 0;
    double t_load_keepusing = 0;

    /* For merge_by_n, we don't need key variables - just keepusing */
    if (merge_by_n) {
        /* Initialize empty keys structure */
        stata_data_init(&g_using_cache.keys);
        g_using_cache.keys.nobs = SF_nobs();  /* Get nobs from Stata */

        /* Load keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = cmerge_get_time_ms();
            rc = ctools_data_load_selective(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0);
            if (rc != STATA_OK) {
                SF_error("cmerge: Failed to load keepusing variables\n");
                return 459;
            }
            if (g_using_cache.keepusing.vars == NULL) {
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = cmerge_get_time_ms() - t_load_keepusing_start;
            g_using_cache.keys.nobs = g_using_cache.keepusing.nobs;
        } else {
            stata_data_init(&g_using_cache.keepusing);
        }
    } else {
        /* Standard merge: Load ONLY key variables */
        double t_load_keys_start = cmerge_get_time_ms();
        rc = ctools_data_load_selective(&g_using_cache.keys, key_indices, nkeys, 0, 0);
        if (rc != STATA_OK) {
            SF_error("cmerge: Failed to load using keys\n");
            return 459;
        }
        /* Critical null check - prevents crash if load failed silently */
        if (g_using_cache.keys.vars == NULL) {
            SF_error("cmerge: FATAL - keys.vars is NULL after load\n");
            return 920;
        }
        t_load_keys = cmerge_get_time_ms() - t_load_keys_start;

        /* Load ONLY keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = cmerge_get_time_ms();
            rc = ctools_data_load_selective(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0);
            if (rc != STATA_OK) {
                stata_data_free(&g_using_cache.keys);
                SF_error("cmerge: Failed to load keepusing variables\n");
                return 459;
            }
            /* Critical null check */
            if (g_using_cache.keepusing.vars == NULL) {
                stata_data_free(&g_using_cache.keys);
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = cmerge_get_time_ms() - t_load_keepusing_start;
        } else {
            stata_data_init(&g_using_cache.keepusing);
        }
    }

    g_using_cache.nobs = g_using_cache.keys.nobs;
    g_using_cache.nkeys = nkeys;
    g_using_cache.n_keepusing = n_keepusing;
    g_using_cache.merge_by_n = merge_by_n;

    double t_sort = 0;
    double t_apply_perm = 0;

    /* Sort using data on keys (skip if sorted option specified or merge_by_n) */
    if (!sorted && !merge_by_n) {
        double t_sort_start = cmerge_get_time_ms();

        int *sort_vars = malloc(nkeys * sizeof(int));
        if (!sort_vars) {
            stata_data_free(&g_using_cache.keys);
            stata_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Memory allocation failed for sort_vars\n");
            return 920;
        }
        for (int i = 0; i < nkeys; i++) {
            sort_vars[i] = i + 1;  /* 1-based within keys structure */
        }

        /* Allocate permutation array to capture ordering before it's applied (overflow-safe) */
        size_t nobs = g_using_cache.nobs;
        size_t *perm = ctools_safe_malloc2(nobs, sizeof(size_t));
        if (!perm) {
            free(sort_vars);
            stata_data_free(&g_using_cache.keys);
            stata_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Memory allocation failed for permutation\n");
            return 920;
        }

        /* Use _with_perm version to get permutation BEFORE it's reset to identity */
        rc = ctools_sort_ips4o_with_perm(&g_using_cache.keys, sort_vars, nkeys, perm);
        if (rc != STATA_OK) {
            free(perm);
            free(sort_vars);
            stata_data_free(&g_using_cache.keys);
            stata_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Failed to sort using data\n");
            return 459;
        }

        t_sort = cmerge_get_time_ms() - t_sort_start;

        /* Apply same permutation to keepusing */
        double t_apply_perm_start = cmerge_get_time_ms();
        if (n_keepusing > 0 && g_using_cache.keepusing.vars != NULL) {
            for (int v = 0; v < n_keepusing; v++) {
                stata_variable *var = &g_using_cache.keepusing.vars[v];
                if (var->type == STATA_TYPE_DOUBLE) {
                    double *new_data = ctools_safe_malloc2(nobs, sizeof(double));
                    if (!new_data) {
                        free(perm);
                        free(sort_vars);
                        stata_data_free(&g_using_cache.keys);
                        stata_data_free(&g_using_cache.keepusing);
                        SF_error("cmerge: Memory allocation failed in permutation\n");
                        return 920;
                    }
                    for (size_t i = 0; i < nobs; i++) {
                        /* Bounds check on permutation index */
                        if (perm[i] >= nobs) {
                            free(new_data);
                            free(perm);
                            free(sort_vars);
                            stata_data_free(&g_using_cache.keys);
                            stata_data_free(&g_using_cache.keepusing);
                            SF_error("cmerge: Invalid permutation index\n");
                            return 459;
                        }
                        new_data[i] = var->data.dbl[perm[i]];
                    }
                    free(var->data.dbl);
                    var->data.dbl = new_data;
                } else {
                    char **new_data = ctools_safe_malloc2(nobs, sizeof(char *));
                    if (!new_data) {
                        free(perm);
                        free(sort_vars);
                        stata_data_free(&g_using_cache.keys);
                        stata_data_free(&g_using_cache.keepusing);
                        SF_error("cmerge: Memory allocation failed in permutation\n");
                        return 920;
                    }
                    for (size_t i = 0; i < nobs; i++) {
                        /* Bounds check on permutation index */
                        if (perm[i] >= nobs) {
                            free(new_data);
                            free(perm);
                            free(sort_vars);
                            stata_data_free(&g_using_cache.keys);
                            stata_data_free(&g_using_cache.keepusing);
                            SF_error("cmerge: Invalid permutation index\n");
                            return 459;
                        }
                        new_data[i] = var->data.str[perm[i]];
                    }
                    free(var->data.str);
                    var->data.str = new_data;
                }
            }
        }
        t_apply_perm = cmerge_get_time_ms() - t_apply_perm_start;

        free(perm);
        free(sort_vars);
    }

    g_using_cache.loaded = 1;

    double t_phase1_total = cmerge_get_time_ms() - t_phase1_start;

    /* Save Phase 1 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_using_nobs", (double)g_using_cache.nobs);
    SF_scal_save("_cmerge_p1_load_keys", t_load_keys / 1000.0);
    SF_scal_save("_cmerge_p1_load_keepusing", t_load_keepusing / 1000.0);
    SF_scal_save("_cmerge_p1_sort", t_sort / 1000.0);
    SF_scal_save("_cmerge_p1_apply_perm", t_apply_perm / 1000.0);
    SF_scal_save("_cmerge_p1_total", t_phase1_total / 1000.0);
    CTOOLS_SAVE_THREAD_INFO("_cmerge");

    return 0;
}

/* ============================================================================
 * Phase 2: Execute merge with full load/permute/store
 *
 * This implementation loads ALL master variables into C memory, applies the
 * merge permutation in C, then stores everything back using parallel I/O.
 * This is more memory-intensive than streaming but enables parallel I/O.
 * ============================================================================ */

static ST_retcode cmerge_execute(const char *args)
{
    ST_retcode rc;
    double t_start, t_total_start;

    t_total_start = cmerge_get_time_ms();

    /* Acquire lock and mark cache as in-use to prevent concurrent clear */
    cache_lock_acquire();

    if (!g_using_cache.loaded) {
        cache_lock_release();
        SF_error("cmerge: Using data not loaded. Call load_using first.\n");
        return 459;
    }

    /* Critical null check - prevents crash if cache was corrupted
       Note: For merge_by_n mode, keys.vars is intentionally NULL */
    if (g_using_cache.keys.vars == NULL && !g_using_cache.merge_by_n) {
        cache_lock_release();
        SF_error("cmerge: FATAL - g_using_cache.keys.vars is NULL!\n");
        return 459;
    }

    /* Mark cache as in-use before releasing lock */
    g_cache_in_use = 1;
    cache_lock_release();

    /* Parse arguments - use static buffers to avoid stack overflow.
     * Thread safety: Stata plugins are single-threaded, so these statics
     * cannot race. They're reinitialized on each call. */
    static char args_copy[8192];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    cmerge_type_t merge_type = MERGE_1_1;
    int nkeys = 0;
    static int master_key_indices[CMERGE_MAX_KEYVARS];
    int orig_row_idx = 0;
    size_t master_nobs = 0;
    size_t master_nvars = 0;
    int n_keepusing = 0;
    static int keepusing_placeholder_indices[CMERGE_MAX_VARS];
    int merge_var_idx = 0;
    int preserve_order = 0;
    int verbose = 0;
    int sorted = 0;
    int update_mode = 0;
    int replace_mode = 0;
    static int shared_flags[CMERGE_MAX_VARS];

    int arg_idx = 0;
    int in_keepusing_placeholders = 0;
    int keepusing_idx = 0;

    char *token = strtok(args_copy, " ");
    while (token != NULL) {
        if (arg_idx == 0) {
            merge_type = (cmerge_type_t)atoi(token);
        }
        else if (arg_idx == 1) {
            nkeys = atoi(token);
            if (nkeys > CMERGE_MAX_KEYVARS) nkeys = CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx >= 2 && arg_idx < 2 + nkeys) {
            master_key_indices[arg_idx - 2] = atoi(token);
        }
        else if (strcmp(token, "orig_row_idx") == 0) {
            token = strtok(NULL, " ");
            if (token) orig_row_idx = atoi(token);
        }
        else if (strcmp(token, "master_nobs") == 0) {
            token = strtok(NULL, " ");
            if (token) master_nobs = (size_t)atol(token);
        }
        else if (strcmp(token, "master_nvars") == 0) {
            token = strtok(NULL, " ");
            if (token) master_nvars = (size_t)atol(token);
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            if (token) n_keepusing = atoi(token);
        }
        else if (strcmp(token, "keepusing_placeholders") == 0) {
            in_keepusing_placeholders = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_placeholders && keepusing_idx < n_keepusing) {
            keepusing_placeholder_indices[keepusing_idx++] = atoi(token);
        }
        else if (strcmp(token, "merge_var_idx") == 0) {
            token = strtok(NULL, " ");
            if (token) merge_var_idx = atoi(token);
            in_keepusing_placeholders = 0;
        }
        else if (strcmp(token, "preserve_order") == 0) {
            token = strtok(NULL, " ");
            if (token) preserve_order = atoi(token);
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "update") == 0) {
            update_mode = 1;
        }
        else if (strcmp(token, "replace") == 0) {
            replace_mode = 1;
        }
        else if (strcmp(token, "shared_flags") == 0) {
            for (int i = 0; i < n_keepusing; i++) {
                token = strtok(NULL, " ");
                shared_flags[i] = token ? atoi(token) : 0;
            }
        }
        else if (strcmp(token, "merge_by_n") == 0) {
            /* Parsed but value comes from cache */
        }
        arg_idx++;
        token = strtok(NULL, " ");
    }

    /* Get merge_by_n from cache (set during load_using phase) */
    int merge_by_n = g_using_cache.merge_by_n;

    (void)verbose;  /* Timing output moved to ado file */

    /* ===================================================================
     * Step 1: Load ALL master variables (parallel I/O)
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    /* Build array of all variable indices to load */
    int *all_var_indices = malloc(master_nvars * sizeof(int));
    if (!all_var_indices) {
        cache_clear_in_use();
        SF_error("cmerge: Failed to allocate var indices\n");
        return 920;
    }
    for (size_t v = 0; v < master_nvars; v++) {
        all_var_indices[v] = (int)(v + 1);
    }

    /* Load all master variables using parallel I/O */
    stata_data master_data;
    rc = ctools_data_load_selective(&master_data, all_var_indices, master_nvars, 1, master_nobs);
    if (rc != STATA_OK) {
        free(all_var_indices);
        cache_clear_in_use();
        SF_error("cmerge: Failed to load master data\n");
        return rc;
    }

    /* Critical null check - prevents crash if load failed silently */
    if (master_data.vars == NULL) {
        free(all_var_indices);
        cache_clear_in_use();
        SF_error("cmerge: FATAL - master_data.vars is NULL after load!\n");
        return 920;
    }
    double t_load = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 2: Sort master data on keys (skip if sorted option specified or merge_by_n)
     * =================================================================== */

    double t_sort = 0;
    size_t *sort_perm = NULL;

    if (!sorted && !merge_by_n) {
        t_start = cmerge_get_time_ms();

        /* Build sort variable indices (within master_data, 1-based) */
        int *sort_vars = malloc(nkeys * sizeof(int));
        if (!sort_vars) {
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            SF_error("cmerge: Failed to allocate sort_vars\n");
            return 920;
        }

        for (int i = 0; i < nkeys; i++) {
            sort_vars[i] = master_key_indices[i];  /* Already 1-based Stata index */
        }

        /* Allocate permutation array to capture sort order (overflow-safe) */
        sort_perm = ctools_safe_malloc2(master_nobs, sizeof(size_t));
        if (!sort_perm) {
            free(sort_vars);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            return 920;
        }

        rc = ctools_sort_ips4o_with_perm(&master_data, sort_vars, nkeys, sort_perm);
        free(sort_vars);

        if (rc != STATA_OK) {
            free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            return rc;
        }

        t_sort = cmerge_get_time_ms() - t_start;
    }

    /* ===================================================================
     * Step 3: Perform merge join
     * For merge_by_n: simple row-by-row matching
     * For standard merge: sorted join using keys
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    cmerge_output_spec_t *output_specs = NULL;
    size_t output_nobs = 0;

    if (merge_by_n) {
        /* Merge by observation number (_n): simple row-by-row matching */
        size_t m_nobs = master_nobs;
        size_t u_nobs = g_using_cache.nobs;
        size_t max_nobs = (m_nobs > u_nobs) ? m_nobs : u_nobs;

        output_specs = ctools_safe_malloc2(max_nobs, sizeof(cmerge_output_spec_t));
        if (!output_specs) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            return 920;
        }

        /* Build output specs: row i matches row i */
        for (size_t i = 0; i < max_nobs; i++) {
            if (i < m_nobs && i < u_nobs) {
                /* Both master and using have this row - matched */
                output_specs[i].master_sorted_row = (int64_t)i;
                output_specs[i].using_sorted_row = (int64_t)i;
                output_specs[i].merge_result = MERGE_RESULT_BOTH;
            } else if (i < m_nobs) {
                /* Only master has this row */
                output_specs[i].master_sorted_row = (int64_t)i;
                output_specs[i].using_sorted_row = -1;
                output_specs[i].merge_result = MERGE_RESULT_MASTER_ONLY;
            } else {
                /* Only using has this row */
                output_specs[i].master_sorted_row = -1;
                output_specs[i].using_sorted_row = (int64_t)i;
                output_specs[i].merge_result = MERGE_RESULT_USING_ONLY;
            }
        }
        output_nobs = max_nobs;
    } else {
        /* Standard merge: create keys view and perform sorted join */
        stata_data master_keys_view;
        stata_data_init(&master_keys_view);
        master_keys_view.nobs = master_data.nobs;
        master_keys_view.nvars = nkeys;
        master_keys_view.vars = malloc(nkeys * sizeof(stata_variable));
        if (!master_keys_view.vars) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            return 920;
        }

        /* Point to the key variables in master_data (0-indexed in vars array) */
        for (int k = 0; k < nkeys; k++) {
            int var_idx = master_key_indices[k] - 1;  /* Convert to 0-based */

            /* Critical bounds check - prevents crash from bad indices */
            if (var_idx < 0 || (size_t)var_idx >= master_data.nvars) {
                SF_error("cmerge: key variable index out of bounds\n");
                free(master_keys_view.vars);
                if (sort_perm) free(sort_perm);
                stata_data_free(&master_data);
                free(all_var_indices);
                cache_clear_in_use();
                return 459;
            }

            master_keys_view.vars[k] = master_data.vars[var_idx];
        }

        int64_t output_nobs_signed = cmerge_sorted_join(
            &master_keys_view, &g_using_cache.keys,
            nkeys, merge_type, &output_specs);

        /* Free the view (but not the underlying data which belongs to master_data) */
        free(master_keys_view.vars);

        if (output_nobs_signed < 0) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            SF_error("cmerge: Merge join failed\n");
            return 920;
        }

        /* Stata max obs is 2^31-1; reject if exceeded */
        if (output_nobs_signed > INT32_MAX) {
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            SF_error("cmerge: Merge result exceeds Stata observation limit\n");
            return 920;
        }

        output_nobs = (size_t)output_nobs_signed;
    }

    double t_merge = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 4: Build master_orig_row mapping using _orig_row variable
     * =================================================================== */

    int64_t *master_orig_rows = ctools_safe_malloc2(output_nobs, sizeof(int64_t));
    if (!master_orig_rows) {
        free(output_specs);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        return 920;
    }

    /* _orig_row is at index orig_row_idx-1 in master_data.vars */
    /* Bounds check for orig_row_idx */
    if (orig_row_idx <= 0 || (size_t)(orig_row_idx - 1) >= master_data.nvars) {
        free(output_specs);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        free(master_orig_rows);
        cache_clear_in_use();
        SF_error("cmerge: orig_row_idx out of bounds\n");
        return 459;
    }

    double *orig_row_data = master_data.vars[orig_row_idx - 1].data.dbl;

    for (size_t i = 0; i < output_nobs; i++) {
        if (output_specs[i].master_sorted_row >= 0) {
            size_t sorted_idx = (size_t)output_specs[i].master_sorted_row;
            master_orig_rows[i] = (int64_t)orig_row_data[sorted_idx] - 1;  /* 0-based */
        } else {
            master_orig_rows[i] = -1;  /* Using-only */
        }
    }

    /* ===================================================================
     * Step 5: Apply preserve_order if requested
     * =================================================================== */

    double t_reorder = 0;
    if (preserve_order && output_nobs > 0) {
        double t_reorder_start = cmerge_get_time_ms();

        cmerge_order_pair_t *pairs = ctools_safe_malloc2(output_nobs, sizeof(cmerge_order_pair_t));
        if (!pairs) {
            free(output_specs);
            free(master_orig_rows);
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            return 920;
        }

        for (size_t i = 0; i < output_nobs; i++) {
            pairs[i].orig_idx = i;
            if (master_orig_rows[i] >= 0) {
                pairs[i].order_key = (size_t)master_orig_rows[i];
            } else {
                pairs[i].order_key = master_nobs + (size_t)output_specs[i].using_sorted_row;
            }
        }

        /* Use radix sort for O(n) performance instead of O(n log n) qsort */
        cmerge_radix_sort_order_pairs(pairs, output_nobs);

        cmerge_output_spec_t *new_specs = ctools_safe_malloc2(output_nobs, sizeof(cmerge_output_spec_t));
        int64_t *new_orig_rows = ctools_safe_malloc2(output_nobs, sizeof(int64_t));
        if (!new_specs || !new_orig_rows) {
            free(pairs);
            if (new_specs) free(new_specs);
            if (new_orig_rows) free(new_orig_rows);
            free(output_specs);
            free(master_orig_rows);
            if (sort_perm) free(sort_perm);
            stata_data_free(&master_data);
            free(all_var_indices);
            cache_clear_in_use();
            return 920;
        }

        for (size_t i = 0; i < output_nobs; i++) {
            size_t src_idx = pairs[i].orig_idx;
            new_specs[i] = output_specs[src_idx];
            new_orig_rows[i] = master_orig_rows[src_idx];
        }

        free(output_specs);
        free(master_orig_rows);
        output_specs = new_specs;
        master_orig_rows = new_orig_rows;
        free(pairs);

        t_reorder = cmerge_get_time_ms() - t_reorder_start;
    }

    /* Count merge results */
    size_t n_master_only = 0, n_using_only = 0, n_matched = 0;
    for (size_t i = 0; i < output_nobs; i++) {
        int8_t m = output_specs[i].merge_result;
        n_master_only += (m == 1);
        n_using_only += (m == 2);
        n_matched += (m == 3);
    }

    /* ===================================================================
     * Step 6: Create output data structure with permuted data
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    /* Determine which variables to include in output (exclude _orig_row and shared keepusing) */
    size_t n_output_vars = 0;
    int *output_var_indices = malloc(master_nvars * sizeof(int));
    int *output_var_stata_idx = malloc(master_nvars * sizeof(int));
    int *output_var_is_key = malloc(master_nvars * sizeof(int));
    int *output_var_key_idx = malloc(master_nvars * sizeof(int));

    if (!output_var_indices || !output_var_stata_idx || !output_var_is_key || !output_var_key_idx) {
        if (output_var_indices) free(output_var_indices);
        if (output_var_stata_idx) free(output_var_stata_idx);
        if (output_var_is_key) free(output_var_is_key);
        if (output_var_key_idx) free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        return 920;
    }

    for (size_t v = 0; v < master_nvars; v++) {
        int stata_idx = (int)(v + 1);
        if (stata_idx == orig_row_idx) continue;  /* Skip _orig_row */

        /* Check if shared keepusing (handled separately) */
        int is_shared_keepusing = 0;
        for (int kv = 0; kv < n_keepusing; kv++) {
            if (keepusing_placeholder_indices[kv] == stata_idx) {
                is_shared_keepusing = 1;
                break;
            }
        }
        if (is_shared_keepusing) continue;

        output_var_indices[n_output_vars] = (int)v;  /* Index in master_data.vars */
        output_var_stata_idx[n_output_vars] = stata_idx;

        /* Check if key variable */
        output_var_is_key[n_output_vars] = 0;
        output_var_key_idx[n_output_vars] = -1;
        for (int k = 0; k < nkeys; k++) {
            if (master_key_indices[k] == stata_idx) {
                output_var_is_key[n_output_vars] = 1;
                output_var_key_idx[n_output_vars] = k;
                break;
            }
        }

        n_output_vars++;
    }

    /* Create output data structure (overflow-safe allocation) */
    stata_data output_data;
    stata_data_init(&output_data);
    output_data.nobs = output_nobs;
    output_data.nvars = n_output_vars;
    size_t vars_alloc_size;
    if (ctools_safe_mul_size(n_output_vars, sizeof(stata_variable), &vars_alloc_size) != 0) {
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        return 920;
    }
    output_data.vars = (stata_variable *)cmerge_aligned_alloc(vars_alloc_size);
    if (!output_data.vars) {
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        return 920;
    }
    memset(output_data.vars, 0, n_output_vars * sizeof(stata_variable));

    /* Count string variables and estimate arena size */
    size_t n_string_vars = 0;
    for (size_t vi = 0; vi < n_output_vars; vi++) {
        int src_var_idx = output_var_indices[vi];
        if (master_data.vars[src_var_idx].type == STATA_TYPE_STRING) {
            n_string_vars++;
        }
    }

    /* Create string arena for merge output: estimate 64 bytes per string cell.
     * Arena enables O(1) bulk free instead of O(n*m) individual frees. */
    cmerge_string_arena *str_arena = NULL;
    if (n_string_vars > 0) {
        /* Overflow check for arena_size = n_string_vars * output_nobs * 64 */
        size_t arena_size = 0;
        if (n_string_vars > 0 && output_nobs > 0) {
            /* Check n_string_vars * output_nobs */
            if (n_string_vars > SIZE_MAX / output_nobs) {
                /* Overflow - use a smaller fallback size */
                arena_size = SIZE_MAX / 64;  /* Will likely fail, triggering strdup fallback */
            } else {
                size_t cell_count = n_string_vars * output_nobs;
                /* Check cell_count * 64 */
                if (cell_count > SIZE_MAX / 64) {
                    arena_size = SIZE_MAX;  /* Will likely fail */
                } else {
                    arena_size = cell_count * 64;
                }
            }
        }
        str_arena = cmerge_arena_create(arena_size);
        /* If arena creation fails, we fall back to strdup (no error) */
    }

    /* Track allocation failures in parallel loop (must be int for OpenMP atomic) */
    int alloc_failed = 0;

    /* Apply permutation to each variable (can be parallelized with OpenMP) */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t vi = 0; vi < n_output_vars; vi++) {
        /* Skip remaining work if allocation already failed - use atomic read */
        int local_alloc_failed;
        #ifdef _OPENMP
        #pragma omp atomic read
        #endif
        local_alloc_failed = alloc_failed;
        if (local_alloc_failed) continue;

        int src_var_idx = output_var_indices[vi];
        stata_variable *src_var = &master_data.vars[src_var_idx];
        stata_variable *dst_var = &output_data.vars[vi];

        dst_var->nobs = output_nobs;
        dst_var->type = src_var->type;
        dst_var->_arena = NULL;

        if (src_var->type == STATA_TYPE_DOUBLE) {
            size_t dbl_alloc_size;
            if (ctools_safe_mul_size(output_nobs, sizeof(double), &dbl_alloc_size) != 0) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }
            dst_var->data.dbl = (double *)cmerge_aligned_alloc(dbl_alloc_size);
            if (!dst_var->data.dbl) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }

            int is_key = output_var_is_key[vi];
            int key_idx = output_var_key_idx[vi];
            double *using_key_data = (is_key && g_using_cache.loaded && key_idx >= 0) ?
                                      g_using_cache.keys.vars[key_idx].data.dbl : NULL;

            /* Use master_sorted_row to index into sorted master_data */
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t sorted_row = output_specs[i].master_sorted_row;
                if (sorted_row >= 0) {
                    dst_var->data.dbl[i] = src_var->data.dbl[sorted_row];
                } else if (using_key_data != NULL) {
                    int64_t using_row = output_specs[i].using_sorted_row;
                    dst_var->data.dbl[i] = (using_row >= 0) ? using_key_data[using_row] : SV_missval;
                } else {
                    dst_var->data.dbl[i] = SV_missval;
                }
            }
        } else {
            /* String variable - use arena allocator to reduce malloc overhead */
            dst_var->data.str = (char **)ctools_safe_calloc2(output_nobs, sizeof(char *));
            if (!dst_var->data.str) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }

            /* Mark this variable as using the shared arena */
            dst_var->_arena = str_arena;

            int is_key = output_var_is_key[vi];
            int key_idx = output_var_key_idx[vi];
            char **using_key_data = (is_key && g_using_cache.loaded && key_idx >= 0) ?
                                     g_using_cache.keys.vars[key_idx].data.str : NULL;

            /* Use master_sorted_row to index into sorted master_data */
            for (size_t i = 0; i < output_nobs; i++) {
                int64_t sorted_row = output_specs[i].master_sorted_row;
                if (sorted_row >= 0 && src_var->data.str[sorted_row]) {
                    dst_var->data.str[i] = cmerge_arena_strdup(str_arena, src_var->data.str[sorted_row]);
                } else if (using_key_data != NULL) {
                    int64_t using_row = output_specs[i].using_sorted_row;
                    if (using_row >= 0 && using_key_data[using_row]) {
                        dst_var->data.str[i] = cmerge_arena_strdup(str_arena, using_key_data[using_row]);
                    } else {
                        dst_var->data.str[i] = cmerge_arena_strdup(str_arena, "");
                    }
                } else {
                    dst_var->data.str[i] = cmerge_arena_strdup(str_arena, "");
                }
            }
        }
    }

    /* Check for allocation failures in the parallel loop */
    if (alloc_failed) {
        SF_error("cmerge: memory allocation failed during output assembly\n");
        /* Cleanup - free any successfully allocated buffers */
        for (size_t vi = 0; vi < n_output_vars; vi++) {
            if (output_data.vars[vi].type == STATA_TYPE_STRING && output_data.vars[vi].data.str) {
                /* Only free strings not owned by arena */
                for (size_t i = 0; i < output_nobs; i++) {
                    if (output_data.vars[vi].data.str[i] &&
                        !cmerge_arena_owns(str_arena, output_data.vars[vi].data.str[i])) {
                        free(output_data.vars[vi].data.str[i]);
                    }
                }
                free(output_data.vars[vi].data.str);
            } else if (output_data.vars[vi].data.dbl) {
                cmerge_aligned_free(output_data.vars[vi].data.dbl);
            }
        }
        cmerge_arena_free(str_arena);
        cmerge_aligned_free(output_data.vars);
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        return 920;
    }

    /* Check for allocation failures during string copying */
    if (str_arena != NULL && str_arena->alloc_failed) {
        SF_error("cmerge: memory allocation failed during string copy\n");
        /* Cleanup - free all allocated resources */
        for (size_t vi = 0; vi < n_output_vars; vi++) {
            if (output_data.vars[vi].type == STATA_TYPE_STRING && output_data.vars[vi].data.str) {
                /* Free fallback strings not owned by arena */
                for (size_t i = 0; i < output_nobs; i++) {
                    if (output_data.vars[vi].data.str[i] &&
                        !cmerge_arena_owns(str_arena, output_data.vars[vi].data.str[i])) {
                        free(output_data.vars[vi].data.str[i]);
                    }
                }
                /* Free the string pointer array */
                free(output_data.vars[vi].data.str);
            } else if (output_data.vars[vi].data.dbl) {
                /* Free numeric buffer */
                cmerge_aligned_free(output_data.vars[vi].data.dbl);
            }
        }
        cmerge_arena_free(str_arena);
        cmerge_aligned_free(output_data.vars);
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        return 920;
    }

    double t_permute = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 7: Store output data using parallel I/O
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    rc = ctools_data_store_selective(&output_data, output_var_stata_idx, n_output_vars, 1);
    if (rc != STATA_OK) {
        /* Cleanup on error - free only non-arena strings */
        for (size_t vi = 0; vi < n_output_vars; vi++) {
            if (output_data.vars[vi].type == STATA_TYPE_STRING && output_data.vars[vi].data.str) {
                /* Only free strings not owned by arena */
                for (size_t i = 0; i < output_nobs; i++) {
                    if (output_data.vars[vi].data.str[i] &&
                        !cmerge_arena_owns(str_arena, output_data.vars[vi].data.str[i])) {
                        free(output_data.vars[vi].data.str[i]);
                    }
                }
                free(output_data.vars[vi].data.str);
            } else if (output_data.vars[vi].data.dbl) {
                cmerge_aligned_free(output_data.vars[vi].data.dbl);
            }
        }
        cmerge_aligned_free(output_data.vars);
        cmerge_arena_free(str_arena);  /* Free arena in one operation */
        free(output_var_indices);
        free(output_var_stata_idx);
        free(output_var_is_key);
        free(output_var_key_idx);
        free(output_specs);
        free(master_orig_rows);
        if (sort_perm) free(sort_perm);
        stata_data_free(&master_data);
        free(all_var_indices);
        cache_clear_in_use();
        SF_error("cmerge: Failed to store output data\n");
        return rc;
    }

    double t_store = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Step 8: Write _merge variable
     * =================================================================== */

    t_start = cmerge_get_time_ms();

    if (merge_var_idx > 0) {
        ST_int mvar = (ST_int)merge_var_idx;
        for (size_t i = 0; i < output_nobs; i++) {
            SF_vstore(mvar, (ST_int)(i + 1), (double)output_specs[i].merge_result);
        }
    }

    /* ===================================================================
     * Step 9: Write keepusing variables from cache (sequential)
     * =================================================================== */

    if (n_keepusing > 0) {
        for (int kv = 0; kv < n_keepusing; kv++) {
            cmerge_keepusing_write_args_t ku_args;
            ku_args.keepusing_idx = kv;
            ku_args.dest_idx = (ST_int)keepusing_placeholder_indices[kv];
            ku_args.is_shared = shared_flags[kv];
            ku_args.update_mode = update_mode;
            ku_args.replace_mode = replace_mode;
            ku_args.specs = output_specs;
            ku_args.output_nobs = output_nobs;
            ku_args.success = 0;
            cmerge_write_keepusing_var_thread(&ku_args);
        }
    }

    double t_write_meta = cmerge_get_time_ms() - t_start;

    /* ===================================================================
     * Cleanup
     * =================================================================== */

    double t_cleanup_start = cmerge_get_time_ms();

    /* Free output data - use arena for O(1) string cleanup */
    for (size_t vi = 0; vi < n_output_vars; vi++) {
        if (output_data.vars[vi].type == STATA_TYPE_STRING && output_data.vars[vi].data.str) {
            /* Only free strings that fell back to strdup (not owned by arena) */
            for (size_t i = 0; i < output_nobs; i++) {
                if (output_data.vars[vi].data.str[i] &&
                    !cmerge_arena_owns(str_arena, output_data.vars[vi].data.str[i])) {
                    free(output_data.vars[vi].data.str[i]);
                }
            }
            free(output_data.vars[vi].data.str);
        } else if (output_data.vars[vi].data.dbl) {
            cmerge_aligned_free(output_data.vars[vi].data.dbl);
        }
    }
    cmerge_aligned_free(output_data.vars);
    cmerge_arena_free(str_arena);  /* O(1) cleanup of all arena strings */

    free(output_var_indices);
    free(output_var_stata_idx);
    free(output_var_is_key);
    free(output_var_key_idx);
    free(output_specs);
    free(master_orig_rows);
    if (sort_perm) free(sort_perm);
    stata_data_free(&master_data);
    free(all_var_indices);

    /* Free using cache (protected by lock) */
    cache_lock_acquire();
    stata_data_free(&g_using_cache.keys);
    stata_data_free(&g_using_cache.keepusing);
    g_using_cache.loaded = 0;
    g_cache_in_use = 0;
    cache_lock_release();

    double t_cleanup = cmerge_get_time_ms() - t_cleanup_start;

    /* Calculate total time */
    double t_total = cmerge_get_time_ms() - t_total_start;

    /* Save merge results to Stata scalars */
    SF_scal_save("_cmerge_N", (double)output_nobs);
    SF_scal_save("_cmerge_N_1", (double)n_master_only);
    SF_scal_save("_cmerge_N_2", (double)n_using_only);
    SF_scal_save("_cmerge_N_3", (double)n_matched);

    /* Save Phase 2 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_p2_load_master", t_load / 1000.0);
    SF_scal_save("_cmerge_p2_sort_master", t_sort / 1000.0);
    SF_scal_save("_cmerge_p2_merge_join", t_merge / 1000.0);
    SF_scal_save("_cmerge_p2_reorder", t_reorder / 1000.0);
    SF_scal_save("_cmerge_p2_permute", t_permute / 1000.0);
    SF_scal_save("_cmerge_p2_store", t_store / 1000.0);
    SF_scal_save("_cmerge_p2_write_meta", t_write_meta / 1000.0);
    SF_scal_save("_cmerge_p2_cleanup", t_cleanup / 1000.0);
    SF_scal_save("_cmerge_p2_total", t_total / 1000.0);
    SF_scal_save("_cmerge_n_output_vars", (double)n_output_vars);
    CTOOLS_SAVE_THREAD_INFO("_cmerge");

    return STATA_OK;
}

/* ============================================================================
 * Clear cached using data
 * ============================================================================ */

static ST_retcode cmerge_clear(void)
{
    cache_lock_acquire();

    /* Check if cache is currently in use by execute */
    if (g_cache_in_use) {
        cache_lock_release();
        SF_error("cmerge: Cannot clear cache while merge is in progress\n");
        return 459;
    }

    if (g_using_cache.loaded) {
        stata_data_free(&g_using_cache.keys);
        stata_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }

    cache_lock_release();
    return 0;
}

/* ============================================================================
 * Plugin entry point
 * ============================================================================ */

ST_retcode cmerge_main(const char *args)
{

    if (args == NULL || strlen(args) == 0) {
        SF_error("cmerge: no arguments specified\n");
        return 198;
    }

    /* Use static buffer to avoid potential stack issues */
    static char args_copy[8192];
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *phase = strtok(args_copy, " ");
    const char *rest = args + strlen(phase);
    while (*rest == ' ') rest++;

    if (strcmp(phase, "load_using") == 0) {
        return cmerge_load_using(rest);
    }
    else if (strcmp(phase, "execute") == 0) {
        return cmerge_execute(rest);
    }
    else if (strcmp(phase, "clear") == 0) {
        return cmerge_clear();
    }
    else {
        char msg[256];
        snprintf(msg, sizeof(msg), "cmerge: unknown phase '%s'\n", phase);
        SF_error(msg);
        return 198;
    }
}
