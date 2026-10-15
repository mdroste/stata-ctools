/*
    cmerge_impl.c
    Optimized C-Accelerated Merge for Stata

    Plugin calls, in order:
    - load_using: load keys and selected payload columns, sort, and cache them.
    - plan: in the untouched master, load the master keys alone, sort and join
      them. Join errors therefore leave the master unchanged.
    - execute: write the output. An identity output reuses the master already
      in Stata; otherwise whole master observations are staged row by row and
      written once in output order, and the staged rows are kept as an undo
      copy. Datasets whose strings are too wide to stage load, assemble and
      keep full columns instead.
    - rollback or commit: restore the master from the undo copy, or release it.
    - On the column path, output string arrays borrow contents from the loaded
      inputs, which remain alive until output storage and cleanup finish.

    Merge types: 1:1, m:1, 1:m, m:m, and 1:1 by observation number.

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

/* C11 atomics for platforms without GCC/Clang/MSVC intrinsics */
#if !defined(_WIN32) && !defined(__GNUC__) && !defined(__clang__)
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_ATOMICS__)
#define CMERGE_USE_C11_ATOMICS 1
#include <stdatomic.h>
#endif
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_runtime.h"
#include "ctools_threads.h"
#include "ctools_spi.h"
#include "ctools_config.h"
#include "ctools_rowstage.h"
#include "cmerge_impl.h"

/* Helper modules */
#include "ctools_sort_pairs.h"
#include "cmerge_keys.h"
#include "cmerge_group_search.h"
#include "cmerge_join.h"
#include "cmerge_io.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#define CMERGE_MAX_KEYVARS 32
#define CMERGE_MAX_VARS 1024

/* Timer: use ctools_timer_ms() directly from ctools_timer.h */

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
#elif defined(CMERGE_USE_C11_ATOMICS)
    /* C11 atomics fallback */
    int expected = 0;
    while (!atomic_compare_exchange_weak((atomic_int *)&g_cache_lock, &expected, 1)) {
        expected = 0;
        /* Spin */
    }
#else
    /* Last resort fallback: volatile spinlock. Not perfectly atomic but provides
     * basic protection. The race window is small (between read and write). */
    while (g_cache_lock) {
        /* Spin */
    }
    g_cache_lock = 1;
#endif
}

static void cache_lock_release(void)
{
#if defined(_WIN32)
    InterlockedExchange((volatile long *)&g_cache_lock, 0);
#elif defined(__GNUC__) || defined(__clang__)
    __sync_lock_release(&g_cache_lock);
#elif defined(CMERGE_USE_C11_ATOMICS)
    /* C11 atomics fallback */
    atomic_store((atomic_int *)&g_cache_lock, 0);
#else
    /* Last resort fallback */
    g_cache_lock = 0;
#endif
}

/* Helper to cleanup cache and clear in-use flag on error paths.
 * This frees the cached using data since the merge operation failed
 * and the cache is no longer useful. Prevents memory leaks in OOM conditions. */
/* The join assumes key order, so a "sorted" hint (the option, or a Stata sort
 * marker) is checked in one pass: 0 if some row sorts before its predecessor.
 * key_cols are 0-based positions in data->vars. Numeric order matches the
 * join (missing codes are ordered doubles above all numbers); strings compare
 * bytewise as in cmerge_compare_keys. */
static int cmerge_keys_in_order(const stata_data *data, const int *key_cols, int nkeys)
{
    for (size_t i = 1; i < data->nobs; i++) {
        int c = 0;
        for (int k = 0; k < nkeys && c == 0; k++) {
            const stata_variable *v = &data->vars[key_cols[k]];
            if (v->type == STATA_TYPE_DOUBLE) {
                double a = v->data.dbl[i - 1], b = v->data.dbl[i];
                c = (a < b) ? -1 : (a > b);
            } else {
                const char *a = v->data.str[i - 1] ? v->data.str[i - 1] : "";
                const char *b = v->data.str[i] ? v->data.str[i] : "";
                c = strcmp(a, b);
            }
        }
        if (c > 0) return 0;
    }
    return 1;
}

static void cache_cleanup_on_error(void)
{
    cache_lock_acquire();
    if (g_using_cache.loaded) {
        ctools_filtered_data_free(&g_using_cache.keys);
        ctools_filtered_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }
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
    if (strlen(args) >= sizeof(args_copy)) return 198;
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    int nkeys = 0;
    int key_indices[CMERGE_MAX_KEYVARS];
    int n_keepusing = 0;
    int keepusing_indices[CMERGE_MAX_VARS];
    int verbose = 0;
    int sorted = 0;
    int sorted_hint = 0;   /* Stata sort marker: trusted only if the check passes */
    int merge_by_n = 0;

    char *token = strtok(args_copy, " ");
    int arg_idx = 0;
    int in_keepusing_indices = 0;
    int keepusing_idx = 0;

    while (token != NULL) {
        if (arg_idx == 0) {
            if (!ctools_safe_atoi(token, &nkeys)) {
                SF_error("cmerge: invalid nkeys value\n");
                return 198;
            }
            if (nkeys > CMERGE_MAX_KEYVARS) {
                SF_error("cmerge: too many key variables (max 32)\n");
                return 198;
            }
        }
        else if (arg_idx <= nkeys) {
            if (!ctools_safe_atoi(token, &key_indices[arg_idx - 1])) {
                SF_error("cmerge: invalid key index value\n");
                return 198;
            }
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            if (token) {
                if (!ctools_safe_atoi(token, &n_keepusing)) {
                    SF_error("cmerge: invalid n_keepusing value\n");
                    return 198;
                }
                if (n_keepusing > CMERGE_MAX_VARS) {
                    SF_error("cmerge: too many keepusing variables (max 1024)\n");
                    return 198;
                }
            }
        }
        else if (strcmp(token, "keepusing_indices") == 0) {
            in_keepusing_indices = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_indices && keepusing_idx < n_keepusing && keepusing_idx < CMERGE_MAX_VARS) {
            if (!ctools_safe_atoi(token, &keepusing_indices[keepusing_idx])) {
                SF_error("cmerge: invalid keepusing index value\n");
                return 198;
            }
            keepusing_idx++;
        }
        else if (strcmp(token, "verbose") == 0) {
            verbose = 1;
        }
        else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        }
        else if (strcmp(token, "sorted_hint") == 0) {
            sorted_hint = 1;
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
        ctools_filtered_data_free(&g_using_cache.keys);
        ctools_filtered_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }

    cache_lock_release();

    double t_phase1_start = ctools_timer_ms();

    double t_load_keys = 0;
    double t_load_keepusing = 0;

    /* For merge_by_n, we don't need key variables - just keepusing */
    if (merge_by_n) {
        /* Initialize empty keys structure */
        ctools_filtered_data_init(&g_using_cache.keys);
        g_using_cache.keys.data.nobs = SF_nobs();  /* Get nobs from Stata */

        /* Load keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = ctools_timer_ms();
            rc = ctools_data_load(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0, CTOOLS_LOAD_SKIP_IF);
            if (rc != STATA_OK) {
                ctools_filtered_data_free(&g_using_cache.keys);
                ctools_filtered_data_free(&g_using_cache.keepusing);
                SF_error("cmerge: Failed to load keepusing variables\n");
                return ctools_stata_rc(rc);
            }
            if (g_using_cache.keepusing.data.vars == NULL) {
                ctools_filtered_data_free(&g_using_cache.keys);
                ctools_filtered_data_free(&g_using_cache.keepusing);
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = ctools_timer_ms() - t_load_keepusing_start;
            g_using_cache.keys.data.nobs = g_using_cache.keepusing.data.nobs;
        } else {
            ctools_filtered_data_init(&g_using_cache.keepusing);
        }
    } else {
        /* Standard merge: Load ONLY key variables */
        double t_load_keys_start = ctools_timer_ms();
        rc = ctools_data_load(&g_using_cache.keys, key_indices, nkeys, 0, 0, CTOOLS_LOAD_SKIP_IF);
        if (rc != STATA_OK) {
            ctools_filtered_data_free(&g_using_cache.keys);
            SF_error("cmerge: Failed to load using keys\n");
            return ctools_stata_rc(rc);
        }
        /* Critical null check - prevents crash if load failed silently */
        if (g_using_cache.keys.data.vars == NULL) {
            ctools_filtered_data_free(&g_using_cache.keys);
            SF_error("cmerge: FATAL - keys.vars is NULL after load\n");
            return 920;
        }
        t_load_keys = ctools_timer_ms() - t_load_keys_start;

        /* Load ONLY keepusing variables */
        if (n_keepusing > 0) {
            double t_load_keepusing_start = ctools_timer_ms();
            rc = ctools_data_load(&g_using_cache.keepusing, keepusing_indices, n_keepusing, 0, 0, CTOOLS_LOAD_SKIP_IF);
            if (rc != STATA_OK) {
                ctools_filtered_data_free(&g_using_cache.keys);
                ctools_filtered_data_free(&g_using_cache.keepusing);
                SF_error("cmerge: Failed to load keepusing variables\n");
                return ctools_stata_rc(rc);
            }
            /* Critical null check */
            if (g_using_cache.keepusing.data.vars == NULL) {
                ctools_filtered_data_free(&g_using_cache.keys);
                ctools_filtered_data_free(&g_using_cache.keepusing);
                SF_error("cmerge: FATAL - keepusing.vars is NULL after load\n");
                return 920;
            }
            t_load_keepusing = ctools_timer_ms() - t_load_keepusing_start;
        } else {
            ctools_filtered_data_init(&g_using_cache.keepusing);
        }
    }

    g_using_cache.nobs = g_using_cache.keys.data.nobs;
    g_using_cache.nkeys = nkeys;
    g_using_cache.n_keepusing = n_keepusing;
    g_using_cache.merge_by_n = merge_by_n;

    double t_sort = 0;
    double t_apply_perm = 0;

    /* merge errors when data promised sorted (the option) are not; a Stata
     * sort marker is trusted only after the same check, else the data are
     * sorted below */
    if ((sorted || sorted_hint) && !merge_by_n && nkeys > 0) {
        int key_cols[CMERGE_MAX_KEYVARS];
        int ok = (nkeys <= CMERGE_MAX_KEYVARS);
        for (int k = 0; ok && k < nkeys; k++) key_cols[k] = k;
        ok = ok && cmerge_keys_in_order(&g_using_cache.keys.data, key_cols, nkeys);
        if (!ok && sorted) {
            ctools_filtered_data_free(&g_using_cache.keys);
            ctools_filtered_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: using data not sorted\n");
            return 5;
        }
        sorted = ok;
    }

    /* Sort using data on keys (skip if sorted option specified or merge_by_n) */
    if (!sorted && !merge_by_n) {
        double t_sort_start = ctools_timer_ms();

        int *sort_vars = malloc(nkeys * sizeof(int));
        if (!sort_vars) {
            ctools_filtered_data_free(&g_using_cache.keys);
            ctools_filtered_data_free(&g_using_cache.keepusing);
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
            ctools_filtered_data_free(&g_using_cache.keys);
            ctools_filtered_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Memory allocation failed for permutation\n");
            return 920;
        }

        /* Use _with_perm version to get permutation BEFORE it's reset to identity */
        rc = ctools_sort_ips4o_with_perm(&g_using_cache.keys.data, sort_vars, nkeys, perm);
        if (rc != STATA_OK) {
            free(perm);
            free(sort_vars);
            ctools_filtered_data_free(&g_using_cache.keys);
            ctools_filtered_data_free(&g_using_cache.keepusing);
            SF_error("cmerge: Failed to sort using data\n");
            return ctools_stata_rc(rc);
        }

        t_sort = ctools_timer_ms() - t_sort_start;

        /* Validate permutation indices once (debug builds only) */
#ifndef NDEBUG
        for (size_t i = 0; i < nobs; i++) {
            if (perm[i] >= nobs) {
                free(perm);
                free(sort_vars);
                ctools_filtered_data_free(&g_using_cache.keys);
                ctools_filtered_data_free(&g_using_cache.keepusing);
                SF_error("cmerge: Invalid permutation index\n");
                return 459;
            }
        }
#endif

        /* Apply same permutation to keepusing (parallel across variables) */
        double t_apply_perm_start = ctools_timer_ms();
        if (n_keepusing > 0 && g_using_cache.keepusing.data.vars != NULL) {
            int perm_alloc_failed = 0;

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (int v = 0; v < n_keepusing; v++) {
                int local_failed;
                #ifdef _OPENMP
                #pragma omp atomic read
                #endif
                local_failed = perm_alloc_failed;
                if (local_failed) continue;

                stata_variable *var = &g_using_cache.keepusing.data.vars[v];
                if (var->type == STATA_TYPE_DOUBLE) {
                    double *new_data = ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(double));
                    if (!new_data) {
                        #ifdef _OPENMP
                        #pragma omp atomic write
                        #endif
                        perm_alloc_failed = 1;
                        continue;
                    }
                    for (size_t i = 0; i < nobs; i++) {
                        new_data[i] = var->data.dbl[perm[i]];
                    }
                    ctools_aligned_free(var->data.dbl);
                    var->data.dbl = new_data;
                } else {
                    char **new_data = ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, nobs, sizeof(char *));
                    if (!new_data) {
                        #ifdef _OPENMP
                        #pragma omp atomic write
                        #endif
                        perm_alloc_failed = 1;
                        continue;
                    }
                    for (size_t i = 0; i < nobs; i++) {
                        new_data[i] = var->data.str[perm[i]];
                    }
                    ctools_aligned_free(var->data.str);
                    var->data.str = new_data;
                }
            }

            if (perm_alloc_failed) {
                free(perm);
                free(sort_vars);
                ctools_filtered_data_free(&g_using_cache.keys);
                ctools_filtered_data_free(&g_using_cache.keepusing);
                SF_error("cmerge: Memory allocation failed in permutation\n");
                return 920;
            }
        }
        t_apply_perm = ctools_timer_ms() - t_apply_perm_start;

        free(perm);
        free(sort_vars);
    }

    g_using_cache.loaded = 1;

    double t_phase1_total = ctools_timer_ms() - t_phase1_start;

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
 * Phase 2: plan, execute, then commit or roll back
 *
 * plan:     runs in the untouched master. Loads the master keys alone, checks a
 *           promised sort order, sorts and joins the keys, and keeps the output
 *           plan, so key, sortedness and uniqueness errors leave Stata unchanged.
 * execute:  runs after the ado has added the new variables and observations.
 *           Before the first write to an original master cell it copies every
 *           original master variable: staged rows, or full columns when the
 *           strings are too wide to stage. The output is written from that
 *           copy, and a failed write restores the original observations from it
 *           before returning the error. An identity output reuses the master
 *           already in Stata. With nocopy (the ado's snapshot or fast modes) no
 *           copy is kept.
 * rollback: after a successful execute (Break, or a failure in the ado's
 *           remaining steps), rewrites the original observations from the copy.
 * commit:   releases the undo copy once the ado has finished.
 * ============================================================================ */

/* The join, computed by plan and consumed by execute. */
typedef struct {
    int ready;
    cmerge_type_t merge_type;
    int merge_by_n;
    int nkeys;
    int key_idx[CMERGE_MAX_KEYVARS]; /* one-based plugin indices of the master keys */
    size_t master_nobs;
    size_t orig_nvars;               /* master variables when planned (the plan call lists all) */
    size_t output_nobs;
    cmerge_output_spec_t *specs;     /* owned [output_nobs] */
    int32_t *orig_rows;              /* owned [output_nobs]: zero-based master row, -1 using-only */
    size_t n_master_only, n_using_only, n_matched;
    int identity;                    /* the output rows are the master rows, in order */
    int hold_undo;                   /* keep an undo copy even if no original cell is written */
    int nocopy;                      /* keep no undo copy (the ado's snapshot or fast mode) */
    double t_load, t_sort, t_join, t_reorder, t_total; /* milliseconds */
} cmerge_plan_t;

/* Every original master variable, observations 1..master_nobs in their original
 * order, held from before the first write to an original master cell until
 * the ado commits or rolls back: staged rows, or loaded columns when the
 * master cannot be staged. */
enum { CMERGE_UNDO_NONE = 0, CMERGE_UNDO_ROWS, CMERGE_UNDO_COLUMNS };
typedef struct {
    int held;                  /* CMERGE_UNDO_* */
    size_t master_nobs;
    ctools_rowstage rows;
    ctools_filtered_data cols;
    int *col_vars;             /* owned [n_cols]: one-based plugin index of each column */
    size_t n_cols;
} cmerge_undo_t;

static cmerge_plan_t g_plan;
static cmerge_undo_t g_undo;

static void cmerge_plan_free(void)
{
    free(g_plan.specs);
    free(g_plan.orig_rows);
    memset(&g_plan, 0, sizeof(g_plan));
}

static void cmerge_undo_free(void)
{
    ctools_rowstage_free(&g_undo.rows);
    ctools_filtered_data_free(&g_undo.cols);
    free(g_undo.col_vars);
    g_undo.col_vars = NULL;
    g_undo.n_cols = 0;
    g_undo.held = CMERGE_UNDO_NONE;
    g_undo.master_nobs = 0;
}

/* Write the original observations back from the undo copy. The in-order row
 * store allocates nothing. */
static stata_retcode cmerge_undo_restore(void)
{
    if ((size_t)SF_nobs() < g_undo.master_nobs) return STATA_ERR_INVALID_INPUT;
    if (g_undo.held == CMERGE_UNDO_ROWS)
        return ctools_rowstage_store(&g_undo.rows, NULL, g_undo.master_nobs, 1);
    if (g_undo.held == CMERGE_UNDO_COLUMNS)
        return ctools_data_store_selective(&g_undo.cols.data, g_undo.col_vars, g_undo.n_cols, 1);
    return STATA_OK;
}

/* Everything one merge holds between plugin calls. */
static void cmerge_release_all(void)
{
    cmerge_undo_free();
    cmerge_plan_free();
    cache_cleanup_on_error();
}

/* Arguments: "<merge_type> <nkeys> <key_idx>... master_nobs <N>
 * preserve_order <0|1> [sorted] [sorted_hint] [hold_undo] [nocopy]
 * [merge_by_n] [verbose]".
 * The plugin varlist lists every master variable in its current order; key
 * indices are one-based within it, and execute uses the same indices.
 * hold_undo: the ado will drop rows (keep()), so an undo copy is needed even
 * when no original cell is written. nocopy: keep no undo copy at all. */
static ST_retcode cmerge_plan(const char *args)
{
    double t_total_start = ctools_timer_ms();
    cmerge_plan_free();
    cmerge_undo_free();

    char args_copy[4096];
    if (strlen(args) >= sizeof(args_copy)) {
        cache_cleanup_on_error();
        return 198;
    }
    memcpy(args_copy, args, strlen(args) + 1);

    cmerge_plan_t plan;
    memset(&plan, 0, sizeof(plan));
    int preserve_order = 0, sorted = 0, sorted_hint = 0, value = 0, bad = 0, arg_idx = 0;
    for (char *token = strtok(args_copy, " "); token && !bad;
         token = strtok(NULL, " "), arg_idx++) {
        if (arg_idx == 0) {
            bad = !ctools_safe_atoi(token, &value) || value < MERGE_1_1 || value > MERGE_M_M;
            plan.merge_type = (cmerge_type_t)value;
        } else if (arg_idx == 1) {
            bad = !ctools_safe_atoi(token, &plan.nkeys) || plan.nkeys < 0 ||
                  plan.nkeys > CMERGE_MAX_KEYVARS;
        } else if (arg_idx < 2 + plan.nkeys) {
            int *key = &plan.key_idx[arg_idx - 2];
            bad = !ctools_safe_atoi(token, key) || *key < 1 || *key > SF_nvars();
        } else if (strcmp(token, "master_nobs") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atozu(token, &plan.master_nobs);
        } else if (strcmp(token, "preserve_order") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &preserve_order);
        } else if (strcmp(token, "sorted") == 0) {
            sorted = 1;
        } else if (strcmp(token, "sorted_hint") == 0) {
            sorted_hint = 1;
        } else if (strcmp(token, "hold_undo") == 0) {
            plan.hold_undo = 1;
        } else if (strcmp(token, "nocopy") == 0) {
            plan.nocopy = 1;
        }
        /* merge_by_n is taken from the using cache; verbose needs nothing here */
    }

    cache_lock_acquire();
    int loaded = g_using_cache.loaded && !g_cache_in_use &&
                 (g_using_cache.merge_by_n || g_using_cache.keys.data.vars != NULL);
    cache_lock_release();
    if (!loaded) {
        cache_cleanup_on_error();
        SF_error("cmerge: Using data not loaded. Call load_using first.\n");
        return 459;
    }
    plan.merge_by_n = g_using_cache.merge_by_n;
    if (bad || plan.master_nobs == 0 || plan.master_nobs > (size_t)SF_nobs() ||
        (!plan.merge_by_n && (plan.nkeys < 1 || plan.nkeys != g_using_cache.nkeys))) {
        cache_cleanup_on_error();
        SF_error("cmerge: invalid plan arguments\n");
        return 198;
    }

    const size_t master_nobs = plan.master_nobs;
    ST_retcode rc = 0;
    ctools_filtered_data keys;
    ctools_filtered_data_init(&keys);
    size_t *sort_perm = NULL;
    cmerge_output_spec_t *specs = NULL;
    int32_t *orig_rows = NULL;
    size_t output_nobs = 0;
    double t_start = ctools_timer_ms();

    if (plan.merge_by_n) {
        /* Observation i of the master meets observation i of the using data. */
        size_t using_nobs = g_using_cache.nobs;
        output_nobs = master_nobs > using_nobs ? master_nobs : using_nobs;
        if (output_nobs > INT32_MAX) {
            SF_error("cmerge: Merge result exceeds Stata observation limit\n");
            rc = 920;
            goto done;
        }
        specs = ctools_safe_malloc2(output_nobs, sizeof(cmerge_output_spec_t));
        if (!specs) {
            rc = 920;
            goto done;
        }
        for (size_t i = 0; i < output_nobs; i++) {
            int in_master = i < master_nobs, in_using = i < using_nobs;
            specs[i].master_sorted_row = in_master ? (int32_t)i : -1;
            specs[i].using_sorted_row = in_using ? (int32_t)i : -1;
            specs[i].merge_result = in_master && in_using ? MERGE_RESULT_BOTH
                                    : in_master          ? MERGE_RESULT_MASTER_ONLY
                                                         : MERGE_RESULT_USING_ONLY;
        }
        plan.t_join = ctools_timer_ms() - t_start;
    } else {
        stata_retcode load_rc = ctools_data_load(&keys, plan.key_idx, (size_t)plan.nkeys, 1,
                                                 master_nobs, CTOOLS_LOAD_SKIP_IF);
        if (load_rc != STATA_OK || keys.data.vars == NULL) {
            SF_error("cmerge: Failed to load master data\n");
            rc = load_rc != STATA_OK ? ctools_stata_rc(load_rc) : 920;
            goto done;
        }
        plan.t_load = ctools_timer_ms() - t_start;

        /* merge errors when data promised sorted (the option) are not; a Stata
         * sort marker is trusted only after the same check, else the data are
         * sorted below */
        int key_cols[CMERGE_MAX_KEYVARS];
        for (int k = 0; k < plan.nkeys; k++) key_cols[k] = k;
        if (sorted || sorted_hint) {
            int ok = cmerge_keys_in_order(&keys.data, key_cols, plan.nkeys);
            if (!ok && sorted) {
                SF_error("cmerge: master data not sorted\n");
                rc = 5;
                goto done;
            }
            sorted = ok;
        }

        if (!sorted) {
            t_start = ctools_timer_ms();
            /* Every algorithm here is stable, so the dispatcher's pick for the
             * key types does not change the order. One-based key positions. */
            int sort_vars[CMERGE_MAX_KEYVARS];
            for (int k = 0; k < plan.nkeys; k++) sort_vars[k] = k + 1;
            stata_retcode sort_rc =
                ctools_sort_dispatch(&keys.data, sort_vars, (size_t)plan.nkeys, SORT_ALG_AUTO);
            if (sort_rc != STATA_OK) {
                rc = ctools_stata_rc(sort_rc);
                goto done;
            }
            sort_perm = ctools_safe_malloc2(master_nobs, sizeof(size_t));
            if (!sort_perm) {
                rc = 920;
                goto done;
            }
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (size_t i = 0; i < master_nobs; i++) {
                sort_perm[i] = (size_t)keys.data.sort_order[i];
            }

            /* Put the keys in sorted order for the join. */
            for (int k = 0; k < plan.nkeys; k++) {
                stata_variable *var = &keys.data.vars[k];
                if (var->type == STATA_TYPE_DOUBLE) {
                    double *sorted_data =
                        ctools_safe_aligned_alloc2(64, master_nobs, sizeof(double));
                    if (!sorted_data) {
                        rc = 920;
                        goto done;
                    }
                    const double *old_data = var->data.dbl;
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static) if(master_nobs >= MIN_OBS_PER_THREAD * 2)
                    #endif
                    for (size_t i = 0; i < master_nobs; i++) {
                        sorted_data[i] = old_data[sort_perm[i]];
                    }
                    ctools_aligned_free(var->data.dbl);
                    var->data.dbl = sorted_data;
                } else {
                    char **sorted_data =
                        ctools_safe_aligned_alloc2(CACHE_LINE_SIZE, master_nobs, sizeof(char *));
                    if (!sorted_data) {
                        rc = 920;
                        goto done;
                    }
                    char *const *old_data = var->data.str;
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(static) if(master_nobs >= MIN_OBS_PER_THREAD * 2)
                    #endif
                    for (size_t i = 0; i < master_nobs; i++) {
                        sorted_data[i] = old_data[sort_perm[i]];
                    }
                    ctools_aligned_free(var->data.str);
                    var->data.str = sorted_data;
                }
            }
            for (size_t i = 0; i < master_nobs; i++) {
                keys.data.sort_order[i] = (perm_idx_t)i;
            }
            plan.t_sort = ctools_timer_ms() - t_start;
        }

        t_start = ctools_timer_ms();
        int64_t joined = cmerge_sorted_join(&keys.data, &g_using_cache.keys.data, plan.nkeys,
                                            plan.merge_type, &specs);
        if (joined < 0) {
            specs = NULL;
            if (joined == -4) {
                SF_error("cmerge: key variable types must match\n");
                rc = 106;
            } else if (joined == -2 || joined == -3) {
                SF_error(joined == -2 ?
                    "cmerge: keys do not uniquely identify master observations\n" :
                    "cmerge: keys do not uniquely identify using observations\n");
                rc = 459;
            } else {
                SF_error("cmerge: Merge join failed\n");
                rc = 920;
            }
            goto done;
        }
        /* Stata max obs is 2^31-1; reject if exceeded */
        if (joined > INT32_MAX) {
            SF_error("cmerge: Merge result exceeds Stata observation limit\n");
            rc = 920;
            goto done;
        }
        output_nobs = (size_t)joined;
        plan.t_join = ctools_timer_ms() - t_start;
    }

    /* Original (zero-based) master row of each output row. */
    orig_rows = ctools_safe_malloc2(output_nobs ? output_nobs : 1, sizeof(int32_t));
    if (!orig_rows) {
        rc = 920;
        goto done;
    }
    for (size_t i = 0; i < output_nobs; i++) {
        int32_t sorted_row = specs[i].master_sorted_row;
        orig_rows[i] = sorted_row < 0 ? -1
                       : sort_perm   ? (int32_t)sort_perm[sorted_row]
                                     : sorted_row;
    }

    /* preserve_order: master rows in their original order, using-only rows after. */
    if (preserve_order && output_nobs > 0) {
        t_start = ctools_timer_ms();
        ctools_order_pair_t *pairs = ctools_safe_malloc2(output_nobs, sizeof(ctools_order_pair_t));
        cmerge_output_spec_t *new_specs =
            ctools_safe_malloc2(output_nobs, sizeof(cmerge_output_spec_t));
        int32_t *new_orig_rows = ctools_safe_malloc2(output_nobs, sizeof(int32_t));
        if (!pairs || !new_specs || !new_orig_rows) {
            free(pairs);
            free(new_specs);
            free(new_orig_rows);
            rc = 920;
            goto done;
        }
        for (size_t i = 0; i < output_nobs; i++) {
            pairs[i].orig_idx = i;
            pairs[i].order_key = orig_rows[i] >= 0
                ? (size_t)orig_rows[i]
                : master_nobs + (size_t)specs[i].using_sorted_row;
        }
        ctools_sort_order_pairs(pairs, output_nobs);
        for (size_t i = 0; i < output_nobs; i++) {
            new_specs[i] = specs[pairs[i].orig_idx];
            new_orig_rows[i] = orig_rows[pairs[i].orig_idx];
        }
        free(pairs);
        free(specs);
        free(orig_rows);
        specs = new_specs;
        orig_rows = new_orig_rows;
        plan.t_reorder = ctools_timer_ms() - t_start;
    }

    /* Count the results; identity: every output row is master row i, in order. */
    plan.identity = 1;
    for (size_t i = 0; i < output_nobs; i++) {
        int8_t m = specs[i].merge_result;
        plan.n_master_only += (m == MERGE_RESULT_MASTER_ONLY);
        plan.n_using_only += (m == MERGE_RESULT_USING_ONLY);
        plan.n_matched += (m == MERGE_RESULT_BOTH);
        if (specs[i].master_sorted_row != (int32_t)i || orig_rows[i] != (int32_t)i)
            plan.identity = 0;
    }


done:
    free(sort_perm);
    ctools_filtered_data_free(&keys);
    if (rc) {
        free(specs);
        free(orig_rows);
        cache_cleanup_on_error();
        return rc;
    }
    plan.output_nobs = output_nobs;
    plan.specs = specs;
    plan.orig_rows = orig_rows;
    plan.orig_nvars = (size_t)SF_nvars();
    plan.t_total = ctools_timer_ms() - t_total_start;
    plan.ready = 1;
    g_plan = plan;
    SF_scal_save("_cmerge_plan_N", (double)output_nobs);
    return 0;
}

/* Load every original master variable (one-based plugin indices) as full
 * columns in their original order, for data that cannot be staged. */
static stata_retcode cmerge_load_columns(ctools_filtered_data *master, int *vars, size_t nvars,
                                         size_t master_nobs)
{
    ctools_filtered_data_init(master);
    stata_retcode rc = ctools_data_load(master, vars, nvars, 1, master_nobs,
                                        CTOOLS_LOAD_SKIP_IF);
    if (rc == STATA_OK && master->data.vars == NULL) rc = STATA_ERR_MEMORY;
    if (rc != STATA_OK) ctools_filtered_data_free(master);
    return rc;
}

/* Column path for master data that cannot be staged: gather the output rows
 * from the loaded columns of every original master variable (vars, one-based
 * plugin indices, in the order of master's columns) and store them. */
static stata_retcode cmerge_store_columns(const ctools_filtered_data *master, int *vars,
                                          size_t nvars, double *t_permute)
{
    const size_t output_nobs = g_plan.output_nobs;
    const int32_t *orig_rows = g_plan.orig_rows;
    const cmerge_output_spec_t *specs = g_plan.specs;
    const size_t using_nobs = g_using_cache.nobs;
    stata_retcode rc;

    double t_start = ctools_timer_ms();
    stata_data output;
    stata_data_init(&output);
    size_t vars_bytes;
    if (ctools_safe_mul_size(nvars, sizeof(stata_variable), &vars_bytes) == 0) {
        output.vars = (stata_variable *)ctools_cacheline_alloc(vars_bytes);
    }
    if (!output.vars) {
        return STATA_ERR_MEMORY;
    }
    memset(output.vars, 0, vars_bytes);
    output.nobs = output_nobs;
    output.nvars = nvars;

    /* Output string arrays borrow contents from the loaded master columns and
     * the using keys, which both outlive the store. */
    int alloc_failed = 0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (size_t j = 0; j < nvars; j++) {
        int failed;
        #ifdef _OPENMP
        #pragma omp atomic read
        #endif
        failed = alloc_failed;
        if (failed) continue;

        const stata_variable *src = &master->data.vars[j];
        stata_variable *dst = &output.vars[j];
        const stata_variable *using_key = NULL;
        for (int k = 0; k < g_plan.nkeys; k++) {
            if (g_plan.key_idx[k] == vars[j]) using_key = &g_using_cache.keys.data.vars[k];
        }
        dst->nobs = output_nobs;
        dst->type = src->type;
        if (src->type == STATA_TYPE_DOUBLE) {
            double *values = (double *)ctools_safe_cacheline_alloc2(output_nobs ? output_nobs : 1,
                                                                    sizeof(double));
            if (!values) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }
            for (size_t i = 0; i < output_nobs; i++) {
                int32_t row = orig_rows[i], using_row = specs[i].using_sorted_row;
                values[i] = row >= 0 ? src->data.dbl[row]
                            : (using_key && using_row >= 0 && (size_t)using_row < using_nobs)
                                ? using_key->data.dbl[using_row]
                                : SV_missval;
            }
            dst->data.dbl = values;
        } else {
            dst->buffer_storage = CTOOLS_BUFFER_HEAP;
            dst->string_storage = CTOOLS_STRINGS_BORROWED;
            char **values = (char **)ctools_safe_calloc2(output_nobs ? output_nobs : 1,
                                                         sizeof(char *));
            if (!values) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                alloc_failed = 1;
                continue;
            }
            for (size_t i = 0; i < output_nobs; i++) {
                int32_t row = orig_rows[i], using_row = specs[i].using_sorted_row;
                char *value = row >= 0 ? src->data.str[row]
                              : (using_key && using_row >= 0 && (size_t)using_row < using_nobs)
                                  ? using_key->data.str[using_row]
                                  : NULL;
                values[i] = value ? value : "";
            }
            dst->data.str = values;
        }
    }
    *t_permute = ctools_timer_ms() - t_start;

    rc = alloc_failed ? STATA_ERR_MEMORY
                      : ctools_data_store_selective(&output, vars, nvars, 1);
    /* Output buffers only borrow strings; the loaded columns stay with the caller. */
    stata_data_free(&output);
    return rc;
}

/* An execute that fails after writing restores the original observations from
 * the undo copy before reporting the error: the runtime releases this
 * command's state after any failed call, so a later rollback would find no
 * copy. The ado then undoes its own changes. */
static ST_retcode cmerge_execute_failed(stata_retcode rc)
{
    if (g_undo.held && cmerge_undo_restore() != STATA_OK) {
        SF_error("cmerge: could not restore the master observations\n");
    }
    cmerge_release_all();
    return ctools_stata_rc(rc);
}

/* Arguments: "<merge_type> <nkeys> <key_idx>... master_nobs <N> master_nvars <K>
 * n_keepusing <n> [keepusing_placeholders <idx>...] [merge_var_idx <idx>]
 * preserve_order <0|1> [update] [replace] [shared_flags <0|1>...] [fault_store]
 * [sorted|sorted_hint] [merge_by_n] [verbose]".
 * The varlist lists every current variable; the original master variables come
 * first, in the order plan saw them. fault_store is a validation hook: the
 * store fails after writing half of the output rows. */
static ST_retcode cmerge_execute(const char *args)
{
    double t_start, t_total_start;

    t_total_start = ctools_timer_ms();

    /* Acquire lock and mark cache as in-use to prevent concurrent clear */
    cache_lock_acquire();

    if (!g_using_cache.loaded || !g_plan.ready) {
        cache_lock_release();
        cmerge_release_all();
        SF_error("cmerge: Using data not loaded or join not planned.\n");
        return 459;
    }

    /* Mark cache as in-use before releasing lock */
    g_cache_in_use = 1;
    cache_lock_release();

    /* Parse arguments - use static buffers to avoid stack overflow.
     * Thread safety: Stata plugins are single-threaded, so these statics
     * cannot race. They're reinitialized on each call. */
    static char args_copy[8192];
    if (strlen(args) >= sizeof(args_copy)) {
        cmerge_release_all();
        return 198;
    }
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    cmerge_type_t merge_type = MERGE_1_1;
    int nkeys = 0;
    static int master_key_indices[CMERGE_MAX_KEYVARS];
    size_t master_nobs = 0;
    size_t master_nvars = 0;
    int n_keepusing = 0;
    static int keepusing_placeholder_indices[CMERGE_MAX_VARS];
    int merge_var_idx = 0;
    int update_mode = 0;
    int replace_mode = 0;
    int fault_store = 0;
    static int shared_flags[CMERGE_MAX_VARS];

    int arg_idx = 0;
    int in_keepusing_placeholders = 0;
    int keepusing_idx = 0;
    int bad = 0;

    char *token = strtok(args_copy, " ");
    int parse_tmp;  /* Temporary for safe parsing */
    while (token != NULL && !bad) {
        if (arg_idx == 0) {
            bad = !ctools_safe_atoi(token, &parse_tmp);
            merge_type = (cmerge_type_t)parse_tmp;
        }
        else if (arg_idx == 1) {
            bad = !ctools_safe_atoi(token, &nkeys) || nkeys < 0 || nkeys > CMERGE_MAX_KEYVARS;
        }
        else if (arg_idx >= 2 && arg_idx < 2 + nkeys) {
            bad = !ctools_safe_atoi(token, &master_key_indices[arg_idx - 2]);
        }
        else if (strcmp(token, "master_nobs") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atozu(token, &master_nobs);
        }
        else if (strcmp(token, "master_nvars") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atozu(token, &master_nvars) ||
                  master_nvars > (size_t)SF_nvars();
        }
        else if (strcmp(token, "n_keepusing") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &n_keepusing) || n_keepusing < 0 ||
                  n_keepusing > CMERGE_MAX_VARS;
        }
        else if (strcmp(token, "keepusing_placeholders") == 0) {
            in_keepusing_placeholders = 1;
            keepusing_idx = 0;
        }
        else if (in_keepusing_placeholders && keepusing_idx < n_keepusing && keepusing_idx < CMERGE_MAX_VARS) {
            bad = !ctools_safe_atoi(token, &keepusing_placeholder_indices[keepusing_idx]);
            keepusing_idx++;
        }
        else if (strcmp(token, "merge_var_idx") == 0) {
            token = strtok(NULL, " ");
            bad = !token || !ctools_safe_atoi(token, &merge_var_idx);
            in_keepusing_placeholders = 0;
        }
        else if (strcmp(token, "preserve_order") == 0) {
            /* Applied by plan */
            token = strtok(NULL, " ");
            bad = !token;
        }
        else if (strcmp(token, "update") == 0) {
            update_mode = 1;
        }
        else if (strcmp(token, "replace") == 0) {
            replace_mode = 1;
        }
        else if (strcmp(token, "fault_store") == 0) {
            fault_store = 1;
        }
        else if (strcmp(token, "shared_flags") == 0) {
            for (int i = 0; i < n_keepusing && i < CMERGE_MAX_VARS && !bad; i++) {
                token = strtok(NULL, " ");
                if (token) {
                    bad = !ctools_safe_atoi(token, &shared_flags[i]);
                } else {
                    shared_flags[i] = 0;
                }
            }
        }
        /* sorted, sorted_hint, merge_by_n and verbose were applied by plan */
        arg_idx++;
        token = strtok(NULL, " ");
    }

    /* The plan must describe this master. */
    int matches = !bad && merge_type == g_plan.merge_type && nkeys == g_plan.nkeys &&
                  master_nobs == g_plan.master_nobs && master_nobs <= (size_t)SF_nobs() &&
                  g_plan.output_nobs <= (size_t)SF_nobs();
    for (int k = 0; matches && k < nkeys; k++) {
        matches = master_key_indices[k] == g_plan.key_idx[k];
    }
    /* _merge must be one of the variables appended after the master's own,
     * or the merge codes would overwrite a master variable. */
    if (merge_var_idx != 0 &&
        (merge_var_idx < 0 || (size_t)merge_var_idx <= g_plan.orig_nvars ||
         merge_var_idx > SF_nvars() || SF_var_is_string(merge_var_idx))) {
        matches = 0;
    }
    if (!matches) {
        cmerge_release_all();
        SF_error("cmerge: invalid execute arguments\n");
        return 198;
    }

    const size_t output_nobs = g_plan.output_nobs;
    cmerge_output_spec_t *output_specs = g_plan.specs;
    const int32_t *master_orig_rows = g_plan.orig_rows;
    const size_t n_using_only = g_plan.n_using_only;

    /* Original master variables: everything except the new keepusing
     * placeholders and _merge, which are all missing here and are filled by
     * the writers below. One-based plugin variable indices. */
    int *stage_vars = NULL;
    size_t n_stage_vars = 0;
    if (master_nvars > 0) {
        stage_vars = ctools_safe_malloc2(master_nvars, sizeof(int));
        if (!stage_vars) {
            cmerge_release_all();
            SF_error("cmerge: Failed to allocate var indices\n");
            return 920;
        }
        for (size_t v = 0; v < master_nvars; v++) {
            int stata_idx = (int)(v + 1);
            int placeholder = (stata_idx == merge_var_idx);
            for (int kv = 0; kv < n_keepusing && !placeholder; kv++) {
                placeholder = (keepusing_placeholder_indices[kv] == stata_idx && !shared_flags[kv]);
            }
            for (int k = 0; k < nkeys; k++) {
                if (master_key_indices[k] == stata_idx) placeholder = 0;
            }
            if (!placeholder) stage_vars[n_stage_vars++] = stata_idx;
        }
    }
    int stage_ok = n_stage_vars > 0 &&
                   ctools_rowstage_supported(stage_vars, n_stage_vars, NULL, master_nobs);

    /* Only a true identity output can reuse the master already in Stata:
     * using-only keys can insert rows anywhere in the output. */
    int skip_master_store = g_plan.identity && n_using_only == 0 && output_nobs == master_nobs;
    /* update and replace rewrite matched master values of shared variables. */
    int shared_rewritten = 0;
    for (int kv = 0; kv < n_keepusing && (update_mode || replace_mode); kv++) {
        shared_rewritten |= shared_flags[kv] != 0;
    }
    int writes_master = n_stage_vars > 0 && (!skip_master_store || shared_rewritten);
    int need_copy = !g_plan.nocopy && n_stage_vars > 0 && (writes_master || g_plan.hold_undo);
    int write_rows = n_stage_vars > 0 && !skip_master_store;

    double t_load = g_plan.t_load;
    double t_permute = 0.0;
    double t_store = 0.0;
    /* The master data the output is written from: whole master observations in
     * their original order, or full columns when they cannot be staged. Kept
     * as the undo copy unless nocopy. Stata has not been modified by this call
     * yet, so a failed load leaves it intact. */
    ctools_rowstage rows;
    ctools_filtered_data cols;
    ctools_rowstage_init(&rows);
    ctools_filtered_data_init(&cols);
    int have = CMERGE_UNDO_NONE;
    if (need_copy || write_rows) {
        t_start = ctools_timer_ms();
        stata_retcode load_rc = stage_ok
            ? ctools_rowstage_load(&rows, stage_vars, n_stage_vars, NULL, 1, master_nobs)
            : cmerge_load_columns(&cols, stage_vars, n_stage_vars, master_nobs);
        t_load += ctools_timer_ms() - t_start;
        if (load_rc != STATA_OK) {
            free(stage_vars);
            cmerge_release_all();
            SF_error("cmerge: Failed to load master data\n");
            return ctools_stata_rc(load_rc);
        }
        have = stage_ok ? CMERGE_UNDO_ROWS : CMERGE_UNDO_COLUMNS;
    }
    if (need_copy) {
        g_undo.held = have;
        g_undo.master_nobs = master_nobs;
        g_undo.rows = rows;
        g_undo.cols = cols;
        g_undo.col_vars = have == CMERGE_UNDO_COLUMNS ? stage_vars : NULL;
        g_undo.n_cols = have == CMERGE_UNDO_COLUMNS ? n_stage_vars : 0;
        if (g_undo.col_vars) stage_vars = NULL;   /* now owned by the undo copy */
        ctools_rowstage_init(&rows);
        ctools_filtered_data_init(&cols);
    }
    const ctools_rowstage *src_rows = need_copy ? &g_undo.rows : &rows;
    const ctools_filtered_data *src_cols = need_copy ? &g_undo.cols : &cols;
    int *col_vars = need_copy && g_undo.col_vars ? g_undo.col_vars : stage_vars;

    if (write_rows) {
        stata_retcode store_rc;
        t_start = ctools_timer_ms();
        if (have == CMERGE_UNDO_ROWS) {
            /* Write each output row once from the staged rows. master_orig_rows
             * holds zero-based master rows, with -1 (CTOOLS_ROW_NONE as
             * perm_idx_t) for using-only rows. */
            size_t n_store = fault_store ? output_nobs / 2 : output_nobs;
            store_rc = ctools_rowstage_store(src_rows, (const perm_idx_t *)master_orig_rows,
                                             n_store, 1);

            /* Using-only rows take their key values from the using data. */
            for (size_t i = 0; i < n_store && store_rc == STATA_OK && n_using_only > 0; i++) {
                if (master_orig_rows[i] >= 0) continue;
                int32_t using_row = output_specs[i].using_sorted_row;
                if (using_row < 0 || (size_t)using_row >= g_using_cache.nobs) continue;
                for (int k = 0; k < nkeys && store_rc == STATA_OK; k++) {
                    const stata_variable *key = &g_using_cache.keys.data.vars[k];
                    ST_int key_idx = (ST_int)master_key_indices[k];
                    int io_rc;
                    if (key->type == STATA_TYPE_DOUBLE) {
                        io_rc = (_stata_)->safestore(key_idx, (ST_int)(i + 1),
                                                     key->data.dbl[using_row]);
                    } else {
                        char *value = key->data.str[using_row];
                        io_rc = SF_sstore(key_idx, (ST_int)(i + 1), value ? value : "");
                    }
                    store_rc = ctools_spi_status(io_rc, STATA_ERR_STATA_WRITE);
                }
            }
        } else {
            store_rc = cmerge_store_columns(src_cols, col_vars, n_stage_vars, &t_permute);
        }
        if (store_rc == STATA_OK && fault_store) store_rc = STATA_ERR_MEMORY;
        t_store = ctools_timer_ms() - t_start;

        if (store_rc != STATA_OK) {
            free(stage_vars);
            ctools_rowstage_free(&rows);
            ctools_filtered_data_free(&cols);
            SF_error("cmerge: Failed to store output data\n");
            return cmerge_execute_failed(store_rc);
        }
    }
    /* Without an undo copy, the source of the output goes now. */
    ctools_rowstage_free(&rows);
    ctools_filtered_data_free(&cols);
    free(stage_vars);

    t_start = ctools_timer_ms();

    /* update: per-row "missing updated" / "nonmissing conflict" flags, filled
     * by the keepusing writers and folded into _merge codes 4 and 5 below. */
    uint8_t *row_flags = NULL;
    if (update_mode && n_keepusing > 0 && output_nobs > 0) {
        row_flags = calloc(output_nobs, 1);
        if (!row_flags) {
            SF_error("cmerge: memory allocation failed\n");
            return cmerge_execute_failed(STATA_ERR_MEMORY);
        }
    }


    /* ===================================================================
     * Step 8: Write keepusing variables from cache (parallel)
     * =================================================================== */

    if (n_keepusing > 0) {
        /* One task per variable and row range: every task owns distinct cells,
         * and a few variables still occupy every worker. */
        size_t ku_chunks = 1;
        if (output_nobs >= (size_t)MIN_OBS_PER_THREAD * 2) {
            int workers = ctools_get_max_threads();
            ku_chunks = output_nobs / MIN_OBS_PER_THREAD;
            if (workers > 0 && ku_chunks > (size_t)workers) ku_chunks = (size_t)workers;
        }
        size_t ku_tasks = (size_t)n_keepusing * ku_chunks;
        cmerge_keepusing_write_args_t *ku_args_array = malloc(ku_tasks * sizeof(cmerge_keepusing_write_args_t));
        if (!ku_args_array) {
            /* Fallback to sequential on allocation failure */
            for (int kv = 0; kv < n_keepusing; kv++) {
                cmerge_keepusing_write_args_t ku_args;
                ku_args.keepusing_idx = kv;
                ku_args.dest_idx = (ST_int)keepusing_placeholder_indices[kv];
                ku_args.is_shared = shared_flags[kv];
                ku_args.update_mode = update_mode;
                ku_args.replace_mode = replace_mode;
                ku_args.specs = output_specs;
                ku_args.row_start = 0;
                ku_args.output_nobs = output_nobs;
                ku_args.row_flags = row_flags;
                ku_args.success = 0;
                cmerge_write_keepusing_var_thread(&ku_args);
            }
        } else {
            /* Prepare all args */
            for (size_t task = 0; task < ku_tasks; task++) {
                int kv = (int)(task / ku_chunks);
                size_t chunk = task % ku_chunks;
                size_t share = output_nobs / ku_chunks, extra = output_nobs % ku_chunks;
                ku_args_array[task].keepusing_idx = kv;
                ku_args_array[task].dest_idx = (ST_int)keepusing_placeholder_indices[kv];
                ku_args_array[task].is_shared = shared_flags[kv];
                ku_args_array[task].update_mode = update_mode;
                ku_args_array[task].replace_mode = replace_mode;
                ku_args_array[task].specs = output_specs;
                ku_args_array[task].row_start = chunk * share + (chunk < extra ? chunk : extra);
                ku_args_array[task].output_nobs = ku_args_array[task].row_start + share + (chunk < extra);
                ku_args_array[task].row_flags = row_flags;
                ku_args_array[task].success = 0;
            }

            /* Submit all keepusing writes to thread pool */
            ctools_persistent_pool *pool = ctools_get_global_pool();
            int pool_ok = 0;
            if (pool != NULL && ku_tasks >= 1) {
                if (ctools_persistent_pool_submit_batch(pool, cmerge_write_keepusing_var_thread,
                                                         ku_args_array, ku_tasks,
                                                         sizeof(cmerge_keepusing_write_args_t)) == 0) {
                    ctools_persistent_pool_wait(pool);
                    pool_ok = 1;
                }
            }
            if (!pool_ok) {
                /* Sequential fallback if no pool or submission failed */
                for (size_t task = 0; task < ku_tasks; task++) {
                    cmerge_write_keepusing_var_thread(&ku_args_array[task]);
                }
            }
            free(ku_args_array);
        }
    }

    /* ===================================================================
     * Step 9: Write _merge variable (after the update flags are known)
     * Note: SF_vstore is the only Stata API for writing values.
     * Workers write disjoint, contiguous row ranges of the one variable.
     * With update, matched rows become 4 (missing updated) or 5 (nonmissing
     * conflict, which takes precedence), as with merge.
     * =================================================================== */

    size_t n_updated = 0, n_conflict = 0;
    int values_changed = 0;   /* update/replace rewrote some master values */
    if (row_flags) {
        for (size_t i = 0; i < output_nobs; i++) {
            if (output_specs[i].merge_result != MERGE_RESULT_BOTH) continue;
            if (row_flags[i] & CMERGE_ROW_CONFLICT) n_conflict++;
            else if (row_flags[i] & CMERGE_ROW_UPDATED) n_updated++;
            if ((row_flags[i] & CMERGE_ROW_UPDATED) ||
                (replace_mode && (row_flags[i] & CMERGE_ROW_CONFLICT)))
                values_changed = 1;
        }
    }
    if (merge_var_idx > 0) {
        ST_int mvar = (ST_int)merge_var_idx;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(output_nobs >= MIN_OBS_PER_THREAD * 2)
        #endif
        for (size_t i = 0; i < output_nobs; i++) {
            double code = (double)output_specs[i].merge_result;
            if (row_flags && output_specs[i].merge_result == MERGE_RESULT_BOTH) {
                if (row_flags[i] & CMERGE_ROW_CONFLICT) code = 5.0;
                else if (row_flags[i] & CMERGE_ROW_UPDATED) code = 4.0;
            }
            (_stata_)->safestore(mvar, (ST_int)(i + 1), code);  /* checked: concurrent */
        }
    }
    free(row_flags);

    double t_write_meta = ctools_timer_ms() - t_start;

    /* ===================================================================
     * Cleanup: the undo copy stays until commit or rollback
     * =================================================================== */

    double t_cleanup_start = ctools_timer_ms();
    const size_t n_master_only = g_plan.n_master_only;
    const size_t n_matched = g_plan.n_matched;
    const int is_identity_permutation = g_plan.identity;
    const double t_sort = g_plan.t_sort, t_merge = g_plan.t_join;
    const double t_reorder = g_plan.t_reorder, t_plan = g_plan.t_total;
    cmerge_plan_free();

    /* Free using cache (protected by lock) */
    cache_lock_acquire();
    ctools_filtered_data_free(&g_using_cache.keys);
    ctools_filtered_data_free(&g_using_cache.keepusing);
    g_using_cache.loaded = 0;
    g_cache_in_use = 0;
    cache_lock_release();

    double t_cleanup = ctools_timer_ms() - t_cleanup_start;

    /* Calculate total time (plan and execute) */
    double t_total = t_plan + ctools_timer_ms() - t_total_start;

    /* Save merge results to Stata scalars */
    SF_scal_save("_cmerge_N", (double)output_nobs);
    SF_scal_save("_cmerge_N_1", (double)n_master_only);
    SF_scal_save("_cmerge_N_2", (double)n_using_only);
    SF_scal_save("_cmerge_N_3", (double)(n_matched - n_updated - n_conflict));
    SF_scal_save("_cmerge_N_4", (double)n_updated);
    SF_scal_save("_cmerge_N_5", (double)n_conflict);
    /* 1 when the output rows are the master rows in their original order
     * and no master value was rewritten (a sort variable may have changed) */
    SF_scal_save("_cmerge_same_order",
                 (double)(is_identity_permutation && n_using_only == 0 &&
                          output_nobs == master_nobs && !values_changed));

    /* Save Phase 2 timing to Stata scalars (in seconds) */
    SF_scal_save("_cmerge_p2_load_master", t_load / 1000.0);
    SF_scal_save("_cmerge_p2_sort_master", t_sort / 1000.0);
    SF_scal_save("_cmerge_p2_merge_join", t_merge / 1000.0);
    SF_scal_save("_cmerge_p2_reorder", t_reorder / 1000.0);
    SF_scal_save("_cmerge_p2_permute", t_permute / 1000.0);
    SF_scal_save("_cmerge_p2_stream", 0.0);
    SF_scal_save("_cmerge_p2_store", t_store / 1000.0);
    SF_scal_save("_cmerge_p2_write_meta", t_write_meta / 1000.0);
    SF_scal_save("_cmerge_p2_cleanup", t_cleanup / 1000.0);
    SF_scal_save("_cmerge_p2_total", t_total / 1000.0);
    SF_scal_save("_cmerge_n_output_vars", (double)master_nvars);
    CTOOLS_SAVE_THREAD_INFO("_cmerge");

    return STATA_OK;
}

/* Write the original observations back from the undo copy, then release
 * everything the merge held. The varlist must list the original master
 * variables first, in their original order and with the storage types they
 * had during execute; the dataset must have at least the original number of
 * observations. Without an undo copy only the release happens. */
static ST_retcode cmerge_rollback(void)
{
    stata_retcode rc = g_undo.held ? cmerge_undo_restore() : STATA_OK;
    cmerge_release_all();
    if (rc != STATA_OK) {
        SF_error("cmerge: could not restore the master observations\n");
        return ctools_stata_rc(rc);
    }
    return 0;
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
        ctools_filtered_data_free(&g_using_cache.keys);
        ctools_filtered_data_free(&g_using_cache.keepusing);
        g_using_cache.loaded = 0;
    }

    cache_lock_release();
    cmerge_plan_free();
    cmerge_undo_free();
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
    if (strlen(args) >= sizeof(args_copy)) return 198;
    strncpy(args_copy, args, sizeof(args_copy) - 1);
    args_copy[sizeof(args_copy) - 1] = '\0';

    char *phase = strtok(args_copy, " ");
    const char *rest = args + strlen(phase);
    while (*rest == ' ') rest++;

    if (strcmp(phase, "load_using") == 0) {
        return cmerge_load_using(rest);
    }
    else if (strcmp(phase, "plan") == 0) {
        return cmerge_plan(rest);
    }
    else if (strcmp(phase, "execute") == 0) {
        return cmerge_execute(rest);
    }
    else if (strcmp(phase, "rollback") == 0) {
        return cmerge_rollback();
    }
    else if (strcmp(phase, "commit") == 0) {
        cmerge_release_all();
        return 0;
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

/* ============================================================================
 * Cleanup function for ctools_cleanup system
 * ============================================================================ */

void cmerge_cleanup_cache(void)
{
    /* The dispatcher calls this only after workers have joined. A previous
     * interrupted execute may have left the in-use marker set. */
    cache_lock_acquire();
    g_cache_in_use = 0;
    ctools_filtered_data_free(&g_using_cache.keys);
    ctools_filtered_data_free(&g_using_cache.keepusing);
    memset(&g_using_cache, 0, sizeof(g_using_cache));
    cache_lock_release();
    cmerge_plan_free();
    cmerge_undo_free();
}
