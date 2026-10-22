/*
    cwinsor_impl.c
    High-performance winsorization for Stata

    Performance approach:
    - Uses ctools_data_load() for parallel data loading with if/in filtering
    - Uses ctools_sort_dispatch() for by-variable sorting
    - O(n) quickselect with single extraction per group
    - Parallel over (variable, group) pairs for by() operations

    Author: ctools project
    License: MIT
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stplugin.h"
#include "ctools_types.h"
#include "ctools_config.h"
#include "ctools_runtime.h"
#include "ctools_select.h"
#include "ctools_parse.h"
#include "ctools_simd.h"
#include "ctools_order.h"
#include "cwinsor_impl.h"

#define CWINSOR_MODULE "cwinsor"
#define DEFAULT_LOWER_PCTL 1.0
#define DEFAULT_UPPER_PCTL 99.0
#define MIN_OBS_WINSOR 1

/* Note: By-variable sorting now uses ctools_sort_dispatch() instead of qsort */

/* Apply permutation to reorder an array */
static void apply_permutation(double *data, const perm_idx_t *perm, size_t n, double *temp)
{
    /* Copy to temp in new order */
    for (size_t i = 0; i < n; i++) {
        temp[i] = data[perm[i]];
    }
    /* Copy back */
    memcpy(data, temp, n * sizeof(double));
}

/* ===========================================================================
   Quickselect - O(n) percentile computation
   Uses shared ctools_quickselect_double from ctools_select.h
   =========================================================================== */

/* Percentile via O(n) quickselect.
 *   h = n * p / 100,  j = floor(h),  g = h - j
 *   g > 0          → sorted[j]
 *   g == 0, j > 0  → (sorted[j-1] + sorted[j]) / 2
 *   g == 0, j == 0 → sorted[0]
 */
static double compute_percentile(double *work, size_t n, double pctl)
{
    double h = (double)n * pctl / 100.0;
    size_t j = (size_t)floor(h);
    double g = h - (double)j;

    if (j >= n)
        return ctools_quickselect_double(work, n, n - 1);

    if (g > 0.0)
        return ctools_quickselect_double(work, n, j);

    if (j == 0)
        return ctools_quickselect_double(work, n, 0);

    /* Exact integer position: average sorted[j-1] and sorted[j] */
    double v1 = ctools_quickselect_double(work, n, j - 1);
    /* After quickselect for j-1, work[j..n-1] >= v1; min of that range is sorted[j] */
    double v2 = work[j];
    for (size_t i = j + 1; i < n; i++) {
        if (work[i] < v2) v2 = work[i];
    }
    return (v1 + v2) / 2.0;
}

/* Compute two percentiles from a work array of non-missing values.
 * Leverages partial ordering: after computing lo, the array is partitioned
 * around position k_lo, so hi only searches work[split..n-1]. */
static void compute_two_percentiles(double *work, size_t n,
                                    double pctl_lo, double pctl_hi,
                                    double *result_lo, double *result_hi)
{
    if (n == 0) {
        *result_lo = SV_missval;
        *result_hi = SV_missval;
        return;
    }
    if (n == 1) {
        *result_lo = work[0];
        *result_hi = work[0];
        return;
    }

    /* Compute lo — partially sorts array around position k_lo */
    *result_lo = compute_percentile(work, n, pctl_lo);

    /* Determine the partition point from the lo computation.
     * After quickselect at position p: work[0..p] <= work[p+1..n-1].
     * For g>0 or j==0: quickselect at j, partition at j.
     * For g==0 && j>0: quickselect at j-1, partition at j-1. */
    double h_lo = (double)n * pctl_lo / 100.0;
    size_t k_lo = (size_t)floor(h_lo);
    if (k_lo >= n) k_lo = n - 1;
    double g_lo = h_lo - (double)k_lo;
    size_t split = (g_lo == 0.0 && k_lo > 0) ? k_lo : (k_lo + 1 < n ? k_lo + 1 : n);

    /* Compute hi percentile position */
    double h_hi = (double)n * pctl_hi / 100.0;
    size_t j_hi = (size_t)floor(h_hi);
    double g_hi = h_hi - (double)j_hi;

    /* If hi position is beyond the partition, search only work[split..n-1] */
    if (j_hi >= split && split < n) {
        size_t off = split;
        size_t len = n - off;

        if (j_hi >= n) {
            /* Maximum of subarray */
            double mx = work[off];
            for (size_t i = off + 1; i < n; i++) {
                if (work[i] > mx) mx = work[i];
            }
            *result_hi = mx;
        } else if (g_hi > 0.0 || j_hi == 0) {
            *result_hi = ctools_quickselect_double(work + off, len, j_hi - off);
        } else {
            /* Exact integer: average sorted[j_hi-1] and sorted[j_hi] */
            double v1;
            if (j_hi - 1 >= off) {
                v1 = ctools_quickselect_double(work + off, len, j_hi - 1 - off);
            } else {
                /* sorted[j_hi-1] is max of work[0..j_hi-1] (left partition) */
                v1 = work[0];
                for (size_t i = 1; i < j_hi; i++) {
                    if (work[i] > v1) v1 = work[i];
                }
            }
            double v2 = work[j_hi];
            for (size_t i = j_hi + 1; i < n; i++) {
                if (work[i] < v2) v2 = work[i];
            }
            *result_hi = (v1 + v2) / 2.0;
        }
    } else {
        /* Rare: hi position within partition, use full computation */
        *result_hi = compute_percentile(work, n, pctl_hi);
    }
}

/* ===========================================================================
   Configuration
   =========================================================================== */

typedef struct {
    double lower_pctl;
    double upper_pctl;
    int trim;
    int return_bounds;
    int verbose;
    int undo;          /* keep the loaded target values for cwinsor rollback */
    int fault_store;   /* validation hook: fail after storing half the rows */
} cwinsor_config;

/* ===========================================================================
   Undo copy for in-place calls
   With option undo the main call keeps the target values it loaded, and
   their observation numbers, until the ado ends the command with
   "cwinsor commit" or "cwinsor rollback" (or the runtime resets it). A store
   that fails part way is undone at once, before the call returns.
   =========================================================================== */

static struct {
    double **values;   /* nvars columns of the loaded (original) values */
    perm_idx_t *obs;   /* 1-based Stata observation of each loaded row */
    size_t nvars, nobs;
} undo;

void cwinsor_cleanup_state(void)
{
    if (undo.values) {
        for (size_t v = 0; v < undo.nvars; v++) free(undo.values[v]);
        free(undo.values);
    }
    free(undo.obs);
    memset(&undo, 0, sizeof(undo));
}

/* Copy the original values before anything changes them; 920 if there is no
 * memory, which fails the call before any change. */
static int cwinsor_undo_save(double **var_data, size_t nvars, const perm_idx_t *obs_map, size_t nobs)
{
    cwinsor_cleanup_state();
    undo.values = (double **)calloc(nvars, sizeof(double *));
    undo.obs = (perm_idx_t *)ctools_safe_malloc2(nobs, sizeof(perm_idx_t));
    if (!undo.values || !undo.obs) {
        cwinsor_cleanup_state();
        return 920;
    }
    undo.nvars = nvars;
    undo.nobs = nobs;
    for (size_t v = 0; v < nvars; v++) {
        undo.values[v] = (double *)ctools_safe_malloc2(nobs, sizeof(double));
        if (!undo.values[v]) {
            cwinsor_cleanup_state();
            return 920;
        }
        memcpy(undo.values[v], var_data[v], nobs * sizeof(double));
    }
    memcpy(undo.obs, obs_map, nobs * sizeof(perm_idx_t));
    return 0;
}

/* Write the original values back into the targets, which are the first
 * nvars variables of the plugin varlist, then free the copy. */
static int cwinsor_undo_restore(void)
{
    int rc = STATA_OK;
    if (undo.nvars == 1) {
        rc = ctools_store_filtered_rowpar(undo.values[0], undo.nobs, 1, undo.obs);
    } else {
        #pragma omp parallel for schedule(static) reduction(max:rc)
        for (size_t v = 0; v < undo.nvars; v++) {
            int vrc = ctools_store_filtered(undo.values[v], undo.nobs, (int)v + 1, undo.obs);
            if (vrc > rc) rc = vrc;
        }
    }
    cwinsor_cleanup_state();
    return rc;
}

/* ===========================================================================
   Core Functions
   =========================================================================== */

static size_t extract_nonmissing(const double *data, size_t nobs, double *work)
{
    const double miss = SV_missval;
    size_t n = 0;
    for (size_t i = 0; i < nobs; i++) {
        double val = data[i];
        if (val < miss) {
            work[n++] = val;
        }
    }
    return n;
}

static void apply_winsor_bounds(double *data, size_t nobs, double lower, double upper)
{
    ctools_simd_clamp(data, nobs, lower, upper, SV_missval);
}

static void apply_trim_bounds(double *data, size_t nobs, double lower, double upper)
{
    ctools_simd_replace_oob(data, nobs, lower, upper, SV_missval);
}

/* ===========================================================================
   Group Detection
   =========================================================================== */

typedef struct {
    size_t start;
    size_t count;
} group_info;

/* ===========================================================================
   Main Entry Point
   =========================================================================== */

ST_retcode cwinsor_main(const char *args)
{
    double t_start, t_load, t_sort, t_groups, t_winsor, t_store;
    int *store_indices = NULL;   /* Plugin-local positions for checked stores */
    int *load_indices = NULL;    /* Sequential plugin-local indices for SF_vdata */
    double **var_data = NULL;
    double **work_arrays = NULL;
    group_info *groups = NULL;
    const perm_idx_t *sort_indices = NULL; /* Borrowed from loaded data */
    double *temp_buffer = NULL;
    perm_idx_t *obs_map = NULL;  /* Maps filtered index to 1-based Stata obs */
    int num_threads = 1;

    t_start = ctools_timer_seconds();
    t_sort = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CWINSOR_MODULE, "no arguments specified");
        return 198;
    }

    /* End an in-place call: put the original values back, or drop the copy. */
    if (!strcmp(args, "commit")) {
        cwinsor_cleanup_state();
        return 0;
    }
    if (!strcmp(args, "rollback")) {
        if (!undo.values) {
            ctools_error(CWINSOR_MODULE, "no original values to restore");
            return 459;
        }
        if (SF_nvars() < (int)undo.nvars) {
            cwinsor_cleanup_state();
            return 198;
        }
        return ctools_stata_rc(cwinsor_undo_restore());
    }
    cwinsor_cleanup_state();

    /* Parse options */
    cwinsor_config config;
    config.lower_pctl = ctools_parse_double_option(args, "p", DEFAULT_LOWER_PCTL);
    config.upper_pctl = ctools_parse_double_option(args, "q", DEFAULT_UPPER_PCTL);
    config.trim = ctools_parse_bool_option(args, "trim");
    config.return_bounds = ctools_parse_bool_option(args, "returnbounds");
    config.verbose = ctools_parse_bool_option(args, "verbose");
    config.undo = ctools_parse_bool_option(args, "undo");
    config.fault_store = ctools_parse_bool_option(args, "fault_store");

    if (!isfinite(config.lower_pctl) || !isfinite(config.upper_pctl) ||
        config.lower_pctl < 0 || config.lower_pctl > 100 ||
        config.upper_pctl < 0 || config.upper_pctl > 100) {
        ctools_error(CWINSOR_MODULE, "percentiles must be between 0 and 100");
        return 198;
    }
    if (config.lower_pctl >= config.upper_pctl) {
        ctools_error(CWINSOR_MODULE, "lower percentile must be less than upper");
        return 198;
    }

    /* Parse: nvars nby store_indices... options... */
    const char *p = args;
    while (*p == ' ' || *p == '\t') p++;

    int nvars_value;
    if (ctools_parse_next_int(&p, &nvars_value) || nvars_value < 0 ||
        nvars_value > SF_nvars() || nvars_value == 0) {
        ctools_error(CWINSOR_MODULE, "invalid nvars");
        return 198;
    }
    size_t nvars = (size_t)nvars_value;

    while (*p == ' ' || *p == '\t') p++;
    int nby_value;
    if (ctools_parse_next_int(&p, &nby_value) || nby_value < 0 ||
        nby_value > SF_nvars() - nvars_value) {
        ctools_error(CWINSOR_MODULE, "invalid nby");
        return 198;
    }
    size_t nby = (size_t)nby_value;
    size_t nload = nvars + nby;
    char bounds_name[256] = {0};
    if (config.return_bounds) {
        if (nby != 0) return 198;
        ST_retcode rc = SF_macro_use("_cwinsor_bounds", bounds_name, sizeof(bounds_name));
        if (rc) return rc;
        if (!bounds_name[0]) return 198;
    }


    /* Parse plugin-local store indices. */
    store_indices = (int *)malloc(nvars * sizeof(int));
    if (!store_indices) {
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }
    if (ctools_parse_int_array(store_indices, nvars, &p) != 0) {
        free(store_indices);
        ctools_error(CWINSOR_MODULE, "failed to parse store indices");
        return 198;
    }

    /* Load sources and by-variables together, sharing one selection/map.
     * Separate output columns, when present, follow these input columns and
     * are addressed by store_indices; in-place sources are also the targets. */
    load_indices = (int *)malloc(nload * sizeof(int));
    if (!load_indices) {
        free(store_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }
    for (size_t i = 0; i < nload; i++) {
        load_indices[i] = (int)(i + 1);
    }

    #ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    #endif

    /* === Load Phase === */
    double load_start = ctools_timer_seconds();

    /* The loaded dataset owns all columns, the sort order, and the row map. */
    ctools_filtered_data target_filtered;
    ctools_filtered_data_init(&target_filtered);
    int load_flags = nby ? 0 : CTOOLS_LOAD_NO_SORT_ORDER;
    if (ctools_parse_bool_option(args, "skipif")) load_flags |= CTOOLS_LOAD_SKIP_IF;
    stata_retcode load_rc = ctools_data_load(&target_filtered, load_indices,
                                              nload, 0, 0, load_flags);
    if (load_rc != STATA_OK) {
        free(store_indices);
        free(load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return ctools_stata_rc(load_rc);
    }

    /* Get filtered observation count and obs_map */
    size_t nobs = target_filtered.data.nobs;
    obs_map = target_filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        ctools_error(CWINSOR_MODULE, "no observations");
        return 2000;
    }

    /* Create var_data pointer array to point to filtered data arrays */
    var_data = (double **)malloc(nvars * sizeof(double *));
    if (!var_data) {
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }
    for (size_t v = 0; v < nvars; v++) {
        var_data[v] = target_filtered.data.vars[v].data.dbl;
    }
    if (ctools_parse_bool_option(args, "newoutputs")) {
        /* gen new = source if !missing(source) initializes all missing codes
         * as ordinary missing. Match it while keeping the inputs read-only. */
        for (size_t v = 0; v < nvars; v++)
            for (size_t i = 0; i < nobs; i++)
                if (var_data[v][i] >= SV_missval) var_data[v][i] = SV_missval;
    }

    /* The undo copy is taken in load order, matching obs_map. */
    if (config.undo && cwinsor_undo_save(var_data, nvars, obs_map, nobs)) {
        free(var_data);
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort Phase using ctools_sort_dispatch() === */
    double sort_start = ctools_timer_seconds();

    if (nby > 0) {
        /* Build 1-based sort key indices for by-variables (heap-allocated to avoid VLA) */
        int *sort_keys = (int *)malloc(nby * sizeof(int));
        if (!sort_keys) {
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
        for (size_t b = 0; b < nby; b++) {
            sort_keys[b] = (int)(nvars + b + 1);  /* 1-based loaded column */
        }

        /* Sort by-variables using ctools infrastructure (order only) */
        stata_retcode sort_rc = ctools_sort_dispatch(&target_filtered.data, sort_keys, nby, SORT_ALG_AUTO);

        free(sort_keys);

        if (sort_rc != STATA_OK) {
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
            ctools_error(CWINSOR_MODULE, "sort failed");
            return 920;
        }

        /* Reuse the shared permutation for gather and write-back. */
        sort_indices = target_filtered.data.sort_order;
        temp_buffer = (double *)malloc(nobs * sizeof(double));
        if (!temp_buffer) {
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }

        /* Apply permutation to target variable arrays. The by-variables stay
         * in load order; group detection reads them through sort_order. */
        for (size_t v = 0; v < nvars; v++) {
            apply_permutation(var_data[v], sort_indices, nobs, temp_buffer);
        }
    }

    t_sort = ctools_timer_seconds() - sort_start;

    /* === Group Detection Phase === */
    double groups_start = ctools_timer_seconds();

    /* Groups are Stata by-groups: exact by-values in sorted order, with
     * strings compared byte-wise and each missing code its own group. */
    size_t ngroups;
    groups = (group_info *)malloc(nobs * sizeof(group_info));
    size_t *group_starts = (size_t *)malloc((nobs + 1) * sizeof(size_t));
    int *group_keys = (int *)malloc((nby > 0 ? nby : 1) * sizeof(int));
    if (!groups || !group_starts || !group_keys) {
        free(groups);
        free(group_starts);
        free(group_keys);
        if (temp_buffer) free(temp_buffer);
        free(var_data);
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    if (nby > 0) {
        for (size_t b = 0; b < nby; b++) group_keys[b] = (int)(nvars + b);
        ngroups = ctools_group_starts(&target_filtered.data, target_filtered.data.sort_order,
                                      group_keys, nby, group_starts);
        for (size_t g = 0; g < ngroups; g++) {
            groups[g].start = group_starts[g];
            groups[g].count = group_starts[g + 1] - group_starts[g];
        }
        /* Trim to actual size */
        if (ngroups < nobs) {
            group_info *trimmed = (group_info *)realloc(groups, ngroups * sizeof(group_info));
            if (trimmed) groups = trimmed;
        }
    } else {
        groups[0].start = 0;
        groups[0].count = nobs;
        ngroups = 1;
    }
    free(group_starts);
    free(group_keys);

    t_groups = ctools_timer_seconds() - groups_start;

    /* === Winsorization Phase === */
    double winsor_start = ctools_timer_seconds();

    size_t max_group_size = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (groups[g].count > max_group_size)
            max_group_size = groups[g].count;
    }

    /* One work buffer per thread that can get a task: without by() there is
     * one task per variable, so few variables need few buffers. The parallel
     * loops below run exactly nwork threads, which index these buffers. */
    size_t total_tasks = ngroups == 1 ? nvars : nvars * ngroups;
    int nwork = (int)(total_tasks < (size_t)num_threads ? total_tasks : (size_t)num_threads);
    if (nwork < 1) nwork = 1;

    work_arrays = (double **)malloc((size_t)nwork * sizeof(double *));
    if (!work_arrays) {
        free(groups);
        if (temp_buffer) free(temp_buffer);
        free(var_data);
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    for (int t = 0; t < nwork; t++) {
        work_arrays[t] = (double *)malloc(max_group_size * sizeof(double));
        if (!work_arrays[t]) {
            for (int j = 0; j < t; j++) free(work_arrays[j]);
            free(work_arrays);
            free(groups);
            if (temp_buffer) free(temp_buffer);
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
    }

    /* Only the caller thread may use the host matrix API. Workers collect
     * exact quantiles here for the ado's winsor2-compatible macro expansion. */
    int store_rc = STATA_OK;
    double *bounds = NULL;
    if (config.return_bounds) {
        bounds = ctools_safe_malloc2(nvars, 2 * sizeof(double));
        if (!bounds) {
            store_rc = STATA_ERR_MEMORY;
            goto cleanup;
        }
        for (size_t v = 0; v < 2 * nvars; v++) bounds[v] = SV_missval;
    }

    if (ngroups == 1) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(nwork)
        #endif
        for (size_t v = 0; v < nvars; v++) {
            #ifdef _OPENMP
            int tid = omp_get_thread_num();
            #else
            int tid = 0;
            #endif

            double *data = var_data[v];
            double *work = work_arrays[tid];
            size_t count = groups[0].count;

            if (count < MIN_OBS_WINSOR) continue;

            size_t n = extract_nonmissing(data, count, work);
            if (n < MIN_OBS_WINSOR) continue;

            double lower, upper;
            compute_two_percentiles(work, n, config.lower_pctl, config.upper_pctl,
                                    &lower, &upper);

            if (lower > upper) {
                double tmp = lower;
                lower = upper;
                upper = tmp;
            }

            if (config.return_bounds) {
                bounds[2 * v] = lower;
                bounds[2 * v + 1] = upper;
            } else if (config.trim) {
                apply_trim_bounds(data, count, lower, upper);
            } else {
                apply_winsor_bounds(data, count, lower, upper);
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nwork)
        #endif
        for (size_t task = 0; task < total_tasks; task++) {
            size_t v = task / ngroups;
            size_t g = task % ngroups;

            #ifdef _OPENMP
            int tid = omp_get_thread_num();
            #else
            int tid = 0;
            #endif

            double *data = var_data[v];
            double *work = work_arrays[tid];
            size_t start = groups[g].start;
            size_t count = groups[g].count;

            if (count < MIN_OBS_WINSOR) continue;

            size_t n = extract_nonmissing(&data[start], count, work);
            if (n < MIN_OBS_WINSOR) continue;

            double lower, upper;
            compute_two_percentiles(work, n, config.lower_pctl, config.upper_pctl,
                                    &lower, &upper);

            if (lower > upper) {
                double tmp = lower;
                lower = upper;
                upper = tmp;
            }

            if (config.trim) {
                apply_trim_bounds(&data[start], count, lower, upper);
            } else {
                apply_winsor_bounds(&data[start], count, lower, upper);
            }
        }
    }

    t_winsor = ctools_timer_seconds() - winsor_start;

    /* === Store Phase using obs_map === */
    double store_start = ctools_timer_seconds();

    if (config.return_bounds) {
        for (size_t v = 0; v < nvars && !store_rc; v++) {
            store_rc = SF_mat_store(bounds_name, v + 1, 1, bounds[2 * v]);
            if (!store_rc)
                store_rc = SF_mat_store(bounds_name, v + 1, 2, bounds[2 * v + 1]);
        }
    }
    perm_idx_t *write_map = obs_map;
    perm_idx_t *composed_map = NULL;
    if (nby > 0 && sort_indices != NULL) {
        composed_map = (perm_idx_t *)ctools_safe_malloc2(nobs, sizeof(perm_idx_t));
        if (!composed_map) store_rc = STATA_ERR_MEMORY;
        else {
            for (size_t j = 0; j < nobs; j++) composed_map[j] = obs_map[sort_indices[j]];
            write_map = composed_map;
        }
    }
    if (!store_rc && !config.return_bounds) {
        /* fault_store writes half the rows and then fails as memory would */
        size_t n_store = config.fault_store ? nobs / 2 : nobs;
        if (nvars == 1) {
            store_rc = ctools_store_filtered_rowpar(var_data[0], n_store, store_indices[0], write_map);
        } else {
            #pragma omp parallel for schedule(static) reduction(max:store_rc)
            for (size_t v = 0; v < nvars; v++) {
                int rc = ctools_store_filtered(var_data[v], n_store, store_indices[v], write_map);
                if (rc > store_rc) store_rc = rc;
            }
        }
        if (!store_rc && config.fault_store) store_rc = STATA_ERR_MEMORY;
        /* A store that failed part way is undone before returning: the
         * runtime frees the copy after any failed call. */
        if (store_rc && undo.values) {
            int restore_rc = cwinsor_undo_restore();
            if (restore_rc) ctools_error(CWINSOR_MODULE, "the original values could not be restored");
        }
    }
    free(composed_map);

    t_store = ctools_timer_seconds() - store_start;

    /* === Cleanup === */
cleanup:
    free(bounds);
    for (int t = 0; t < nwork; t++) {
        free(work_arrays[t]);
    }
    free(work_arrays);
    free(groups);
    if (temp_buffer) free(temp_buffer);
    /* Free pointer arrays (not the data - owned by stata_data structs) */
    free(var_data);
    /* Free stata_data structures (which own the actual data arrays) */
    ctools_filtered_data_free(&target_filtered);
    free(store_indices);
    free(load_indices);

    if (store_rc) {
        cwinsor_cleanup_state();
        return ctools_stata_rc(store_rc);
    }

    /* Store timing scalars */
    double total = ctools_timer_seconds() - t_start;
    SF_scal_save("_cwinsor_time_load", t_load);
    SF_scal_save("_cwinsor_time_sort", t_sort);
    SF_scal_save("_cwinsor_time_groups", t_groups);
    SF_scal_save("_cwinsor_time_winsor", t_winsor);
    SF_scal_save("_cwinsor_time_store", t_store);
    SF_scal_save("_cwinsor_time_total", total);
    SF_scal_save("_cwinsor_nobs", (double)nobs);
    SF_scal_save("_cwinsor_nvars", (double)nvars);
    SF_scal_save("_cwinsor_ngroups", (double)ngroups);

    #ifdef _OPENMP
    SF_scal_save("_cwinsor_threads", (double)num_threads);
    #else
    SF_scal_save("_cwinsor_threads", 1.0);
    #endif

    /* Verbose output handled by .ado file using stored scalars */

    return 0;
}
