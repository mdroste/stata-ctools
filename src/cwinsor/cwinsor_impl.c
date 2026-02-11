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
#include "cwinsor_impl.h"

#define CWINSOR_MODULE "cwinsor"
#define DEFAULT_LOWER_PCTL 1.0
#define DEFAULT_UPPER_PCTL 99.0
#define MIN_OBS_WINSOR 3

/* Note: By-variable sorting now uses ctools_sort_dispatch() instead of qsort */

/* Apply permutation to reorder an array */
static void apply_permutation(double *data, const size_t *perm, size_t n, double *temp)
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
    int verbose;
} cwinsor_config;

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

/* Detect groups in a single pass (speculatively allocated at max size).
 * Replaces the previous two-pass approach (count_groups + detect_groups). */
static void detect_groups(double **by_data, size_t nobs, size_t nby,
                          group_info *groups, size_t *ngroups)
{
    if (nobs == 0) {
        *ngroups = 0;
        return;
    }

    const double miss = SV_missval;
    size_t g = 0;
    groups[0].start = 0;
    groups[0].count = 1;

    for (size_t i = 1; i < nobs; i++) {
        int same_group = 1;

        for (size_t b = 0; b < nby && same_group; b++) {
            double prev = by_data[b][i - 1];
            double curr = by_data[b][i];

            int prev_miss = (prev >= miss);
            int curr_miss = (curr >= miss);

            if (prev_miss && curr_miss) continue;
            if (prev_miss || curr_miss || prev != curr) {
                same_group = 0;
            }
        }

        if (same_group) {
            groups[g].count++;
        } else {
            g++;
            groups[g].start = i;
            groups[g].count = 1;
        }
    }

    *ngroups = g + 1;
}

/* ===========================================================================
   Main Entry Point
   =========================================================================== */

ST_retcode cwinsor_main(const char *args)
{
    double t_start, t_load, t_sort, t_groups, t_winsor, t_store;
    int *store_indices = NULL;   /* Global dataset positions for SF_vstore */
    int *load_indices = NULL;    /* Sequential plugin-local indices for SF_vdata */
    int *by_load_indices = NULL; /* Sequential plugin-local indices for by-vars */
    double **var_data = NULL;
    double **by_data = NULL;
    double **work_arrays = NULL;
    group_info *groups = NULL;
    size_t *sort_indices = NULL;
    double *temp_buffer = NULL;
    perm_idx_t *obs_map = NULL;  /* Maps filtered index to 1-based Stata obs */
    int num_threads = 1;

    t_start = ctools_timer_seconds();
    t_sort = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CWINSOR_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    cwinsor_config config;
    config.lower_pctl = ctools_parse_double_option(args, "p", DEFAULT_LOWER_PCTL);
    config.upper_pctl = ctools_parse_double_option(args, "q", DEFAULT_UPPER_PCTL);
    config.trim = ctools_parse_bool_option(args, "trim");
    config.verbose = ctools_parse_bool_option(args, "verbose");

    if (config.lower_pctl < 0 || config.lower_pctl > 100 ||
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

    char *end;
    size_t nvars = (size_t)strtol(p, &end, 10);
    if (end == p || nvars == 0) {
        ctools_error(CWINSOR_MODULE, "invalid nvars");
        return 198;
    }
    p = end;

    while (*p == ' ' || *p == '\t') p++;
    size_t nby = (size_t)strtol(p, &end, 10);
    if (end == p) {
        ctools_error(CWINSOR_MODULE, "invalid nby");
        return 198;
    }
    p = end;

    /* Parse global store indices (for SF_vstore which uses global dataset positions) */
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

    /* Generate sequential plugin-local load indices (stack-allocated).
     * Plugin call receives: target_var1 target_var2 ... by_var1 by_var2 ...
     * So plugin-local positions are: targets = 1..nvars, by-vars = nvars+1..nvars+nby */
    /* Heap-allocate index arrays (avoids VLA stack overflow for large nvars/nby) */
    load_indices = (int *)malloc(nvars * sizeof(int));
    if (!load_indices) {
        free(store_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }
    for (size_t i = 0; i < nvars; i++) {
        load_indices[i] = (int)(i + 1);
    }

    if (nby > 0) {
        by_load_indices = (int *)malloc(nby * sizeof(int));
        if (!by_load_indices) {
            free(load_indices);
            free(store_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
        for (size_t i = 0; i < nby; i++) {
            by_load_indices[i] = (int)(nvars + i + 1);
        }
    }

    #ifdef _OPENMP
    num_threads = ctools_get_max_threads();
    #endif

    /* === Load Phase === */
    double load_start = ctools_timer_seconds();

    /* Load target variables using plugin-local indices (handles if/in at load time).
     * Single-variable: row-parallel load (parallelizes across observations).
     * Multi-variable: column-parallel load (one thread per variable). */
    ctools_filtered_data target_filtered;
    ctools_filtered_data_init(&target_filtered);
    stata_retcode load_rc = ctools_data_load(&target_filtered, load_indices,
                                              nvars, 0, 0, 0);
    if (load_rc != STATA_OK) {
        free(store_indices);
        free(load_indices);
        free(by_load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    /* Get filtered observation count and obs_map */
    size_t nobs = target_filtered.data.nobs;
    obs_map = target_filtered.obs_map;

    if (nobs == 0) {
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        free(by_load_indices);
        ctools_error(CWINSOR_MODULE, "no observations");
        return 2000;
    }

    /* Create var_data pointer array to point to filtered data arrays */
    var_data = (double **)malloc(nvars * sizeof(double *));
    if (!var_data) {
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        free(by_load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }
    for (size_t v = 0; v < nvars; v++) {
        var_data[v] = target_filtered.data.vars[v].data.dbl;
    }

    /* Load by-variables if specified */
    ctools_filtered_data by_filtered;
    ctools_filtered_data_init(&by_filtered);

    if (nby > 0) {
        load_rc = ctools_data_load(&by_filtered, by_load_indices,
                                             nby, 0, 0, 0);
        if (load_rc != STATA_OK) {
            ctools_filtered_data_free(&target_filtered);
            free(var_data);
            free(store_indices);
            free(load_indices);
        free(by_load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }

        /* Create by_data pointer array */
        by_data = (double **)malloc(nby * sizeof(double *));
        if (!by_data) {
            ctools_filtered_data_free(&by_filtered);
            ctools_filtered_data_free(&target_filtered);
            free(var_data);
            free(store_indices);
            free(load_indices);
        free(by_load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
        for (size_t b = 0; b < nby; b++) {
            by_data[b] = by_filtered.data.vars[b].data.dbl;
        }
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort Phase using ctools_sort_dispatch() === */
    double sort_start = ctools_timer_seconds();

    if (nby > 0) {
        /* Build 1-based sort key indices for by-variables (heap-allocated to avoid VLA) */
        int *sort_keys = (int *)malloc(nby * sizeof(int));
        if (!sort_keys) {
            if (by_data) free(by_data);
            ctools_filtered_data_free(&by_filtered);
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
            free(by_load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
        for (size_t b = 0; b < nby; b++) {
            sort_keys[b] = (int)(b + 1);  /* 1-based indices into by_filtered.data */
        }

        /* Sort by-variables using ctools infrastructure (order only) */
        stata_retcode sort_rc = ctools_sort_dispatch(&by_filtered.data, sort_keys, nby, SORT_ALG_AUTO);

        free(sort_keys);

        if (sort_rc != STATA_OK) {
            if (by_data) free(by_data);
            ctools_filtered_data_free(&by_filtered);
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
            free(by_load_indices);
            ctools_error(CWINSOR_MODULE, "sort failed");
            return 920;
        }

        /* Allocate sort_indices for sorting and later un-sorting during store */
        sort_indices = (size_t *)malloc(nobs * sizeof(size_t));
        temp_buffer = (double *)malloc(nobs * sizeof(double));

        if (!sort_indices || !temp_buffer) {
            if (sort_indices) free(sort_indices);
            if (temp_buffer) free(temp_buffer);
            if (by_data) free(by_data);
            ctools_filtered_data_free(&by_filtered);
            free(var_data);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
        free(by_load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }

        /* Copy sort_order from by_filtered.data (converting perm_idx_t to size_t) */
        for (size_t i = 0; i < nobs; i++) {
            sort_indices[i] = (size_t)by_filtered.data.sort_order[i];
        }

        /* Apply permutation to target variable arrays */
        for (size_t v = 0; v < nvars; v++) {
            apply_permutation(var_data[v], sort_indices, nobs, temp_buffer);
        }

        /* Apply permutation to by-variable arrays */
        for (size_t b = 0; b < nby; b++) {
            apply_permutation(by_data[b], sort_indices, nobs, temp_buffer);
        }
    }

    t_sort = ctools_timer_seconds() - sort_start;

    /* === Group Detection Phase === */
    double groups_start = ctools_timer_seconds();

    /* Single-pass group detection with speculative allocation at max size */
    size_t ngroups;
    groups = (group_info *)malloc(nobs * sizeof(group_info));
    if (!groups) {
        if (sort_indices) free(sort_indices);
        if (temp_buffer) free(temp_buffer);
        if (by_data) free(by_data);
        free(var_data);
        ctools_filtered_data_free(&by_filtered);
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        free(by_load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    if (nby > 0) {
        detect_groups(by_data, nobs, nby, groups, &ngroups);
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

    t_groups = ctools_timer_seconds() - groups_start;

    /* === Winsorization Phase === */
    double winsor_start = ctools_timer_seconds();

    size_t max_group_size = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (groups[g].count > max_group_size)
            max_group_size = groups[g].count;
    }

    work_arrays = (double **)malloc(num_threads * sizeof(double *));
    if (!work_arrays) {
        free(groups);
        if (sort_indices) free(sort_indices);
        if (temp_buffer) free(temp_buffer);
        if (by_data) free(by_data);
        free(var_data);
        ctools_filtered_data_free(&by_filtered);
        ctools_filtered_data_free(&target_filtered);
        free(store_indices);
        free(load_indices);
        free(by_load_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    for (int t = 0; t < num_threads; t++) {
        work_arrays[t] = (double *)malloc(max_group_size * sizeof(double));
        if (!work_arrays[t]) {
            for (int j = 0; j < t; j++) free(work_arrays[j]);
            free(work_arrays);
            free(groups);
            if (sort_indices) free(sort_indices);
            if (temp_buffer) free(temp_buffer);
                if (by_data) free(by_data);
            free(var_data);
            ctools_filtered_data_free(&by_filtered);
            ctools_filtered_data_free(&target_filtered);
            free(store_indices);
            free(load_indices);
        free(by_load_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
    }

    if (ngroups == 1) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
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

            if (config.trim) {
                apply_trim_bounds(data, count, lower, upper);
            } else {
                apply_winsor_bounds(data, count, lower, upper);
            }
        }
    } else {
        size_t total_tasks = nvars * ngroups;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1)
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

    if (nby > 0 && sort_indices != NULL) {
        /* Store in original order using sort_indices directly.
         * sort_indices[j] = original (pre-sort) index of the j-th sorted element.
         * data[j] = winsorized value at sorted position j.
         * obs_map[sort_indices[j]] = Stata observation for that original index.
         * This avoids computing and storing an inverse permutation. */
        if (nvars == 1) {
            /* Single variable: compose obs_map for row-parallel store */
            perm_idx_t *composed_map = (perm_idx_t *)malloc(nobs * sizeof(perm_idx_t));
            if (composed_map) {
                for (size_t j = 0; j < nobs; j++) {
                    composed_map[j] = obs_map[sort_indices[j]];
                }
                ctools_store_filtered_rowpar(var_data[0], nobs, store_indices[0], composed_map);
                free(composed_map);
            } else {
                /* Fallback: sequential store */
                for (size_t j = 0; j < nobs; j++) {
                    SF_vstore(store_indices[0], (ST_int)obs_map[sort_indices[j]], var_data[0][j]);
                }
            }
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (size_t v = 0; v < nvars; v++) {
                int idx = store_indices[v];
                double *data = var_data[v];
                for (size_t j = 0; j < nobs; j++) {
                    SF_vstore(idx, (ST_int)obs_map[sort_indices[j]], data[j]);
                }
            }
        }
    } else {
        if (nvars == 1) {
            /* Single variable: row-parallel store */
            ctools_store_filtered_rowpar(var_data[0], nobs, store_indices[0], obs_map);
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (size_t v = 0; v < nvars; v++) {
                int idx = store_indices[v];
                double *data = var_data[v];
                for (size_t i = 0; i < nobs; i++) {
                    SF_vstore(idx, (ST_int)obs_map[i], data[i]);
                }
            }
        }
    }

    t_store = ctools_timer_seconds() - store_start;

    /* === Cleanup === */
    for (int t = 0; t < num_threads; t++) {
        free(work_arrays[t]);
    }
    free(work_arrays);
    free(groups);
    if (sort_indices) free(sort_indices);
    if (temp_buffer) free(temp_buffer);
    /* Free pointer arrays (not the data - owned by stata_data structs) */
    if (by_data) free(by_data);
    free(var_data);
    /* Free stata_data structures (which own the actual data arrays) */
    ctools_filtered_data_free(&by_filtered);
    ctools_filtered_data_free(&target_filtered);
    free(store_indices);
    free(load_indices);
    free(by_load_indices);

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
