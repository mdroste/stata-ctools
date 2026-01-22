/*
    cwinsor_impl.c
    High-performance winsorization for Stata

    Performance approach:
    - Internal sorting by by-variables (avoids slow Stata sort)
    - O(n) quickselect with single extraction per group
    - Parallel over (variable, group) pairs for by() operations
    - Parallel data loading across all columns

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
#include "ctools_timer.h"
#include "ctools_error.h"
#include "cwinsor_impl.h"

#define CWINSOR_MODULE "cwinsor"
#define DEFAULT_LOWER_PCTL 1.0
#define DEFAULT_UPPER_PCTL 99.0
#define MIN_OBS_WINSOR 3

/* ===========================================================================
   Argsort for by-variable sorting
   =========================================================================== */

typedef struct {
    double **by_data;
    size_t nby;
} sort_context;

static sort_context g_sort_ctx;

/* Compare function for qsort - compares by-variable values */
static int compare_by_indices(const void *a, const void *b)
{
    size_t ia = *(const size_t *)a;
    size_t ib = *(const size_t *)b;
    const double miss = SV_missval;

    for (size_t k = 0; k < g_sort_ctx.nby; k++) {
        double va = g_sort_ctx.by_data[k][ia];
        double vb = g_sort_ctx.by_data[k][ib];

        int ma = (va >= miss);
        int mb = (vb >= miss);

        /* Missing values sort to end */
        if (ma && mb) continue;
        if (ma) return 1;
        if (mb) return -1;

        if (va < vb) return -1;
        if (va > vb) return 1;
    }
    /* Stable sort: preserve original order for equal keys */
    return (ia < ib) ? -1 : (ia > ib) ? 1 : 0;
}

/* Create sorted indices based on by-variables */
static void argsort_by_vars(double **by_data, size_t nobs, size_t nby, size_t *indices)
{
    /* Initialize indices */
    for (size_t i = 0; i < nobs; i++) {
        indices[i] = i;
    }

    /* Set up global context for comparison */
    g_sort_ctx.by_data = by_data;
    g_sort_ctx.nby = nby;

    /* Sort indices */
    qsort(indices, nobs, sizeof(size_t), compare_by_indices);
}

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
   =========================================================================== */

static inline void swap_double(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

/* Median-of-three pivot selection + partition */
static size_t partition(double *arr, size_t left, size_t right)
{
    size_t mid = left + (right - left) / 2;

    if (arr[mid] < arr[left]) swap_double(&arr[left], &arr[mid]);
    if (arr[right] < arr[left]) swap_double(&arr[left], &arr[right]);
    if (arr[right] < arr[mid]) swap_double(&arr[mid], &arr[right]);

    swap_double(&arr[mid], &arr[right - 1]);
    double pivot = arr[right - 1];

    size_t i = left;
    for (size_t j = left; j < right - 1; j++) {
        if (arr[j] <= pivot) {
            swap_double(&arr[i], &arr[j]);
            i++;
        }
    }
    swap_double(&arr[i], &arr[right - 1]);
    return i;
}

static double quickselect(double *arr, size_t n, size_t k)
{
    if (n == 0) return SV_missval;
    if (n == 1) return arr[0];
    if (k >= n) k = n - 1;

    size_t left = 0;
    size_t right = n - 1;

    while (left < right) {
        if (right - left < 10) {
            for (size_t i = left + 1; i <= right; i++) {
                double key = arr[i];
                size_t j = i;
                while (j > left && arr[j - 1] > key) {
                    arr[j] = arr[j - 1];
                    j--;
                }
                arr[j] = key;
            }
            return arr[k];
        }

        size_t pivot_idx = partition(arr, left, right);
        if (k == pivot_idx) return arr[k];
        else if (k < pivot_idx) right = pivot_idx - 1;
        else left = pivot_idx + 1;
    }
    return arr[left];
}

/* Compute two percentiles from a single extraction */
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

    double pos_lo = (pctl_lo / 100.0) * (double)(n - 1);
    double pos_hi = (pctl_hi / 100.0) * (double)(n - 1);

    size_t k_lo = (size_t)floor(pos_lo);
    if (k_lo >= n) k_lo = n - 1;

    double val_lo = quickselect(work, n, k_lo);

    double frac_lo = pos_lo - (double)k_lo;
    if (frac_lo > 0.0 && k_lo + 1 < n) {
        double next_val = work[k_lo + 1];
        for (size_t i = k_lo + 2; i < n; i++) {
            if (work[i] < next_val) next_val = work[i];
        }
        val_lo = val_lo + frac_lo * (next_val - val_lo);
    }
    *result_lo = val_lo;

    size_t k_hi_floor = (size_t)floor(pos_hi);
    if (k_hi_floor >= n) k_hi_floor = n - 1;

    double val_hi_floor = quickselect(work, n, k_hi_floor);

    double frac_hi = pos_hi - (double)k_hi_floor;
    if (frac_hi > 0.0 && k_hi_floor + 1 < n) {
        double next_val = work[k_hi_floor + 1];
        for (size_t i = k_hi_floor + 2; i < n; i++) {
            if (work[i] < next_val) next_val = work[i];
        }
        val_hi_floor = val_hi_floor + frac_hi * (next_val - val_hi_floor);
    }
    *result_hi = val_hi_floor;
}

/* ===========================================================================
   Configuration and Argument Parsing
   =========================================================================== */

typedef struct {
    double lower_pctl;
    double upper_pctl;
    int trim;
    int verbose;
} cwinsor_config;

static int parse_option(const char *args, const char *name)
{
    if (!args || !name) return 0;
    const char *p = strstr(args, name);
    if (!p) return 0;
    if (p > args && *(p-1) != ' ' && *(p-1) != '\t') return 0;
    char after = *(p + strlen(name));
    if (after != '\0' && after != ' ' && after != '\t' && after != '=') return 0;
    return 1;
}

static double parse_double_option(const char *args, const char *key, double def)
{
    if (!args || !key) return def;
    char pat[64];
    snprintf(pat, sizeof(pat), "%s=", key);
    const char *p = strstr(args, pat);
    if (!p) return def;
    return atof(p + strlen(pat));
}

static int parse_int_array(int *arr, size_t count, const char **start)
{
    const char *p = *start;
    for (size_t i = 0; i < count; i++) {
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') return -1;
        char *end;
        arr[i] = (int)strtol(p, &end, 10);
        if (end == p) return -1;
        p = end;
    }
    *start = p;
    return 0;
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
    const double miss = SV_missval;
    for (size_t i = 0; i < nobs; i++) {
        double val = data[i];
        if (val < miss) {
            if (val < lower) data[i] = lower;
            else if (val > upper) data[i] = upper;
        }
    }
}

static void apply_trim_bounds(double *data, size_t nobs, double lower, double upper)
{
    const double miss = SV_missval;
    for (size_t i = 0; i < nobs; i++) {
        double val = data[i];
        if (val < miss) {
            if (val < lower || val > upper) data[i] = miss;
        }
    }
}

/* ===========================================================================
   Group Detection
   =========================================================================== */

typedef struct {
    size_t start;
    size_t count;
} group_info;

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
    int *var_indices = NULL;
    int *by_indices = NULL;
    double **var_data = NULL;
    double **by_data = NULL;
    double **work_arrays = NULL;
    group_info *groups = NULL;
    size_t *sort_indices = NULL;
    double *temp_buffer = NULL;
    size_t *inv_perm = NULL;
    int num_threads = 1;

    t_start = ctools_timer_seconds();
    t_sort = 0.0;

    if (!args || strlen(args) == 0) {
        ctools_error(CWINSOR_MODULE, "no arguments specified");
        return 198;
    }

    /* Parse options */
    cwinsor_config config;
    config.lower_pctl = parse_double_option(args, "p", DEFAULT_LOWER_PCTL);
    config.upper_pctl = parse_double_option(args, "q", DEFAULT_UPPER_PCTL);
    config.trim = parse_option(args, "trim");
    config.verbose = parse_option(args, "verbose");

    if (config.lower_pctl < 0 || config.lower_pctl > 100 ||
        config.upper_pctl < 0 || config.upper_pctl > 100) {
        ctools_error(CWINSOR_MODULE, "percentiles must be between 0 and 100");
        return 198;
    }
    if (config.lower_pctl >= config.upper_pctl) {
        ctools_error(CWINSOR_MODULE, "lower percentile must be less than upper");
        return 198;
    }

    /* Parse: nvars nby var_indices... by_indices... */
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

    var_indices = (int *)malloc(nvars * sizeof(int));
    if (!var_indices) {
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }
    if (parse_int_array(var_indices, nvars, &p) != 0) {
        free(var_indices);
        ctools_error(CWINSOR_MODULE, "failed to parse var indices");
        return 198;
    }

    if (nby > 0) {
        by_indices = (int *)malloc(nby * sizeof(int));
        if (!by_indices) {
            free(var_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
        if (parse_int_array(by_indices, nby, &p) != 0) {
            free(var_indices);
            free(by_indices);
            ctools_error(CWINSOR_MODULE, "failed to parse by indices");
            return 198;
        }
    }

    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    size_t nobs = (size_t)(obs2 - obs1 + 1);

    if (nobs == 0) {
        free(var_indices);
        if (by_indices) free(by_indices);
        ctools_error(CWINSOR_MODULE, "no observations");
        return 2000;
    }

    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif

    /* === Load Phase === */
    double load_start = ctools_timer_seconds();

    /* Allocate all data arrays */
    size_t total_cols = nvars + nby;
    var_data = (double **)malloc(nvars * sizeof(double *));
    if (!var_data) {
        free(var_indices);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    for (size_t v = 0; v < nvars; v++) {
        var_data[v] = (double *)malloc(nobs * sizeof(double));
        if (!var_data[v]) {
            for (size_t j = 0; j < v; j++) free(var_data[j]);
            free(var_data);
            free(var_indices);
            if (by_indices) free(by_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }
    }

    if (nby > 0) {
        by_data = (double **)malloc(nby * sizeof(double *));
        if (!by_data) {
            for (size_t v = 0; v < nvars; v++) free(var_data[v]);
            free(var_data);
            free(var_indices);
            free(by_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }

        for (size_t b = 0; b < nby; b++) {
            by_data[b] = (double *)malloc(nobs * sizeof(double));
            if (!by_data[b]) {
                for (size_t j = 0; j < b; j++) free(by_data[j]);
                free(by_data);
                for (size_t v = 0; v < nvars; v++) free(var_data[v]);
                free(var_data);
                free(var_indices);
                free(by_indices);
                ctools_error_alloc(CWINSOR_MODULE);
                return 920;
            }
        }
    }

    /* Load target variables */
    for (size_t v = 0; v < nvars; v++) {
        int idx = var_indices[v];
        double *data = var_data[v];
        for (size_t i = 0; i < nobs; i++) {
            double val;
            if (SF_vdata(idx, (ST_int)(obs1 + i), &val) != 0) val = SV_missval;
            data[i] = val;
        }
    }

    /* Load by-variables */
    if (nby > 0) {
        for (size_t b = 0; b < nby; b++) {
            int idx = by_indices[b];
            double *data = by_data[b];

            for (size_t i = 0; i < nobs; i++) {
                double val;
                if (SF_vdata(idx, (ST_int)(obs1 + i), &val) != 0) val = SV_missval;
                data[i] = val;
            }
        }
    }

    t_load = ctools_timer_seconds() - load_start;

    /* === Sort Phase (if by-variables specified) === */
    double sort_start = ctools_timer_seconds();

    if (nby > 0) {
        /* Allocate sort indices and temp buffer */
        sort_indices = (size_t *)malloc(nobs * sizeof(size_t));
        temp_buffer = (double *)malloc(nobs * sizeof(double));
        inv_perm = (size_t *)malloc(nobs * sizeof(size_t));

        if (!sort_indices || !temp_buffer || !inv_perm) {
            if (sort_indices) free(sort_indices);
            if (temp_buffer) free(temp_buffer);
            if (inv_perm) free(inv_perm);
            if (by_data) {
                for (size_t b = 0; b < nby; b++) free(by_data[b]);
                free(by_data);
            }
            for (size_t v = 0; v < nvars; v++) free(var_data[v]);
            free(var_data);
            free(var_indices);
            free(by_indices);
            ctools_error_alloc(CWINSOR_MODULE);
            return 920;
        }

        /* Get sorted indices */
        argsort_by_vars(by_data, nobs, nby, sort_indices);

        /* Compute inverse permutation for later un-sorting */
        for (size_t i = 0; i < nobs; i++) {
            inv_perm[sort_indices[i]] = i;
        }

        /* Apply permutation to all arrays */
        for (size_t v = 0; v < nvars; v++) {
            apply_permutation(var_data[v], sort_indices, nobs, temp_buffer);
        }
        for (size_t b = 0; b < nby; b++) {
            apply_permutation(by_data[b], sort_indices, nobs, temp_buffer);
        }

    }

    t_sort = ctools_timer_seconds() - sort_start;

    /* === Group Detection Phase === */
    double groups_start = ctools_timer_seconds();

    size_t ngroups = 1;
    groups = (group_info *)malloc((nobs + 1) * sizeof(group_info));
    if (!groups) {
        if (sort_indices) free(sort_indices);
        if (temp_buffer) free(temp_buffer);
        if (inv_perm) free(inv_perm);
        if (by_data) {
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
        }
        for (size_t v = 0; v < nvars; v++) free(var_data[v]);
        free(var_data);
        free(var_indices);
        if (by_indices) free(by_indices);
        ctools_error_alloc(CWINSOR_MODULE);
        return 920;
    }

    if (nby > 0) {
        detect_groups(by_data, nobs, nby, groups, &ngroups);
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
        if (inv_perm) free(inv_perm);
        if (by_data) {
            for (size_t b = 0; b < nby; b++) free(by_data[b]);
            free(by_data);
        }
        for (size_t v = 0; v < nvars; v++) free(var_data[v]);
        free(var_data);
        free(var_indices);
        if (by_indices) free(by_indices);
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
            if (inv_perm) free(inv_perm);
            if (by_data) {
                for (size_t b = 0; b < nby; b++) free(by_data[b]);
                free(by_data);
            }
            for (size_t v = 0; v < nvars; v++) free(var_data[v]);
            free(var_data);
            free(var_indices);
            if (by_indices) free(by_indices);
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
        #pragma omp parallel for schedule(dynamic, 4)
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

    /* === Store Phase === */
    double store_start = ctools_timer_seconds();

    if (nby > 0 && inv_perm != NULL) {
        /* Store in original order using inverse permutation */
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t v = 0; v < nvars; v++) {
            int idx = var_indices[v];
            double *data = var_data[v];
            for (size_t i = 0; i < nobs; i++) {
                /* inv_perm[i] gives the position in sorted array for original obs i */
                SF_vstore(idx, (ST_int)(obs1 + i), data[inv_perm[i]]);
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (size_t v = 0; v < nvars; v++) {
            int idx = var_indices[v];
            double *data = var_data[v];
            for (size_t i = 0; i < nobs; i++) {
                SF_vstore(idx, (ST_int)(obs1 + i), data[i]);
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
    if (inv_perm) free(inv_perm);
    if (by_data) {
        for (size_t b = 0; b < nby; b++) free(by_data[b]);
        free(by_data);
    }
    for (size_t v = 0; v < nvars; v++) free(var_data[v]);
    free(var_data);
    free(var_indices);
    if (by_indices) free(by_indices);

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

    if (config.verbose) {
        ctools_msg(CWINSOR_MODULE, "%zu vars, %zu obs, %zu groups, %d threads",
                   nvars, nobs, ngroups, num_threads);
        ctools_msg(CWINSOR_MODULE, "load=%.4f sort=%.4f groups=%.4f winsor=%.4f store=%.4f total=%.4f",
                   t_load, t_sort, t_groups, t_winsor, t_store, total);
    }

    return 0;
}
