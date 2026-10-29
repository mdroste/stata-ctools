/*
 * cbinscatter_impl.c
 *
 * Main implementation for cbinscatter module
 * Part of the ctools Stata plugin suite
 *
 * Performance notes:
 * - Direct SF_vdata loading (no extra copies)
 * - OpenMP only for embarrassingly parallel operations
 * - No VLA array reductions (too slow)
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
#include "ctools_runtime.h"
#include "../ctools_config.h"
#include "../ctools_hdfe_utils.h"
#include "hdfe/hdfe_remap.h"
#include "cbinscatter_impl.h"
#include "cbinscatter_types.h"
#include "cbinscatter_bins.h"
#include "cbinscatter_binsreg.h"
#include "cbinscatter_resid.h"
#include "cbinscatter_fit.h"

/* Use SF_is_missing() from stplugin.h for missing value checks */

/* ========================================================================
 * Type Initialization and Cleanup
 * ======================================================================== */

static void cbinscatter_config_init(BinscatterConfig *config) {
    config->nquantiles = 20;
    config->linetype = 0;
    config->compute_se = 0;
    config->discrete = 0;
    config->has_controls = 0;
    config->num_controls = 0;
    config->has_absorb = 0;
    config->num_absorb = 0;
    config->has_by = 0;
    config->has_weights = 0;
    config->weight_type = 0;
    config->verbose = 0;
    config->reportreg = 0;
    config->method = 0;  /* 0=classic, 1=binsreg */
    config->maxiter = 10000;
    config->tolerance = 1e-8;
    config->touse_mode = 0;
    config->check_weights = 0;
}

static void cbinscatter_results_init(BinscatterResults *results) {
    results->num_by_groups = 0;
    results->groups = NULL;
    results->nquantiles = 0;
    results->total_obs = 0;
    results->obs_dropped = 0;
}

static void cbinscatter_results_free(BinscatterResults *results) {
    if (results->groups) {
        for (ST_int g = 0; g < results->num_by_groups; g++) {
            if (results->groups[g].bins) {
                free(results->groups[g].bins);
            }
            if (results->groups[g].fit_coefs) {
                free(results->groups[g].fit_coefs);
            }
        }
        free(results->groups);
    }
    cbinscatter_results_init(results);
}

/* ========================================================================
 * Read Configuration from Stata Scalars
 * ======================================================================== */

static ST_retcode read_config(BinscatterConfig *config) {
    ST_double val;

    cbinscatter_config_init(config);

    if (SF_scal_use("__cbinscatter_nquantiles", &val) == 0) {
        config->nquantiles = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_linetype", &val) == 0) {
        config->linetype = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_compute_se", &val) == 0) {
        config->compute_se = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_discrete", &val) == 0) {
        config->discrete = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_has_controls", &val) == 0) {
        config->has_controls = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_num_controls", &val) == 0) {
        config->num_controls = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_has_absorb", &val) == 0) {
        config->has_absorb = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_num_absorb", &val) == 0) {
        config->num_absorb = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_has_by", &val) == 0) {
        config->has_by = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_has_weights", &val) == 0) {
        config->has_weights = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_weight_type", &val) == 0) {
        config->weight_type = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_verbose", &val) == 0) {
        config->verbose = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_reportreg", &val) == 0) {
        config->reportreg = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_method", &val) == 0) {
        config->method = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_maxiter", &val) == 0) {
        config->maxiter = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_tolerance", &val) == 0) {
        config->tolerance = val;
    }
    if (SF_scal_use("__cbinscatter_touse_mode", &val) == 0) {
        config->touse_mode = (ST_int)val;
    }
    if (SF_scal_use("__cbinscatter_check_weights", &val) == 0) {
        config->check_weights = (ST_int)val;
    }

    return CBINSCATTER_OK;
}

/* ========================================================================
 * Load Data from Stata using ctools_data_load()
 * ======================================================================== */

static ST_retcode load_data(
    BinscatterConfig *config,
    ST_double **y_out,
    ST_double **x_out,
    ST_double **controls_out,
    ST_int **fe_vars_out,
    ST_int **by_groups_out,
    ST_double **weights_out,
    ST_int *N_out,
    ST_int *N_valid_out
) {
    ST_int i, j;
    ST_double val;

    ST_double *y = NULL;
    ST_double *x = NULL;
    ST_double *controls = NULL;
    ST_int *fe_vars = NULL;
    ST_double *fe_codes = NULL;
    ST_int *by_groups = NULL;
    ST_double *weights = NULL;
    ST_int *valid_idx = NULL;
    ST_int N_valid = 0;

    /* Count total variables to load via ctools infrastructure
     * Variables are ordered: y, x, controls..., absorb..., by, weights */
    int num_vars_to_load = 2;  /* y and x */
    num_vars_to_load += config->has_controls ? config->num_controls : 0;
    num_vars_to_load += config->has_absorb ? config->num_absorb : 0;
    num_vars_to_load += config->has_by ? 1 : 0;
    num_vars_to_load += config->has_weights ? 1 : 0;
    num_vars_to_load += (config->touse_mode == 1) ? 1 : 0;

    /* Build array of variable indices for all variables */
    int *var_indices = (int *)ctools_safe_malloc2(num_vars_to_load, sizeof(int));
    if (!var_indices) {
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Track indices into loaded_data for each variable type */
    int vidx = 0;
    int y_data_idx = vidx;
    var_indices[vidx++] = 1;  /* y */
    int x_data_idx = vidx;
    var_indices[vidx++] = 2;  /* x */

    int controls_data_idx = vidx;  /* Track where controls start in loaded data */
    for (j = 0; j < (config->has_controls ? config->num_controls : 0); j++) {
        var_indices[vidx++] = 3 + j;  /* controls */
    }

    int absorb_data_idx = vidx;  /* Track where absorb starts in loaded data */
    int absorb_var_start = 3 + (config->has_controls ? config->num_controls : 0);
    for (j = 0; j < (config->has_absorb ? config->num_absorb : 0); j++) {
        var_indices[vidx++] = absorb_var_start + j;  /* absorb vars */
    }

    int by_data_idx = vidx;  /* Track where by-var is in loaded data */
    if (config->has_by) {
        int by_var = absorb_var_start + (config->has_absorb ? config->num_absorb : 0);
        var_indices[vidx++] = by_var;  /* by variable */
    }

    int weights_data_idx = -1;
    if (config->has_weights) {
        /* weights is the last variable */
        weights_data_idx = vidx;
        int w_var = 3;
        if (config->has_controls) w_var += config->num_controls;
        if (config->has_absorb) w_var += config->num_absorb;
        if (config->has_by) w_var += 1;
        var_indices[vidx++] = w_var;
    }

    /* The sample marker follows every other variable */
    int touse_data_idx = -1;
    if (config->touse_mode == 1) {
        touse_data_idx = vidx;
        var_indices[vidx] = vidx + 1;
        vidx++;
    }

    /* Load all variables. With a marker variable (or no if/in) the rows are
     * selected below rather than by SF_ifobs(), which costs one serial SPI
     * call per observation. */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    int load_flags = CTOOLS_LOAD_NO_SORT_ORDER |
                     (config->touse_mode > 0 ? CTOOLS_LOAD_SKIP_IF : CTOOLS_LOAD_CHECK_IF);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices,
                                                       num_vars_to_load, 0, 0, load_flags);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        return ctools_stata_rc(load_rc);
    }

    ST_int N = (ST_int)filtered.data.nobs;  /* Number of if/in filtered observations */

    if (N == 0) {
        ctools_filtered_data_free(&filtered);
        if (config->touse_mode > 0) {
            SF_error("no observations\n");
            SF_scal_save("__cbinscatter_reported", 1.0);
            return 3301;
        }
        *y_out = NULL;
        *x_out = NULL;
        *controls_out = NULL;
        *fe_vars_out = NULL;
        *by_groups_out = NULL;
        *weights_out = NULL;
        *N_out = 0;
        *N_valid_out = 0;
        return CBINSCATTER_OK;
    }

    /* Allocate and copy y and x */
    y = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    x = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    if (!y || !x) {
        free(y); free(x);
        ctools_filtered_data_free(&filtered);
        return CBINSCATTER_ERR_MEMORY;
    }
    memcpy(y, filtered.data.vars[y_data_idx].data.dbl, N * sizeof(ST_double));
    memcpy(x, filtered.data.vars[x_data_idx].data.dbl, N * sizeof(ST_double));

    /* Handle controls - need to reorganize to column-major layout expected by cbinscatter */
    if (config->has_controls && config->num_controls > 0) {
        controls = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)config->num_controls, sizeof(ST_double));
        if (!controls) {
            free(y); free(x);
            ctools_filtered_data_free(&filtered);
            return CBINSCATTER_ERR_MEMORY;
        }
        /* Copy from filtered data (contiguous per var) to column-major layout */
        for (j = 0; j < config->num_controls; j++) {
            ST_double *src = filtered.data.vars[controls_data_idx + j].data.dbl;
            ST_double *dst = &controls[(size_t)j * N];
            memcpy(dst, src, N * sizeof(ST_double));
        }
    }

    /* Handle absorb variables. FE values are only labels: keep them as numeric
     * codes (string variables become group codes) until the estimation sample
     * is known, then map each variable to dense levels 1..L below. */
    if (config->has_absorb && config->num_absorb > 0) {
        fe_vars = (ST_int *)ctools_safe_malloc3((size_t)N, (size_t)config->num_absorb, sizeof(ST_int));
        fe_codes = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)config->num_absorb, sizeof(ST_double));
        if (!fe_vars || !fe_codes) {
            free(y); free(x); free(controls); free(fe_vars); free(fe_codes);
            ctools_filtered_data_free(&filtered);
            return CBINSCATTER_ERR_MEMORY;
        }
        for (j = 0; j < config->num_absorb; j++) {
            stata_variable *var = &filtered.data.vars[absorb_data_idx + j];
            ST_double *dst = &fe_codes[(size_t)j * N];
            if (var->type == STATA_TYPE_STRING) {
                /* Group ids 0..L-1 ("" = missing -> -1), fe_vars as scratch */
                ST_int *ids = &fe_vars[(size_t)j * N], n_groups = 0;
                if (ctools_strings_to_cluster_ids(var->data.str, N, ids, &n_groups) != 0) {
                    free(y); free(x); free(controls); free(fe_vars); free(fe_codes);
                    ctools_filtered_data_free(&filtered);
                    return CBINSCATTER_ERR_MEMORY;
                }
                for (i = 0; i < N; i++) dst[i] = ids[i] < 0 ? SV_missval : (ST_double)ids[i];
            } else {
                memcpy(dst, var->data.dbl, (size_t)N * sizeof(ST_double));
            }
        }
    }

    /* Handle by variable - convert double to integer with missing handling */
    if (config->has_by) {
        by_groups = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
        if (!by_groups) {
            free(y); free(x); free(controls); free(fe_vars); free(fe_codes);
            ctools_filtered_data_free(&filtered);
            return CBINSCATTER_ERR_MEMORY;
        }
        ST_double *src = filtered.data.vars[by_data_idx].data.dbl;
        for (i = 0; i < N; i++) {
            val = src[i];
            /* Set missings to -1 so they're filtered by the < 0 check later */
            by_groups[i] = SF_is_missing(val) ? -1 : (ST_int)val;
        }
    }

    /* Get weights pointer */
    if (config->has_weights && weights_data_idx >= 0) {
        weights = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
        if (!weights) {
            free(y); free(x); free(controls); free(fe_vars); free(fe_codes); free(by_groups);
            ctools_filtered_data_free(&filtered);
            return CBINSCATTER_ERR_MEMORY;
        }
        memcpy(weights, filtered.data.vars[weights_data_idx].data.dbl, N * sizeof(ST_double));
    }

    /* Allocate valid_idx array */
    valid_idx = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
    if (!valid_idx) {
        free(y); free(x); free(controls); free(fe_vars); free(fe_codes); free(by_groups); free(weights);
        ctools_filtered_data_free(&filtered);
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Build valid index array - check for missing values in variables.
     * The requested sample is the loaded rows (if condition), the rows
     * marked by the marker variable, or every row with a usable weight. */
    const ST_double *touse = touse_data_idx >= 0 ? filtered.data.vars[touse_data_idx].data.dbl : NULL;
    int64_t N_requested = 0;
    int negative_weight = 0, fractional_fweight = 0;
    N_valid = 0;
    {
        /* Flags first, one column at a time over blocks of rows (sequential
         * reads; blocks in parallel), then the index list in row order. A
         * value is missing when it is not below Stata's missing code, which
         * for stored data is what SF_is_missing() tests. */
        const ST_double missval = SV_missval;
        const ST_int block = 16384;
        const ST_int num_blocks = N / block + (N % block != 0); /* no N + block overflow */
        ST_int blk;
        #pragma omp parallel for schedule(static) private(j) \
            reduction(| : negative_weight, fractional_fweight) reduction(+ : N_requested) if (N > 100000)
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_int i0 = blk * block;
            const ST_int i1 = (N - i0 < block) ? N : i0 + block;
            ST_int r;
            /* y and x */
            for (r = i0; r < i1; r++) valid_idx[r] = (y[r] < missval) && (x[r] < missval);
            if (controls) {
                for (j = 0; j < config->num_controls; j++) {
                    const ST_double *col = controls + (size_t)j * N;
                    for (r = i0; r < i1; r++) if (!(col[r] < missval)) valid_idx[r] = 0;
                }
            }
            /* absorb variables (any non-missing code is a valid level) */
            if (fe_codes) {
                for (j = 0; j < config->num_absorb; j++) {
                    const ST_double *col = fe_codes + (size_t)j * N;
                    for (r = i0; r < i1; r++) if (!(col[r] < missval)) valid_idx[r] = 0;
                }
            }
            if (by_groups) {
                for (r = i0; r < i1; r++) if (by_groups[r] < 0) valid_idx[r] = 0;
            }
            if (weights) {
                for (r = i0; r < i1; r++) if (!(weights[r] < missval) || weights[r] <= 0) valid_idx[r] = 0;
            }
            /* Requested rows; their weights must be positive (missing and zero
             * weights are excluded, as marksample does) */
            int64_t requested = 0;
            for (r = i0; r < i1; r++) {
                int req = touse ? (touse[r] != 0 && touse[r] < missval)
                        : (config->touse_mode == 2 && weights) ? (weights[r] < missval && weights[r] != 0)
                        : 1;
                if (!req) valid_idx[r] = 0;
                requested += req;
                if (req && config->check_weights && weights[r] < missval) {
                    if (weights[r] < 0) negative_weight |= 1;
                    if (config->weight_type == 2 && weights[r] != floor(weights[r])) fractional_fweight |= 1;
                }
            }
            N_requested += requested;
        }
        for (i = 0; i < N; i++) {
            if (valid_idx[i]) valid_idx[N_valid++] = i;
        }
    }

    /* Done with filtered data */
    ctools_filtered_data_free(&filtered);

    if (negative_weight || fractional_fweight || N_requested == 0) {
        free(y); free(x); free(controls); free(fe_vars); free(fe_codes); free(by_groups); free(weights);
        free(valid_idx);
        SF_scal_save("__cbinscatter_reported", 1.0);
        if (negative_weight) {
            SF_error("ctools: weights must be positive\n");
            return 198;
        }
        if (fractional_fweight) {
            SF_error("frequency weights must be integers\n");
            return 401;
        }
        SF_error("no observations\n");
        return 3301;
    }

    /* Compact data using valid indices (streaming compaction).
     * IMPORTANT: Multi-column arrays (controls, fe_vars) must be compacted
     * column-by-column to avoid cross-column overwrites. When compacting
     * column j+1 with stride N_valid, the destination range can overlap with
     * column j's source range (stride N), corrupting unread source data. */
    if (N_valid < N && N_valid > 0) {
        /* Compact y, x, by_groups, weights (single-column, always safe) */
        for (ST_int k = 0; k < N_valid; k++) {
            ST_int src = valid_idx[k];
            if (src != k) {
                y[k] = y[src];
                x[k] = x[src];
                if (by_groups) by_groups[k] = by_groups[src];
                if (weights) weights[k] = weights[src];
            }
        }
        /* Compact controls column-by-column.
         * Cannot skip when src == k because destination stride (N_valid) differs
         * from source stride (N), so j*N_valid+k != j*N+k for j > 0. */
        if (controls) {
            for (j = 0; j < config->num_controls; j++) {
                for (ST_int k = 0; k < N_valid; k++) {
                    controls[(size_t)j * N_valid + k] = controls[(size_t)j * N + valid_idx[k]];
                }
            }
        }
        /* Compact FE codes column-by-column */
        if (fe_codes) {
            for (j = 0; j < config->num_absorb; j++) {
                for (ST_int k = 0; k < N_valid; k++) {
                    fe_codes[(size_t)j * N_valid + k] = fe_codes[(size_t)j * N + valid_idx[k]];
                }
            }
        }
    }

    free(valid_idx);

    /* Map each FE variable to dense levels 1..L over the estimation sample
     * (negative, fractional and large codes and strings are all ordinary
     * labels; the residualizers index their level arrays by these). */
    if (fe_codes && N_valid > 0) {
        for (j = 0; j < config->num_absorb; j++) {
            ST_int num_levels = 0, *counts = NULL;
            if (remap_and_count(&fe_codes[(size_t)j * N_valid], N_valid,
                                &fe_vars[(size_t)j * N_valid], &num_levels,
                                &counts, NULL, NULL) != 0) {
                free(y); free(x); free(controls); free(fe_vars); free(fe_codes);
                free(by_groups); free(weights);
                return CBINSCATTER_ERR_MEMORY;
            }
            free(counts);
        }
    }
    free(fe_codes);

    *y_out = y;
    *x_out = x;
    *controls_out = controls;
    *fe_vars_out = fe_vars;
    *by_groups_out = by_groups;
    *weights_out = weights;
    *N_out = (ST_int)N_requested;
    *N_valid_out = N_valid;

    return CBINSCATTER_OK;
}

/* ========================================================================
 * Store Results to Stata
 * ======================================================================== */

static ST_retcode store_results(
    BinscatterResults *results,
    BinscatterConfig *config
) {
    ST_int g, b, row;

    /* Store scalars */
    { ST_retcode rc = SF_scal_save("__cbinscatter_N", (ST_double)results->total_obs); if (rc) return rc; }
    { ST_retcode rc = SF_scal_save("__cbinscatter_N_dropped", (ST_double)results->obs_dropped); if (rc) return rc; }
    { ST_retcode rc = SF_scal_save("__cbinscatter_num_groups", (ST_double)results->num_by_groups); if (rc) return rc; }

    /* Store bin data to matrix __cbinscatter_bins */
    row = 1;
    for (g = 0; g < results->num_by_groups; g++) {
        ByGroupResult *group = &results->groups[g];
        for (b = 0; b < group->num_bins; b++) {
            BinStats *bin = &group->bins[b];
            { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_bins", row, 1, (ST_double)group->by_group_id); if (rc) return rc; }
            { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_bins", row, 2, (ST_double)bin->bin_id); if (rc) return rc; }
            { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_bins", row, 3, bin->x_mean); if (rc) return rc; }
            { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_bins", row, 4, bin->y_mean); if (rc) return rc; }
            { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_bins", row, 5, (ST_double)bin->n_obs); if (rc) return rc; }
            row++;
        }
    }

    /* Store fit coefficients if line fitting was done */
    if (config->linetype > 0) {
        for (g = 0; g < results->num_by_groups; g++) {
            ByGroupResult *group = &results->groups[g];
            if (group->fit_coefs && group->fit_order > 0) {
                for (ST_int j = 0; j <= group->fit_order && j < 4; j++) {
                    { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_coefs", g + 1, j + 1, group->fit_coefs[j]); if (rc) return rc; }
                }
                { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_fit_stats", g + 1, 1, group->fit_r2); if (rc) return rc; }
                { ST_retcode rc = (_stata_)->safematstore("__cbinscatter_fit_stats", g + 1, 2, (ST_double)group->fit_n); if (rc) return rc; }
            }
        }
    }

    return CBINSCATTER_OK;
}

/* ========================================================================
 * Rows in by-group order
 * ======================================================================== */

/*
 * Reorder the N rows of y, x, *w, the K control columns, the F FE columns
 * and by (each column N long) by group: groups 1..G in order, rows with any
 * other id last, row order kept within each group. Blocks of rows are counted
 * and then scattered in parallel. start[g] (g = 0..G) is group g + 1's first
 * row and start[G] the end of the groups.
 */
static ST_retcode partition_by_group(ST_double **y, ST_double **x, ST_double **w,
                                     ST_double **ctrl, ST_int K, ST_int **fe, ST_int F,
                                     ST_int **by, ST_int N, ST_int G, ST_int *start)
{
    const ST_int num_blocks = cbs_num_blocks(N);
    const ST_int nb = G + 1;  /* bucket G: ids outside 1..G */
    ST_int *offsets = (ST_int *)calloc((size_t)num_blocks * nb, sizeof(ST_int));
    ST_double *ny = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    ST_double *nx = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));
    ST_double *nw = *w ? (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double)) : NULL;
    ST_double *nc = (*ctrl && K > 0) ? (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K, sizeof(ST_double)) : NULL;
    ST_int *nf = (*fe && F > 0) ? (ST_int *)ctools_safe_malloc3((size_t)N, (size_t)F, sizeof(ST_int)) : NULL;
    ST_int *nby = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
    if (!offsets || !ny || !nx || (*w && !nw) || (*ctrl && K > 0 && !nc) ||
        (*fe && F > 0 && !nf) || !nby) {
        free(offsets); free(ny); free(nx); free(nw); free(nc); free(nf); free(nby);
        return CBINSCATTER_ERR_MEMORY;
    }
    const ST_int *b_in = *by;
    ST_int blk;

    #pragma omp parallel for schedule(static) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        ST_int *cnt = offsets + (size_t)blk * nb;
        for (ST_int i = i0; i < i1; i++) {
            ST_int b = b_in[i] - 1;
            cnt[(b >= 0 && b < G) ? b : G]++;
        }
    }
    /* Bucket by bucket, blocks in order: each block's first row per bucket */
    ST_int run = 0;
    for (ST_int b = 0; b < nb; b++) {
        start[b] = run;
        for (blk = 0; blk < num_blocks; blk++) {
            ST_int c = offsets[(size_t)blk * nb + b];
            offsets[(size_t)blk * nb + b] = run;
            run += c;
        }
    }

    const ST_double *y_in = *y, *x_in = *x, *w_in = *w, *c_in = *ctrl;
    const ST_int *f_in = *fe;
    #pragma omp parallel for schedule(static) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        ST_int *pos = offsets + (size_t)blk * nb;
        for (ST_int i = i0; i < i1; i++) {
            ST_int b = b_in[i] - 1;
            ST_int at = pos[(b >= 0 && b < G) ? b : G]++;
            ny[at] = y_in[i];
            nx[at] = x_in[i];
            if (nw) nw[at] = w_in[i];
            for (ST_int k = 0; nc && k < K; k++) nc[(size_t)k * N + at] = c_in[(size_t)k * N + i];
            for (ST_int k = 0; nf && k < F; k++) nf[(size_t)k * N + at] = f_in[(size_t)k * N + i];
            nby[at] = b_in[i];
        }
    }

    free(offsets);
    free(*y); *y = ny;
    free(*x); *x = nx;
    if (nw) { free(*w); *w = nw; }
    if (nc) { free(*ctrl); *ctrl = nc; }
    if (nf) { free(*fe); *fe = nf; }
    free(*by); *by = nby;
    return CBINSCATTER_OK;
}

/* ========================================================================
 * Main Computation Function
 * ======================================================================== */

static ST_retcode do_compute_bins(void) {
    BinscatterConfig config;
    BinscatterResults results;
    ST_retcode rc;
    double t_start, t_load, t_resid, t_bins, t_fit, t_total;

    ST_double *y = NULL, *x = NULL;
    ST_double *controls = NULL;
    ST_int *fe_vars = NULL;
    ST_int *by_groups = NULL;
    ST_double *weights = NULL;
    ST_int *group_sizes = NULL;
    ST_int *group_start = NULL;
    ST_double *pooled_cutpoints = NULL;
    ST_int pooled_nq = 0;
    ST_int N, N_valid;

    cbinscatter_results_init(&results);

    t_start = ctools_timer_seconds();

    /* Read configuration */
    rc = read_config(&config);
    if (rc != CBINSCATTER_OK) return rc;

    /* Verbose output handled by .ado file using stored scalars */

    /* Load data */
    rc = load_data(&config, &y, &x, &controls, &fe_vars, &by_groups, &weights,
                   &N, &N_valid);
    if (rc != CBINSCATTER_OK) goto cleanup;

    t_load = ctools_timer_seconds();

    /* The requested sample, before incomplete observations are excluded */
    SF_scal_save("__cbinscatter_N_touse", (ST_double)N);

    if (N_valid < 2) {
        SF_error("cbinscatter: insufficient valid observations\n");
        rc = CBINSCATTER_ERR_NOOBS;
        goto cleanup;
    }

    results.total_obs = N_valid;
    results.obs_dropped = N - N_valid;

    /* Verbose output handled by .ado file using stored scalars */

    /* Save raw means for addback after residualization (classic method) */
    ST_double y_raw_mean = 0.0, x_raw_mean = 0.0;

    /* Residualize if needed */
    t_resid = t_load;
    if (config.has_controls || config.has_absorb) {
        ST_int dropped = 0;

        /* Classic method: save raw means before in-place residualization.
         * binscatter adds sample means back to residuals so plot axes show
         * values in the original scale; we do the same for compatibility. */
        if (config.method == 0) {
            if (weights != NULL) {
                ST_double total_w = 0.0;
                for (ST_int i = 0; i < N_valid; i++) {
                    y_raw_mean += weights[i] * y[i];
                    x_raw_mean += weights[i] * x[i];
                    total_w += weights[i];
                }
                if (total_w > 0) {
                    y_raw_mean /= total_w;
                    x_raw_mean /= total_w;
                }
            } else {
                for (ST_int i = 0; i < N_valid; i++) {
                    y_raw_mean += y[i];
                    x_raw_mean += x[i];
                }
                y_raw_mean /= N_valid;
                x_raw_mean /= N_valid;
            }
        }

        if (config.verbose) {
            if (config.method == 1) {
                SF_display("cbinscatter: using binsreg method (bin on raw X, regression adjustment)\n");
            } else {
                SF_display("cbinscatter: using classic method (residualize both X and Y)\n");
            }
        }

        /*
         * Method 0 (classic): Residualize both Y and X, bin on residualized X
         * Method 1 (binsreg): Bin on raw X, use regression Y ~ bins + controls
         *   (Cattaneo et al. "On Binscatter" approach)
         *
         * For binsreg method:
         *   - Controls only: Skip residualization here; binsreg function handles it
         *   - Absorb only: DO NOT residualize here - handled by adjust_bins_binsreg_hdfe
         *   - Absorb + Controls: DO NOT residualize here - handled by adjust_bins_binsreg_hdfe
         *
         * The binsreg HDFE adjustment uses the FWL theorem correctly by demeaning
         * BOTH Y and bin indicators, which is different from simple HDFE demeaning.
         */
        if (config.method == 1) {
            /* Binsreg method: skip residualization here.
             * For absorb cases, adjust_bins_binsreg_hdfe will handle the FWL approach.
             * For controls-only, adjust_bins_binsreg handles it via regression. */
        } else {
            /* Classic method: Residualize both Y and X */
            if (config.has_absorb && config.has_controls) {
                rc = combined_residualize(y, x, controls, fe_vars, N_valid,
                                          config.num_controls, config.num_absorb,
                                          weights, config.weight_type,
                                          config.maxiter, config.tolerance,
                                          config.verbose, &dropped);
            } else if (config.has_absorb) {
                rc = hdfe_residualize(y, x, fe_vars, N_valid, config.num_absorb,
                                      weights, config.weight_type,
                                      config.maxiter, config.tolerance,
                                      config.verbose, &dropped);
            } else {
                rc = ols_residualize(y, x, controls, N_valid, config.num_controls,
                                     weights, config.weight_type);
            }

            if (rc != CBINSCATTER_OK) {
                if (rc != CBINSCATTER_ERR_HDFE)
                    SF_error("cbinscatter: residualization failed\n");
                goto cleanup;
            }

            /* Like binscatter, add the means back to the residuals themselves,
             * so bins, bin means and the fit line share the original scale */
            for (ST_int i = 0; i < N_valid; i++) {
                y[i] += y_raw_mean;
                x[i] += x_raw_mean;
            }
        }

        results.obs_dropped += dropped;
        t_resid = ctools_timer_seconds();

        if (config.verbose) {
            char msg[256];
            snprintf(msg, sizeof(msg), "cbinscatter: residualization complete (%.3f sec)\n",
                     t_resid - t_load);
            SF_display(msg);
        }
    }

    /* Count by-groups */
    ST_int num_groups = 1;
    ST_int min_group = 0;
    if (config.has_by && by_groups != NULL) {
        ST_int max_group = by_groups[0];
        min_group = by_groups[0];
        for (ST_int i = 1; i < N_valid; i++) {
            if (by_groups[i] > max_group) max_group = by_groups[i];
            if (by_groups[i] < min_group) min_group = by_groups[i];
        }
        num_groups = max_group - min_group + 1;

        /* Remap by_groups to be 1-based if they start from 0 */
        if (min_group == 0) {
            for (ST_int i = 0; i < N_valid; i++) {
                by_groups[i] += 1;
            }
        }
    }

    results.num_by_groups = num_groups;
    results.nquantiles = config.nquantiles;
    results.groups = (ByGroupResult *)calloc(num_groups, sizeof(ByGroupResult));
    if (!results.groups) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Pre-compute group sizes in one O(N) pass instead of O(N*G) */
    group_sizes = (ST_int *)calloc(num_groups, sizeof(ST_int));
    if (!group_sizes) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }
    if (config.has_by) {
        for (ST_int i = 0; i < N_valid; i++) {
            ST_int g_id = by_groups[i] - 1;
            if (g_id >= 0 && g_id < num_groups) {
                group_sizes[g_id]++;
            }
        }
    } else {
        group_sizes[0] = N_valid;
    }

    /* With by(), each group becomes one contiguous range of rows */
    if (config.has_by) {
        group_start = (ST_int *)calloc((size_t)num_groups + 1, sizeof(ST_int));
        if (!group_start) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }
        rc = partition_by_group(&y, &x, &weights, &controls, config.num_controls,
                                &fe_vars, config.num_absorb, &by_groups, N_valid,
                                num_groups, group_start);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* Binsreg method with controls/absorb fits each group's line jointly with
     * them (inside the loop); otherwise lines are fit after binning */
    int binsreg_joint_fit = (config.method == 1 &&
                             ((config.has_controls && config.num_controls > 0) ||
                              config.has_absorb));
    double t_group_fit = 0.0;

    /* binscatter's by() bins every group on the quantiles of the pooled
     * (residualized) x; binsreg's method bins each group separately */
    if (config.has_by && config.method == 0 && !config.discrete) {
        rc = compute_pooled_cutpoints(x, weights, N_valid, config.nquantiles,
                                      &pooled_cutpoints, &pooled_nq);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* store_results writes every group's bins into __cbinscatter_bins, which
     * the ado sized (500 rows per group for discrete). Bins that cannot be
     * stored fail here, before the binsreg adjustment does N x bins work. */
    const ST_int bin_capacity = SF_row("__cbinscatter_bins");
    ST_int bins_used = 0;

    /* Process each by-group */
    for (ST_int g = 0; g < num_groups; g++) {
        ByGroupResult *group = &results.groups[g];
        group->by_group_id = g + 1;

        ST_int n_group = group_sizes[g];

        if (n_group == 1 && config.method == 0 && config.linetype > 0) {
            rc = 2001; /* binscatter's group regression needs two observations */
            goto cleanup;
        }
        if (n_group == 0) {
            group->num_bins = 0;
            group->bins = NULL;
            continue;
        }

        /* The group's rows (contiguous after partition_by_group) */
        const ST_int first = config.has_by ? group_start[g] : 0;
        ST_double *y_group = y + first;
        ST_double *x_group = x + first;
        ST_double *w_group = weights ? weights + first : NULL;
        ST_double *c_group = NULL;  /* Controls for this group */
        ST_int *fe_group = NULL;    /* FE vars for this group */

        /* For binsreg method with controls or absorb, we need group controls/FE */
        int use_binsreg_ctrl = (config.method == 1 && config.has_controls && config.num_controls > 0);
        int use_binsreg_hdfe = (config.method == 1 && config.has_absorb);

        if (use_binsreg_ctrl || use_binsreg_hdfe) {
            if (config.has_controls && config.num_controls > 0) {
                c_group = (ST_double *)ctools_safe_malloc3((size_t)n_group, (size_t)config.num_controls, sizeof(ST_double));
            }
        }
        if (use_binsreg_hdfe && config.num_absorb > 0) {
            fe_group = (ST_int *)ctools_safe_malloc3((size_t)n_group, (size_t)config.num_absorb, sizeof(ST_int));
        }

        if ((use_binsreg_ctrl && config.has_controls && !c_group) ||
            (use_binsreg_hdfe && config.num_absorb > 0 && !fe_group)) {
            free(c_group); free(fe_group);
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        /* The binsreg adjustments take the group's controls and FE columns
         * contiguously */
        if (c_group) {
            for (ST_int k = 0; k < config.num_controls; k++) {
                memcpy(&c_group[(size_t)k * n_group], &controls[(size_t)k * N_valid + first],
                       (size_t)n_group * sizeof(ST_double));
            }
        }
        if (fe_group) {
            for (ST_int k = 0; k < config.num_absorb; k++) {
                memcpy(&fe_group[(size_t)k * n_group], &fe_vars[(size_t)k * N_valid + first],
                       (size_t)n_group * sizeof(ST_int));
            }
        }

        /* The binsreg adjustments need each observation's bin */
        ST_int *bin_ids = NULL;
        if (use_binsreg_ctrl || use_binsreg_hdfe) {
            bin_ids = (ST_int *)malloc((size_t)n_group * sizeof(ST_int));
            if (!bin_ids) {
                free(c_group); free(fe_group);
                rc = CBINSCATTER_ERR_MEMORY;
                goto cleanup;
            }
        }

        /* Compute bins (bin assignments and X means) */
        if (config.discrete) {
            rc = compute_bins_discrete(y_group, x_group, w_group, n_group,
                                       config.compute_se, group, bin_ids);
        } else {
            rc = compute_bins_single_group(y_group, x_group, w_group, n_group,
                                           &config, pooled_cutpoints, pooled_nq,
                                           group, bin_ids);
        }

        if (rc == CBINSCATTER_OK) {
            if (group->num_bins > bin_capacity - bins_used) {
                SF_error("cbinscatter: x has too many distinct values for option discrete "
                         "(at most 500 per group)\n");
                rc = 503; /* as the store of a row past the matrix would fail */
            } else {
                bins_used += group->num_bins;
            }
        }

        /* Binsreg method: regression-adjust the Y-means of the bins assigned
         * above, discrete or quantile (renumbered to index group->bins) */
        if (rc == CBINSCATTER_OK && bin_ids && group->bins != NULL && group->num_bins > 0) {
            rc = densify_bin_ids(bin_ids, n_group, group);
            if (rc == CBINSCATTER_OK && use_binsreg_hdfe) {
                /* FWL-based HDFE adjustment (controls, if any, enter jointly) */
                rc = adjust_bins_binsreg_hdfe(y_group, c_group, fe_group, bin_ids, w_group,
                                               n_group,
                                               config.has_controls ? config.num_controls : 0,
                                               config.num_absorb,
                                               config.weight_type,
                                               config.maxiter, config.tolerance,
                                               group);
            } else if (rc == CBINSCATTER_OK) {
                rc = adjust_bins_binsreg(y_group, c_group, bin_ids, w_group, n_group,
                                         config.num_controls, group);
            }
        }
        free(bin_ids);

        /* Binsreg method: fit the line jointly with the controls/FE and
         * evaluate it at the controls' means, as binsreg's polyreg() does */
        if (rc == CBINSCATTER_OK && binsreg_joint_fit && config.linetype > 0) {
            double t_fit_start = ctools_timer_seconds();
            rc = fit_group_polynomial(y_group, x_group, w_group, n_group,
                                      c_group, c_group ? config.num_controls : 0,
                                      fe_group, fe_group ? config.num_absorb : 0,
                                      &config, group);
            t_group_fit += ctools_timer_seconds() - t_fit_start;
            /* As in fit_all_groups: a group without a line is not an error */
            if (rc != CBINSCATTER_ERR_MEMORY && rc != CBINSCATTER_ERR_HDFE && num_groups > 1)
                rc = CBINSCATTER_OK;
        }

        free(c_group);
        free(fe_group);

        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    free(group_sizes);
    group_sizes = NULL;

    t_bins = ctools_timer_seconds();

    /* Verbose output handled by .ado file using stored scalars */

    /* Fit lines if requested (on the mean-restored residuals for the classic
     * method with controls/absorb, as binscatter does) */
    t_fit = t_bins;
    if (config.linetype > 0 && !binsreg_joint_fit) {
        rc = fit_all_groups(y, x, group_start, weights, N_valid, &config, &results);
        if (rc != CBINSCATTER_OK) goto cleanup;

        t_fit = ctools_timer_seconds();

        /* Verbose output handled by .ado file using stored scalars */
    }

    /* Store results to Stata */
    rc = store_results(&results, &config);

    t_total = ctools_timer_seconds() - t_start;

    /* Store timing info */
    SF_scal_save("__cbinscatter_time_load", t_load - t_start);
    SF_scal_save("__cbinscatter_time_resid", t_resid - t_load);
    SF_scal_save("__cbinscatter_time_bins", t_bins - t_resid - t_group_fit);
    SF_scal_save("__cbinscatter_time_fit", t_fit - t_bins + t_group_fit);
    SF_scal_save("__cbinscatter_time_total", t_total);
    CTOOLS_SAVE_THREAD_INFO("__cbinscatter");

    /* Verbose output handled by .ado file using stored scalars */

cleanup:
    if (rc == CBINSCATTER_ERR_HDFE)
        SF_error("cbinscatter: the absorbed effects did not converge\n");
    free(y);
    free(x);
    free(controls);
    free(fe_vars);
    free(by_groups);
    free(weights);
    free(group_sizes);
    free(group_start);
    free(pooled_cutpoints);
    cbinscatter_results_free(&results);

    return rc;
}

/* ========================================================================
 * Main Entry Point
 * ======================================================================== */

ST_retcode cbinscatter_main(const char *args) {
    if (args == NULL || strlen(args) == 0) {
        SF_error("cbinscatter: no subcommand specified\n");
        return CBINSCATTER_ERR_SYNTAX;
    }

    if (strcmp(args, "compute_bins") == 0) {
        return do_compute_bins();
    }

    char msg[256];
    snprintf(msg, sizeof(msg), "cbinscatter: unknown subcommand '%s'\n", args);
    SF_error(msg);
    return CBINSCATTER_ERR_SYNTAX;
}
