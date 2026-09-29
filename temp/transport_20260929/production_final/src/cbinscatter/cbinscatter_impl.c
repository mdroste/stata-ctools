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
#include "../creghdfe/creghdfe_utils.h"
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

    /* Load all variables with if/in filtering using ctools_data_load().
     * This eliminates the need to call SF_ifobs() during the validation loop. */
    ctools_filtered_data filtered;
    ctools_filtered_data_init(&filtered);
    stata_retcode load_rc = ctools_data_load(&filtered, var_indices,
                                                       num_vars_to_load, 0, 0, CTOOLS_LOAD_CHECK_IF | CTOOLS_LOAD_NO_SORT_ORDER);
    free(var_indices);

    if (load_rc != STATA_OK) {
        ctools_filtered_data_free(&filtered);
        return ctools_stata_rc(load_rc);
    }

    ST_int N = (ST_int)filtered.data.nobs;  /* Number of if/in filtered observations */

    if (N == 0) {
        ctools_filtered_data_free(&filtered);
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

    /* Done with filtered data */
    ctools_filtered_data_free(&filtered);

    /* Allocate valid_idx array */
    valid_idx = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
    if (!valid_idx) {
        free(y); free(x); free(controls); free(fe_vars); free(fe_codes); free(by_groups); free(weights);
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Build valid index array - check for missing values in variables.
     * Note: SF_ifobs is already handled by ctools_data_load(). */
    N_valid = 0;
    for (i = 0; i < N; i++) {
        int valid = 1;

        /* Check y and x for missing values */
        if (SF_is_missing(y[i]) || SF_is_missing(x[i])) {
            valid = 0;
        }

        /* Check controls */
        if (valid && controls) {
            for (j = 0; j < config->num_controls; j++) {
                if (SF_is_missing(controls[(size_t)j * N + i])) {
                    valid = 0;
                    break;
                }
            }
        }

        /* Check absorb variables (any non-missing code is a valid level) */
        if (valid && fe_codes) {
            for (j = 0; j < config->num_absorb; j++) {
                if (SF_is_missing(fe_codes[(size_t)j * N + i])) {
                    valid = 0;
                    break;
                }
            }
        }

        /* Check by variable */
        if (valid && by_groups) {
            if (by_groups[i] < 0) {
                valid = 0;
            }
        }

        /* Check weights */
        if (valid && weights) {
            if (SF_is_missing(weights[i]) || weights[i] <= 0) {
                valid = 0;
            }
        }

        if (valid) {
            valid_idx[N_valid++] = i;
        }
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
    *N_out = N;
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

    /* Binsreg method with controls/absorb fits each group's line jointly with
     * them (inside the loop); otherwise lines are fit after binning */
    int binsreg_joint_fit = (config.method == 1 &&
                             ((config.has_controls && config.num_controls > 0) ||
                              config.has_absorb));
    double t_group_fit = 0.0;

    /* Process each by-group */
    for (ST_int g = 0; g < num_groups; g++) {
        ByGroupResult *group = &results.groups[g];
        group->by_group_id = g + 1;

        ST_int n_group = group_sizes[g];

        if (n_group < 2) {
            group->num_bins = 0;
            group->bins = NULL;
            continue;
        }

        /* Allocate and extract group data */
        ST_double *y_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        ST_double *x_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        ST_double *w_group = NULL;
        ST_double *c_group = NULL;  /* Controls for this group */
        ST_int *fe_group = NULL;    /* FE vars for this group */

        if (weights) {
            w_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        }

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

        if (!y_group || !x_group || (weights && !w_group) ||
            (use_binsreg_ctrl && config.has_controls && !c_group) ||
            (use_binsreg_hdfe && config.num_absorb > 0 && !fe_group)) {
            free(y_group); free(x_group); free(w_group); free(c_group); free(fe_group);
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        ST_int j = 0;
        for (ST_int i = 0; i < N_valid; i++) {
            if (!config.has_by || by_groups[i] == g + 1) {
                y_group[j] = y[i];
                x_group[j] = x[i];
                if (w_group) w_group[j] = weights[i];
                if (c_group) {
                    for (ST_int k = 0; k < config.num_controls; k++) {
                        c_group[(size_t)k * n_group + j] = controls[(size_t)k * N_valid + i];
                    }
                }
                if (fe_group) {
                    for (ST_int k = 0; k < config.num_absorb; k++) {
                        fe_group[(size_t)k * n_group + j] = fe_vars[(size_t)k * N_valid + i];
                    }
                }
                j++;
            }
        }

        /* The binsreg adjustments need each observation's bin */
        ST_int *bin_ids = NULL;
        if (use_binsreg_ctrl || use_binsreg_hdfe) {
            bin_ids = (ST_int *)malloc((size_t)n_group * sizeof(ST_int));
            if (!bin_ids) {
                free(y_group); free(x_group); free(w_group); free(c_group); free(fe_group);
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
                                           &config, group, bin_ids);
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
            if (rc != CBINSCATTER_ERR_MEMORY && num_groups > 1) rc = CBINSCATTER_OK;
        }

        free(y_group);
        free(x_group);
        free(w_group);
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
        rc = fit_all_groups(y, x, by_groups, weights, N_valid, &config, &results);
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
    free(y);
    free(x);
    free(controls);
    free(fe_vars);
    free(by_groups);
    free(weights);
    free(group_sizes);
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
