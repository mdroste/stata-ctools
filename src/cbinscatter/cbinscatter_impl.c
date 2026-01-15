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
#include "ctools_timer.h"
#include "cbinscatter_impl.h"
#include "cbinscatter_types.h"
#include "cbinscatter_bins.h"
#include "cbinscatter_resid.h"
#include "cbinscatter_fit.h"

/* Stata missing value check */
#define STATA_MISSING 8.988465674311579e+307
#define IS_MISSING(x) ((x) >= STATA_MISSING)

/* ========================================================================
 * Type Initialization and Cleanup
 * ======================================================================== */

void cbinscatter_config_init(BinscatterConfig *config) {
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

void cbinscatter_results_init(BinscatterResults *results) {
    results->num_by_groups = 0;
    results->groups = NULL;
    results->nquantiles = 0;
    results->total_obs = 0;
    results->obs_dropped = 0;
}

void cbinscatter_results_free(BinscatterResults *results) {
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

void cbinscatter_workdata_init(BinscatterWorkData *work) {
    work->N = 0;
    work->y = NULL;
    work->x = NULL;
    work->weights = NULL;
    work->by_groups = NULL;
    work->sort_idx = NULL;
    work->y_owned = 0;
    work->x_owned = 0;
}

void cbinscatter_workdata_free(BinscatterWorkData *work) {
    if (work->y_owned && work->y) free(work->y);
    if (work->x_owned && work->x) free(work->x);
    cbinscatter_workdata_init(work);
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
 * Load Data from Stata (Optimized with 8x unrolling)
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
    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();
    ST_int N = in2 - in1 + 1;
    ST_int i, j, var_idx;
    ST_double val;

    ST_double *y = NULL;
    ST_double *x = NULL;
    ST_double *controls = NULL;
    ST_int *fe_vars = NULL;
    ST_int *by_groups = NULL;
    ST_double *weights = NULL;
    ST_int *valid_idx = NULL;
    ST_int N_valid = 0;

    /* Allocate arrays */
    y = (ST_double *)malloc(N * sizeof(ST_double));
    x = (ST_double *)malloc(N * sizeof(ST_double));
    valid_idx = (ST_int *)malloc(N * sizeof(ST_int));

    if (!y || !x || !valid_idx) {
        free(y); free(x); free(valid_idx);
        return CBINSCATTER_ERR_MEMORY;
    }

    if (config->has_controls && config->num_controls > 0) {
        controls = (ST_double *)malloc(N * config->num_controls * sizeof(ST_double));
        if (!controls) {
            free(y); free(x); free(valid_idx);
            return CBINSCATTER_ERR_MEMORY;
        }
    }

    if (config->has_absorb && config->num_absorb > 0) {
        fe_vars = (ST_int *)malloc(N * config->num_absorb * sizeof(ST_int));
        if (!fe_vars) {
            free(y); free(x); free(valid_idx); free(controls);
            return CBINSCATTER_ERR_MEMORY;
        }
    }

    if (config->has_by) {
        by_groups = (ST_int *)malloc(N * sizeof(ST_int));
        if (!by_groups) {
            free(y); free(x); free(valid_idx); free(controls); free(fe_vars);
            return CBINSCATTER_ERR_MEMORY;
        }
    }

    if (config->has_weights) {
        weights = (ST_double *)malloc(N * sizeof(ST_double));
        if (!weights) {
            free(y); free(x); free(valid_idx); free(controls);
            free(fe_vars); free(by_groups);
            return CBINSCATTER_ERR_MEMORY;
        }
    }

    /* Read data from Stata with 8x unrolling for main variables */
    var_idx = 1;

    /* Read y (var 1) - 8x unrolled */
    {
        ST_int i_end = N - (N % 8);
        ST_double v0, v1, v2, v3, v4, v5, v6, v7;
        for (i = 0; i < i_end; i += 8) {
            SF_vdata(var_idx, i + in1, &v0);
            SF_vdata(var_idx, i + 1 + in1, &v1);
            SF_vdata(var_idx, i + 2 + in1, &v2);
            SF_vdata(var_idx, i + 3 + in1, &v3);
            SF_vdata(var_idx, i + 4 + in1, &v4);
            SF_vdata(var_idx, i + 5 + in1, &v5);
            SF_vdata(var_idx, i + 6 + in1, &v6);
            SF_vdata(var_idx, i + 7 + in1, &v7);
            y[i] = v0; y[i+1] = v1; y[i+2] = v2; y[i+3] = v3;
            y[i+4] = v4; y[i+5] = v5; y[i+6] = v6; y[i+7] = v7;
        }
        for (; i < N; i++) {
            SF_vdata(var_idx, i + in1, &y[i]);
        }
    }
    var_idx++;

    /* Read x (var 2) - 8x unrolled */
    {
        ST_int i_end = N - (N % 8);
        ST_double v0, v1, v2, v3, v4, v5, v6, v7;
        for (i = 0; i < i_end; i += 8) {
            SF_vdata(var_idx, i + in1, &v0);
            SF_vdata(var_idx, i + 1 + in1, &v1);
            SF_vdata(var_idx, i + 2 + in1, &v2);
            SF_vdata(var_idx, i + 3 + in1, &v3);
            SF_vdata(var_idx, i + 4 + in1, &v4);
            SF_vdata(var_idx, i + 5 + in1, &v5);
            SF_vdata(var_idx, i + 6 + in1, &v6);
            SF_vdata(var_idx, i + 7 + in1, &v7);
            x[i] = v0; x[i+1] = v1; x[i+2] = v2; x[i+3] = v3;
            x[i+4] = v4; x[i+5] = v5; x[i+6] = v6; x[i+7] = v7;
        }
        for (; i < N; i++) {
            SF_vdata(var_idx, i + in1, &x[i]);
        }
    }
    var_idx++;

    /* Read controls - 8x unrolled */
    if (controls) {
        for (j = 0; j < config->num_controls; j++) {
            ST_double *col = &controls[j * N];
            ST_int i_end = N - (N % 8);
            ST_double v0, v1, v2, v3, v4, v5, v6, v7;
            for (i = 0; i < i_end; i += 8) {
                SF_vdata(var_idx, i + in1, &v0);
                SF_vdata(var_idx, i + 1 + in1, &v1);
                SF_vdata(var_idx, i + 2 + in1, &v2);
                SF_vdata(var_idx, i + 3 + in1, &v3);
                SF_vdata(var_idx, i + 4 + in1, &v4);
                SF_vdata(var_idx, i + 5 + in1, &v5);
                SF_vdata(var_idx, i + 6 + in1, &v6);
                SF_vdata(var_idx, i + 7 + in1, &v7);
                col[i] = v0; col[i+1] = v1; col[i+2] = v2; col[i+3] = v3;
                col[i+4] = v4; col[i+5] = v5; col[i+6] = v6; col[i+7] = v7;
            }
            for (; i < N; i++) {
                SF_vdata(var_idx, i + in1, &col[i]);
            }
            var_idx++;
        }
    }

    /* Read absorb variables (as integers) */
    if (fe_vars) {
        for (j = 0; j < config->num_absorb; j++) {
            for (i = 0; i < N; i++) {
                SF_vdata(var_idx, i + in1, &val);
                fe_vars[j * N + i] = (ST_int)val;
            }
            var_idx++;
        }
    }

    /* Read by variable */
    if (by_groups) {
        for (i = 0; i < N; i++) {
            SF_vdata(var_idx, i + in1, &val);
            by_groups[i] = (ST_int)val;
        }
        var_idx++;
    }

    /* Read weights - 8x unrolled */
    if (weights) {
        ST_int i_end = N - (N % 8);
        ST_double v0, v1, v2, v3, v4, v5, v6, v7;
        for (i = 0; i < i_end; i += 8) {
            SF_vdata(var_idx, i + in1, &v0);
            SF_vdata(var_idx, i + 1 + in1, &v1);
            SF_vdata(var_idx, i + 2 + in1, &v2);
            SF_vdata(var_idx, i + 3 + in1, &v3);
            SF_vdata(var_idx, i + 4 + in1, &v4);
            SF_vdata(var_idx, i + 5 + in1, &v5);
            SF_vdata(var_idx, i + 6 + in1, &v6);
            SF_vdata(var_idx, i + 7 + in1, &v7);
            weights[i] = v0; weights[i+1] = v1; weights[i+2] = v2; weights[i+3] = v3;
            weights[i+4] = v4; weights[i+5] = v5; weights[i+6] = v6; weights[i+7] = v7;
        }
        for (; i < N; i++) {
            SF_vdata(var_idx, i + in1, &weights[i]);
        }
        var_idx++;
    }

    /* Build valid index array (streaming - single pass) */
    N_valid = 0;
    for (i = 0; i < N; i++) {
        int valid = 1;

        /* Check y and x */
        if (IS_MISSING(y[i]) || IS_MISSING(x[i])) {
            valid = 0;
        }

        /* Check controls */
        if (valid && controls) {
            for (j = 0; j < config->num_controls; j++) {
                if (IS_MISSING(controls[j * N + i])) {
                    valid = 0;
                    break;
                }
            }
        }

        /* Check absorb variables */
        if (valid && fe_vars) {
            for (j = 0; j < config->num_absorb; j++) {
                if (fe_vars[j * N + i] <= 0) {
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
            if (IS_MISSING(weights[i]) || weights[i] <= 0) {
                valid = 0;
            }
        }

        if (valid) {
            valid_idx[N_valid++] = i;
        }
    }

    /* Compact data using valid indices (streaming compaction) */
    if (N_valid < N && N_valid > 0) {
        for (ST_int k = 0; k < N_valid; k++) {
            ST_int src = valid_idx[k];
            if (src != k) {
                y[k] = y[src];
                x[k] = x[src];
                if (controls) {
                    for (j = 0; j < config->num_controls; j++) {
                        controls[j * N_valid + k] = controls[j * N + src];
                    }
                }
                if (fe_vars) {
                    for (j = 0; j < config->num_absorb; j++) {
                        fe_vars[j * N_valid + k] = fe_vars[j * N + src];
                    }
                }
                if (by_groups) {
                    by_groups[k] = by_groups[src];
                }
                if (weights) {
                    weights[k] = weights[src];
                }
            }
        }
    }

    free(valid_idx);

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
    SF_scal_save("__cbinscatter_N", (ST_double)results->total_obs);
    SF_scal_save("__cbinscatter_N_dropped", (ST_double)results->obs_dropped);
    SF_scal_save("__cbinscatter_num_groups", (ST_double)results->num_by_groups);

    /* Store bin data to matrix __cbinscatter_bins */
    row = 1;
    for (g = 0; g < results->num_by_groups; g++) {
        ByGroupResult *group = &results->groups[g];
        for (b = 0; b < group->num_bins; b++) {
            BinStats *bin = &group->bins[b];
            SF_mat_store("__cbinscatter_bins", row, 1, (ST_double)group->by_group_id);
            SF_mat_store("__cbinscatter_bins", row, 2, (ST_double)bin->bin_id);
            SF_mat_store("__cbinscatter_bins", row, 3, bin->x_mean);
            SF_mat_store("__cbinscatter_bins", row, 4, bin->y_mean);
            SF_mat_store("__cbinscatter_bins", row, 5, bin->x_se);
            SF_mat_store("__cbinscatter_bins", row, 6, bin->y_se);
            SF_mat_store("__cbinscatter_bins", row, 7, (ST_double)bin->n_obs);
            row++;
        }
    }

    /* Store fit coefficients if line fitting was done */
    if (config->linetype > 0) {
        for (g = 0; g < results->num_by_groups; g++) {
            ByGroupResult *group = &results->groups[g];
            if (group->fit_coefs && group->fit_order > 0) {
                for (ST_int j = 0; j <= group->fit_order && j < 4; j++) {
                    SF_mat_store("__cbinscatter_coefs", g + 1, j + 1, group->fit_coefs[j]);
                }
                SF_mat_store("__cbinscatter_fit_stats", g + 1, 1, group->fit_r2);
                SF_mat_store("__cbinscatter_fit_stats", g + 1, 2, (ST_double)group->fit_n);
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
    ST_int N, N_valid;

    cbinscatter_results_init(&results);

    t_start = ctools_timer_seconds();

    /* Read configuration */
    rc = read_config(&config);
    if (rc != CBINSCATTER_OK) return rc;

    if (config.verbose) {
        SF_display("cbinscatter: starting computation\n");
    }

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

    if (config.verbose) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cbinscatter: loaded %d observations (%d dropped)\n",
                 N_valid, N - N_valid);
        SF_display(msg);
    }

    /* Residualize if needed */
    t_resid = t_load;
    if (config.has_controls || config.has_absorb) {
        ST_int dropped = 0;

        if (config.verbose) {
            if (config.method == 1) {
                SF_display("cbinscatter: using binsreg method (bin on raw X, residualize Y only)\n");
            } else {
                SF_display("cbinscatter: using classic method (residualize both X and Y)\n");
            }
        }

        /*
         * Method 0 (classic): Residualize both Y and X, bin on residualized X
         * Method 1 (binsreg): Only residualize Y, bin on raw X
         *   (Cattaneo et al. "On Binscatter" approach)
         */
        if (config.method == 1) {
            /* Binsreg method: Y-only residualization */
            if (config.has_absorb && config.has_controls) {
                rc = combined_residualize_y_only(y, controls, fe_vars, N_valid,
                                                  config.num_controls, config.num_absorb,
                                                  weights, config.weight_type,
                                                  config.maxiter, config.tolerance,
                                                  config.verbose, &dropped);
            } else if (config.has_absorb) {
                rc = hdfe_residualize_y_only(y, fe_vars, N_valid, config.num_absorb,
                                              weights, config.weight_type,
                                              config.maxiter, config.tolerance,
                                              config.verbose, &dropped);
            } else {
                rc = ols_residualize_y_only(y, controls, N_valid, config.num_controls,
                                             weights, config.weight_type);
            }
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
        }

        if (rc != CBINSCATTER_OK) {
            SF_error("cbinscatter: residualization failed\n");
            goto cleanup;
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

    /* Process each by-group */
    for (ST_int g = 0; g < num_groups; g++) {
        ByGroupResult *group = &results.groups[g];
        group->by_group_id = g + 1;

        /* Count observations in this group */
        ST_int n_group = 0;
        if (config.has_by) {
            for (ST_int i = 0; i < N_valid; i++) {
                if (by_groups[i] == g + 1) n_group++;
            }
        } else {
            n_group = N_valid;
        }

        if (n_group < 2) {
            group->num_bins = 0;
            group->bins = NULL;
            continue;
        }

        /* Allocate and extract group data */
        ST_double *y_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        ST_double *x_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        ST_double *w_group = NULL;
        if (weights) {
            w_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        }

        if (!y_group || !x_group || (weights && !w_group)) {
            free(y_group); free(x_group); free(w_group);
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        ST_int j = 0;
        for (ST_int i = 0; i < N_valid; i++) {
            if (!config.has_by || by_groups[i] == g + 1) {
                y_group[j] = y[i];
                x_group[j] = x[i];
                if (w_group) w_group[j] = weights[i];
                j++;
            }
        }

        /* Compute bins */
        if (config.discrete) {
            rc = compute_bins_discrete(y_group, x_group, w_group, n_group,
                                       config.compute_se, group);
        } else {
            rc = compute_bins_single_group(y_group, x_group, w_group, n_group,
                                           &config, group);
        }

        free(y_group);
        free(x_group);
        free(w_group);

        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    t_bins = ctools_timer_seconds();

    if (config.verbose) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cbinscatter: bin computation complete (%.3f sec)\n",
                 t_bins - t_resid);
        SF_display(msg);
    }

    /* Fit lines if requested */
    t_fit = t_bins;
    if (config.linetype > 0) {
        rc = fit_all_groups(y, x, by_groups, weights, N_valid, &config, &results);
        if (rc != CBINSCATTER_OK) goto cleanup;

        t_fit = ctools_timer_seconds();

        if (config.verbose) {
            char msg[256];
            snprintf(msg, sizeof(msg), "cbinscatter: line fitting complete (%.3f sec)\n",
                     t_fit - t_bins);
            SF_display(msg);
        }
    }

    /* Store results to Stata */
    rc = store_results(&results, &config);

    t_total = ctools_timer_seconds() - t_start;

    /* Store timing info */
    SF_scal_save("__cbinscatter_time_load", t_load - t_start);
    SF_scal_save("__cbinscatter_time_resid", t_resid - t_load);
    SF_scal_save("__cbinscatter_time_bins", t_bins - t_resid);
    SF_scal_save("__cbinscatter_time_fit", t_fit - t_bins);
    SF_scal_save("__cbinscatter_time_total", t_total);

    if (config.verbose) {
        char msg[256];
        snprintf(msg, sizeof(msg), "cbinscatter: total time %.3f sec\n", t_total);
        SF_display(msg);
    }

cleanup:
    free(y);
    free(x);
    free(controls);
    free(fe_vars);
    free(by_groups);
    free(weights);
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
