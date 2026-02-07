/*
 * cbinscatter_binsreg.c
 *
 * Binsreg method implementation for cbinscatter
 * Implements the Cattaneo et al. "On Binscatter" approach
 * Part of the ctools Stata plugin suite
 *
 * The key insight: binsreg runs reg Y ibn.xcat W, nocon
 * which gives different results than first residualizing Y on W
 * then computing bin means. This file implements the correct approach.
 *
 * For HDFE case with absorb(): we must run FWL correctly by
 * demeaning BOTH Y AND the bin indicators within FE groups.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"
#include "cbinscatter_binsreg.h"
#include "cbinscatter_resid.h"
#include "../ctools_config.h"
#include "../ctools_ols.h"

/* ========================================================================
 * Helper: Compute within-bin means
 * ======================================================================== */

static ST_retcode compute_bin_means(
    const ST_double *vals,      /* Values to average (N x 1) */
    const ST_int *bin_ids,      /* Bin assignments (1-based) */
    const ST_double *weights,   /* Weights (NULL if unweighted) */
    ST_int N,
    ST_int num_bins,
    ST_double *bin_means        /* Output: mean for each bin (num_bins x 1) */
) {
    ST_double *sum_w = (ST_double *)calloc(num_bins, sizeof(ST_double));
    ST_double *sum_v = (ST_double *)calloc(num_bins, sizeof(ST_double));

    if (!sum_w || !sum_v) {
        free(sum_w);
        free(sum_v);
        return CBINSCATTER_ERR_MEMORY;
    }

    /* Accumulate weighted sums */
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i] - 1;  /* Convert to 0-based */
        if (b < 0 || b >= num_bins) continue;

        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        sum_w[b] += w;
        sum_v[b] += w * vals[i];
    }

    /* Compute means */
    for (ST_int b = 0; b < num_bins; b++) {
        if (sum_w[b] > 0) {
            bin_means[b] = sum_v[b] / sum_w[b];
        } else {
            bin_means[b] = 0.0;
        }
    }

    free(sum_w);
    free(sum_v);
    return CBINSCATTER_OK;
}

/* ========================================================================
 * Adjust bin Y-means using binsreg regression method (no FE)
 * ======================================================================== */

ST_retcode adjust_bins_binsreg(
    const ST_double *y,
    const ST_double *controls,
    const ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int K,
    ByGroupResult *result
) {
    ST_retcode rc = CBINSCATTER_OK;
    ST_int num_bins = result->num_bins;

    /* Allocate all working arrays in a single block for cache locality.
     * Layout: [y_bin_means | w_bin_means | gamma | XtX | Xty | w_overall_means | y_demeaned | w_demeaned] */
    ST_double *y_bin_means = NULL;      /* E[Y|bin] for each bin */
    ST_double *w_bin_means = NULL;      /* E[W_k|bin] for each bin, all K controls */
    ST_double *y_demeaned = NULL;       /* Y - E[Y|bin] */
    ST_double *w_demeaned = NULL;       /* W - E[W|bin] for all controls */
    ST_double *gamma = NULL;            /* Regression coefficients on controls */
    ST_double *XtX = NULL;              /* K x K normal equations matrix */
    ST_double *Xty = NULL;              /* K x 1 right-hand side */
    ST_double *w_overall_means = NULL;  /* Overall mean of each control */
    void *binsreg_block = NULL;         /* Single allocation block */

    {
        /* Compute total size for consolidated allocation */
        size_t sz_small = (size_t)num_bins              /* y_bin_means */
                        + (size_t)num_bins * K          /* w_bin_means */
                        + (size_t)K                     /* gamma */
                        + (size_t)K * K                 /* XtX */
                        + (size_t)K                     /* Xty */
                        + (size_t)K;                    /* w_overall_means */
        size_t sz_large = (size_t)N + (size_t)N * K;    /* y_demeaned + w_demeaned */
        size_t total_doubles = sz_small + sz_large;

        binsreg_block = calloc(total_doubles, sizeof(ST_double));
        if (!binsreg_block) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        /* Partition the block */
        ST_double *p = (ST_double *)binsreg_block;
        y_bin_means = p;      p += num_bins;
        w_bin_means = p;      p += (size_t)num_bins * K;
        gamma = p;            p += K;
        XtX = p;              p += (size_t)K * K;
        Xty = p;              p += K;
        w_overall_means = p;  p += K;
        y_demeaned = p;       p += N;
        w_demeaned = p;       /* p += (size_t)N * K; */
    }

    /* Step 1: Compute within-bin means of Y */
    rc = compute_bin_means(y, bin_ids, weights, N, num_bins, y_bin_means);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Step 2: Compute within-bin means of each control variable */
    for (ST_int k = 0; k < K; k++) {
        rc = compute_bin_means(&controls[k * N], bin_ids, weights, N, num_bins,
                               &w_bin_means[k * num_bins]);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* Step 3: Create demeaned versions */
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i] - 1;
        if (b < 0 || b >= num_bins) {
            y_demeaned[i] = 0.0;
            for (ST_int k = 0; k < K; k++) {
                w_demeaned[k * N + i] = 0.0;
            }
            continue;
        }

        y_demeaned[i] = y[i] - y_bin_means[b];
        for (ST_int k = 0; k < K; k++) {
            w_demeaned[k * N + i] = controls[k * N + i] - w_bin_means[k * num_bins + b];
        }
    }

    /* Step 4: Run OLS: Y_demeaned ~ W_demeaned to get gamma */
    /* Build X'X and X'y */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights != NULL) ? weights[i] : 1.0;

        for (ST_int j = 0; j < K; j++) {
            ST_double wj = w_demeaned[j * N + i];
            Xty[j] += w * wj * y_demeaned[i];

            for (ST_int k = j; k < K; k++) {
                XtX[j * K + k] += w * wj * w_demeaned[k * N + i];
            }
        }
    }

    /* Mirror to lower triangle */
    for (ST_int j = 0; j < K; j++) {
        for (ST_int k = j + 1; k < K; k++) {
            XtX[k * K + j] = XtX[j * K + k];
        }
    }

    /* Solve (X'X) * gamma = X'y using shared Cholesky */
    if (ctools_solve_cholesky(XtX, Xty, K, gamma) != 0) {
        /* If singular, gamma remains zero - bin means are just raw Y means */
        memset(gamma, 0, K * sizeof(ST_double));
    }

    /* Step 5: Compute overall mean of each control (for "at mean" evaluation) */
    ST_double total_weight = 0.0;
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i] - 1;
        if (b < 0 || b >= num_bins) continue;
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        total_weight += w;
        for (ST_int k = 0; k < K; k++) {
            w_overall_means[k] += w * controls[k * N + i];
        }
    }
    if (total_weight > 0) {
        for (ST_int k = 0; k < K; k++) {
            w_overall_means[k] /= total_weight;
        }
    }

    /* Step 6: Update bin y_means
     *
     * binsreg evaluates predictions at the mean of controls:
     *   y_plot = beta_j + sum_k(mean(W_k) * gamma_k)
     *
     * where beta_j = E[Y|bin=j] - E[W_k|bin=j] * gamma_k
     *
     * So: y_plot = E[Y|bin=j] - sum_k((E[W_k|bin=j] - mean(W_k)) * gamma_k)
     */
    for (ST_int b = 0; b < num_bins; b++) {
        ST_double adjustment = 0.0;
        for (ST_int k = 0; k < K; k++) {
            /* Adjustment is (E[W_k|bin] - mean(W_k)) * gamma_k */
            adjustment += (w_bin_means[k * num_bins + b] - w_overall_means[k]) * gamma[k];
        }
        result->bins[b].y_mean = y_bin_means[b] - adjustment;
    }

cleanup:
    free(binsreg_block);

    return rc;
}

/* ========================================================================
 * Adjust bin Y-means using binsreg method with fixed effects
 *
 * This implements the correct FWL approach:
 * 1. Save mean(Y) before HDFE demeaning
 * 2. HDFE demean Y
 * 3. Create bin indicators and HDFE demean each one
 * 4. Regress demeaned Y on demeaned bin indicators
 * 5. Compute constant: const = mean(Y) - sum_k(N_k/N * beta_k)
 * 6. Final values: beta_k + const
 * ======================================================================== */

ST_retcode adjust_bins_binsreg_hdfe(
    ST_double *y,
    ST_double *controls,
    ST_int *fe_vars,
    const ST_int *bin_ids,
    const ST_double *weights,
    ST_int N,
    ST_int K_ctrl,
    ST_int K_fe,
    ST_int weight_type,
    ST_int maxiter,
    ST_double tolerance,
    ByGroupResult *result
) {
    ST_retcode rc = CBINSCATTER_OK;
    ST_int num_bins = result->num_bins;
    ST_int dropped = 0;

    /* Working arrays */
    ST_double *y_hdfe = NULL;
    ST_double *controls_hdfe = NULL;
    ST_double *bin_indicators = NULL;      /* N x num_bins matrix */
    ST_double *bin_indicators_dm = NULL;   /* HDFE-demeaned bin indicators */
    ST_double *beta = NULL;                /* Bin coefficients from FWL regression */
    ST_double *bin_weights = NULL;         /* N_k/N for each bin */
    ST_double *XtX = NULL;                 /* num_bins x num_bins */
    ST_double *Xty = NULL;                 /* num_bins x 1 */

    /* Step 1: Compute mean(Y) before HDFE demeaning */
    ST_double y_mean = 0.0;
    ST_double total_weight = 0.0;
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i] - 1;
        if (b < 0 || b >= num_bins) continue;
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        y_mean += w * y[i];
        total_weight += w;
    }
    if (total_weight > 0) {
        y_mean /= total_weight;
    }

    /* Allocate arrays */
    y_hdfe = (ST_double *)malloc(N * sizeof(ST_double));
    bin_indicators = (ST_double *)malloc((size_t)N * num_bins * sizeof(ST_double));
    bin_indicators_dm = (ST_double *)malloc((size_t)N * num_bins * sizeof(ST_double));
    beta = (ST_double *)calloc(num_bins, sizeof(ST_double));
    bin_weights = (ST_double *)calloc(num_bins, sizeof(ST_double));
    XtX = (ST_double *)calloc((size_t)num_bins * num_bins, sizeof(ST_double));
    Xty = (ST_double *)calloc(num_bins, sizeof(ST_double));

    if (!y_hdfe || !bin_indicators || !bin_indicators_dm ||
        !beta || !bin_weights || !XtX || !Xty) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    memcpy(y_hdfe, y, N * sizeof(ST_double));

    if (K_ctrl > 0 && controls != NULL) {
        controls_hdfe = (ST_double *)malloc((size_t)N * K_ctrl * sizeof(ST_double));
        if (!controls_hdfe) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }
        memcpy(controls_hdfe, controls, (size_t)N * K_ctrl * sizeof(ST_double));
    }

    /* Step 2: Create bin indicators and compute bin weights */
    memset(bin_indicators, 0, (size_t)N * num_bins * sizeof(ST_double));
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i] - 1;
        if (b >= 0 && b < num_bins) {
            bin_indicators[b * N + i] = 1.0;
            ST_double w = (weights != NULL) ? weights[i] : 1.0;
            bin_weights[b] += w;
        }
    }
    /* Normalize bin weights to fractions */
    for (ST_int b = 0; b < num_bins; b++) {
        bin_weights[b] /= total_weight;
    }

    /* Steps 3-5: HDFE demean Y, bin indicators, and controls simultaneously.
     * Using batch residualization reduces FE level array traversals from
     * (1 + num_bins + K_ctrl) * iters * 2G to iters * 2G. */
    memcpy(bin_indicators_dm, bin_indicators, (size_t)N * num_bins * sizeof(ST_double));

    {
        ST_int K_batch = 1 + num_bins + K_ctrl;
        ST_double **batch_vars = (ST_double **)malloc(K_batch * sizeof(ST_double *));
        if (!batch_vars) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        /* First variable: Y */
        batch_vars[0] = y_hdfe;
        /* Next num_bins: bin indicators */
        for (ST_int b = 0; b < num_bins; b++) {
            batch_vars[1 + b] = &bin_indicators_dm[b * N];
        }
        /* Last K_ctrl: control variables */
        if (controls_hdfe != NULL) {
            for (ST_int k = 0; k < K_ctrl; k++) {
                batch_vars[1 + num_bins + k] = &controls_hdfe[k * N];
            }
        }

        rc = hdfe_residualize_batch(batch_vars, K_batch, fe_vars, N, K_fe,
                                     weights, weight_type,
                                     maxiter, tolerance, 0, &dropped);
        free(batch_vars);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* Step 6: Run OLS: y_dm ~ bin_indicators_dm + controls_dm
     * This runs a SINGLE regression with both bins and controls.
     * For controls + absorb, this is the correct FWL approach.
     */

    ST_int num_regressors = num_bins + K_ctrl;
    ST_double *XtX_full = (ST_double *)calloc((size_t)num_regressors * num_regressors, sizeof(ST_double));
    ST_double *Xty_full = (ST_double *)calloc(num_regressors, sizeof(ST_double));
    ST_double *coef = (ST_double *)calloc(num_regressors, sizeof(ST_double));

    if (!XtX_full || !Xty_full || !coef) {
        free(XtX_full);
        free(Xty_full);
        free(coef);
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Build X'X and X'y for [bin_indicators_dm | controls_hdfe] */
    for (ST_int i = 0; i < N; i++) {
        ST_int b = bin_ids[i] - 1;
        if (b < 0 || b >= num_bins) continue;

        ST_double w = (weights != NULL) ? weights[i] : 1.0;

        /* Bin indicators part */
        for (ST_int j = 0; j < num_bins; j++) {
            ST_double xj = bin_indicators_dm[j * N + i];
            Xty_full[j] += w * xj * y_hdfe[i];

            /* Bin x Bin */
            for (ST_int k = j; k < num_bins; k++) {
                XtX_full[j * num_regressors + k] += w * xj * bin_indicators_dm[k * N + i];
            }

            /* Bin x Control */
            if (controls_hdfe != NULL) {
                for (ST_int k = 0; k < K_ctrl; k++) {
                    ST_int col_idx = num_bins + k;
                    XtX_full[j * num_regressors + col_idx] += w * xj * controls_hdfe[k * N + i];
                }
            }
        }

        /* Controls part */
        if (controls_hdfe != NULL) {
            for (ST_int j = 0; j < K_ctrl; j++) {
                ST_int row_idx = num_bins + j;
                ST_double wj = controls_hdfe[j * N + i];
                Xty_full[row_idx] += w * wj * y_hdfe[i];

                /* Control x Control */
                for (ST_int k = j; k < K_ctrl; k++) {
                    ST_int col_idx = num_bins + k;
                    XtX_full[row_idx * num_regressors + col_idx] += w * wj * controls_hdfe[k * N + i];
                }
            }
        }
    }

    /* Mirror to lower triangle */
    for (ST_int j = 0; j < num_regressors; j++) {
        for (ST_int k = j + 1; k < num_regressors; k++) {
            XtX_full[k * num_regressors + j] = XtX_full[j * num_regressors + k];
        }
    }

    /* Solve using regularized Cholesky (add small ridge for stability) */
    ST_double ridge = 1e-12;
    for (ST_int j = 0; j < num_regressors; j++) {
        XtX_full[j * num_regressors + j] += ridge;
    }

    if (ctools_solve_cholesky(XtX_full, Xty_full, num_regressors, coef) != 0) {
        /* Fallback: use simple bin means adjusted by y_mean */
        rc = compute_bin_means(y_hdfe, bin_ids, weights, N, num_bins, beta);
        if (rc == CBINSCATTER_OK) {
            for (ST_int b = 0; b < num_bins; b++) {
                result->bins[b].y_mean = beta[b] + y_mean;
            }
        }
        free(XtX_full);
        free(Xty_full);
        free(coef);
        goto cleanup;
    }

    /* Extract bin coefficients (first num_bins elements) */
    for (ST_int b = 0; b < num_bins; b++) {
        beta[b] = coef[b];
    }

    /* Step 7: Compute constant = mean(Y) - sum_k(N_k/N * beta_k)
     * Use Kahan summation for better numerical precision */
    ST_double sum_weighted_beta = 0.0;
    ST_double kahan_c = 0.0;  /* Compensation for lost low-order bits */
    for (ST_int b = 0; b < num_bins; b++) {
        ST_double y_val = bin_weights[b] * beta[b] - kahan_c;
        ST_double t = sum_weighted_beta + y_val;
        kahan_c = (t - sum_weighted_beta) - y_val;
        sum_weighted_beta = t;
    }
    ST_double const_adj = y_mean - sum_weighted_beta;

    /* Final values: beta[b] + const_adj */
    for (ST_int b = 0; b < num_bins; b++) {
        result->bins[b].y_mean = beta[b] + const_adj;
    }

    free(XtX_full);
    free(Xty_full);
    free(coef);

cleanup:
    free(y_hdfe);
    free(controls_hdfe);
    free(bin_indicators);
    free(bin_indicators_dm);
    free(beta);
    free(bin_weights);
    free(XtX);
    free(Xty);

    return rc;
}
