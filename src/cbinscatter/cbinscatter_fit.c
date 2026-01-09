/*
 * cbinscatter_fit.c
 *
 * Line fitting implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 *
 * Optimizations:
 * - Single-pass linear fitting with closed-form solution
 * - No design matrix allocation for any polynomial order
 * - Merged computation of coefficients and R²
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"
#include "cbinscatter_fit.h"
#include "cbinscatter_resid.h"

/* ========================================================================
 * Fast Linear Fit (Single Pass)
 *
 * For y = a + b*x, computes coefficients and R² in one pass using:
 *   b = (Σw*x*y - Σw*x*Σw*y/Σw) / (Σw*x² - (Σw*x)²/Σw)
 *   a = ȳ - b*x̄
 *   R² = 1 - RSS/TSS
 * ======================================================================== */

static ST_retcode fit_linear_fast(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_double *coefs,
    ST_double *r2
) {
    ST_double sum_w = 0.0, sum_x = 0.0, sum_y = 0.0;
    ST_double sum_xx = 0.0, sum_xy = 0.0, sum_yy = 0.0;
    ST_int i;

    /* Single pass: accumulate all needed sums */
    if (weights == NULL) {
        for (i = 0; i < N; i++) {
            ST_double xi = x[i];
            ST_double yi = y[i];
            sum_x += xi;
            sum_y += yi;
            sum_xx += xi * xi;
            sum_xy += xi * yi;
            sum_yy += yi * yi;
        }
        sum_w = (ST_double)N;
    } else {
        for (i = 0; i < N; i++) {
            ST_double w = weights[i];
            ST_double xi = x[i];
            ST_double yi = y[i];
            ST_double wxi = w * xi;
            ST_double wyi = w * yi;
            sum_w += w;
            sum_x += wxi;
            sum_y += wyi;
            sum_xx += wxi * xi;
            sum_xy += wxi * yi;
            sum_yy += wyi * yi;
        }
    }

    /* Compute means */
    ST_double x_mean = sum_x / sum_w;
    ST_double y_mean = sum_y / sum_w;

    /* Compute slope: Cov(x,y) / Var(x) */
    ST_double var_x = sum_xx / sum_w - x_mean * x_mean;
    ST_double cov_xy = sum_xy / sum_w - x_mean * y_mean;

    if (var_x <= 0.0) {
        /* Degenerate case: constant x */
        coefs[0] = y_mean;
        coefs[1] = 0.0;
        *r2 = 0.0;
        return CBINSCATTER_OK;
    }

    ST_double slope = cov_xy / var_x;
    ST_double intercept = y_mean - slope * x_mean;

    coefs[0] = intercept;
    coefs[1] = slope;

    /* Compute R² = 1 - RSS/TSS
     * TSS = Σw(y - ȳ)² = Σwy² - (Σwy)²/Σw
     * RSS = Σw(y - ŷ)² = Σw(y - a - bx)²
     *
     * For efficiency, use: R² = (slope * Cov(x,y)) / Var(y)
     * which equals correlation² for linear regression
     */
    ST_double var_y = sum_yy / sum_w - y_mean * y_mean;
    if (var_y > 0.0) {
        *r2 = (slope * cov_xy) / var_y;
        if (*r2 < 0.0) *r2 = 0.0;
        if (*r2 > 1.0) *r2 = 1.0;
    } else {
        *r2 = 0.0;
    }

    return CBINSCATTER_OK;
}

/* ========================================================================
 * Fast Quadratic Fit (Two Pass, No Design Matrix)
 *
 * For y = a + b*x + c*x², we need X'X (3x3) and X'y (3x1)
 * Computed directly without storing the design matrix.
 * ======================================================================== */

static ST_retcode fit_quadratic_fast(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_double *coefs,
    ST_double *r2
) {
    /* Accumulators for X'X (symmetric, so store 6 unique elements) */
    ST_double s_1 = 0.0;    /* Σw */
    ST_double s_x = 0.0;    /* Σwx */
    ST_double s_x2 = 0.0;   /* Σwx² */
    ST_double s_x3 = 0.0;   /* Σwx³ */
    ST_double s_x4 = 0.0;   /* Σwx⁴ */
    /* Accumulators for X'y */
    ST_double s_y = 0.0;    /* Σwy */
    ST_double s_xy = 0.0;   /* Σwxy */
    ST_double s_x2y = 0.0;  /* Σwx²y */
    /* For R² */
    ST_double s_yy = 0.0;   /* Σwy² */
    ST_int i;

    /* Single pass to accumulate */
    if (weights == NULL) {
        for (i = 0; i < N; i++) {
            ST_double xi = x[i];
            ST_double yi = y[i];
            ST_double x2 = xi * xi;
            ST_double x3 = x2 * xi;
            ST_double x4 = x2 * x2;

            s_1 += 1.0;
            s_x += xi;
            s_x2 += x2;
            s_x3 += x3;
            s_x4 += x4;
            s_y += yi;
            s_xy += xi * yi;
            s_x2y += x2 * yi;
            s_yy += yi * yi;
        }
    } else {
        for (i = 0; i < N; i++) {
            ST_double w = weights[i];
            ST_double xi = x[i];
            ST_double yi = y[i];
            ST_double x2 = xi * xi;

            s_1 += w;
            s_x += w * xi;
            s_x2 += w * x2;
            s_x3 += w * x2 * xi;
            s_x4 += w * x2 * x2;
            s_y += w * yi;
            s_xy += w * xi * yi;
            s_x2y += w * x2 * yi;
            s_yy += w * yi * yi;
        }
    }

    /* Build X'X matrix (3x3 symmetric) */
    ST_double XtX[9] = {
        s_1,  s_x,  s_x2,
        s_x,  s_x2, s_x3,
        s_x2, s_x3, s_x4
    };

    /* Build X'y vector */
    ST_double Xty[3] = { s_y, s_xy, s_x2y };

    /* Cholesky solve */
    ST_retcode rc = cholesky_decompose(XtX, 3);
    if (rc != CBINSCATTER_OK) return rc;

    coefs[0] = Xty[0];
    coefs[1] = Xty[1];
    coefs[2] = Xty[2];
    rc = cholesky_solve(XtX, 3, coefs);
    if (rc != CBINSCATTER_OK) return rc;

    /* Compute R² */
    ST_double y_mean = s_y / s_1;
    ST_double tss = s_yy - s_y * s_y / s_1;

    if (tss > 0.0) {
        /* RSS = Σw(y - ŷ)² = Σwy² - 2*coefs'*X'y + coefs'*X'X*coefs
         * But simpler: RSS = TSS - ESS where ESS = coefs'*X'y - n*ȳ²  (for centered)
         * Actually compute RSS directly in second pass for accuracy */
        ST_double rss = 0.0;
        if (weights == NULL) {
            for (i = 0; i < N; i++) {
                ST_double xi = x[i];
                ST_double fitted = coefs[0] + coefs[1] * xi + coefs[2] * xi * xi;
                ST_double resid = y[i] - fitted;
                rss += resid * resid;
            }
        } else {
            for (i = 0; i < N; i++) {
                ST_double xi = x[i];
                ST_double fitted = coefs[0] + coefs[1] * xi + coefs[2] * xi * xi;
                ST_double resid = y[i] - fitted;
                rss += weights[i] * resid * resid;
            }
        }
        *r2 = 1.0 - rss / tss;
        if (*r2 < 0.0) *r2 = 0.0;
    } else {
        *r2 = 0.0;
    }

    return CBINSCATTER_OK;
}

/* ========================================================================
 * General Polynomial Fitting (No Design Matrix Storage)
 *
 * For higher-order polynomials, compute X'X and X'y directly.
 * ======================================================================== */

ST_retcode fit_polynomial(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_int order,
    ST_double *coefs,
    ST_double *r2
) {
    /* Use specialized fast paths for common cases */
    if (order == 1) {
        return fit_linear_fast(y, x, weights, N, coefs, r2);
    }
    if (order == 2) {
        return fit_quadratic_fast(y, x, weights, N, coefs, r2);
    }

    /* General case for order >= 3 */
    ST_int K = order + 1;
    ST_int K2 = K * K;
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_double *xpows = NULL;  /* Powers of x for current observation */
    ST_retcode rc = CBINSCATTER_OK;
    ST_int i, j, k;
    ST_double sum_w = 0.0, sum_y = 0.0, sum_yy = 0.0;

    /* Allocate only what we need - no design matrix */
    XtX = (ST_double *)calloc(K2, sizeof(ST_double));
    Xty = (ST_double *)calloc(K, sizeof(ST_double));
    xpows = (ST_double *)malloc(K * sizeof(ST_double));

    if (!XtX || !Xty || !xpows) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Single pass: compute X'X and X'y directly */
    for (i = 0; i < N; i++) {
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        ST_double xi = x[i];
        ST_double yi = y[i];
        ST_double wyi = w * yi;

        /* Compute powers of x */
        xpows[0] = 1.0;
        for (j = 1; j < K; j++) {
            xpows[j] = xpows[j-1] * xi;
        }

        /* Accumulate X'y */
        for (j = 0; j < K; j++) {
            Xty[j] += xpows[j] * wyi;
        }

        /* Accumulate upper triangle of X'X */
        for (j = 0; j < K; j++) {
            ST_double wxj = w * xpows[j];
            for (k = j; k < K; k++) {
                XtX[j * K + k] += wxj * xpows[k];
            }
        }

        /* For R² computation */
        sum_w += w;
        sum_y += wyi;
        sum_yy += wyi * yi;
    }

    /* Mirror to lower triangle */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            XtX[k * K + j] = XtX[j * K + k];
        }
    }

    /* Cholesky solve */
    rc = cholesky_decompose(XtX, K);
    if (rc != CBINSCATTER_OK) goto cleanup;

    memcpy(coefs, Xty, K * sizeof(ST_double));
    rc = cholesky_solve(XtX, K, coefs);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Compute R² with second pass for RSS */
    ST_double y_mean = sum_y / sum_w;
    ST_double tss = sum_yy - sum_y * sum_y / sum_w;

    if (tss > 0.0) {
        ST_double rss = 0.0;
        for (i = 0; i < N; i++) {
            ST_double w = (weights != NULL) ? weights[i] : 1.0;
            ST_double xi = x[i];

            /* Horner's method for polynomial evaluation */
            ST_double fitted = coefs[K - 1];
            for (j = K - 2; j >= 0; j--) {
                fitted = fitted * xi + coefs[j];
            }

            ST_double resid = y[i] - fitted;
            rss += w * resid * resid;
        }
        *r2 = 1.0 - rss / tss;
        if (*r2 < 0.0) *r2 = 0.0;
    } else {
        *r2 = 0.0;
    }

cleanup:
    free(XtX);
    free(Xty);
    free(xpows);
    return rc;
}

/* ========================================================================
 * Fit All Groups (Optimized)
 * ======================================================================== */

ST_retcode fit_all_groups(
    const ST_double *y,
    const ST_double *x,
    const ST_int *by_groups,
    const ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    BinscatterResults *results
) {
    ST_int g, i;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int order = config->linetype;
    ST_int max_coefs = order + 1;

    if (order < 1) return CBINSCATTER_OK;

    /* For single group (no by()), fit directly without copying */
    if (results->num_by_groups == 1 && by_groups == NULL) {
        ByGroupResult *group = &results->groups[0];

        if (N < order + 1) {
            group->fit_order = 0;
            group->fit_coefs = NULL;
            return CBINSCATTER_OK;
        }

        group->fit_coefs = (ST_double *)calloc(max_coefs, sizeof(ST_double));
        if (!group->fit_coefs) return CBINSCATTER_ERR_MEMORY;

        rc = fit_polynomial(y, x, weights, N, order,
                            group->fit_coefs, &group->fit_r2);

        if (rc == CBINSCATTER_OK) {
            group->fit_order = order;
            group->fit_n = N;
        } else {
            free(group->fit_coefs);
            group->fit_coefs = NULL;
            group->fit_order = 0;
        }
        return rc;
    }

    /* Multiple groups: need to extract data per group */
    ST_double *y_group = NULL, *x_group = NULL, *w_group = NULL;
    ST_int *group_counts = NULL;
    ST_int max_group_size = 0;

    /* First pass: count observations per group */
    group_counts = (ST_int *)calloc(results->num_by_groups, sizeof(ST_int));
    if (!group_counts) return CBINSCATTER_ERR_MEMORY;

    for (i = 0; i < N; i++) {
        ST_int g_id = by_groups[i] - 1;
        if (g_id >= 0 && g_id < results->num_by_groups) {
            group_counts[g_id]++;
            if (group_counts[g_id] > max_group_size) {
                max_group_size = group_counts[g_id];
            }
        }
    }

    /* Allocate buffers for largest group */
    y_group = (ST_double *)malloc(max_group_size * sizeof(ST_double));
    x_group = (ST_double *)malloc(max_group_size * sizeof(ST_double));
    if (weights != NULL) {
        w_group = (ST_double *)malloc(max_group_size * sizeof(ST_double));
    }

    if (!y_group || !x_group || (weights != NULL && !w_group)) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Process each group */
    for (g = 0; g < results->num_by_groups; g++) {
        ByGroupResult *group = &results->groups[g];
        ST_int n_group = group_counts[g];

        if (n_group < order + 1) {
            group->fit_order = 0;
            group->fit_coefs = NULL;
            continue;
        }

        /* Extract group data */
        ST_int j = 0;
        for (i = 0; i < N; i++) {
            if (by_groups[i] == g + 1) {
                y_group[j] = y[i];
                x_group[j] = x[i];
                if (w_group != NULL) {
                    w_group[j] = weights[i];
                }
                j++;
            }
        }

        /* Allocate and fit */
        group->fit_coefs = (ST_double *)calloc(max_coefs, sizeof(ST_double));
        if (!group->fit_coefs) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto cleanup;
        }

        rc = fit_polynomial(y_group, x_group, w_group, n_group, order,
                            group->fit_coefs, &group->fit_r2);

        if (rc == CBINSCATTER_OK) {
            group->fit_order = order;
            group->fit_n = n_group;
        } else {
            free(group->fit_coefs);
            group->fit_coefs = NULL;
            group->fit_order = 0;
            rc = CBINSCATTER_OK;  /* Don't fail entire operation */
        }
    }

cleanup:
    free(y_group);
    free(x_group);
    free(w_group);
    free(group_counts);
    return rc;
}

/* ========================================================================
 * Evaluate Polynomial
 * ======================================================================== */

void eval_polynomial(
    const ST_double *x,
    ST_int N,
    const ST_double *coefs,
    ST_int order,
    ST_double *fitted
) {
    ST_int i, j;

    /* Use Horner's method for numerical stability */
    for (i = 0; i < N; i++) {
        ST_double xi = x[i];
        ST_double val = coefs[order];
        for (j = order - 1; j >= 0; j--) {
            val = val * xi + coefs[j];
        }
        fitted[i] = val;
    }
}
