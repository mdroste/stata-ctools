/*
 * cbinscatter_fit.c
 *
 * Line fitting implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"
#include "cbinscatter_fit.h"
#include "cbinscatter_resid.h"

/* ========================================================================
 * Polynomial Fitting
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
    ST_int K = order + 1;  /* Number of coefficients */
    ST_double *X = NULL;   /* Design matrix (N x K) */
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int i, j, k;
    ST_double tss = 0.0, rss = 0.0, y_mean = 0.0, sum_w = 0.0;

    /* Allocate design matrix */
    X = (ST_double *)malloc(N * K * sizeof(ST_double));
    XtX = (ST_double *)malloc(K * K * sizeof(ST_double));
    Xty = (ST_double *)calloc(K, sizeof(ST_double));

    if (!X || !XtX || !Xty) {
        rc = CBINSCATTER_ERR_MEMORY;
        goto cleanup;
    }

    /* Build design matrix: [1, x, x^2, ...] */
    for (i = 0; i < N; i++) {
        ST_double xi = x[i];
        ST_double xpow = 1.0;
        for (j = 0; j < K; j++) {
            X[j * N + i] = xpow;  /* Column-major */
            xpow *= xi;
        }
    }

    /* Initialize X'X to zero */
    memset(XtX, 0, K * K * sizeof(ST_double));

    /* Compute X'X and X'y sequentially */
    for (i = 0; i < N; i++) {
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        ST_double wyi = w * y[i];

        for (j = 0; j < K; j++) {
            ST_double xj = X[j * N + i];
            ST_double wxj = w * xj;

            /* X'y */
            Xty[j] += xj * wyi;

            /* Upper triangle of X'X */
            for (k = j; k < K; k++) {
                XtX[j * K + k] += wxj * X[k * N + i];
            }
        }
    }

    /* Mirror to lower triangle */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            XtX[k * K + j] = XtX[j * K + k];
        }
    }

    /* Cholesky decomposition of X'X */
    rc = cholesky_decompose(XtX, K);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Solve for coefficients */
    memcpy(coefs, Xty, K * sizeof(ST_double));
    rc = cholesky_solve(XtX, K, coefs);
    if (rc != CBINSCATTER_OK) goto cleanup;

    /* Compute R-squared - first compute y mean */
    for (i = 0; i < N; i++) {
        ST_double w = (weights != NULL) ? weights[i] : 1.0;
        y_mean += w * y[i];
        sum_w += w;
    }
    y_mean /= sum_w;

    /* Compute TSS and RSS */
    for (i = 0; i < N; i++) {
        ST_double w = (weights != NULL) ? weights[i] : 1.0;

        /* Compute fitted value using Horner's method for better numerical stability */
        ST_double xi = x[i];
        ST_double fitted = coefs[K - 1];
        for (j = K - 2; j >= 0; j--) {
            fitted = fitted * xi + coefs[j];
        }

        ST_double residual = y[i] - fitted;
        ST_double dev = y[i] - y_mean;
        rss += w * residual * residual;
        tss += w * dev * dev;
    }

    *r2 = (tss > 0) ? (1.0 - rss / tss) : 0.0;

cleanup:
    free(X);
    free(XtX);
    free(Xty);

    return rc;
}

/* ========================================================================
 * Fit All Groups
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
    ST_int g, i, n_group;
    ST_retcode rc = CBINSCATTER_OK;
    ST_int order = config->linetype;  /* 1=linear, 2=quadratic, etc. */
    ST_double *y_group = NULL, *x_group = NULL, *w_group = NULL;
    ST_int max_coefs = order + 1;

    if (order < 1) return CBINSCATTER_OK;  /* No fitting requested */

    for (g = 0; g < results->num_by_groups; g++) {
        ByGroupResult *group = &results->groups[g];

        /* Count observations in this group */
        if (by_groups != NULL) {
            n_group = 0;
            for (i = 0; i < N; i++) {
                if (by_groups[i] == g + 1) n_group++;
            }
        } else {
            n_group = N;
        }

        if (n_group < order + 1) {
            /* Not enough observations for this polynomial order */
            group->fit_order = 0;
            group->fit_coefs = NULL;
            continue;
        }

        /* Allocate group data */
        y_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        x_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        if (weights != NULL) {
            w_group = (ST_double *)malloc(n_group * sizeof(ST_double));
        }

        if (!y_group || !x_group || (weights != NULL && !w_group)) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto group_cleanup;
        }

        /* Extract group data */
        ST_int j = 0;
        for (i = 0; i < N; i++) {
            if (by_groups == NULL || by_groups[i] == g + 1) {
                y_group[j] = y[i];
                x_group[j] = x[i];
                if (w_group != NULL) {
                    w_group[j] = weights[i];
                }
                j++;
            }
        }

        /* Allocate coefficients */
        group->fit_coefs = (ST_double *)calloc(max_coefs, sizeof(ST_double));
        if (!group->fit_coefs) {
            rc = CBINSCATTER_ERR_MEMORY;
            goto group_cleanup;
        }

        /* Fit polynomial */
        rc = fit_polynomial(y_group, x_group, w_group, n_group, order,
                            group->fit_coefs, &group->fit_r2);

        if (rc == CBINSCATTER_OK) {
            group->fit_order = order;
            group->fit_n = n_group;
        } else {
            /* Fitting failed - clear coefficients */
            free(group->fit_coefs);
            group->fit_coefs = NULL;
            group->fit_order = 0;
            rc = CBINSCATTER_OK;  /* Don't fail entire operation */
        }

group_cleanup:
        free(y_group);
        free(x_group);
        free(w_group);
        y_group = x_group = w_group = NULL;

        if (rc != CBINSCATTER_OK) break;
    }

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
