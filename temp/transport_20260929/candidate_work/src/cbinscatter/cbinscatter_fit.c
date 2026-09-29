/*
 * cbinscatter_fit.c
 *
 * Line fitting implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 *
 * Optimizations:
 * - No design matrix allocation unless fixed effects must be absorbed
 * - Powers of centered, scaled x keep the normal equations well conditioned
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stplugin.h"
#include "cbinscatter_fit.h"
#include "cbinscatter_resid.h"
#include "../ctools_config.h"
#include "../ctools_ols.h"

/* ========================================================================
 * Polynomial Fit, Optionally Jointly with Controls and Absorbed Effects
 *
 * Weighted least squares of y on 1, d, ..., d^p with d = (x - x̄)/s, plus
 * the controls centered at their weighted means. Redundant terms are
 * omitted (coefficient 0), as regress does. Centering the controls
 * evaluates the line at their means, binsreg's polyreg() convention. With
 * absorbed effects, all columns are first demeaned within FE levels (FWL)
 * and the intercept is ȳ - Σ b_j mean(d^j), matching reghdfe's _cons.
 * The coefficients are returned in powers of x.
 * ======================================================================== */

/* Rewrite Σ_j a[j] ((x - center)/scale)^j as Σ_k coefs[k] x^k */
static void poly_to_x_basis(const ST_double *a, ST_int p, ST_double center,
                            ST_double scale, ST_double *coefs)
{
    for (ST_int k = 0; k <= p; k++) coefs[k] = 0.0;
    for (ST_int j = 0; j <= p; j++) {
        ST_double aj = a[j], binom = 1.0, cpow = 1.0;
        for (ST_int t = 0; t < j; t++) aj /= scale;
        /* (x - c)^j = Σ_k C(j,k) x^k (-c)^(j-k), from k = j down */
        for (ST_int k = j; k >= 0; k--) {
            coefs[k] += aj * binom * cpow;
            binom = binom * k / (j - k + 1);
            cpow *= -center;
        }
    }
}

static ST_retcode fit_polynomial_joint(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    ST_int order,
    const ST_double *controls,   /* K x N column-major, or NULL */
    ST_int K,
    const ST_int *fe_vars,       /* G x N column-major levels 1..L, or NULL */
    ST_int G,
    const BinscatterConfig *config,
    ST_double *coefs,            /* order + 1 */
    ST_double *r2
) {
    ST_retcode rc = CBINSCATTER_OK;
    ST_double sw = 0.0, sx = 0.0, sy = 0.0;
    ST_double xmin = HUGE_VAL, xmax = -HUGE_VAL;
    ST_double *wbar = NULL, *XtX = NULL, *Xty = NULL, *b = NULL, *v = NULL;
    ST_double *cols = NULL, **vars = NULL, *dmean = NULL;
    ST_double a[101];
    ST_int i, j, k;

    if (order < 1 || order > 100 || K < 0 || (G > 0 && (!fe_vars || !config))) {
        return CBINSCATTER_ERR_INVALID;
    }
    if (!controls) K = 0;
    for (j = 0; j <= order; j++) coefs[j] = 0.0;
    *r2 = 0.0;

    wbar = (ST_double *)calloc((size_t)K + 1, sizeof(ST_double));
    if (!wbar) return CBINSCATTER_ERR_MEMORY;

    /* Pass 1: weighted means and the range of x */
    for (i = 0; i < N; i++) {
        ST_double w = weights ? weights[i] : 1.0;
        sw += w;
        sx += w * x[i];
        sy += w * y[i];
        for (k = 0; k < K; k++) wbar[k] += w * controls[(size_t)k * N + i];
        if (x[i] < xmin) xmin = x[i];
        if (x[i] > xmax) xmax = x[i];
    }
    if (!(sw > 0.0)) goto cleanup;
    ST_double xbar = sx / sw, ybar = sy / sw;
    for (k = 0; k < K; k++) wbar[k] /= sw;

    /* A constant x leaves only the intercept (regress omits x) */
    ST_int p = (xmax > xmin) ? order : 0;
    ST_double scale = p ? fmax(xmax - xbar, xbar - xmin) : 1.0;

    ST_int m = 1 + p + K;                 /* constant, d..d^p, controls */
    ST_int q = (G > 0) ? m - 1 : m;       /* absorbed FE replace the constant */
    ST_int off = (G > 0) ? 1 : 0;         /* column of d^1 within the regressors */
    XtX = (ST_double *)calloc((size_t)q * q + 1, sizeof(ST_double));
    Xty = (ST_double *)calloc((size_t)q + 1, sizeof(ST_double));
    b = (ST_double *)calloc((size_t)q + 1, sizeof(ST_double));
    v = (ST_double *)calloc((size_t)m, sizeof(ST_double));
    if (!XtX || !Xty || !b || !v) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }

    ST_double tss = 0.0, rss = 0.0;
    if (G > 0) {
        /* Columns [y, d^1..d^p, controls], demeaned within FE levels together */
        cols = (ST_double *)ctools_safe_malloc3((size_t)m, (size_t)N, sizeof(ST_double));
        vars = (ST_double **)malloc((size_t)m * sizeof(ST_double *));
        dmean = (ST_double *)calloc((size_t)p + 1, sizeof(ST_double));
        if (!cols || !vars || !dmean) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }
        for (j = 0; j < m; j++) vars[j] = &cols[(size_t)j * N];
        for (i = 0; i < N; i++) {
            ST_double w = weights ? weights[i] : 1.0;
            ST_double d = (x[i] - xbar) / scale, pw = 1.0;
            vars[0][i] = y[i];
            tss += w * (y[i] - ybar) * (y[i] - ybar);
            for (j = 1; j <= p; j++) {
                pw *= d;
                vars[j][i] = pw;
                dmean[j] += w * pw;
            }
            for (k = 0; k < K; k++) vars[1 + p + k][i] = controls[(size_t)k * N + i];
        }
        for (j = 1; j <= p; j++) dmean[j] /= sw;

        ST_int dropped = 0;
        rc = hdfe_residualize_batch(vars, m, fe_vars, N, G, weights, config->weight_type,
                                    config->maxiter, config->tolerance, 0, &dropped);
        if (rc != CBINSCATTER_OK) goto cleanup;

        for (i = 0; i < N; i++) {
            ST_double w = weights ? weights[i] : 1.0;
            for (j = 0; j < q; j++) {
                ST_double wv = w * vars[1 + j][i];
                Xty[j] += wv * vars[0][i];
                for (k = j; k < q; k++) XtX[(size_t)j * q + k] += wv * vars[1 + k][i];
            }
        }
    } else {
        /* Streamed: [1, d..d^p, W - w̄] on y - ȳ */
        for (i = 0; i < N; i++) {
            ST_double w = weights ? weights[i] : 1.0;
            ST_double d = (x[i] - xbar) / scale, r = y[i] - ybar;
            v[0] = 1.0;
            for (j = 1; j <= p; j++) v[j] = v[j - 1] * d;
            for (k = 0; k < K; k++) v[1 + p + k] = controls[(size_t)k * N + i] - wbar[k];
            tss += w * r * r;
            for (j = 0; j < q; j++) {
                ST_double wv = w * v[j];
                Xty[j] += wv * r;
                for (k = j; k < q; k++) XtX[(size_t)j * q + k] += wv * v[k];
            }
        }
    }

    if (q > 0) {
        for (j = 0; j < q; j++)
            for (k = j + 1; k < q; k++) XtX[(size_t)k * q + j] = XtX[(size_t)j * q + k];
        rc = cbinscatter_solve_with_collinearity(XtX, Xty, q, b);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* RSS from the same (demeaned or centered) columns */
    for (i = 0; i < N; i++) {
        ST_double w = weights ? weights[i] : 1.0, e;
        if (G > 0) {
            e = vars[0][i];
            for (j = 0; j < q; j++) e -= b[j] * vars[1 + j][i];
        } else {
            ST_double d = (x[i] - xbar) / scale;
            v[0] = 1.0;
            for (j = 1; j <= p; j++) v[j] = v[j - 1] * d;
            for (k = 0; k < K; k++) v[1 + p + k] = controls[(size_t)k * N + i] - wbar[k];
            e = y[i] - ybar;
            for (j = 0; j < q; j++) e -= b[j] * v[j];
        }
        rss += w * e * e;
    }
    if (tss > 0.0) {
        *r2 = 1.0 - rss / tss;
        if (*r2 < 0.0) *r2 = 0.0;
        if (*r2 > 1.0) *r2 = 1.0;
    }

    /* Polynomial in d with the controls at their means, then in powers of x */
    a[0] = ybar + (G > 0 ? 0.0 : b[0]);
    for (j = 1; j <= p; j++) {
        a[j] = b[j - off];
        if (G > 0) a[0] -= a[j] * dmean[j];
    }
    poly_to_x_basis(a, p, xbar, scale, coefs);

cleanup:
    free(wbar); free(XtX); free(Xty); free(b); free(v);
    free(cols); free(vars); free(dmean);
    return rc;
}

/* ========================================================================
 * Fit One Group
 * ======================================================================== */

ST_retcode fit_group_polynomial(
    const ST_double *y,
    const ST_double *x,
    const ST_double *weights,
    ST_int N,
    const ST_double *controls,
    ST_int K,
    const ST_int *fe_vars,
    ST_int G,
    const BinscatterConfig *config,
    ByGroupResult *group
) {
    ST_int order = config->linetype;

    group->fit_order = 0;
    group->fit_coefs = NULL;
    if (order < 1 || N < order + 1) return CBINSCATTER_OK;

    group->fit_coefs = (ST_double *)calloc((size_t)order + 1, sizeof(ST_double));
    if (!group->fit_coefs) return CBINSCATTER_ERR_MEMORY;

    ST_retcode rc = fit_polynomial_joint(y, x, weights, N, order, controls, K,
                                         fe_vars, G, config, group->fit_coefs,
                                         &group->fit_r2);
    if (rc == CBINSCATTER_OK) {
        group->fit_order = order;
        group->fit_n = N;
    } else {
        free(group->fit_coefs);
        group->fit_coefs = NULL;
    }
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

    if (order < 1) return CBINSCATTER_OK;

    /* For single group (no by()), fit directly without copying */
    if (results->num_by_groups == 1 && by_groups == NULL) {
        return fit_group_polynomial(y, x, weights, N, NULL, 0, NULL, 0,
                                    config, &results->groups[0]);
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

        rc = fit_group_polynomial(y_group, x_group, w_group, n_group, NULL, 0,
                                  NULL, 0, config, group);
        if (rc == CBINSCATTER_ERR_MEMORY) goto cleanup;
        rc = CBINSCATTER_OK;  /* Don't fail entire operation */
    }

cleanup:
    free(y_group);
    free(x_group);
    free(w_group);
    free(group_counts);
    return rc;
}
