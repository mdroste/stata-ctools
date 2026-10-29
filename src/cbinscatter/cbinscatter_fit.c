/*
 * cbinscatter_fit.c
 *
 * Line fitting implementation for cbinscatter
 * Part of the ctools Stata plugin suite
 *
 * Optimizations:
 * - No design matrix allocation unless fixed effects must be absorbed
 * - Powers of centered, scaled x keep the normal equations well conditioned
 * - Passes over the rows run in blocks in parallel; each block's sums are
 *   added in block order, so results do not depend on the thread count
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
    ST_int j, k;

    if (order < 1 || order > 100 || K < 0 || (G > 0 && (!fe_vars || !config))) {
        return CBINSCATTER_ERR_INVALID;
    }
    if (!controls) K = 0;
    for (j = 0; j <= order; j++) coefs[j] = 0.0;
    *r2 = 0.0;

    const ST_int num_blocks = cbs_num_blocks(N);
    ST_double *blk_sums = NULL;
    ST_int blk;

    wbar = (ST_double *)calloc((size_t)K + 1, sizeof(ST_double));
    if (!wbar) return CBINSCATTER_ERR_MEMORY;

    /* Pass 1: weighted means and the range of x. Per block: sw, sx, sy,
     * xmin, xmax, then the K control sums. */
    {
        const size_t per = 5 + (size_t)K;
        blk_sums = (ST_double *)ctools_safe_calloc2((size_t)num_blocks * per + 1, sizeof(ST_double));
        if (!blk_sums) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }
        #pragma omp parallel for schedule(static) if (num_blocks > 1)
        for (blk = 0; blk < num_blocks; blk++) {
            ST_int i0, i1;
            cbs_block_range(blk, N, &i0, &i1);
            ST_double *acc = blk_sums + (size_t)blk * per;
            ST_double bsw = 0.0, bsx = 0.0, bsy = 0.0, bmin = HUGE_VAL, bmax = -HUGE_VAL;
            for (ST_int ii = i0; ii < i1; ii++) {
                ST_double w = weights ? weights[ii] : 1.0;
                bsw += w;
                bsx += w * x[ii];
                bsy += w * y[ii];
                for (ST_int kk = 0; kk < K; kk++) acc[5 + kk] += w * controls[(size_t)kk * N + ii];
                if (x[ii] < bmin) bmin = x[ii];
                if (x[ii] > bmax) bmax = x[ii];
            }
            acc[0] = bsw; acc[1] = bsx; acc[2] = bsy; acc[3] = bmin; acc[4] = bmax;
        }
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_double *acc = blk_sums + (size_t)blk * per;
            sw += acc[0];
            sx += acc[1];
            sy += acc[2];
            if (acc[3] < xmin) xmin = acc[3];
            if (acc[4] > xmax) xmax = acc[4];
            for (k = 0; k < K; k++) wbar[k] += acc[5 + k];
        }
        free(blk_sums);
        blk_sums = NULL;
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
    v = (ST_double *)calloc((size_t)num_blocks * m, sizeof(ST_double));
    if (!XtX || !Xty || !b || !v) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }

    ST_double tss = 0.0, rss = 0.0;
    if (G > 0) {
        /* Columns [y, d^1..d^p, controls], demeaned within FE levels together.
         * Per block: tss, then the sums of d^1..d^p. */
        const size_t per = 1 + (size_t)p;
        cols = (ST_double *)ctools_safe_malloc3((size_t)m, (size_t)N, sizeof(ST_double));
        vars = (ST_double **)malloc((size_t)m * sizeof(ST_double *));
        dmean = (ST_double *)calloc((size_t)p + 1, sizeof(ST_double));
        blk_sums = (ST_double *)ctools_safe_calloc2((size_t)num_blocks * per + 1, sizeof(ST_double));
        if (!cols || !vars || !dmean || !blk_sums) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }
        for (j = 0; j < m; j++) vars[j] = &cols[(size_t)j * N];
        #pragma omp parallel for schedule(static) if (num_blocks > 1)
        for (blk = 0; blk < num_blocks; blk++) {
            ST_int i0, i1;
            cbs_block_range(blk, N, &i0, &i1);
            ST_double *acc = blk_sums + (size_t)blk * per;
            for (ST_int ii = i0; ii < i1; ii++) {
                ST_double w = weights ? weights[ii] : 1.0;
                ST_double d = (x[ii] - xbar) / scale, pw = 1.0;
                vars[0][ii] = y[ii];
                acc[0] += w * (y[ii] - ybar) * (y[ii] - ybar);
                for (ST_int jj = 1; jj <= p; jj++) {
                    pw *= d;
                    vars[jj][ii] = pw;
                    acc[jj] += w * pw;
                }
                for (ST_int kk = 0; kk < K; kk++) vars[1 + p + kk][ii] = controls[(size_t)kk * N + ii];
            }
        }
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_double *acc = blk_sums + (size_t)blk * per;
            tss += acc[0];
            for (j = 1; j <= p; j++) dmean[j] += acc[j];
        }
        free(blk_sums);
        blk_sums = NULL;
        for (j = 1; j <= p; j++) dmean[j] /= sw;

        ST_int dropped = 0;
        rc = hdfe_residualize_batch(vars, m, fe_vars, N, G, weights, config->weight_type,
                                    config->maxiter, config->tolerance, 0, &dropped);
        if (rc != CBINSCATTER_OK) goto cleanup;

        cbinscatter_normal_equations((const ST_double *const *)(vars + 1), q, vars[0], weights,
                                     NULL, 0, N, XtX, Xty);
    } else {
        /* Streamed: [1, d..d^p, W - w̄] on y - ȳ. Per block: tss, Xty, then
         * the upper triangle of XtX. */
        const size_t per = 1 + (size_t)q + (size_t)q * q;
        blk_sums = (ST_double *)ctools_safe_calloc2((size_t)num_blocks * per + 1, sizeof(ST_double));
        if (!blk_sums) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }
        #pragma omp parallel for schedule(static) if (num_blocks > 1)
        for (blk = 0; blk < num_blocks; blk++) {
            ST_int i0, i1;
            cbs_block_range(blk, N, &i0, &i1);
            ST_double *acc = blk_sums + (size_t)blk * per, *bXty = acc + 1, *bXtX = bXty + q;
            ST_double *vb = v + (size_t)blk * m;
            for (ST_int ii = i0; ii < i1; ii++) {
                ST_double w = weights ? weights[ii] : 1.0;
                ST_double d = (x[ii] - xbar) / scale, r = y[ii] - ybar;
                vb[0] = 1.0;
                for (ST_int jj = 1; jj <= p; jj++) vb[jj] = vb[jj - 1] * d;
                for (ST_int kk = 0; kk < K; kk++) vb[1 + p + kk] = controls[(size_t)kk * N + ii] - wbar[kk];
                acc[0] += w * r * r;
                for (ST_int jj = 0; jj < q; jj++) {
                    ST_double wv = w * vb[jj];
                    bXty[jj] += wv * r;
                    for (ST_int kk = jj; kk < q; kk++) bXtX[(size_t)jj * q + kk] += wv * vb[kk];
                }
            }
        }
        for (blk = 0; blk < num_blocks; blk++) {
            const ST_double *acc = blk_sums + (size_t)blk * per, *bXty = acc + 1, *bXtX = bXty + q;
            tss += acc[0];
            for (j = 0; j < q; j++) {
                Xty[j] += bXty[j];
                for (k = j; k < q; k++) XtX[(size_t)j * q + k] += bXtX[(size_t)j * q + k];
            }
        }
        free(blk_sums);
        blk_sums = NULL;
    }

    if (q > 0) {
        for (j = 0; j < q; j++)
            for (k = j + 1; k < q; k++) XtX[(size_t)k * q + j] = XtX[(size_t)j * q + k];
        rc = cbinscatter_solve_with_collinearity(XtX, Xty, q, b);
        if (rc != CBINSCATTER_OK) goto cleanup;
    }

    /* RSS from the same (demeaned or centered) columns, per block */
    blk_sums = (ST_double *)ctools_safe_calloc2((size_t)num_blocks + 1, sizeof(ST_double));
    if (!blk_sums) { rc = CBINSCATTER_ERR_MEMORY; goto cleanup; }
    #pragma omp parallel for schedule(static) if (num_blocks > 1)
    for (blk = 0; blk < num_blocks; blk++) {
        ST_int i0, i1;
        cbs_block_range(blk, N, &i0, &i1);
        ST_double *vb = v + (size_t)blk * m, brss = 0.0;
        for (ST_int ii = i0; ii < i1; ii++) {
            ST_double w = weights ? weights[ii] : 1.0, e;
            if (G > 0) {
                e = vars[0][ii];
                for (ST_int jj = 0; jj < q; jj++) e -= b[jj] * vars[1 + jj][ii];
            } else {
                ST_double d = (x[ii] - xbar) / scale;
                vb[0] = 1.0;
                for (ST_int jj = 1; jj <= p; jj++) vb[jj] = vb[jj - 1] * d;
                for (ST_int kk = 0; kk < K; kk++) vb[1 + p + kk] = controls[(size_t)kk * N + ii] - wbar[kk];
                e = y[ii] - ybar;
                for (ST_int jj = 0; jj < q; jj++) e -= b[jj] * vb[jj];
            }
            brss += w * e * e;
        }
        blk_sums[blk] = brss;
    }
    for (blk = 0; blk < num_blocks; blk++) rss += blk_sums[blk];
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
    free(cols); free(vars); free(dmean); free(blk_sums);
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
    const ST_int *group_start,
    const ST_double *weights,
    ST_int N,
    const BinscatterConfig *config,
    BinscatterResults *results
) {
    ST_int order = config->linetype;

    if (order < 1) return CBINSCATTER_OK;

    /* For single group (no by()), fit directly without copying */
    if (results->num_by_groups == 1 && group_start == NULL) {
        return fit_group_polynomial(y, x, weights, N, NULL, 0, NULL, 0,
                                    config, &results->groups[0]);
    }

    /* Each group is a contiguous range of rows */
    for (ST_int g = 0; g < results->num_by_groups; g++) {
        ByGroupResult *group = &results->groups[g];
        ST_int first = group_start[g], n_group = group_start[g + 1] - first;

        if (n_group < order + 1) {
            group->fit_order = 0;
            group->fit_coefs = NULL;
            continue;
        }

        ST_retcode rc = fit_group_polynomial(y + first, x + first, weights ? weights + first : NULL,
                                             n_group, NULL, 0, NULL, 0, config, group);
        if (rc == CBINSCATTER_ERR_MEMORY) return rc;
        /* Don't fail entire operation */
    }
    return CBINSCATTER_OK;
}
