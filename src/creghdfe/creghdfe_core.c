/*
 * creghdfe_core.c
 *
 * Collinearity rules and OLS on partialled data (see creghdfe_core.h). Moved
 * from do_full_regression's STEP 8 without changing the arithmetic.
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_core.h"

#include "../ctools_ols.h"
#include "../ctools_types.h"

ST_retcode creghdfe_collinear(ST_double *data, ST_int N, ST_int K, const ST_double *norm2,
                              const ST_double *count_w, const ST_double *weights, ST_int weight_type,
                              ST_double collinear_tol, ST_int *fe_flagged,
                              ST_int *is_collinear, ST_int *depvar_explained, ST_double *xtx,
                              ST_double *xty_all, ST_double *scratch)
{
    const ST_int K_x = K - 1;
    ST_int k, j;
    if (K_x > 0)
    {
        memset(is_collinear, 0, (size_t)K_x * sizeof(ST_int));
        memset(fe_flagged, 0, (size_t)K_x * sizeof(ST_int));
    }

    /* Columns collinear with the fixed effects, by reghdfe's rule
     * (check_collinear_with_fe): the partialled sum of squares is at most
     * collinear_tol = min(1e-6, tolerance()/10) times the uncentered sum of
     * squares before partialling (both unweighted, so the rule does not
     * depend on the units or on standardization). Such a regressor is
     * omitted. A dependent variable that passes the rule is perfectly
     * explained by the effects; it is zeroed, as reghdfe does, so its
     * coefficients, rss and tss_within are exactly zero. */
    *depvar_explained = 0;
    ST_int num_fe_collinear = 0;
    for (k = 0; k < K; k++)
    {
        ST_double *col = &data[(size_t)k * N];
        ST_double ss;
        if (count_w)
        {
            ss = 0.0;
            for (ST_int i = 0; i < N; i++)
                ss += count_w[i] * col[i] * col[i];
        }
        else
            ss = fast_dot(col, col, N);
        ST_double ratio = ss / norm2[k];
        ST_int flagged = (ratio <= collinear_tol);
        if (k == 0)
        {
            if (flagged)
            {
                memset(col, 0, (size_t)N * sizeof(ST_double));
                *depvar_explained = 1;
            }
            continue;
        }
        if (flagged)
        {
            is_collinear[k - 1] = 1;
            num_fe_collinear++;
        }
        fe_flagged[k - 1] = flagged;
    }

    /* The remaining collinearity as reghdfe_rmcoll: Mata's invsym() of the
     * (weighted) cross products of the partialled regressors, sweeping them
     * in order, and invsym() pivoting on the largest diagonal instead when
     * that drops a different number of them. A regressor is collinear when
     * its diagonal of the generalized inverse is zero. */
    if (K_x > 0)
    {
        if (weights)
            compute_xtx_xty_weighted(data, weights, weight_type, N, K, xtx, xty_all);
        else
            compute_xtx_xty(data, N, K, xtx, xty_all);
        ST_double *xx_rm = xtx;
        if (num_fe_collinear > 0)
        {
            /* Those columns leave the cross products, as reghdfe trims them */
            xx_rm = scratch;
            memcpy(xx_rm, xtx, (size_t)K_x * K_x * sizeof(ST_double));
            for (k = 0; k < K_x; k++)
            {
                if (!is_collinear[k])
                    continue;
                for (j = 0; j < K_x; j++)
                {
                    xx_rm[(size_t)k * K_x + j] = 0.0;
                    xx_rm[(size_t)j * K_x + k] = 0.0;
                }
            }
        }
        ST_double *ginv = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
        if (!ginv)
            return 920;
        ST_int rank_ordered = 0, rank_pivoted = 0;
        ST_int inv_rc = ctools_invsym_ordered(xx_rm, K_x, ginv, &rank_ordered);
        if (!inv_rc && rank_ordered < K_x - num_fe_collinear)
        {
            ST_double *ginv_alt = (ST_double *)malloc((size_t)K_x * K_x * sizeof(ST_double));
            if (!ginv_alt)
            {
                free(ginv);
                return 920;
            }
            inv_rc = ctools_invsym(xx_rm, K_x, ginv_alt, &rank_pivoted);
            if (!inv_rc && rank_pivoted != rank_ordered)
            {
                free(ginv);
                ginv = ginv_alt;
            }
            else
            {
                free(ginv_alt);
            }
        }
        if (inv_rc)
        {
            free(ginv);
            return 920;
        }
        for (k = 0; k < K_x; k++)
        {
            if (ginv[(size_t)k * K_x + k] == 0.0)
                is_collinear[k] = 1;
        }
        free(ginv);
    }
    return 0;
}

ST_retcode creghdfe_ols(const creghdfe_ols_in *in, creghdfe_ols_out *out)
{
    const ST_int N = in->N, K_x = in->K - 1;
    const ST_double *weights = in->weights;
    ST_int i, j, k, idx;
    memset(out, 0, sizeof(*out));

    ST_int num_collinear = 0;
    for (k = 0; k < K_x; k++)
    {
        if (in->is_collinear[k])
            num_collinear++;
    }
    const ST_int K_keep = K_x - num_collinear;
    out->K_keep = K_keep;

    /* K_keep == 0 (every regressor omitted) continues through the same path,
     * as reghdfe does: the model is then the fixed effects and the constant.
     * K_alloc keeps the K_keep-sized allocations nonempty in that case. */
    const size_t K_alloc = (K_keep > 0) ? (size_t)K_keep : 1;

    /* Build index of non-collinear variables */
    out->keep_idx = (ST_int *)malloc(K_alloc * sizeof(ST_int));
    if (!out->keep_idx)
        return 920;
    idx = 0;
    for (k = 0; k < K_x; k++)
    {
        if (!in->is_collinear[k])
            out->keep_idx[idx++] = k + 1; /* +1 to skip y */
    }

    const ST_int K_with_cons = K_keep + 1;
    /* y and the kept regressors. When nothing is omitted this is the data
     * matrix itself and no copy is made. */
    const int keep_is_data = (K_keep == K_x);
    out->data_keep = keep_is_data
                         ? in->data
                         : (ST_double *)malloc((size_t)N * (K_keep + 1) * sizeof(ST_double));
    out->xtx_keep = (ST_double *)malloc(K_alloc * K_alloc * sizeof(ST_double));
    out->xty_keep = (ST_double *)malloc(K_alloc * sizeof(ST_double));
    out->beta_keep = (ST_double *)malloc(K_with_cons * sizeof(ST_double));
    out->inv_xx_keep = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));
    out->means_x = (ST_double *)malloc(K_alloc * sizeof(ST_double));
    if (!out->data_keep || !out->xtx_keep || !out->xty_keep || !out->beta_keep ||
        !out->inv_xx_keep || !out->means_x)
        return 920;
    ST_double *data_keep = out->data_keep, *xtx_keep = out->xtx_keep;
    ST_double *xty_keep = out->xty_keep, *beta_keep = out->beta_keep;
    ST_double *means_x = out->means_x;
    const ST_int *keep_idx = out->keep_idx;

    /* Copy y and non-collinear X */
    if (!keep_is_data)
    {
        memcpy(data_keep, in->data, (size_t)N * sizeof(ST_double));
        for (k = 0; k < K_keep; k++)
        {
            memcpy(data_keep + (size_t)(k + 1) * N, in->data + (size_t)keep_idx[k] * N,
                   (size_t)N * sizeof(ST_double));
        }
    }
    for (k = 0; k < K_keep; k++)
    {
        if (in->standardize)
            means_x[k] = in->means[keep_idx[k]] / in->stdevs[keep_idx[k]];
        else
            means_x[k] = in->means[keep_idx[k]];
    }

    /* Compute X'X and X'y on non-collinear data (weighted if using weights);
     * when nothing is omitted they are the cross products computed before */
    ST_double tss_within;
    const int reuse_xx = keep_is_data && K_keep > 0;
    if (reuse_xx)
    {
        memcpy(xtx_keep, in->xtx, (size_t)K_keep * K_keep * sizeof(ST_double));
        memcpy(xty_keep, in->xty_all, (size_t)K_keep * sizeof(ST_double));
    }
    if (weights)
    {
        if (!reuse_xx)
            compute_xtx_xty_weighted(data_keep, weights, in->weight_type, N, K_keep + 1, xtx_keep,
                                     xty_keep);
        /* Weighted TSS_within = sum(w_i * y_i^2)
         * Use quad-precision or Kahan summation based on option. */
        if (in->use_quad)
        {
            tss_within = dd_sum_sq_weighted(data_keep, weights, N);
        }
        else
        {
            tss_within = 0.0;
            for (idx = 0; idx < N; idx++)
            {
                volatile ST_double sq = data_keep[idx] * data_keep[idx];
                tss_within += weights[idx] * sq;
            }
        }
    }
    else
    {
        if (!reuse_xx)
            compute_xtx_xty(data_keep, N, K_keep + 1, xtx_keep, xty_keep);
        /* Compute tss_within: use quad-precision or Kahan summation. */
        if (in->use_quad)
        {
            tss_within = dd_sum_sq(data_keep, N);
        }
        else
        {
            tss_within = 0.0;
            for (idx = 0; idx < N; idx++)
            {
                volatile ST_double sq = data_keep[idx] * data_keep[idx];
                tss_within += sq;
            }
        }
    }
    out->tss_within = tss_within;

    /* Cholesky and invert */
    ST_double *inv_xx_x = (ST_double *)malloc(K_alloc * K_alloc * sizeof(ST_double));
    if (!inv_xx_x)
        return 920;
    if (K_keep > 0)
    {
        memcpy(inv_xx_x, xtx_keep, (size_t)K_keep * K_keep * sizeof(ST_double));
        if (ctools_cholesky(inv_xx_x, K_keep) != 0)
        {
            free(inv_xx_x);
            return 198;
        }
        stata_retcode inverse_status_0 = ctools_invert_from_cholesky(inv_xx_x, K_keep, inv_xx_x);
        if (inverse_status_0 != STATA_OK)
        {
            free(inv_xx_x);
            return ctools_stata_rc(inverse_status_0);
        }
    }

    /* Compute beta */
    for (i = 0; i < K_keep; i++)
    {
        beta_keep[i] = 0.0;
        for (j = 0; j < K_keep; j++)
            beta_keep[i] += inv_xx_x[i * K_keep + j] * xty_keep[j];
    }

    /* Compute RSS */
    ST_double rss = tss_within;
    for (k = 0; k < K_keep; k++)
        rss -= beta_keep[k] * xty_keep[k];
    out->rss = rss;

    /* Extend inv_xx with the constant (block partition formula). For
     * fweight the 1'W1 term is sum(weights). */
    ST_double N_corner = (in->weight_type == 2 && weights) ? in->sum_weights : (ST_double)N;
    ctools_inverse_add_constant(inv_xx_x, K_keep, means_x, N_corner, out->inv_xx_keep);

    /* Compute constant */
    ST_double y_mean = in->means[0];
    ST_double xb_mean = 0.0;
    for (k = 0; k < K_keep; k++)
        xb_mean += means_x[k] * beta_keep[k];
    beta_keep[K_keep] = y_mean - xb_mean;

    free(inv_xx_x);
    return 0;
}
