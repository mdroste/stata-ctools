/*
 * cqreg_irls.c
 *
 * Iteratively Reweighted Least Squares solver for quantile regression.
 * This is a simpler, more robust alternative to IPM.
 * Part of the ctools suite.
 */

#include "cqreg_ipm.h"  /* Reuse the state structure */
#include "cqreg_linalg.h"
#include "../stplugin.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* Minimum weight to avoid division by zero */
#define IRLS_MIN_WEIGHT 1e-8

/* Default epsilon for residuals near zero */
#define IRLS_EPSILON 1e-6

/*
 * Solve weighted least squares: (X'WX)^(-1) X'Wy
 * Returns 0 on success, -1 on failure.
 */
static ST_int solve_weighted_ls(ST_int N, ST_int K,
                                const ST_double *X,   /* N x K column-major */
                                const ST_double *y,
                                const ST_double *w,   /* N weights */
                                ST_double *beta,
                                ST_double *XWX,       /* K x K workspace */
                                ST_double *XWy,       /* K workspace */
                                ST_double *L)         /* K x K workspace for Cholesky */
{
    ST_int i, j, k;

    /* Compute X'WX */
    memset(XWX, 0, K * K * sizeof(ST_double));
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        /* Diagonal */
        ST_double diag = 0.0;
        for (i = 0; i < N; i++) {
            diag += w[i] * Xj[i] * Xj[i];
        }
        XWX[j * K + j] = diag;
        /* Off-diagonal (upper triangle) */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += w[i] * Xj[i] * Xk[i];
            }
            XWX[j * K + k] = sum;
            XWX[k * K + j] = sum;
        }
    }

    /* Compute X'Wy */
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        ST_double sum = 0.0;
        for (i = 0; i < N; i++) {
            sum += w[i] * Xj[i] * y[i];
        }
        XWy[j] = sum;
    }

    /* Copy XWX to L for Cholesky */
    memcpy(L, XWX, K * K * sizeof(ST_double));

    /* Add small regularization for numerical stability */
    for (j = 0; j < K; j++) {
        L[j * K + j] += 1e-10;
    }

    /* Cholesky factorization */
    if (cqreg_cholesky(L, K) != 0) {
        return -1;
    }

    /* Solve L*L'*beta = XWy */
    memcpy(beta, XWy, K * sizeof(ST_double));
    cqreg_solve_cholesky(L, beta, K);

    return 0;
}

/*
 * IRLS solver for quantile regression.
 * Much simpler than IPM, works well for small/medium datasets.
 */
ST_int cqreg_irls_solve(cqreg_ipm_state *ipm,
                        const ST_double *y,
                        const ST_double *X,
                        ST_double q,
                        ST_double *beta)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int iter, i;

    ST_double *w = ipm->D;        /* Reuse D as weights */
    ST_double *r = ipm->r_primal; /* Reuse for residuals */
    ST_double *XWX = ipm->XDX;    /* Reuse workspace */
    ST_double *XWy = ipm->work_K; /* Workspace */
    ST_double *L = ipm->L;        /* Cholesky factor */


    /* Initialize beta to zero */
    memset(beta, 0, K * sizeof(ST_double));

    /* First iteration: equal weights (OLS-like) */
    for (i = 0; i < N; i++) {
        w[i] = 1.0;
    }

    /* Initial solve */
    if (solve_weighted_ls(N, K, X, y, w, beta, XWX, XWy, L) != 0) {
        return -1;
    }

    /* IRLS iterations */
    ST_double prev_obj = 1e30;
    ST_double epsilon = IRLS_EPSILON * cqreg_vmax_abs(y, N);  /* Scale epsilon by y */
    if (epsilon < IRLS_EPSILON) epsilon = IRLS_EPSILON;

    for (iter = 0; iter < ipm->config.maxiter; iter++) {
        /* Compute residuals: r = y - X*beta */
        cqreg_matvec_col(r, X, beta, N, K);
        for (i = 0; i < N; i++) {
            r[i] = y[i] - r[i];
        }

        /* Compute objective: sum of check function */
        ST_double obj = 0.0;
        for (i = 0; i < N; i++) {
            if (r[i] >= 0) {
                obj += q * r[i];
            } else {
                obj += (q - 1.0) * r[i];
            }
        }

        if (iter == 0 || iter % 10 == 0) {
        }

        /* Check convergence (use more relaxed tolerance for IRLS) */
        ST_double rel_change = fabs(prev_obj - obj) / (1.0 + fabs(obj));
        if (rel_change < 1e-6) {  /* More relaxed than IPM tolerance */
            ipm->converged = 1;
            ipm->iterations = iter + 1;
            break;
        }
        prev_obj = obj;

        /* Update weights for quantile regression
         * w[i] = q / |r[i]| if r[i] > 0
         * w[i] = (1-q) / |r[i]| if r[i] < 0
         * with smoothing near zero
         */
        for (i = 0; i < N; i++) {
            ST_double abs_r = fabs(r[i]);
            if (abs_r < epsilon) abs_r = epsilon;

            if (r[i] >= 0) {
                w[i] = q / abs_r;
            } else {
                w[i] = (1.0 - q) / abs_r;
            }

            /* Clamp weights */
            if (w[i] < IRLS_MIN_WEIGHT) w[i] = IRLS_MIN_WEIGHT;
            if (w[i] > 1.0 / IRLS_MIN_WEIGHT) w[i] = 1.0 / IRLS_MIN_WEIGHT;
        }

        /* Solve weighted least squares */
        if (solve_weighted_ls(N, K, X, y, w, beta, XWX, XWy, L) != 0) {
            return -(iter + 1);
        }

    }

    ipm->iterations = iter + 1;
    if (!ipm->converged) {
        ipm->converged = 0;
    }

    /* Store final residuals in u (for compatibility) */
    cqreg_matvec_col(r, X, beta, N, K);
    for (i = 0; i < N; i++) {
        r[i] = y[i] - r[i];
        /* Store in u (positive residuals) and v (negative residuals) */
        if (r[i] >= 0) {
            ipm->u[i] = r[i];
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -r[i];
        }
    }

    /* Copy solution to ipm->beta */
    cqreg_vcopy(ipm->beta, beta, K);

    return ipm->converged ? ipm->iterations : -ipm->iterations;
}
