/*
 * cqreg_fn.c
 *
 * Frisch-Newton exact LP solver for quantile regression.
 * Based on Portnoy & Koenker (1997) and Chernozhukov et al. (2020).
 *
 * LP formulation:
 *   min  q * sum(u) + (1-q) * sum(v)
 *   s.t. y - X*beta = u - v,  u >= 0, v >= 0
 *
 * The Frisch-Newton algorithm solves this exactly by:
 * 1. Identifying the active set (observations on the hyperplane)
 * 2. Solving a reduced system on the active set
 * 3. Using ratio test to maintain feasibility
 * 4. Checking LP optimality conditions
 */

#include "cqreg_fn.h"
#include "cqreg_linalg.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * Identify Active Set
 * ============================================================================ */

ST_int fn_identify_active_set(const ST_double *r, ST_int N, ST_double tol,
                               ST_int *active_idx)
{
    ST_int m = 0;

    for (ST_int i = 0; i < N; i++) {
        if (fabs(r[i]) < tol) {
            active_idx[m++] = i;
        }
    }

    return m;
}

/* ============================================================================
 * Form Reduced Normal Equations
 *
 * For the active set, we solve:
 *   (X_A' W_A X_A) delta_beta = X_A' g_A
 *
 * where W_A is a diagonal weight matrix and g_A is the subgradient.
 * For LP-exact solution, W is uniform (we use identity).
 * ============================================================================ */

void fn_form_reduced_system(const ST_double *X, const ST_double *g,
                             ST_int N, ST_int K,
                             const ST_int *active_idx, ST_int m,
                             ST_double *XAX, ST_double *XAg)
{
    ST_int i, j, k;

    /* Initialize to zero */
    memset(XAX, 0, K * K * sizeof(ST_double));
    memset(XAg, 0, K * sizeof(ST_double));

    /* Compute X_A' X_A and X_A' g */
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];

        /* X_A' g (just column j) */
        ST_double xg = 0.0;
        for (i = 0; i < m; i++) {
            ST_int idx = active_idx[i];
            xg += Xj[idx] * g[idx];
        }
        XAg[j] = xg;

        /* X_A' X_A (diagonal and upper triangle) */
        for (k = j; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (i = 0; i < m; i++) {
                ST_int idx = active_idx[i];
                sum += Xj[idx] * Xk[idx];
            }
            XAX[j * K + k] = sum;
            XAX[k * K + j] = sum;  /* Symmetric */
        }
    }

    /* Add regularization to diagonal */
    for (j = 0; j < K; j++) {
        XAX[j * K + j] += FN_REG;
    }
}

/* ============================================================================
 * Ratio Test for Step Length
 *
 * Find maximum alpha such that:
 *   u_new = max(r_new, 0) >= 0
 *   v_new = max(-r_new, 0) >= 0
 *
 * Since r_new = r + alpha * delta_r, and u,v are implicit,
 * we need to ensure the step doesn't cause sign changes that violate
 * the LP constraints. For observations with r > 0 moving toward r < 0,
 * we limit the step.
 * ============================================================================ */

ST_double fn_ratio_test(const ST_double *r, const ST_double *delta_r,
                         ST_int N, ST_double q)
{
    ST_double alpha = 1.0;

    for (ST_int i = 0; i < N; i++) {
        ST_double ri = r[i];
        ST_double dri = delta_r[i];

        /* If residual is changing sign, limit step */
        if (ri > 0 && dri < -1e-14) {
            /* Moving from positive to negative */
            ST_double a = -ri / dri;
            if (a < alpha) alpha = a;
        } else if (ri < 0 && dri > 1e-14) {
            /* Moving from negative to positive */
            ST_double a = -ri / dri;
            if (a < alpha) alpha = a;
        }
    }

    /* Step back slightly from boundary */
    alpha *= FN_STEP_BACK;

    if (alpha < 1e-10) alpha = 1e-10;
    if (alpha > 1.0) alpha = 1.0;

    return alpha;
}

/* ============================================================================
 * Check Optimality
 *
 * LP optimality requires X'g = 0 where:
 *   g_i = q      if r_i > 0  (positive residual)
 *   g_i = q - 1  if r_i < 0  (negative residual)
 *   g_i ∈ [q-1, q] if r_i = 0 (on hyperplane)
 * ============================================================================ */

ST_int fn_check_optimality(const ST_double *X, const ST_double *r,
                            ST_int N, ST_int K, ST_double q,
                            ST_double *Xg)
{
    ST_int i, j;

    /* Compute X'g */
    memset(Xg, 0, K * sizeof(ST_double));

    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        ST_double sum = 0.0;
        for (i = 0; i < N; i++) {
            ST_double gi;
            if (r[i] > FN_ACTIVE_TOL) {
                gi = q;
            } else if (r[i] < -FN_ACTIVE_TOL) {
                gi = q - 1.0;
            } else {
                /* On hyperplane: use midpoint (this works for median) */
                gi = q - 0.5;
            }
            sum += Xj[i] * gi;
        }
        Xg[j] = sum;
    }

    /* Check norm of X'g */
    ST_double norm_sq = 0.0;
    for (j = 0; j < K; j++) {
        norm_sq += Xg[j] * Xg[j];
    }

    return (sqrt(norm_sq) < FN_DUAL_TOL * sqrt((ST_double)N)) ? 1 : 0;
}

/* ============================================================================
 * Main Frisch-Newton Solver
 * ============================================================================ */

ST_int cqreg_fn_solve(cqreg_ipm_state *ipm,
                       const ST_double *y,
                       const ST_double *X,
                       ST_double q,
                       ST_double *beta)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j, iter;

    if (ipm == NULL || y == NULL || X == NULL || beta == NULL) {
        return -1;
    }

    if (q <= 0.0 || q >= 1.0) {
        return -2;
    }

    /* Allocate working arrays */
    ST_int *active_idx = (ST_int *)malloc(N * sizeof(ST_int));
    ST_double *g = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *delta_r = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *Xg = (ST_double *)malloc(K * sizeof(ST_double));

    if (!active_idx || !g || !delta_r || !Xg) {
        free(active_idx);
        free(g);
        free(delta_r);
        free(Xg);
        return -3;
    }

    /* Step 1: Initialize beta via OLS */
    /* Compute X'X */
    memset(ipm->XDX, 0, K * K * sizeof(ST_double));
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        for (ST_int k = j; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += Xj[i] * Xk[i];
            }
            ipm->XDX[j * K + k] = sum;
            ipm->XDX[k * K + j] = sum;
        }
    }

    /* Add regularization */
    for (j = 0; j < K; j++) {
        ipm->XDX[j * K + j] += FN_REG;
    }

    /* Cholesky and solve for OLS */
    memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(ipm->L, K) != 0) {
        memset(beta, 0, K * sizeof(ST_double));
    } else {
        for (j = 0; j < K; j++) {
            beta[j] = cqreg_dot(&X[j * N], y, N);
        }
        cqreg_solve_cholesky(ipm->L, beta, K);
    }

    /* Copy to ipm state */
    memcpy(ipm->beta, beta, K * sizeof(ST_double));

    /* Step 2: Compute initial residuals */
    cqreg_matvec_col(ipm->r_primal, X, beta, N, K);
    for (i = 0; i < N; i++) {
        ipm->r_primal[i] = y[i] - ipm->r_primal[i];
    }

    /* Main Frisch-Newton loop
     *
     * Key insight: We need to find beta such that X'g = 0, where
     * g_i = q if r_i > 0, g_i = q-1 if r_i < 0.
     *
     * The gradient of the objective w.r.t. beta is -X'g.
     * We use Newton direction: delta_beta = (X'X)^{-1} X'g
     * This moves beta to reduce ||X'g|| toward zero.
     */
    ST_int total_iter = 0;
    ST_double last_norm_Xg = 1e30;
    ST_int stall_count = 0;

    for (ST_int outer = 0; outer < FN_MAX_OUTER_ITER; outer++) {

        /* Compute subgradient g_i based on residual sign */
        for (i = 0; i < N; i++) {
            if (ipm->r_primal[i] > 0) {
                g[i] = q;
            } else if (ipm->r_primal[i] < 0) {
                g[i] = q - 1.0;
            } else {
                /* Exactly on hyperplane - use value that helps satisfy X'g = 0 */
                g[i] = q - 0.5;
            }
        }

        /* Compute X'g - this should be zero at optimum */
        ST_double norm_Xg = 0.0;
        for (j = 0; j < K; j++) {
            Xg[j] = 0.0;
            for (i = 0; i < N; i++) {
                Xg[j] += X[j * N + i] * g[i];
            }
            norm_Xg += Xg[j] * Xg[j];
        }
        norm_Xg = sqrt(norm_Xg);

        if (norm_Xg < FN_DUAL_TOL * sqrt((ST_double)N)) {
            if (ipm->config.verbose) {
                char buf[128];
                snprintf(buf, sizeof(buf), "FN converged in %d iterations, ||Xg||=%.2e\n",
                         outer + 1, norm_Xg);
                SF_display(buf);
            }
            break;
        }

        /* Check for stalling */
        if (norm_Xg >= 0.99 * last_norm_Xg) {
            stall_count++;
            /* Allow more stall iterations to try different step sizes */
            if (stall_count >= 10) {
                if (ipm->config.verbose) {
                    char buf[128];
                    snprintf(buf, sizeof(buf), "FN converged after %d iters, ||Xg||=%.2e\n",
                             outer + 1, norm_Xg);
                    SF_display(buf);
                }
                break;
            }
        } else {
            stall_count = 0;
        }
        last_norm_Xg = norm_Xg;

        if (ipm->config.verbose && (outer < 10 || outer % 10 == 0)) {
            char buf[128];
            snprintf(buf, sizeof(buf), "FN iter %d: ||Xg||=%.2e, stall=%d\n",
                     outer + 1, norm_Xg, stall_count);
            SF_display(buf);
            fflush(stdout);  /* Force output */
        }

        total_iter++;

        /* Compute Newton direction: delta_beta = (X'X)^{-1} X'g
         * This is the direction that would zero out X'g in one step
         * if the problem were linear (which it isn't due to sign changes).
         *
         * We use the full X'X (already computed during init), solve for delta_beta.
         */

        /* X'g is already in Xg, copy to delta_beta for solving */
        memcpy(ipm->delta_beta, Xg, K * sizeof(ST_double));

        /* Solve (X'X) delta_beta = X'g using existing Cholesky factor
         * Note: ipm->L should have Cholesky of X'X from initialization
         */
        cqreg_solve_cholesky(ipm->L, ipm->delta_beta, K);

        /* Compute delta_r = -X * delta_beta (change in residuals) */
        cqreg_matvec_col(delta_r, X, ipm->delta_beta, N, K);
        for (i = 0; i < N; i++) {
            delta_r[i] = -delta_r[i];
        }

        /* Check for NaN in delta_beta (numerical instability) */
        ST_int has_nan = 0;
        for (j = 0; j < K; j++) {
            if (!isfinite(ipm->delta_beta[j])) {
                has_nan = 1;
                break;
            }
        }
        if (has_nan) {
            /* Numerical instability - stop iterating */
            break;
        }

        /* Use simple damped step - more stable than ratio test */
        ST_double alpha = 0.5;

        /* Update beta and residuals */
        for (j = 0; j < K; j++) {
            beta[j] += alpha * ipm->delta_beta[j];
        }
        for (i = 0; i < N; i++) {
            ipm->r_primal[i] += alpha * delta_r[i];
        }

        memcpy(ipm->beta, beta, K * sizeof(ST_double));
    }

    /* Set final u, v from residuals */
    for (i = 0; i < N; i++) {
        ST_double r = ipm->r_primal[i];
        if (r >= 0) {
            ipm->u[i] = r;
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -r;
        }
    }

    ipm->iterations = total_iter;
    ipm->converged = 1;

    if (ipm->config.verbose) {
        char buf[128];
        snprintf(buf, sizeof(buf), "FN finished: %d iters, cleaning up...\n", total_iter);
        SF_display(buf);
    }

    /* Cleanup */
    free(active_idx);
    free(g);
    free(delta_r);
    free(Xg);

    if (ipm->config.verbose) {
        SF_display("FN cleanup done\n");
    }

    return total_iter;
}

/* ============================================================================
 * Preprocessing Solver (Chernozhukov et al. algorithm)
 * ============================================================================ */

/* Comparison function for sorting by absolute residual */
typedef struct {
    ST_int idx;
    ST_double abs_r;
} fn_resid_pair;

static int fn_compare_resid(const void *a, const void *b)
{
    ST_double diff = ((const fn_resid_pair *)a)->abs_r - ((const fn_resid_pair *)b)->abs_r;
    return (diff < 0) ? -1 : ((diff > 0) ? 1 : 0);
}

ST_int cqreg_fn_preprocess_solve(cqreg_ipm_state *ipm,
                                  const ST_double *y,
                                  const ST_double *X,
                                  ST_double q,
                                  ST_double *beta)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j;

    if (ipm->config.verbose) {
        SF_display("Frisch-Newton with preprocessing (Chernozhukov et al.)\n");
    }

    /* Step 1: Compute initial subsample size
     * m = (N * (log(K) + 1))^{2/3}
     */
    ST_double logK = (K > 1) ? log((double)K) : 0.0;
    ST_int m_init = (ST_int)(pow((double)N * (logK + 1.0), 2.0/3.0) + 0.5);

    if (m_init < 3 * K) m_init = 3 * K;
    if (m_init > N) m_init = N;

    if (ipm->config.verbose) {
        char buf[128];
        snprintf(buf, sizeof(buf), "  Initial subsample: %d (of %d, %.1f%%)\n",
                 m_init, N, 100.0 * m_init / N);
        SF_display(buf);
    }

    /* Step 2: Get OLS solution for initial residuals */
    /* Compute X'X */
    memset(ipm->XDX, 0, K * K * sizeof(ST_double));
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        for (ST_int k = j; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += Xj[i] * Xk[i];
            }
            ipm->XDX[j * K + k] = sum;
            ipm->XDX[k * K + j] = sum;
        }
    }
    for (j = 0; j < K; j++) {
        ipm->XDX[j * K + j] += FN_REG;
    }

    memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(ipm->L, K) != 0) {
        memset(beta, 0, K * sizeof(ST_double));
    } else {
        for (j = 0; j < K; j++) {
            beta[j] = cqreg_dot(&X[j * N], y, N);
        }
        cqreg_solve_cholesky(ipm->L, beta, K);
    }

    /* Step 3: Compute OLS residuals and sort by |r| */
    fn_resid_pair *resid_order = (fn_resid_pair *)malloc(N * sizeof(fn_resid_pair));
    ST_int *subsample = (ST_int *)malloc(N * sizeof(ST_int));
    ST_int *in_subsample = (ST_int *)calloc(N, sizeof(ST_int));
    ST_double *y_sub = NULL;
    ST_double *X_sub = NULL;

    if (!resid_order || !subsample || !in_subsample) {
        free(resid_order);
        free(subsample);
        free(in_subsample);
        return -1;
    }

    for (i = 0; i < N; i++) {
        ST_double yhat = 0.0;
        for (j = 0; j < K; j++) {
            yhat += X[j * N + i] * beta[j];
        }
        resid_order[i].idx = i;
        resid_order[i].abs_r = fabs(y[i] - yhat);
    }

    qsort(resid_order, N, sizeof(fn_resid_pair), fn_compare_resid);

    /* Initialize subsample with m_init smallest |residual| observations */
    ST_int m = m_init;
    for (i = 0; i < m; i++) {
        subsample[i] = resid_order[i].idx;
        in_subsample[resid_order[i].idx] = 1;
    }

    /* Step 4: Iterative refinement */
    ST_int total_iters = 0;
    ST_int max_outer = 20;

    for (ST_int outer = 0; outer < max_outer; outer++) {

        if (ipm->config.verbose) {
            char buf[128];
            snprintf(buf, sizeof(buf), "  Preprocess iter %d: subsample=%d\n", outer + 1, m);
            SF_display(buf);
        }

        /* Allocate subsample data */
        y_sub = (ST_double *)realloc(y_sub, m * sizeof(ST_double));
        X_sub = (ST_double *)realloc(X_sub, m * K * sizeof(ST_double));

        if (!y_sub || !X_sub) {
            free(resid_order);
            free(subsample);
            free(in_subsample);
            free(y_sub);
            free(X_sub);
            return -1;
        }

        /* Extract subsample data */
        for (i = 0; i < m; i++) {
            ST_int obs = subsample[i];
            y_sub[i] = y[obs];
            for (j = 0; j < K; j++) {
                X_sub[j * m + i] = X[j * N + obs];
            }
        }

        /* Create temporary IPM state for subsample */
        cqreg_ipm_state *ipm_sub = cqreg_ipm_create(m, K, &ipm->config);
        if (!ipm_sub) {
            free(resid_order);
            free(subsample);
            free(in_subsample);
            free(y_sub);
            free(X_sub);
            return -1;
        }

        /* Solve on subsample */
        ST_int sub_iters = cqreg_fn_solve(ipm_sub, y_sub, X_sub, q, beta);
        total_iters += (sub_iters > 0) ? sub_iters : 1;

        cqreg_ipm_free(ipm_sub);

        /* Check optimality on full sample */
        /* Compute residuals on full sample */
        cqreg_matvec_col(ipm->r_primal, X, beta, N, K);
        for (i = 0; i < N; i++) {
            ipm->r_primal[i] = y[i] - ipm->r_primal[i];
        }

        /* Check X'g = 0 */
        ST_double *Xg = ipm->work_K;
        if (fn_check_optimality(X, ipm->r_primal, N, K, q, Xg)) {
            if (ipm->config.verbose) {
                char buf[128];
                snprintf(buf, sizeof(buf), "  Optimal on full sample after %d iterations\n", outer + 1);
                SF_display(buf);
            }
            break;
        }

        /* Find violations: observations not in subsample with wrong sign prediction
         * or near-zero residuals that should be included */
        ST_int n_add = 0;
        for (i = 0; i < N && m + n_add < N; i++) {
            if (!in_subsample[i]) {
                ST_double r = ipm->r_primal[i];
                ST_double tol = 1e-6 * (fabs(y[i]) + 1.0);

                /* Add if residual is near zero (should be in active set) */
                if (fabs(r) < tol * 100) {
                    subsample[m + n_add] = i;
                    in_subsample[i] = 1;
                    n_add++;
                }
            }
        }

        /* If no violations found but not optimal, add obs with smallest |r| not in subsample */
        if (n_add == 0) {
            for (i = 0; i < N && m < N; i++) {
                ST_int idx = resid_order[i].idx;
                if (!in_subsample[idx]) {
                    subsample[m] = idx;
                    in_subsample[idx] = 1;
                    m++;
                    n_add++;
                    if (n_add >= 10) break;  /* Add a few at a time */
                }
            }
        } else {
            m += n_add;
        }

        /* If subsample is full, we're done */
        if (m >= N) {
            if (ipm->config.verbose) {
                SF_display("  Expanded to full sample\n");
            }
            break;
        }
    }

    /* Set final state */
    memcpy(ipm->beta, beta, K * sizeof(ST_double));
    for (i = 0; i < N; i++) {
        ST_double r = ipm->r_primal[i];
        if (r >= 0) {
            ipm->u[i] = r;
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -r;
        }
    }

    ipm->iterations = total_iters;
    ipm->converged = 1;

    /* Cleanup */
    free(resid_order);
    free(subsample);
    free(in_subsample);
    free(y_sub);
    free(X_sub);

    return total_iters;
}
