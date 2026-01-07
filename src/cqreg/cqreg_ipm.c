/*
 * cqreg_ipm.c
 *
 * Interior Point Method solver for quantile regression.
 * Implements primal-dual IPM with Mehrotra predictor-corrector.
 * Part of the ctools suite.
 *
 * LP formulation:
 *   min  q * sum(u) + (1-q) * sum(v)
 *   s.t. y - X*beta = u - v
 *        u, v >= 0
 *
 * Dual variables: lambda_u (for u >= 0), lambda_v (for v >= 0)
 * At optimum: lambda_u = q for positive residuals, lambda_v = 1-q for negative
 */

#include "cqreg_ipm.h"
#include "cqreg_linalg.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Small constant to prevent division by zero */
#define IPM_EPS 1e-14

/* Minimum value for primal/dual variables */
#define IPM_MIN_VAL 1e-12

/* Step-back factor from boundary */
#define IPM_STEP_BACK 0.9995

/* ============================================================================
 * Initialization
 * ============================================================================ */

void cqreg_ipm_initialize(cqreg_ipm_state *ipm,
                          const ST_double *y,
                          const ST_double *X,
                          ST_double q)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j;

    /* Compute initial beta via OLS: beta = (X'X)^{-1} X'y */
    /* This gives a good starting point close to the solution */

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

    /* Add regularization for stability */
    for (j = 0; j < K; j++) {
        ipm->XDX[j * K + j] += 1e-10;
    }

    /* Cholesky factorization */
    memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(ipm->L, K) != 0) {
        /* If Cholesky fails, start with zero beta */
        memset(ipm->beta, 0, K * sizeof(ST_double));
    } else {
        /* Compute X'y */
        for (j = 0; j < K; j++) {
            ipm->beta[j] = cqreg_dot(&X[j * N], y, N);
        }
        /* Solve (X'X)^{-1} X'y */
        cqreg_solve_cholesky(ipm->L, ipm->beta, K);
    }

    /* Compute initial residuals: r = y - X * beta */
    cqreg_matvec_col(ipm->r_primal, X, ipm->beta, N, K);
    for (i = 0; i < N; i++) {
        ipm->r_primal[i] = y[i] - ipm->r_primal[i];
    }

    /* Initialize primal variables u, v based on residual sign
     * We need u - v = r and u, v > 0
     * Start with u, v sufficiently away from zero
     */
    ST_double r_scale = 0.0;
    for (i = 0; i < N; i++) {
        r_scale += fabs(ipm->r_primal[i]);
    }
    r_scale = (r_scale / N) + 1.0;  /* Average absolute residual + 1 */

    for (i = 0; i < N; i++) {
        ST_double r = ipm->r_primal[i];
        if (r >= 0) {
            ipm->u[i] = r + r_scale;
            ipm->v[i] = r_scale;
        } else {
            ipm->u[i] = r_scale;
            ipm->v[i] = -r + r_scale;
        }
    }

    /* Initialize dual variables for complementarity
     * We want u * lambda_u ≈ mu and v * lambda_v ≈ mu
     * where mu is the initial barrier parameter
     */
    ST_double mu_target = r_scale;

    for (i = 0; i < N; i++) {
        ipm->lambda_u[i] = mu_target / ipm->u[i];
        ipm->lambda_v[i] = mu_target / ipm->v[i];

        /* Clamp to prevent extreme values */
        if (ipm->lambda_u[i] < IPM_MIN_VAL) ipm->lambda_u[i] = IPM_MIN_VAL;
        if (ipm->lambda_v[i] < IPM_MIN_VAL) ipm->lambda_v[i] = IPM_MIN_VAL;
        if (ipm->lambda_u[i] > 1.0 - IPM_MIN_VAL) ipm->lambda_u[i] = 1.0 - IPM_MIN_VAL;
        if (ipm->lambda_v[i] > 1.0 - IPM_MIN_VAL) ipm->lambda_v[i] = 1.0 - IPM_MIN_VAL;
    }

    /* Compute initial barrier parameter */
    ipm->mu = cqreg_ipm_complementarity(ipm) / (2.0 * N);

    ipm->converged = 0;
    ipm->iterations = 0;
}

/* ============================================================================
 * Residual Computation
 * ============================================================================ */

void cqreg_ipm_compute_residuals(cqreg_ipm_state *ipm,
                                 const ST_double *y,
                                 const ST_double *X,
                                 ST_double q)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i;

    /* Primal residual: r_p = y - X*beta - u + v
     * At feasibility: r_p = 0
     */
    cqreg_matvec_col(ipm->r_primal, X, ipm->beta, N, K);
    for (i = 0; i < N; i++) {
        ipm->r_primal[i] = y[i] - ipm->r_primal[i] - ipm->u[i] + ipm->v[i];
    }

    /* Dual residual: r_d = q*X'1 - X'*lambda_u
     *
     * The KKT conditions give: lambda_u = q - λ, lambda_v = (1-q) + λ
     * where λ is the dual variable for the equality constraint.
     * The dual constraint is X'λ = 0, i.e., X'(q - lambda_u) = 0.
     * So: r_d = X'λ = q*X'1 - X'*lambda_u
     */
    for (ST_int j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        ST_double col_sum = 0.0;
        ST_double weighted_sum = 0.0;
        for (i = 0; i < N; i++) {
            col_sum += Xj[i];
            weighted_sum += Xj[i] * ipm->lambda_u[i];
        }
        ipm->r_dual[j] = q * col_sum - weighted_sum;
    }

    /* Compute infeasibility measures */
    ipm->primal_infeas = cqreg_vmax_abs(ipm->r_primal, N);
    ipm->dual_infeas = cqreg_vmax_abs(ipm->r_dual, K);

    /* Compute primal objective: q * sum(u) + (1-q) * sum(v) */
    ST_double sum_u = 0.0, sum_v = 0.0;
    for (i = 0; i < N; i++) {
        sum_u += ipm->u[i];
        sum_v += ipm->v[i];
    }
    ipm->primal_obj = q * sum_u + (1.0 - q) * sum_v;

    /* Dual objective: y' * (lambda_u - lambda_v) */
    ST_double dual_obj = 0.0;
    for (i = 0; i < N; i++) {
        dual_obj += y[i] * (ipm->lambda_u[i] - ipm->lambda_v[i]);
    }
    ipm->dual_obj = dual_obj;

    /* Duality gap */
    ipm->gap = fabs(ipm->primal_obj - ipm->dual_obj);
}

/* ============================================================================
 * Convergence Check
 * ============================================================================ */

ST_int cqreg_ipm_check_convergence(cqreg_ipm_state *ipm)
{
    /* Check primal feasibility */
    if (ipm->primal_infeas > ipm->config.tol_primal) return 0;

    /* Check dual feasibility */
    if (ipm->dual_infeas > ipm->config.tol_dual) return 0;

    /* Check complementarity (via barrier parameter) */
    if (ipm->mu > ipm->config.tol_gap) return 0;

    return 1;
}

/* ============================================================================
 * Form Normal Equations and Solve
 * ============================================================================ */

/*
 * Form and solve the reduced system for delta_beta.
 *
 * The KKT system for the perturbed problem is:
 *
 *  [ 0    -X     I    -I   0    0   ] [delta_beta    ]   [ -r_d          ]
 *  [-X'    0     0     0   I   -I   ] [delta_lambda  ]   [ -r_p          ]
 *  [ 0     0    U^-1  0   0    0   ] [delta_u       ] = [ mu*e - u*lu   ]
 *  [ 0     0     0    V^-1 0    0   ] [delta_v       ]   [ mu*e - v*lv   ]
 *  [ Lu    0     0     0   U    0   ] [delta_lu      ]   [               ]
 *  [ 0    Lv     0     0   0    V   ] [delta_lv      ]   [               ]
 *
 * After elimination, we get the reduced system:
 *   X' * D^{-1} * X * delta_beta = rhs
 * where D = diag(lambda_u/u + lambda_v/v)
 */
static ST_int solve_newton_system(cqreg_ipm_state *ipm,
                                   const ST_double *X,
                                   const ST_double *y,
                                   ST_double q,
                                   ST_double sigma,
                                   const ST_double *delta_u_aff,
                                   const ST_double *delta_v_aff,
                                   const ST_double *delta_lu_aff,
                                   const ST_double *delta_lv_aff)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j, k;

    ST_double mu_sigma = sigma * ipm->mu;

    /* Step 1: Compute diagonal scaling D[i] = u[i]/lambda_u[i] + v[i]/lambda_v[i]
     * This is for the reduced system: (X' D^{-1} X) Δβ = rhs
     */
    for (i = 0; i < N; i++) {
        ipm->D[i] = ipm->u[i] / ipm->lambda_u[i] + ipm->v[i] / ipm->lambda_v[i];
    }

    /* Step 2: Compute X' * D^{-1} * X */
    memset(ipm->XDX, 0, K * K * sizeof(ST_double));
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];

        /* Diagonal */
        ST_double diag = 0.0;
        for (i = 0; i < N; i++) {
            diag += Xj[i] * Xj[i] / ipm->D[i];
        }
        ipm->XDX[j * K + j] = diag;

        /* Off-diagonal (upper triangle) */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += Xj[i] * Xk[i] / ipm->D[i];
            }
            ipm->XDX[j * K + k] = sum;
            ipm->XDX[k * K + j] = sum;
        }
    }

    /* Add small regularization */
    for (j = 0; j < K; j++) {
        ipm->XDX[j * K + j] += 1e-12;
    }

    /* Step 3: Cholesky factorization */
    memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(ipm->L, K) != 0) {
        return -1;  /* Cholesky failed */
    }

    /* Step 4: Compute RHS for reduced system
     *
     * From the KKT derivation:
     *   (X' D^{-1} X) Δβ = -r_d + X' D^{-1} ξ
     * where ξ[i] = r_p[i] + comp_u[i]/lu[i] - comp_v[i]/lv[i]
     *
     * r_p = y - Xβ - u + v (primal residual, should be 0)
     * r_d = X'λ = q*X'1 - X'*lambda_u (dual residual, should be 0)
     */
    for (i = 0; i < N; i++) {
        ST_double u = ipm->u[i];
        ST_double v = ipm->v[i];
        ST_double lu = ipm->lambda_u[i];
        ST_double lv = ipm->lambda_v[i];

        /* Centering terms */
        ST_double comp_u = mu_sigma - u * lu;
        ST_double comp_v = mu_sigma - v * lv;

        /* Corrector terms (Mehrotra) */
        if (delta_u_aff != NULL) {
            comp_u -= delta_u_aff[i] * delta_lu_aff[i];
            comp_v -= delta_v_aff[i] * delta_lv_aff[i];
        }

        /* Modified RHS: ξ = r_p + comp_u/lu - comp_v/lv */
        ST_double rhs_i = ipm->r_primal[i] + comp_u / lu - comp_v / lv;

        /* Store D^{-1} * ξ */
        ipm->work_N[i] = rhs_i / ipm->D[i];
    }

    /* Step 5: Compute -r_d + X' * D^{-1} * ξ */
    for (j = 0; j < K; j++) {
        ST_double sum = -ipm->r_dual[j];
        const ST_double *Xj = &X[j * N];
        for (i = 0; i < N; i++) {
            sum += Xj[i] * ipm->work_N[i];
        }
        ipm->delta_beta[j] = sum;
    }

    /* Step 6: Solve for delta_beta */
    cqreg_solve_cholesky(ipm->L, ipm->delta_beta, K);

    /* Step 7: Back-substitute to get delta_u, delta_v, delta_lu, delta_lv */
    cqreg_matvec_col(ipm->work_N, X, ipm->delta_beta, N, K);

    for (i = 0; i < N; i++) {
        ST_double u = ipm->u[i];
        ST_double v = ipm->v[i];
        ST_double lu = ipm->lambda_u[i];
        ST_double lv = ipm->lambda_v[i];

        ST_double comp_u = mu_sigma - u * lu;
        ST_double comp_v = mu_sigma - v * lv;

        if (delta_u_aff != NULL) {
            comp_u -= delta_u_aff[i] * delta_lu_aff[i];
            comp_v -= delta_v_aff[i] * delta_lv_aff[i];
        }

        /* From the derivation:
         *   D*Δλ = X*Δβ - r_p - comp_u/lu + comp_v/lv
         *   Δλ = (X*Δβ - r_p - comp_u/lu + comp_v/lv) / D
         */
        ST_double D_i = u / lu + v / lv;
        ST_double delta_lambda = (ipm->work_N[i] - ipm->r_primal[i] - comp_u / lu + comp_v / lv) / D_i;

        /* Back-substitute for Δu, Δv, Δlambda_u, Δlambda_v */
        ipm->delta_u[i] = (comp_u + u * delta_lambda) / lu;
        ipm->delta_v[i] = (comp_v - v * delta_lambda) / lv;
        ipm->delta_lambda_u[i] = -delta_lambda;
        ipm->delta_lambda_v[i] = delta_lambda;
    }

    return 0;
}

/* ============================================================================
 * Affine Direction (Predictor Step)
 * ============================================================================ */

void cqreg_ipm_affine_direction(cqreg_ipm_state *ipm,
                                const ST_double *X,
                                const ST_double *y,
                                ST_double q)
{
    /* Solve with sigma = 0 (no centering) and no corrector */
    if (solve_newton_system(ipm, X, y, q, 0.0, NULL, NULL, NULL, NULL) != 0) {
        /* If solve fails, use zeros */
        memset(ipm->delta_beta, 0, ipm->K * sizeof(ST_double));
        memset(ipm->delta_u, 0, ipm->N * sizeof(ST_double));
        memset(ipm->delta_v, 0, ipm->N * sizeof(ST_double));
        memset(ipm->delta_lambda_u, 0, ipm->N * sizeof(ST_double));
        memset(ipm->delta_lambda_v, 0, ipm->N * sizeof(ST_double));
    }

    /* Copy to affine directions */
    cqreg_vcopy(ipm->delta_u_aff, ipm->delta_u, ipm->N);
    cqreg_vcopy(ipm->delta_v_aff, ipm->delta_v, ipm->N);
    cqreg_vcopy(ipm->delta_lambda_u_aff, ipm->delta_lambda_u, ipm->N);
    cqreg_vcopy(ipm->delta_lambda_v_aff, ipm->delta_lambda_v, ipm->N);
}

/* ============================================================================
 * Combined Direction (Predictor-Corrector)
 * ============================================================================ */

void cqreg_ipm_combined_direction(cqreg_ipm_state *ipm,
                                  const ST_double *X,
                                  const ST_double *y,
                                  ST_double q,
                                  ST_double sigma)
{
    /* Solve with centering and corrector */
    solve_newton_system(ipm, X, y, q, sigma,
                        ipm->delta_u_aff, ipm->delta_v_aff,
                        ipm->delta_lambda_u_aff, ipm->delta_lambda_v_aff);
}

/* ============================================================================
 * Step Length Computation
 * ============================================================================ */

ST_double cqreg_ipm_step_length(cqreg_ipm_state *ipm, ST_int affine)
{
    ST_int N = ipm->N;
    ST_double alpha = 1.0;

    ST_double *du = affine ? ipm->delta_u_aff : ipm->delta_u;
    ST_double *dv = affine ? ipm->delta_v_aff : ipm->delta_v;
    ST_double *dlu = affine ? ipm->delta_lambda_u_aff : ipm->delta_lambda_u;
    ST_double *dlv = affine ? ipm->delta_lambda_v_aff : ipm->delta_lambda_v;

    /* Find maximum step maintaining positivity */
    for (ST_int i = 0; i < N; i++) {
        /* u + alpha * delta_u > 0 */
        if (du[i] < -IPM_EPS) {
            ST_double a = -ipm->u[i] / du[i];
            if (a < alpha) alpha = a;
        }

        /* v + alpha * delta_v > 0 */
        if (dv[i] < -IPM_EPS) {
            ST_double a = -ipm->v[i] / dv[i];
            if (a < alpha) alpha = a;
        }

        /* lambda_u + alpha * delta_lambda_u > 0 */
        if (dlu[i] < -IPM_EPS) {
            ST_double a = -ipm->lambda_u[i] / dlu[i];
            if (a < alpha) alpha = a;
        }

        /* lambda_v + alpha * delta_lambda_v > 0 */
        if (dlv[i] < -IPM_EPS) {
            ST_double a = -ipm->lambda_v[i] / dlv[i];
            if (a < alpha) alpha = a;
        }
    }

    /* Step back from boundary */
    alpha *= IPM_STEP_BACK;

    /* Ensure minimum step */
    if (alpha < IPM_MIN_VAL) {
        alpha = IPM_MIN_VAL;
    }

    return alpha;
}

/* ============================================================================
 * Variable Update
 * ============================================================================ */

void cqreg_ipm_update_variables(cqreg_ipm_state *ipm, ST_double alpha)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;

    /* Update beta */
    for (ST_int k = 0; k < K; k++) {
        ipm->beta[k] += alpha * ipm->delta_beta[k];
    }

    /* Update primal and dual variables */
    for (ST_int i = 0; i < N; i++) {
        ipm->u[i] += alpha * ipm->delta_u[i];
        ipm->v[i] += alpha * ipm->delta_v[i];
        ipm->lambda_u[i] += alpha * ipm->delta_lambda_u[i];
        ipm->lambda_v[i] += alpha * ipm->delta_lambda_v[i];

        /* Ensure positivity (numerical safety) */
        if (ipm->u[i] < IPM_MIN_VAL) ipm->u[i] = IPM_MIN_VAL;
        if (ipm->v[i] < IPM_MIN_VAL) ipm->v[i] = IPM_MIN_VAL;
        if (ipm->lambda_u[i] < IPM_MIN_VAL) ipm->lambda_u[i] = IPM_MIN_VAL;
        if (ipm->lambda_v[i] < IPM_MIN_VAL) ipm->lambda_v[i] = IPM_MIN_VAL;
    }
}

/* ============================================================================
 * Complementarity
 * ============================================================================ */

ST_double cqreg_ipm_complementarity(const cqreg_ipm_state *ipm)
{
    ST_int N = ipm->N;
    ST_double sum = 0.0;

    for (ST_int i = 0; i < N; i++) {
        sum += ipm->u[i] * ipm->lambda_u[i] + ipm->v[i] * ipm->lambda_v[i];
    }

    return sum;
}

/* ============================================================================
 * Scaling Matrix (for external use)
 * ============================================================================ */

void cqreg_ipm_compute_scaling(cqreg_ipm_state *ipm)
{
    ST_int N = ipm->N;

    for (ST_int i = 0; i < N; i++) {
        ipm->D[i] = ipm->u[i] / ipm->lambda_u[i] + ipm->v[i] / ipm->lambda_v[i];
    }
}

/* ============================================================================
 * Form Normal Equations (for external use)
 * ============================================================================ */

ST_int cqreg_ipm_form_normal_equations(cqreg_ipm_state *ipm, const ST_double *X)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;

    /* X' * D^{-1} * X */
    cqreg_xtdx(ipm->XDX, X, ipm->D, N, K);

    /* Cholesky */
    memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
    for (ST_int j = 0; j < K; j++) {
        ipm->L[j * K + j] += 1e-12;
    }

    return cqreg_cholesky(ipm->L, K);
}

/* ============================================================================
 * Main Solver
 * ============================================================================ */

ST_int cqreg_ipm_solve(cqreg_ipm_state *ipm,
                       const ST_double *y,
                       const ST_double *X,
                       ST_double q,
                       ST_double *beta)
{
    ST_int iter;
    ST_int N = ipm->N;

    if (ipm == NULL || y == NULL || X == NULL || beta == NULL) {
        return -1;
    }

    if (q <= 0.0 || q >= 1.0) {
        return -2;
    }

    /* Initialize variables */
    cqreg_ipm_initialize(ipm, y, X, q);

    for (iter = 0; iter < ipm->config.maxiter; iter++) {
        /* Compute residuals */
        cqreg_ipm_compute_residuals(ipm, y, X, q);

        /* Check convergence */
        if (cqreg_ipm_check_convergence(ipm)) {
            ipm->converged = 1;
            break;
        }

        ST_double sigma;

        if (ipm->config.use_mehrotra) {
            /* Predictor step (affine scaling, sigma = 0) */
            cqreg_ipm_affine_direction(ipm, X, y, q);

            /* Compute step length for affine direction */
            ST_double alpha_aff = cqreg_ipm_step_length(ipm, 1);

            /* Compute mu after affine step */
            ST_double mu_aff = 0.0;
            for (ST_int i = 0; i < N; i++) {
                ST_double u_new = ipm->u[i] + alpha_aff * ipm->delta_u_aff[i];
                ST_double v_new = ipm->v[i] + alpha_aff * ipm->delta_v_aff[i];
                ST_double lu_new = ipm->lambda_u[i] + alpha_aff * ipm->delta_lambda_u_aff[i];
                ST_double lv_new = ipm->lambda_v[i] + alpha_aff * ipm->delta_lambda_v_aff[i];
                mu_aff += u_new * lu_new + v_new * lv_new;
            }
            mu_aff /= (2.0 * N);

            /* Adaptive centering parameter */
            ST_double ratio = mu_aff / ipm->mu;
            sigma = ratio * ratio * ratio;  /* sigma = (mu_aff / mu)^3 */
            if (sigma > 1.0) sigma = 1.0;
            if (sigma < 0.0001) sigma = 0.0001;

            /* Corrector step (combined direction) */
            cqreg_ipm_combined_direction(ipm, X, y, q, sigma);
        } else {
            /* Simple centering (no Mehrotra) */
            sigma = ipm->config.sigma;
            cqreg_ipm_affine_direction(ipm, X, y, q);

            /* Re-solve with centering */
            solve_newton_system(ipm, X, y, q, sigma, NULL, NULL, NULL, NULL);
        }

        /* Compute step length */
        ST_double alpha = cqreg_ipm_step_length(ipm, 0);

        /* Update variables */
        cqreg_ipm_update_variables(ipm, alpha);

        /* Update barrier parameter */
        ipm->mu = cqreg_ipm_complementarity(ipm) / (2.0 * N);
    }

    ipm->iterations = iter + 1;

    /* Copy solution */
    cqreg_vcopy(beta, ipm->beta, ipm->K);

    return ipm->converged ? iter + 1 : -(iter + 1);
}

/* ============================================================================
 * Result Extraction
 * ============================================================================ */

ST_double cqreg_ipm_get_objective(const cqreg_ipm_state *ipm, ST_double q)
{
    ST_int N = ipm->N;
    ST_double obj = 0.0;

    for (ST_int i = 0; i < N; i++) {
        obj += q * ipm->u[i] + (1.0 - q) * ipm->v[i];
    }

    return obj;
}

void cqreg_ipm_get_residuals(const cqreg_ipm_state *ipm, ST_double *residuals)
{
    ST_int N = ipm->N;

    /* Residuals = u - v */
    for (ST_int i = 0; i < N; i++) {
        residuals[i] = ipm->u[i] - ipm->v[i];
    }
}
