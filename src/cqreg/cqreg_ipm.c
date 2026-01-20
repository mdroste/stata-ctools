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
    (void)q;  /* Unused - reserved for future use */
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j;

    /* Compute initial beta via OLS: beta = (X'X)^{-1} X'y */
    /* This gives a good starting point close to the solution */

    /* Compute X'X with 8x loop unrolling */
    ST_int N8 = N - (N & 7);
    memset(ipm->XDX, 0, K * K * sizeof(ST_double));
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        /* Diagonal with 8x unrolling */
        ST_double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
        ST_double d4 = 0.0, d5 = 0.0, d6 = 0.0, d7 = 0.0;
        for (i = 0; i < N8; i += 8) {
            d0 += Xj[i]     * Xj[i];
            d1 += Xj[i + 1] * Xj[i + 1];
            d2 += Xj[i + 2] * Xj[i + 2];
            d3 += Xj[i + 3] * Xj[i + 3];
            d4 += Xj[i + 4] * Xj[i + 4];
            d5 += Xj[i + 5] * Xj[i + 5];
            d6 += Xj[i + 6] * Xj[i + 6];
            d7 += Xj[i + 7] * Xj[i + 7];
        }
        for (; i < N; i++) {
            d0 += Xj[i] * Xj[i];
        }
        ipm->XDX[j * K + j] = ((d0 + d4) + (d1 + d5)) + ((d2 + d6) + (d3 + d7));

        /* Off-diagonal with 8x unrolling */
        for (ST_int k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
            for (i = 0; i < N8; i += 8) {
                s0 += Xj[i]     * Xk[i];
                s1 += Xj[i + 1] * Xk[i + 1];
                s2 += Xj[i + 2] * Xk[i + 2];
                s3 += Xj[i + 3] * Xk[i + 3];
                s4 += Xj[i + 4] * Xk[i + 4];
                s5 += Xj[i + 5] * Xk[i + 5];
                s6 += Xj[i + 6] * Xk[i + 6];
                s7 += Xj[i + 7] * Xk[i + 7];
            }
            for (; i < N; i++) {
                s0 += Xj[i] * Xk[i];
            }
            ST_double total = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
            ipm->XDX[j * K + k] = total;
            ipm->XDX[k * K + j] = total;
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
    (void)y;  /* Unused - residuals computed from current state */
    (void)q;  /* Unused - quantile stored in ipm state */
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

        /* Modified RHS: ξ = r_p - comp_u/lu + comp_v/lv */
        ST_double rhs_i = ipm->r_primal[i] - comp_u / lu + comp_v / lv;

        /* Store D^{-1} * ξ */
        ipm->work_N[i] = rhs_i / ipm->D[i];
    }

    /* Step 5: Compute X' * D^{-1} * ξ + r_d */
    for (j = 0; j < K; j++) {
        ST_double sum = ipm->r_dual[j];
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
         *   Δλ = (ξ - X*Δβ) / D
         *   where ξ = r_p - comp_u/lu + comp_v/lv
         */
        ST_double D_i = u / lu + v / lv;
        ST_double delta_lambda = (ipm->r_primal[i] - comp_u / lu + comp_v / lv - ipm->work_N[i]) / D_i;

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
    ST_int K = ipm->K;
    ST_int i, j;

    if (ipm == NULL || y == NULL || X == NULL || beta == NULL) {
        return -1;
    }

    if (q <= 0.0 || q >= 1.0) {
        return -2;
    }

    /* Initialize β with OLS solution */
    cqreg_ipm_initialize(ipm, y, X, q);

    /*
     * Interior Point Method for Quantile Regression
     *
     * We solve the barrier problem:
     *   min_β  Σ ρ_μ(r_i)   where r = y - Xβ
     *
     * where ρ_μ(r) is a smoothed check function:
     *   ρ_μ(r) = q*r + μ*log(1 + exp(-r/μ))  [Huber-like smoothing]
     *
     * As μ → 0, ρ_μ → ρ_q (the quantile check function).
     *
     * The gradient is:
     *   ∂ρ_μ/∂r = q - 1/(1 + exp(r/μ)) = q - sigmoid(-r/μ)
     *
     * The Hessian is:
     *   ∂²ρ_μ/∂r² = (1/μ) * sigmoid(-r/μ) * sigmoid(r/μ)
     *             = (1/μ) * exp(-r/μ) / (1 + exp(-r/μ))²
     *
     * Newton's method:
     *   (X' W X) Δβ = X' g
     * where W = diag(∂²ρ_μ/∂r²) and g = ∂ρ_μ/∂r
     */

    ST_double mu = 1000.0;  /* Large initial smoothing parameter */
    ST_double mu_min = 0.00001;  /* Stop before numerical instability */
    ST_double mu_factor = 0.5;  /* Reduction factor */

    /* Allocate workspace */
    ST_double *w = ipm->D;  /* Weights */
    ST_double *g = ipm->work_N;  /* Gradient */

    /* Pre-compute constant */
    ST_double inv_mu = 1.0 / mu;

    for (iter = 0; iter < ipm->config.maxiter; iter++) {
        /* Compute residuals r = y - Xβ */
        cqreg_matvec_col(ipm->r_primal, X, ipm->beta, N, K);

        /* FUSED LOOP: Compute residuals, objective, gradient, and weights in single pass
         * This reduces memory bandwidth by only traversing the data once */
        ST_double obj = 0.0;
        inv_mu = 1.0 / mu;

        for (i = 0; i < N; i++) {
            /* Compute residual */
            ST_double r = y[i] - ipm->r_primal[i];
            ipm->r_primal[i] = r;

            /* Compute z = -r/mu using pre-computed inverse */
            ST_double z = -r * inv_mu;

            /* Numerically stable sigmoid and objective */
            ST_double sig, log1pexp;
            if (z > 20.0) {
                sig = 1.0;
                log1pexp = z;
            } else if (z < -20.0) {
                sig = 0.0;
                log1pexp = 0.0;
            } else {
                ST_double exp_z = exp(z);
                sig = exp_z / (1.0 + exp_z);  /* sigmoid(z) = exp(z)/(1+exp(z)) */
                log1pexp = log(1.0 + exp_z);
            }

            /* Accumulate objective */
            obj += q * r + mu * log1pexp;

            /* Gradient and Hessian weight */
            g[i] = q - sig;
            ST_double w_val = sig * (1.0 - sig) * inv_mu;

            /* Clamp weights */
            w[i] = (w_val < 1e-12) ? 1e-12 : ((w_val > 1e12) ? 1e12 : w_val);
        }


        /* Form weighted normal equations: (X'WX) Δβ = X'g
         * Parallelized across columns with 8x loop unrolling for CPU pipelining */
        memset(ipm->XDX, 0, K * K * sizeof(ST_double));
        ST_int N8 = N - (N & 7);  /* N rounded down to multiple of 8 */

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) if(K > 2 && N > 1000)
        #endif
        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];
            ST_int ii;

            /* Diagonal element with 8x unrolling */
            ST_double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
            ST_double d4 = 0.0, d5 = 0.0, d6 = 0.0, d7 = 0.0;
            for (ii = 0; ii < N8; ii += 8) {
                d0 += w[ii]     * Xj[ii]     * Xj[ii];
                d1 += w[ii + 1] * Xj[ii + 1] * Xj[ii + 1];
                d2 += w[ii + 2] * Xj[ii + 2] * Xj[ii + 2];
                d3 += w[ii + 3] * Xj[ii + 3] * Xj[ii + 3];
                d4 += w[ii + 4] * Xj[ii + 4] * Xj[ii + 4];
                d5 += w[ii + 5] * Xj[ii + 5] * Xj[ii + 5];
                d6 += w[ii + 6] * Xj[ii + 6] * Xj[ii + 6];
                d7 += w[ii + 7] * Xj[ii + 7] * Xj[ii + 7];
            }
            for (ii = N8; ii < N; ii++) {
                d0 += w[ii] * Xj[ii] * Xj[ii];
            }
            ipm->XDX[j * K + j] = ((d0 + d4) + (d1 + d5)) + ((d2 + d6) + (d3 + d7)) + 1e-10;

            /* Off-diagonal elements with 8x unrolling */
            for (ST_int kk = j + 1; kk < K; kk++) {
                const ST_double *Xk = &X[kk * N];
                ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
                ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
                for (ii = 0; ii < N8; ii += 8) {
                    s0 += w[ii]     * Xj[ii]     * Xk[ii];
                    s1 += w[ii + 1] * Xj[ii + 1] * Xk[ii + 1];
                    s2 += w[ii + 2] * Xj[ii + 2] * Xk[ii + 2];
                    s3 += w[ii + 3] * Xj[ii + 3] * Xk[ii + 3];
                    s4 += w[ii + 4] * Xj[ii + 4] * Xk[ii + 4];
                    s5 += w[ii + 5] * Xj[ii + 5] * Xk[ii + 5];
                    s6 += w[ii + 6] * Xj[ii + 6] * Xk[ii + 6];
                    s7 += w[ii + 7] * Xj[ii + 7] * Xk[ii + 7];
                }
                for (ii = N8; ii < N; ii++) {
                    s0 += w[ii] * Xj[ii] * Xk[ii];
                }
                ST_double total = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
                ipm->XDX[j * K + kk] = total;
                ipm->XDX[kk * K + j] = total;
            }
        }

        /* Compute X'g with 8x unrolling - parallelized */
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(K > 2 && N > 1000)
        #endif
        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];
            ST_int ii;
            ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
            for (ii = 0; ii < N8; ii += 8) {
                s0 += Xj[ii]     * g[ii];
                s1 += Xj[ii + 1] * g[ii + 1];
                s2 += Xj[ii + 2] * g[ii + 2];
                s3 += Xj[ii + 3] * g[ii + 3];
                s4 += Xj[ii + 4] * g[ii + 4];
                s5 += Xj[ii + 5] * g[ii + 5];
                s6 += Xj[ii + 6] * g[ii + 6];
                s7 += Xj[ii + 7] * g[ii + 7];
            }
            for (ii = N8; ii < N; ii++) {
                s0 += Xj[ii] * g[ii];
            }
            ipm->delta_beta[j] = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
        }

        /* Cholesky factorization */
        memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
        if (cqreg_cholesky(ipm->L, K) != 0) {
            /* Cholesky failed - likely mu is too small */
            break;
        }

        /* Solve for Δβ */
        cqreg_solve_cholesky(ipm->L, ipm->delta_beta, K);

        /* Pre-compute X * delta_beta once (cache for line search)
         * Use delta_v as temporary storage since we don't use it in smoothed IPM */
        ST_double *X_delta_beta = ipm->delta_v;
        cqreg_matvec_col(X_delta_beta, X, ipm->delta_beta, N, K);

        /* Pre-compute delta_beta norm squared for Armijo condition */
        ST_double delta_norm_sq = 0.0;
        for (j = 0; j < K; j++) {
            delta_norm_sq += ipm->delta_beta[j] * ipm->delta_beta[j];
        }

        /* Line search with backtracking
         * New residual = old_residual - alpha * X * delta_beta
         * OPTIMIZED: 8x loop unrolling for vectorization */
        ST_double alpha = 1.0;
        ST_double obj_new;
        for (ST_int ls = 0; ls < 20; ls++) {
            /* Compute new objective using cached X*delta_beta with 8x unrolling */
            ST_double obj0 = 0.0, obj1 = 0.0, obj2 = 0.0, obj3 = 0.0;
            ST_double obj4 = 0.0, obj5 = 0.0, obj6 = 0.0, obj7 = 0.0;
            ST_int ii;

            for (ii = 0; ii < N8; ii += 8) {
                /* Unroll 8 iterations for ILP */
                ST_double r0 = ipm->r_primal[ii]     - alpha * X_delta_beta[ii];
                ST_double r1 = ipm->r_primal[ii + 1] - alpha * X_delta_beta[ii + 1];
                ST_double r2 = ipm->r_primal[ii + 2] - alpha * X_delta_beta[ii + 2];
                ST_double r3 = ipm->r_primal[ii + 3] - alpha * X_delta_beta[ii + 3];
                ST_double r4 = ipm->r_primal[ii + 4] - alpha * X_delta_beta[ii + 4];
                ST_double r5 = ipm->r_primal[ii + 5] - alpha * X_delta_beta[ii + 5];
                ST_double r6 = ipm->r_primal[ii + 6] - alpha * X_delta_beta[ii + 6];
                ST_double r7 = ipm->r_primal[ii + 7] - alpha * X_delta_beta[ii + 7];

                ST_double z0 = -r0 * inv_mu, z1 = -r1 * inv_mu;
                ST_double z2 = -r2 * inv_mu, z3 = -r3 * inv_mu;
                ST_double z4 = -r4 * inv_mu, z5 = -r5 * inv_mu;
                ST_double z6 = -r6 * inv_mu, z7 = -r7 * inv_mu;

                /* Branchless approximation for log(1+exp(z)) in typical range */
                #define LS_CONTRIB(r_val, z_val) \
                    ((z_val > 20.0) ? (q * r_val + mu * z_val) : \
                     ((z_val < -20.0) ? (q * r_val) : \
                      (q * r_val + mu * log(1.0 + exp(z_val)))))

                obj0 += LS_CONTRIB(r0, z0);
                obj1 += LS_CONTRIB(r1, z1);
                obj2 += LS_CONTRIB(r2, z2);
                obj3 += LS_CONTRIB(r3, z3);
                obj4 += LS_CONTRIB(r4, z4);
                obj5 += LS_CONTRIB(r5, z5);
                obj6 += LS_CONTRIB(r6, z6);
                obj7 += LS_CONTRIB(r7, z7);
                #undef LS_CONTRIB
            }
            /* Remainder loop */
            for (; ii < N; ii++) {
                ST_double r = ipm->r_primal[ii] - alpha * X_delta_beta[ii];
                ST_double z = -r * inv_mu;
                if (z > 20.0) {
                    obj0 += q * r + mu * z;
                } else if (z < -20.0) {
                    obj0 += q * r;
                } else {
                    obj0 += q * r + mu * log(1.0 + exp(z));
                }
            }
            obj_new = ((obj0 + obj4) + (obj1 + obj5)) + ((obj2 + obj6) + (obj3 + obj7));

            if (obj_new < obj + 1e-4 * alpha * delta_norm_sq) {
                break;  /* Sufficient decrease */
            }
            alpha *= 0.5;
        }

        /* Update β */
        for (j = 0; j < K; j++) {
            ipm->beta[j] += alpha * ipm->delta_beta[j];
        }

        /* Check convergence based on step size and mu
         * Use pre-computed delta_norm_sq */
        ST_double step_norm = sqrt(delta_norm_sq) * alpha;

        /* Converged when step is small and we've reached minimum mu */
        if (step_norm < 1e-6 && mu <= mu_min) {
            ipm->converged = 1;
            break;
        }

        /* Reduce smoothing parameter when converged for current mu */
        if (step_norm < 1e-5) {
            mu *= mu_factor;
            if (mu < mu_min) {
                mu = mu_min;
            }
        }

        ipm->mu = mu;
    }

    ipm->iterations = iter + 1;

    /* Compute final residuals and u, v for reporting */
    cqreg_matvec_col(ipm->r_primal, X, ipm->beta, N, K);
    for (i = 0; i < N; i++) {
        ST_double r = y[i] - ipm->r_primal[i];
        ipm->r_primal[i] = r;
        if (r >= 0) {
            ipm->u[i] = r;
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -r;
        }
        /* Set dual variables at LP optimum */
        if (r > 0) {
            ipm->lambda_u[i] = 0.0;
            ipm->lambda_v[i] = 1.0 - q;
        } else if (r < 0) {
            ipm->lambda_u[i] = 1.0;
            ipm->lambda_v[i] = 0.0;
        } else {
            ipm->lambda_u[i] = q;
            ipm->lambda_v[i] = 1.0 - q;
        }
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

/* ============================================================================
 * Preprocessing Algorithm (Chernozhukov et al. 2020)
 *
 * The key insight: QR solution interpolates exactly K data points.
 * Observations far from the hyperplane don't affect the solution.
 *
 * Algorithm:
 * 1. OLS to get initial beta
 * 2. Compute residuals, predict signs
 * 3. Select m = O((n*log(k))^{2/3}) observations with smallest |r|
 * 4. Solve QR on subsample
 * 5. Check optimality: X'g = 0 where g_i = q - I(r_i < 0)
 * 6. Add violated observations (wrong sign prediction), repeat
 * ============================================================================ */

/* Comparison function for qsort - sort by absolute residual */
typedef struct {
    ST_int idx;
    ST_double abs_resid;
} resid_pair;

static int compare_resid_pair(const void *a, const void *b)
{
    ST_double diff = ((const resid_pair *)a)->abs_resid - ((const resid_pair *)b)->abs_resid;
    return (diff < 0) ? -1 : ((diff > 0) ? 1 : 0);
}

/* Comparison function for qsort - sort doubles ascending */
static int compare_double(const void *a, const void *b)
{
    ST_double diff = *(const ST_double *)a - *(const ST_double *)b;
    return (diff < 0) ? -1 : ((diff > 0) ? 1 : 0);
}

/*
 * Solve QR on a subsample and return beta.
 * Uses the smoothed IPM on the subsample.
 */
static ST_int solve_subsample_qr(
    const ST_double *y_full, const ST_double *X_full,
    ST_int N_full, ST_int K,
    const ST_int *active_set, ST_int m,
    ST_double q, ST_double *beta,
    const cqreg_ipm_config *config)
{
    ST_int i, j;

    /* Allocate subsample data - cast to size_t to prevent 32-bit overflow */
    ST_double *y_sub = (ST_double *)malloc((size_t)m * sizeof(ST_double));
    ST_double *X_sub = (ST_double *)malloc((size_t)m * K * sizeof(ST_double));
    if (!y_sub || !X_sub) {
        free(y_sub);
        free(X_sub);
        return -1;
    }

    /* Extract subsample */
    for (i = 0; i < m; i++) {
        ST_int obs = active_set[i];
        y_sub[i] = y_full[obs];
        for (j = 0; j < K; j++) {
            X_sub[j * m + i] = X_full[j * N_full + obs];
        }
    }

    /* Create IPM state for subsample */
    cqreg_ipm_state *ipm_sub = cqreg_ipm_create(m, K, config);
    if (!ipm_sub) {
        free(y_sub);
        free(X_sub);
        return -1;
    }

    /* Solve on subsample */
    ST_int result = cqreg_ipm_solve(ipm_sub, y_sub, X_sub, q, beta);

    /* Cleanup */
    cqreg_ipm_free(ipm_sub);
    free(y_sub);
    free(X_sub);

    return result;
}

/*
 * Check optimality on full sample.
 * Returns number of violated observations (0 = optimal).
 * Also identifies which observations are violated.
 *
 * A solution β is optimal iff:
 *   1. The subgradient g is well-defined:
 *      g_i = q - I(r_i < 0)
 *      - For r_i > 0: g_i = q
 *      - For r_i < 0: g_i = q - 1
 *      - For r_i = 0: g_i ∈ [q-1, q] (need to solve for optimal g_i)
 *   2. X'g = 0 (dual feasibility / stationarity)
 *
 * This function checks:
 *   - Observations with sign changes from prediction (clear violations)
 *   - Observations where X'g ≠ 0 would be resolved by including them
 */
static ST_int check_optimality(
    const ST_double *y, const ST_double *X,
    ST_int N, ST_int K,
    const ST_double *beta, ST_double q,
    const ST_int *in_active_set,  /* Boolean: 1 if in active set */
    const ST_int *predicted_sign, /* 1 = positive, -1 = negative, 0 = uncertain */
    ST_int *violations,           /* Output: indices of violated obs */
    ST_int *n_violations)
{
    ST_int i, j;
    *n_violations = 0;

    /* Allocate space for residuals and subgradient */
    ST_double *r = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *g = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *Xg = (ST_double *)calloc(K, sizeof(ST_double));

    if (!r || !g || !Xg) {
        free(r);
        free(g);
        free(Xg);
        return -1;
    }

    /* Compute residuals for all observations and gather stats */
    ST_double sum_r2 = 0.0;
    for (i = 0; i < N; i++) {
        ST_double yhat = 0.0;
        for (j = 0; j < K; j++) {
            yhat += X[j * N + i] * beta[j];
        }
        r[i] = y[i] - yhat;
        sum_r2 += r[i] * r[i];
    }

    /* Compute residual scale for tolerance */
    ST_double sigma_r = sqrt(sum_r2 / N);
    ST_double zero_tol = 1e-6 * sigma_r;  /* Tolerance for "zero" residual */
    if (zero_tol < 1e-12) zero_tol = 1e-12;

    /* Compute subgradient */
    for (i = 0; i < N; i++) {
        if (r[i] > zero_tol) {
            g[i] = q;  /* Positive residual */
        } else if (r[i] < -zero_tol) {
            g[i] = q - 1.0;  /* Negative residual */
        } else {
            g[i] = 0.0;  /* At kink - will adjust later */
        }
    }

    /* Compute X'g */
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        for (i = 0; i < N; i++) {
            Xg[j] += Xj[i] * g[i];
        }
    }

    /* Check dual feasibility: X'g should be zero
     * Compute norm of X'g relative to scale */
    ST_double Xg_norm = 0.0;
    ST_double X_scale = 0.0;
    for (j = 0; j < K; j++) {
        Xg_norm += Xg[j] * Xg[j];
        /* Estimate scale from X'X diagonal */
        ST_double col_norm = 0.0;
        for (i = 0; i < N; i++) {
            col_norm += X[j * N + i] * X[j * N + i];
        }
        X_scale += col_norm;
    }
    Xg_norm = sqrt(Xg_norm);
    X_scale = sqrt(X_scale / K) * sqrt((double)N);

    /* If dual constraint is violated, find observations to add */
    ST_double dual_tol = 1e-6 * X_scale;

    if (Xg_norm > dual_tol) {
        /* Dual constraint violated - find ALL observations that could fix it.
         * Be aggressive: add any observation not in active set that either:
         * 1. Has a different sign than predicted (clear violation)
         * 2. Has uncertain prediction (predicted_sign == 0)
         * 3. Has small residual (near boundary)
         * 4. If still no violations, add batch of smallest |r| observations
         */

        /* Threshold for "small residual" - observations near the quantile boundary */
        ST_double boundary_tol = 0.1 * sigma_r;

        for (i = 0; i < N; i++) {
            /* Check observations not in active set */
            if (!in_active_set[i]) {
                /* Near-zero residual - near the quantile boundary */
                if (fabs(r[i]) < boundary_tol) {
                    violations[(*n_violations)++] = i;
                }
                /* Uncertain prediction - could go either way */
                else if (predicted_sign[i] == 0) {
                    violations[(*n_violations)++] = i;
                }
                /* Sign changed from prediction */
                else if (predicted_sign[i] > 0 && r[i] < -zero_tol) {
                    violations[(*n_violations)++] = i;
                }
                else if (predicted_sign[i] < 0 && r[i] > zero_tol) {
                    violations[(*n_violations)++] = i;
                }
            }
        }

        /* If still no violations found, add a batch of observations
         * with smallest |residual| to make progress */
        if (*n_violations == 0) {
            /* Add up to 10% more observations or at least 10 */
            ST_int batch_size = N / 10;
            if (batch_size < 10) batch_size = 10;
            if (batch_size > N / 2) batch_size = N / 2;

            /* Find observations with smallest |r| not in active set */
            for (ST_int added = 0; added < batch_size; added++) {
                ST_double min_abs_r = 1e30;
                ST_int min_idx = -1;
                for (i = 0; i < N; i++) {
                    if (!in_active_set[i]) {
                        /* Check if already added as violation */
                        ST_int already_added = 0;
                        for (ST_int v = 0; v < *n_violations; v++) {
                            if (violations[v] == i) {
                                already_added = 1;
                                break;
                            }
                        }
                        if (!already_added && fabs(r[i]) < min_abs_r) {
                            min_abs_r = fabs(r[i]);
                            min_idx = i;
                        }
                    }
                }
                if (min_idx >= 0) {
                    violations[(*n_violations)++] = min_idx;
                } else {
                    break;  /* No more observations to add */
                }
            }
        }
    }

    free(r);
    free(g);
    free(Xg);

    return *n_violations;
}

/*
 * Main preprocessing solver.
 */
ST_int cqreg_preprocess_solve(cqreg_ipm_state *ipm,
                               const ST_double *y,
                               const ST_double *X,
                               ST_double q,
                               ST_double *beta)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j;
    ST_int total_iters = 0;


    /* Step 1: Initial subsample size
     * Per Chernozhukov, Fernández-Val, Melly (2020), use:
     *   m = (n * (log(K) + 1))^(2/3)
     * This is the optimal rate for the preprocessing algorithm.
     * We add a safety factor of 1.5 since we use smoothed IPM (not exact LP).
     *
     * IMPORTANT: For quantiles far from 0.5, OLS predictions are less accurate.
     * Increase subsample size based on distance from median.
     * At extreme quantiles (0.1 or 0.9), use 3x the base size.
     */
    ST_double log_term = log((ST_double)K) + 1.0;
    ST_double m_raw = pow((ST_double)N * log_term, 2.0/3.0);

    /* Quantile adjustment factor: increases subsample for non-median quantiles */
    ST_double q_dist = fabs(q - 0.5);  /* Distance from median */
    ST_double q_factor = 1.5 + 3.0 * q_dist;  /* 1.5 at median, 3.0 at extreme */

    ST_int m_init = (ST_int)(q_factor * m_raw + 0.5);

    /* Ensure reasonable bounds:
     * - At least 10*K to have enough obs for stable regression
     * - At least 100 for numerical stability
     * - Cap at 80% of N to get meaningful speedup (otherwise just solve directly)
     */
    if (m_init < 10 * K) m_init = 10 * K;
    if (m_init < 100) m_init = 100;
    if (m_init > (ST_int)(0.8 * N)) m_init = N;  /* Fall through to full solve */
    if (m_init > N) m_init = N;


    /* Step 2: Compute OLS solution for initial residuals */
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
        ipm->XDX[j * K + j] += 1e-10;
    }

    /* Cholesky factorization */
    memcpy(ipm->L, ipm->XDX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(ipm->L, K) != 0) {
        /* Fallback: use zero initial beta */
        memset(beta, 0, K * sizeof(ST_double));
    } else {
        /* Compute X'y */
        for (j = 0; j < K; j++) {
            beta[j] = cqreg_dot(&X[j * N], y, N);
        }
        /* Solve for OLS beta */
        cqreg_solve_cholesky(ipm->L, beta, K);
    }

    /* Step 3: Compute residuals and sort by |r| */
    resid_pair *resid_order = (resid_pair *)malloc(N * sizeof(resid_pair));
    ST_int *active_set = (ST_int *)malloc(N * sizeof(ST_int));
    ST_int *in_active_set = (ST_int *)calloc(N, sizeof(ST_int));  /* Boolean: 1 if in active set */
    ST_int *predicted_sign = (ST_int *)malloc(N * sizeof(ST_int)); /* 1=pos, -1=neg, 0=uncertain */
    ST_int *violations = (ST_int *)malloc(N * sizeof(ST_int));

    if (!resid_order || !active_set || !in_active_set || !predicted_sign || !violations) {
        free(resid_order);
        free(active_set);
        free(in_active_set);
        free(predicted_sign);
        free(violations);
        return -1;
    }

    /* Compute OLS residuals.
     * Key insight from Chernozhukov et al.: The quantile regression hyperplane
     * will have approximately (1-q)*N observations with positive residuals.
     * OLS gives the conditional mean, so we need to adjust predictions.
     *
     * Strategy: Sort OLS residuals and use ranks to predict QR signs.
     * The top (1-q)*N by OLS residual rank should have positive QR residuals.
     */

    /* First pass: compute residuals and gather stats */
    ST_double *ols_resid = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double sum_r2 = 0.0;
    for (i = 0; i < N; i++) {
        ST_double yhat = 0.0;
        for (j = 0; j < K; j++) {
            yhat += X[j * N + i] * beta[j];
        }
        ST_double r = y[i] - yhat;
        ols_resid[i] = r;
        resid_order[i].idx = i;
        resid_order[i].abs_resid = fabs(r);
        sum_r2 += r * r;
    }

    /* Compute scale for uncertainty band */
    ST_double sigma_r = sqrt(sum_r2 / N);
    ST_double sign_tol = 0.05 * sigma_r;  /* 5% of std dev is "uncertain" */
    if (sign_tol < 1e-10) sign_tol = 1e-10;

    /* Find q-th quantile of OLS residuals using selection algorithm.
     * We need the value at rank q*N to determine the cutoff. */
    ST_double *resid_copy = (ST_double *)malloc(N * sizeof(ST_double));
    memcpy(resid_copy, ols_resid, N * sizeof(ST_double));

    /* Sort to find q-th quantile */
    qsort(resid_copy, N, sizeof(ST_double), compare_double);

    /* Find the qth quantile cutoff in OLS residuals.
     * Observations below this should have negative QR residuals,
     * observations above should have positive QR residuals. */
    ST_int q_rank = (ST_int)(q * N);
    if (q_rank >= N) q_rank = N - 1;
    ST_double q_cutoff = resid_copy[q_rank];

    /* Record predicted signs based on OLS residual relative to cutoff */
    for (i = 0; i < N; i++) {
        ST_double r = ols_resid[i];
        ST_double dist_from_cutoff = r - q_cutoff;

        if (dist_from_cutoff > sign_tol) {
            predicted_sign[i] = 1;  /* Predicted positive */
        } else if (dist_from_cutoff < -sign_tol) {
            predicted_sign[i] = -1; /* Predicted negative */
        } else {
            predicted_sign[i] = 0;  /* Uncertain - near cutoff */
        }
    }

    free(resid_copy);
    free(ols_resid);

    /* Sort by absolute residual (ascending) for active set selection */
    qsort(resid_order, N, sizeof(resid_pair), compare_resid_pair);

    /* Step 4: Initialize active set with m_init smallest |residual| observations */
    ST_int m = m_init;
    for (i = 0; i < m; i++) {
        active_set[i] = resid_order[i].idx;
        in_active_set[resid_order[i].idx] = 1;
    }

    /* Step 5: Iterative refinement loop */
    ST_int max_outer_iters = 20;
    for (ST_int outer_iter = 0; outer_iter < max_outer_iters; outer_iter++) {

        /* Solve QR on active set */
        ST_int iters = solve_subsample_qr(y, X, N, K, active_set, m, q, beta, &ipm->config);
        if (iters < 0) {
            /* Solver failed - fall back to full solve */
            free(resid_order);
            free(active_set);
            free(in_active_set);
            free(predicted_sign);
            free(violations);
            return cqreg_ipm_solve(ipm, y, X, q, beta);
        }
        total_iters += (iters > 0) ? iters : -iters;

        /* Check optimality on full sample */
        ST_int n_violations = 0;
        check_optimality(y, X, N, K, beta, q, in_active_set, predicted_sign, violations, &n_violations);

        if (n_violations == 0) {
            /* Optimal! */
            break;
        }

        /* Add violated observations to active set */
        for (i = 0; i < n_violations && m < N; i++) {
            ST_int obs = violations[i];
            if (!in_active_set[obs]) {
                active_set[m++] = obs;
                in_active_set[obs] = 1;
            }
        }

        /* If active set is now the full sample, just solve directly */
        if (m >= N) {
            free(resid_order);
            free(active_set);
            free(in_active_set);
            free(predicted_sign);
            free(violations);
            return cqreg_ipm_solve(ipm, y, X, q, beta);
        }
    }

    /* Update IPM state with final results */
    cqreg_matvec_col(ipm->r_primal, X, beta, N, K);
    for (i = 0; i < N; i++) {
        ST_double r = y[i] - ipm->r_primal[i];
        ipm->r_primal[i] = r;
        if (r >= 0) {
            ipm->u[i] = r;
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -r;
        }
    }

    /* Copy solution */
    cqreg_vcopy(ipm->beta, beta, K);
    ipm->iterations = total_iters;
    ipm->converged = 1;

    /* Cleanup */
    free(resid_order);
    free(active_set);
    free(in_active_set);
    free(predicted_sign);
    free(violations);

    return total_iters;
}
