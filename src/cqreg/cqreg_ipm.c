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
    ST_int i, j, k;

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

        /* Debug output */
        if (ipm->config.verbose && (iter % 10 == 0 || iter < 5)) {
            char buf[256];
            snprintf(buf, sizeof(buf), "IPM iter %3d: mu=%10.4e obj=%10.2f\n", iter, mu, obj);
            SF_display(buf);
        }

        /* Form weighted normal equations: (X'WX) Δβ = X'g
         * Using 8x loop unrolling for better CPU pipelining */
        memset(ipm->XDX, 0, K * K * sizeof(ST_double));
        ST_int N8 = N - (N & 7);  /* N rounded down to multiple of 8 */

        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];

            /* Diagonal element with 8x unrolling */
            ST_double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
            ST_double d4 = 0.0, d5 = 0.0, d6 = 0.0, d7 = 0.0;
            for (i = 0; i < N8; i += 8) {
                d0 += w[i]     * Xj[i]     * Xj[i];
                d1 += w[i + 1] * Xj[i + 1] * Xj[i + 1];
                d2 += w[i + 2] * Xj[i + 2] * Xj[i + 2];
                d3 += w[i + 3] * Xj[i + 3] * Xj[i + 3];
                d4 += w[i + 4] * Xj[i + 4] * Xj[i + 4];
                d5 += w[i + 5] * Xj[i + 5] * Xj[i + 5];
                d6 += w[i + 6] * Xj[i + 6] * Xj[i + 6];
                d7 += w[i + 7] * Xj[i + 7] * Xj[i + 7];
            }
            for (i = N8; i < N; i++) {
                d0 += w[i] * Xj[i] * Xj[i];
            }
            ipm->XDX[j * K + j] = ((d0 + d4) + (d1 + d5)) + ((d2 + d6) + (d3 + d7)) + 1e-10;

            /* Off-diagonal elements with 8x unrolling */
            for (k = j + 1; k < K; k++) {
                const ST_double *Xk = &X[k * N];
                ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
                ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
                for (i = 0; i < N8; i += 8) {
                    s0 += w[i]     * Xj[i]     * Xk[i];
                    s1 += w[i + 1] * Xj[i + 1] * Xk[i + 1];
                    s2 += w[i + 2] * Xj[i + 2] * Xk[i + 2];
                    s3 += w[i + 3] * Xj[i + 3] * Xk[i + 3];
                    s4 += w[i + 4] * Xj[i + 4] * Xk[i + 4];
                    s5 += w[i + 5] * Xj[i + 5] * Xk[i + 5];
                    s6 += w[i + 6] * Xj[i + 6] * Xk[i + 6];
                    s7 += w[i + 7] * Xj[i + 7] * Xk[i + 7];
                }
                for (i = N8; i < N; i++) {
                    s0 += w[i] * Xj[i] * Xk[i];
                }
                ST_double total = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
                ipm->XDX[j * K + k] = total;
                ipm->XDX[k * K + j] = total;
            }
        }

        /* Compute X'g with 8x unrolling */
        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];
            ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
            for (i = 0; i < N8; i += 8) {
                s0 += Xj[i]     * g[i];
                s1 += Xj[i + 1] * g[i + 1];
                s2 += Xj[i + 2] * g[i + 2];
                s3 += Xj[i + 3] * g[i + 3];
                s4 += Xj[i + 4] * g[i + 4];
                s5 += Xj[i + 5] * g[i + 5];
                s6 += Xj[i + 6] * g[i + 6];
                s7 += Xj[i + 7] * g[i + 7];
            }
            for (i = N8; i < N; i++) {
                s0 += Xj[i] * g[i];
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
         * New residual = old_residual - alpha * X * delta_beta */
        ST_double alpha = 1.0;
        ST_double obj_new;
        for (ST_int ls = 0; ls < 20; ls++) {
            /* Compute new objective using cached X*delta_beta */
            obj_new = 0.0;
            for (i = 0; i < N; i++) {
                ST_double r = ipm->r_primal[i] - alpha * X_delta_beta[i];
                ST_double z = -r * inv_mu;
                if (z > 20.0) {
                    obj_new += q * r + mu * z;
                } else if (z < -20.0) {
                    obj_new += q * r;
                } else {
                    obj_new += q * r + mu * log(1.0 + exp(z));
                }
            }

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

    /* Debug: print solution statistics */
    if (ipm->config.verbose) {
        ST_double sum_u = 0.0, sum_v = 0.0, sum_u_lam_u = 0.0, sum_v_lam_v = 0.0;
        ST_double min_u = 1e30, min_v = 1e30, max_u = 0.0, max_v = 0.0;
        ST_double min_lu = 1e30, min_lv = 1e30, max_lu = 0.0, max_lv = 0.0;
        for (ST_int i = 0; i < N; i++) {
            sum_u += ipm->u[i];
            sum_v += ipm->v[i];
            sum_u_lam_u += ipm->u[i] * ipm->lambda_u[i];
            sum_v_lam_v += ipm->v[i] * ipm->lambda_v[i];
            if (ipm->u[i] < min_u) min_u = ipm->u[i];
            if (ipm->v[i] < min_v) min_v = ipm->v[i];
            if (ipm->u[i] > max_u) max_u = ipm->u[i];
            if (ipm->v[i] > max_v) max_v = ipm->v[i];
            if (ipm->lambda_u[i] < min_lu) min_lu = ipm->lambda_u[i];
            if (ipm->lambda_v[i] < min_lv) min_lv = ipm->lambda_v[i];
            if (ipm->lambda_u[i] > max_lu) max_lu = ipm->lambda_u[i];
            if (ipm->lambda_v[i] > max_lv) max_lv = ipm->lambda_v[i];
        }
        char buf[512];
        snprintf(buf, sizeof(buf), "Solution stats: sum_u=%.2f sum_v=%.2f sum_u+v=%.2f\n", sum_u, sum_v, sum_u + sum_v);
        SF_display(buf);
        snprintf(buf, sizeof(buf), "  u: min=%.6e max=%.2f  v: min=%.6e max=%.2f\n", min_u, max_u, min_v, max_v);
        SF_display(buf);
        snprintf(buf, sizeof(buf), "  lu: min=%.6e max=%.6f  lv: min=%.6e max=%.6f\n", min_lu, max_lu, min_lv, max_lv);
        SF_display(buf);
        snprintf(buf, sizeof(buf), "  compl: sum(u*lu)=%.6e sum(v*lv)=%.6e\n", sum_u_lam_u, sum_v_lam_v);
        SF_display(buf);
    }

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

    /* Allocate subsample data */
    ST_double *y_sub = (ST_double *)malloc(m * sizeof(ST_double));
    ST_double *X_sub = (ST_double *)malloc(m * K * sizeof(ST_double));
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

    /* Compute residuals and subgradient for all observations */
    for (i = 0; i < N; i++) {
        ST_double yhat = 0.0;
        for (j = 0; j < K; j++) {
            yhat += X[j * N + i] * beta[j];
        }
        r[i] = y[i] - yhat;

        /* Subgradient: g_i = q - I(r_i < 0) */
        ST_double tol = 1e-8 * (fabs(y[i]) + 1.0);
        if (r[i] > tol) {
            g[i] = q;  /* Positive residual */
        } else if (r[i] < -tol) {
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
        /* Dual constraint violated - find observations whose sign might be wrong
         * or near-zero observations that should be in active set */

        for (i = 0; i < N; i++) {
            ST_double tol = 1e-6 * (fabs(y[i]) + 1.0);

            /* Check observations not in active set */
            if (!in_active_set[i]) {
                /* Near-zero residual - should be in active set */
                if (fabs(r[i]) < tol * 100) {  /* Wider tolerance for boundary obs */
                    violations[(*n_violations)++] = i;
                }
                /* Sign changed from prediction */
                else if (predicted_sign[i] > 0 && r[i] < -tol) {
                    violations[(*n_violations)++] = i;
                }
                else if (predicted_sign[i] < 0 && r[i] > tol) {
                    violations[(*n_violations)++] = i;
                }
            }
        }

        /* If no violations found but X'g ≠ 0, add obs with smallest |residual| not in active set */
        if (*n_violations == 0) {
            /* Find obs with smallest non-zero |residual| not in active set */
            ST_double min_abs_r = 1e30;
            ST_int min_idx = -1;
            for (i = 0; i < N; i++) {
                if (!in_active_set[i]) {
                    ST_double abs_r = fabs(r[i]);
                    if (abs_r < min_abs_r) {
                        min_abs_r = abs_r;
                        min_idx = i;
                    }
                }
            }
            if (min_idx >= 0) {
                violations[(*n_violations)++] = min_idx;
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

    if (ipm->config.verbose) {
        SF_display("Preprocessing QR solver (Chernozhukov et al. 2020)\n");
    }

    /* Step 1: Initial subsample size
     * For smoothed IPM, we need a larger subsample than exact LP methods.
     * Use m = max(N^0.7, 10*K) to ensure good coverage.
     */
    ST_int m_init = (ST_int)(pow((double)N, 0.7) + 0.5);

    /* Ensure reasonable bounds */
    if (m_init < 10 * K) m_init = 10 * K;
    if (m_init > N / 2) m_init = N / 2;  /* Cap at 50% to get speedup */
    if (m_init > N) m_init = N;

    if (ipm->config.verbose) {
        char buf[128];
        snprintf(buf, sizeof(buf), "  Initial subsample size: %d (of %d, %.1f%%)\n",
                 m_init, N, 100.0 * m_init / N);
        SF_display(buf);
    }

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

    /* Compute OLS residuals and record their signs */
    for (i = 0; i < N; i++) {
        ST_double yhat = 0.0;
        for (j = 0; j < K; j++) {
            yhat += X[j * N + i] * beta[j];
        }
        ST_double r = y[i] - yhat;
        resid_order[i].idx = i;
        resid_order[i].abs_resid = fabs(r);

        /* Record predicted sign from OLS residuals */
        ST_double tol = 1e-6 * (fabs(y[i]) + 1.0);
        if (r > tol) {
            predicted_sign[i] = 1;  /* Positive */
        } else if (r < -tol) {
            predicted_sign[i] = -1; /* Negative */
        } else {
            predicted_sign[i] = 0;  /* Uncertain - near zero */
        }
    }

    /* Sort by absolute residual (ascending) */
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

        if (ipm->config.verbose) {
            char buf[128];
            snprintf(buf, sizeof(buf), "  Outer iteration %d: active set size = %d\n",
                     outer_iter + 1, m);
            SF_display(buf);
        }

        /* Solve QR on active set */
        ST_int iters = solve_subsample_qr(y, X, N, K, active_set, m, q, beta, &ipm->config);
        if (iters < 0) {
            /* Solver failed - fall back to full solve */
            if (ipm->config.verbose) {
                SF_display("  Subsample solve failed, falling back to full solve\n");
            }
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
            if (ipm->config.verbose) {
                char buf[128];
                snprintf(buf, sizeof(buf), "  Converged in %d outer iterations, %d total IPM iterations\n",
                         outer_iter + 1, total_iters);
                SF_display(buf);
            }
            break;
        }

        if (ipm->config.verbose) {
            char buf[128];
            snprintf(buf, sizeof(buf), "  Found %d violations\n", n_violations);
            SF_display(buf);
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
            if (ipm->config.verbose) {
                SF_display("  Active set expanded to full sample\n");
            }
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
