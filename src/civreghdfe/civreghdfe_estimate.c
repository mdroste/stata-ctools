/*
    civreghdfe_estimate.c
    Core IV Estimation: 2SLS, LIML, Fuller, GMM2S, CUE

    This module implements the k-class family of IV estimators.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "civreghdfe_estimate.h"
#include "civreghdfe_impl.h"
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"
#include "civreghdfe_tests.h"
#include "civreghdfe_ffirst.h"
#include "../ctools_config.h"
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */

/* OpenMP for parallel residual computation */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Convenience aliases for the matrix functions */
#define matmul_atb  ctools_matmul_atb
#define matmul_ab   ctools_matmul_ab
#define matmul_atdb ctools_matmul_atdb
#define compute_liml_lambda civreghdfe_compute_liml_lambda

/* Shared OLS functions */
#include "../ctools_ols.h"

/*
    Compute GMM2S (two-step efficient GMM) estimator (ivreg2's s_egmm).

    Step 1: the caller builds S, the moment covariance of the chosen VCE,
            from the 2SLS residuals (civ_omega_build)
    Step 2: W = S^-1; beta = (X'ZWZ'X)^-1 X'ZWZ'y
*/
ST_retcode ivest_compute_gmm2s(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *S,
    ST_double *beta,
    ST_double *resid,
    ST_double *XZWZX_inv_out
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    ST_int K_total = ctx->K_total;
    ST_int i, j, k;

    /* Optimal weighting matrix W = S^-1 */
    ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ST_double *ZOmegaZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));

    if (!ZOmegaZ || !ZOmegaZ_inv) {
        free(ZOmegaZ); free(ZOmegaZ_inv);
        return 920;
    }
    memcpy(ZOmegaZ, S, (size_t)K_iv * K_iv * sizeof(ST_double));

    /* Invert Z'ΩZ */
    memcpy(ZOmegaZ_inv, ZOmegaZ, K_iv * K_iv * sizeof(ST_double));
    if (ctools_cholesky(ZOmegaZ_inv, K_iv) != 0) {
        if (ctx->report) SF_error("civreghdfe: GMM2S optimal weighting matrix is singular\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        return 198;
    }
    stata_retcode inverse_status_0 = ctools_invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv);
    if (inverse_status_0 != STATA_OK) {
        if (ctx->report) SF_error("civreghdfe: Failed to invert GMM2S weighting matrix\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        return ctools_stata_rc(inverse_status_0);
    }

    /* Compute GMM estimator: β = (X'ZWZ'X)^-1 X'ZWZ'y */
    ST_double *XZW = (ST_double *)calloc((size_t)K_total * K_iv, sizeof(ST_double));
    ST_double *XZWZX = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
    ST_double *XZWZy = (ST_double *)calloc((size_t)K_total, sizeof(ST_double));
    ST_double *XZWZX_L = (ST_double *)malloc((size_t)K_total * K_total * sizeof(ST_double));
    ST_double *beta_temp = (ST_double *)calloc((size_t)K_total, sizeof(ST_double));
    if (!XZW || !XZWZX || !XZWZy || !XZWZX_L || !beta_temp) {
        free(ZOmegaZ); free(ZOmegaZ_inv);
        free(XZW); free(XZWZX); free(XZWZy);
        free(XZWZX_L); free(beta_temp);
        return 920;
    }

    /* 1. XZW = (Z'X)' * W (K_total x K_iv) */
    for (i = 0; i < K_total; i++) {
        for (j = 0; j < K_iv; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += ctx->ZtX[i * K_iv + k] * ZOmegaZ_inv[j * K_iv + k];
            }
            XZW[j * K_total + i] = sum;
        }
    }

    /* 2. XZWZX = XZW * Z'X (K_total x K_total) */
    for (i = 0; i < K_total; i++) {
        for (j = 0; j < K_total; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += XZW[k * K_total + i] * ctx->ZtX[j * K_iv + k];
            }
            XZWZX[j * K_total + i] = sum;
        }
    }

    /* 3. XZWZy = XZW * Z'y (K_total x 1) */
    for (i = 0; i < K_total; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += XZW[k * K_total + i] * ctx->Zty[k];
        }
        XZWZy[i] = sum;
    }

    /* 4. Solve XZWZX * beta = XZWZy */
    memcpy(XZWZX_L, XZWZX, K_total * K_total * sizeof(ST_double));

    if (ctools_cholesky(XZWZX_L, K_total) != 0) {
        if (ctx->report) SF_error("civreghdfe: GMM2S X'ZWZ'X matrix is singular\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        free(XZW); free(XZWZX); free(XZWZy);
        free(XZWZX_L); free(beta_temp);
        return 198;
    }

    /* Forward substitution */
    for (i = 0; i < K_total; i++) {
        ST_double sum = XZWZy[i];
        for (j = 0; j < i; j++) {
            sum -= XZWZX_L[i * K_total + j] * beta_temp[j];
        }
        beta_temp[i] = sum / XZWZX_L[i * K_total + i];
    }

    /* Backward substitution */
    for (i = K_total - 1; i >= 0; i--) {
        ST_double sum = beta_temp[i];
        for (j = i + 1; j < K_total; j++) {
            sum -= XZWZX_L[j * K_total + i] * beta[j];
        }
        beta[i] = sum / XZWZX_L[i * K_total + i];
    }

    /* Compute residuals */
    if (resid) {
        #pragma omp parallel for schedule(static) if(N > 10000)
        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            /* the index is declared here: one shared across threads would race */
            for (ST_int kk = 0; kk < K_total; kk++) {
                pred += ctx->X_all[(size_t)kk * N + i] * beta[kk];
            }
            resid[i] = y[i] - pred;
        }
    }

    /* Output the GMM Hessian inverse (X'ZWZ'X)^-1 if requested; beta is
       already the GMM estimate, so a failure here is an error (the caller
       must not fall back to the 2SLS VCE) */
    ST_retcode out_rc = STATA_OK;
    if (XZWZX_inv_out) {
        ST_double *XZWZX_inv = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
        if (!XZWZX_inv) {
            out_rc = 920;
        } else {
            /* Re-compute Cholesky and inverse from XZWZX */
            memcpy(XZWZX_inv, XZWZX, K_total * K_total * sizeof(ST_double));
            if (ctools_cholesky(XZWZX_inv, K_total) == 0) {
                stata_retcode inverse_status = ctools_invert_from_cholesky(XZWZX_inv, K_total, XZWZX_inv);
                if (inverse_status != STATA_OK)
                    out_rc = ctools_stata_rc(inverse_status);
                else
                    memcpy(XZWZX_inv_out, XZWZX_inv, K_total * K_total * sizeof(ST_double));
            }
            free(XZWZX_inv);
        }
    }

    free(ZOmegaZ); free(ZOmegaZ_inv);
    free(XZW); free(XZWZX); free(XZWZy);
    free(XZWZX_L); free(beta_temp);

    return out_rc;
}

/*
    Helper: Compute residuals given beta
*/
static void cue_compute_residuals(
    const ST_double *y,
    const ST_double *X_all,
    const ST_double *beta,
    ST_int N,
    ST_int K_total,
    ST_double *resid
)
{
    ST_int i;
    #pragma omp parallel for schedule(static) if(N > 10000)
    for (i = 0; i < N; i++) {
        ST_double pred = 0.0;
        for (ST_int k = 0; k < K_total; k++) {
            pred += X_all[(size_t)k * N + i] * beta[k];
        }
        resid[i] = y[i] - pred;
    }
}

/*
    Helper: CUE objective Q(b) = g' S(b)^-1 g / N with g = Z'We, where S(b) is
    the moment covariance of the chosen VCE (ctx->om: robust, cluster,
    two-way, HAC, AC, Driscoll-Kraay, Kiefer) at the residuals of b.
    N * Q is ivreg2's m_cuecrit J.
*/
static ST_double cue_objective(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *beta,
    ST_double *work_resid  /* N-sized work buffer */
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    const ST_double *Z = ctx->Z;

    /* Compute residuals */
    cue_compute_residuals(y, ctx->X_all, beta, N, ctx->K_total, work_resid);

    /* g = Z'We (weighted moment conditions), then S(b)^-1 g. When memory
       runs out, ctx->alloc_failed is set so the optimizer stops with r(920)
       instead of taking the large value as a poor step. */
    ST_double *g = (ST_double *)calloc((size_t)K_iv * (K_iv + 2), sizeof(ST_double));
    if (!g) {
        #pragma omp atomic write
        ctx->alloc_failed = 1;
        return 1e30;
    }
    ST_double *Sg = g + K_iv, *S = Sg + K_iv;
    for (ST_int j = 0; j < K_iv; j++) {
        const ST_double *zj = Z + (size_t)j * N;
        ST_double s = 0.0;
        for (ST_int i = 0; i < N; i++) {
            ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
            s += zj[i] * w * work_resid[i];
        }
        g[j] = s;
    }

    ST_double Q = 1e30;  /* returned if S cannot be built or solved */
    if (civ_omega_build(ctx->om, work_resid, 1, Z, K_iv, S) == 0 &&
        civ_sym_solve(S, K_iv, g, 1, Sg) == 0) {
        ST_double gWinvg = 0.0;
        for (ST_int j = 0; j < K_iv; j++) gWinvg += g[j] * Sg[j];
        Q = gWinvg / N;
    } else {
        /* both fail only when memory runs out */
        #pragma omp atomic write
        ctx->alloc_failed = 1;
    }

    free(g);
    return Q;
}

/*
    Eigenvalues (into w) and eigenvectors (columns of V) of the symmetric
    K x K matrix A (column-major, overwritten) by cyclic Jacobi sweeps until
    the off-diagonal mass is below 1e-15 of the matrix's. Returns 0, or -1
    when the sweeps do not converge or a value is not finite.
*/
static int cue_sym_eigen(ST_int K, ST_double *A, ST_double *V, ST_double *w)
{
    ST_double total = 0.0;
    for (ST_int j = 0; j < K; j++) {
        for (ST_int i = 0; i < K; i++) {
            V[(size_t)j * K + i] = (i == j) ? 1.0 : 0.0;
            total += A[(size_t)j * K + i] * A[(size_t)j * K + i];
        }
    }
    if (!isfinite(total)) return -1;
    int done = 0;
    for (int sweep = 0; sweep < 100 && !done; sweep++) {
        ST_double off = 0.0;
        for (ST_int j = 0; j < K; j++)
            for (ST_int i = 0; i < j; i++)
                off += A[(size_t)j * K + i] * A[(size_t)j * K + i];
        if (off <= 1e-30 * total) {
            done = 1;
            break;
        }
        for (ST_int p = 0; p < K - 1; p++) {
            for (ST_int q = p + 1; q < K; q++) {
                const ST_double apq = A[(size_t)q * K + p];
                if (apq == 0.0) continue;
                const ST_double app = A[(size_t)p * K + p], aqq = A[(size_t)q * K + q];
                const ST_double theta = 0.5 * atan2(2.0 * apq, aqq - app);
                const ST_double c = cos(theta), s = sin(theta);
                for (ST_int i = 0; i < K; i++) {
                    const ST_double aip = A[(size_t)p * K + i], aiq = A[(size_t)q * K + i];
                    A[(size_t)p * K + i] = c * aip - s * aiq;
                    A[(size_t)q * K + i] = s * aip + c * aiq;
                }
                for (ST_int i = 0; i < K; i++) {
                    const ST_double api = A[(size_t)i * K + p], aqi = A[(size_t)i * K + q];
                    A[(size_t)i * K + p] = c * api - s * aqi;
                    A[(size_t)i * K + q] = s * api + c * aqi;
                }
                for (ST_int i = 0; i < K; i++) {
                    const ST_double vip = V[(size_t)p * K + i], viq = V[(size_t)q * K + i];
                    V[(size_t)p * K + i] = c * vip - s * viq;
                    V[(size_t)q * K + i] = s * vip + c * viq;
                }
            }
        }
    }
    for (ST_int i = 0; i < K; i++) {
        w[i] = A[(size_t)i * K + i];
        if (!isfinite(w[i])) return -1;
    }
    return done ? 0 : -1;
}

/*
    The CUE step where the Hessian H is not positive definite (a nonconvex
    region, typical at the 2SLS start of a weakly identified model): the
    Newton step on H with its eigenvalues replaced by their absolute values
    (floored at 1e-8 of the largest), in the coefficients' units u, so that
    directions of negative curvature are descended at their own scale
    instead of by a normalized steepest-descent step. Returns 1 with the
    direction in d, or 0 (H zero, the eigensystem failed, or no memory).
*/
static int cue_modified_newton(ST_int K, const ST_double *H, const ST_double *g,
                               const ST_double *u, ST_double *d)
{
    ST_double *A = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));
    ST_double *V = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));
    ST_double *w = (ST_double *)malloc((size_t)K * 2 * sizeof(ST_double));
    int ok = 0;
    if (A && V && w) {
        for (ST_int j = 0; j < K; j++)
            for (ST_int i = 0; i < K; i++)
                A[(size_t)j * K + i] = u[i] * H[(size_t)j * K + i] * u[j];
        if (cue_sym_eigen(K, A, V, w) == 0) {
            ST_double lmax = 0.0;
            for (ST_int i = 0; i < K; i++)
                if (fabs(w[i]) > lmax) lmax = fabs(w[i]);
            if (lmax > 0.0 && isfinite(lmax)) {
                ST_double *dt = w + K;
                for (ST_int i = 0; i < K; i++) dt[i] = 0.0;
                for (ST_int e = 0; e < K; e++) {
                    const ST_double *ve = V + (size_t)e * K;
                    ST_double c = 0.0;
                    for (ST_int i = 0; i < K; i++) c += ve[i] * u[i] * g[i];
                    ST_double lam = fabs(w[e]);
                    if (lam < 1e-8 * lmax) lam = 1e-8 * lmax;
                    for (ST_int i = 0; i < K; i++) dt[i] -= ve[i] * c / lam;
                }
                ok = 1;
                for (ST_int i = 0; i < K; i++) {
                    d[i] = u[i] * dt[i];
                    if (!isfinite(d[i])) ok = 0;
                }
            }
        }
    }
    free(A);
    free(V);
    free(w);
    return ok;
}

/*
    Compute CUE (Continuously Updated Estimator).

    CUE minimizes: Q(β) = g(β)' S(β)^-1 g(β)
    where g(β) = Z'W(y - Xβ) and S(β) is the moment covariance of the chosen
    VCE (ctx->om) at the residuals of β: robust, cluster, two-way, HAC, AC,
    Driscoll-Kraay or Kiefer, as ivreg2's m_cuecrit calls m_omega.

    Uses Newton's method with numerical derivatives to minimize Q(β).
    Step sizes scale with max(|β_k|, u_k) (matching Stata's optimize()), and
    convergence is tested on the gradient in units u_k = rms(y)/rms(x_k), so
    the iterations do not depend on the units of y or of the regressors
    (Q itself is invariant to them).
    For homoskedastic CUE (no kernel), the caller should use LIML directly.
*/
ST_retcode ivest_compute_cue(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *initial_beta,
    ST_double *beta,
    ST_double *resid,
    ST_int max_iter,
    ST_double tol,
    ST_double *XZWZX_inv_out
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    ST_int K_total = ctx->K_total;
    const ST_double *Z = ctx->Z;
    ST_int i, j, k;

    /* Initialize beta from input (2SLS estimate) */
    memcpy(beta, initial_beta, K_total * sizeof(ST_double));

    /* Allocate working arrays */
    ST_double *current_resid = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *gradient = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *hessian = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *direction = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *unit = (ST_double *)malloc(K_total * sizeof(ST_double));

    if (!current_resid || !gradient || !hessian || !direction || !unit) {
        free(current_resid); free(gradient);
        free(hessian); free(direction); free(unit);
        return 920;
    }

    /* Natural unit of each coefficient: rms(y)/rms(x_k) */
    {
        ST_double yy = 0.0;
        for (i = 0; i < N; i++) yy += y[i] * y[i];
        for (k = 0; k < K_total; k++) {
            const ST_double *xk = ctx->X_all + (size_t)k * N;
            ST_double xx = 0.0;
            for (i = 0; i < N; i++) xx += xk[i] * xk[i];
            unit[k] = (yy > 0.0 && xx > 0.0) ? sqrt(yy / xx) : 1.0;
            if (!isfinite(unit[k]) || unit[k] <= 0.0) unit[k] = 1.0;
        }
    }

    /* Pre-allocate per-thread buffers for parallel CUE gradient/Hessian.
       Each thread needs its own beta_test and work_resid buffers since
       cue_objective modifies work_resid. */
    int cue_nthreads = 1;
#ifdef _OPENMP
    /* A bootstrap replicate (its moment covariance carries the draw counts)
       runs on a worker of the replicate scheduler: one thread, so that it
       opens no nested team and allocates one set of buffers */
    cue_nthreads = (ctx->om && ctx->om->mult) ? 1 : ctools_get_max_threads();
    if (cue_nthreads > K_total) cue_nthreads = K_total;
    if (cue_nthreads < 1) cue_nthreads = 1;
#endif
    ST_double **thread_beta = (ST_double **)malloc(cue_nthreads * sizeof(ST_double *));
    ST_double **thread_resid = (ST_double **)malloc(cue_nthreads * sizeof(ST_double *));
    ST_double **thread_grad = (ST_double **)malloc(cue_nthreads * sizeof(ST_double *));
    if (!thread_beta || !thread_resid || !thread_grad) {
        free(current_resid); free(gradient);
        free(hessian); free(direction); free(unit);
        free(thread_beta); free(thread_resid); free(thread_grad);
        return 920;
    }
    int alloc_ok = 1;
    for (int t = 0; t < cue_nthreads; t++) {
        thread_beta[t] = (ST_double *)malloc(K_total * sizeof(ST_double));
        thread_resid[t] = (ST_double *)malloc(N * sizeof(ST_double));
        thread_grad[t] = (ST_double *)calloc(K_total, sizeof(ST_double));
        if (!thread_beta[t] || !thread_resid[t] || !thread_grad[t]) alloc_ok = 0;
    }
    if (!alloc_ok) {
        for (int t = 0; t < cue_nthreads; t++) {
            if (thread_beta[t]) free(thread_beta[t]);
            if (thread_resid[t]) free(thread_resid[t]);
            if (thread_grad[t]) free(thread_grad[t]);
        }
        free(thread_beta); free(thread_resid); free(thread_grad);
        free(current_resid); free(gradient);
        free(hessian); free(direction); free(unit);
        return 920;
    }

    if (ctx->verbose) {
        SF_display("civreghdfe: Computing CUE (Newton's method)\n");
    }

    /* cue_objective sets ctx->alloc_failed when memory runs out */
    ST_retcode out_rc = STATA_OK;
    ctx->alloc_failed = 0;

    /* Compute initial objective */
    ST_double Q_current = cue_objective(ctx, y, beta, current_resid);
    if (ctx->alloc_failed) {
        out_rc = 920;
        goto cue_done;
    }

    if (ctx->verbose) {
        char buf[128];
        snprintf(buf, sizeof(buf), "civreghdfe: CUE initial Q = %.8e\n", Q_current);
        SF_display(buf);
    }

    ST_int iter;

    for (iter = 0; iter < max_iter; iter++) {
        /* Adaptive step sizes (matching Stata's optimize() d0 convention) */
        /* eps_g ≈ max(|β_k|, u_k) * eps^(1/3) for gradient */
        /* eps_h ≈ max(|β_k|, u_k) * eps^(1/4) for Hessian */

        /* Compute gradient using central differences with adaptive step.
           Each k is independent — parallelize with per-thread buffers. */
        #pragma omp parallel for schedule(dynamic) num_threads(cue_nthreads) if(K_total > 1)
        for (k = 0; k < K_total; k++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            ST_double *my_beta = thread_beta[tid];
            ST_double *my_resid = thread_resid[tid];

            ST_double scale = fabs(beta[k]);
            if (scale < unit[k]) scale = unit[k];
            ST_double eps_g = scale * 6e-6;  /* eps^(1/3) ≈ 6e-6 */

            memcpy(my_beta, beta, K_total * sizeof(ST_double));

            my_beta[k] = beta[k] + eps_g;
            ST_double Q_plus = cue_objective(ctx, y, my_beta, my_resid);

            my_beta[k] = beta[k] - eps_g;
            ST_double Q_minus = cue_objective(ctx, y, my_beta, my_resid);

            gradient[k] = (Q_plus - Q_minus) / (2.0 * eps_g);
        }
        if (ctx->alloc_failed) {
            out_rc = 920;
            goto cue_done;
        }

        /* Gradient norm in the coefficients' units */
        ST_double grad_norm = 0.0;
        for (k = 0; k < K_total; k++) {
            grad_norm += gradient[k] * unit[k] * gradient[k] * unit[k];
        }
        grad_norm = sqrt(grad_norm);

        if (grad_norm < tol) {
            if (ctx->verbose) {
                char buf[128];
                snprintf(buf, sizeof(buf), "civreghdfe: CUE converged in %d iterations (grad=%.2e)\n", iter, grad_norm);
                SF_display(buf);
            }
            break;
        }

        /* Compute Hessian using forward differences on gradient with adaptive step.
           Each column j is independent — parallelize with per-thread buffers. */
        #pragma omp parallel for schedule(dynamic) num_threads(cue_nthreads) if(K_total > 1)
        for (j = 0; j < K_total; j++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            ST_double *my_beta = thread_beta[tid];
            ST_double *my_resid = thread_resid[tid];
            ST_double *my_grad2 = thread_grad[tid];

            ST_double scale_j = fabs(beta[j]);
            if (scale_j < unit[j]) scale_j = unit[j];
            ST_double eps_h = scale_j * 1.2e-4;  /* eps^(1/4) ≈ 1.2e-4 */

            memcpy(my_beta, beta, K_total * sizeof(ST_double));
            my_beta[j] = beta[j] + eps_h;

            /* Compute gradient at shifted point */
            for (ST_int kk = 0; kk < K_total; kk++) {
                ST_double scale_k = fabs(my_beta[kk]);
                if (scale_k < unit[kk]) scale_k = unit[kk];
                ST_double eps_g = scale_k * 6e-6;

                ST_double saved = my_beta[kk];
                my_beta[kk] = saved + eps_g;
                ST_double Q_plus = cue_objective(ctx, y, my_beta, my_resid);

                my_beta[kk] = saved - eps_g;
                ST_double Q_minus = cue_objective(ctx, y, my_beta, my_resid);

                my_beta[kk] = saved;
                my_grad2[kk] = (Q_plus - Q_minus) / (2.0 * eps_g);
            }

            /* H[:,j] = (gradient2 - gradient) / eps_h */
            for (ST_int kk = 0; kk < K_total; kk++) {
                hessian[j * K_total + kk] = (my_grad2[kk] - gradient[kk]) / eps_h;
            }
        }
        if (ctx->alloc_failed) {
            out_rc = 920;
            goto cue_done;
        }

        /* Symmetrize Hessian */
        for (j = 0; j < K_total; j++) {
            for (k = j + 1; k < K_total; k++) {
                ST_double avg = 0.5 * (hessian[j * K_total + k] + hessian[k * K_total + j]);
                hessian[j * K_total + k] = avg;
                hessian[k * K_total + j] = avg;
            }
        }

        /* Solve H * direction = -gradient for Newton direction */
        ST_double *H_copy = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
        if (!H_copy) {
            for (int t = 0; t < cue_nthreads; t++) {
                free(thread_beta[t]); free(thread_resid[t]); free(thread_grad[t]);
            }
            free(thread_beta); free(thread_resid); free(thread_grad);
            free(current_resid); free(gradient);
            free(hessian); free(direction); free(unit);
            return 920;
        }
        memcpy(H_copy, hessian, K_total * K_total * sizeof(ST_double));

        ST_int newton_ok = 0;
        if (ctools_cholesky(H_copy, K_total) == 0) {
            out_rc = ctools_stata_rc(ctools_invert_from_cholesky(H_copy, K_total, H_copy));
            if (out_rc) {
                free(H_copy);
                goto cue_done;
            }
            /* direction = -H^{-1} * gradient */
            for (j = 0; j < K_total; j++) {
                direction[j] = 0.0;
                for (k = 0; k < K_total; k++) {
                    direction[j] -= H_copy[k * K_total + j] * gradient[k];
                }
            }
            newton_ok = 1;
        }
        free(H_copy);

        ST_int modified = 0;
        if (!newton_ok && cue_modified_newton(K_total, hessian, gradient, unit, direction)) {
            /* Hessian not positive definite: the absolute-eigenvalue
               Newton step */
            newton_ok = 1;
            modified = 1;
        }
        if (!newton_ok) {
            /* No usable curvature at all: steepest descent with scaling */
            ST_double inv_grad_norm = 1.0 / (grad_norm + 1e-30);
            for (k = 0; k < K_total; k++) {
                direction[k] = -gradient[k] * unit[k] * unit[k] * inv_grad_norm;
            }
        }

        /* Backtracking line search along Newton direction */
        ST_double step = 1.0;
        ST_double armijo_c = 1e-4;
        ST_int line_search_iter = 0;
        ST_int max_ls_iter = 40;
        ST_double Q_new;

        /* Compute directional derivative for Armijo: g'd */
        ST_double dirderiv = 0.0;
        for (k = 0; k < K_total; k++) {
            dirderiv += gradient[k] * direction[k];
        }
        /* Ensure descent direction */
        if (dirderiv > 0) {
            for (k = 0; k < K_total; k++) {
                direction[k] = -gradient[k] * unit[k] * unit[k];
            }
            dirderiv = -grad_norm * grad_norm;
        }

        /* Use thread_beta[0] as scratch for line search (single-threaded) */
        ST_double *beta_test = thread_beta[0];
        while (line_search_iter < max_ls_iter) {
            for (k = 0; k < K_total; k++) {
                beta_test[k] = beta[k] + step * direction[k];
            }

            Q_new = cue_objective(ctx, y, beta_test, current_resid);
            if (ctx->alloc_failed) {
                out_rc = 920;
                goto cue_done;
            }

            /* Armijo condition */
            if (Q_new <= Q_current + armijo_c * step * dirderiv) {
                break;
            }

            step *= 0.5;
            line_search_iter++;
        }

        if (line_search_iter >= max_ls_iter) {
            if (ctx->verbose) {
                SF_display("civreghdfe: CUE line search failed, stopping\n");
            }
            break;
        }

        /* No decrease: the gradient is at its rounding floor */
        if (!(Q_new < Q_current)) {
            if (ctx->verbose) {
                SF_display("civreghdfe: CUE objective cannot decrease further, stopping\n");
            }
            break;
        }

        /* Update beta */
        memcpy(beta, beta_test, K_total * sizeof(ST_double));
        Q_current = Q_new;

        if (ctx->verbose) {
            char buf[256];
            snprintf(buf, sizeof(buf), "civreghdfe: CUE iter %d: Q=%.8e grad=%.2e step=%.4f %s\n",
                     iter, Q_current, grad_norm, step,
                     modified ? "modified newton" : newton_ok ? "newton" : "gradient");
            SF_display(buf);
        }
    }

    if (iter >= max_iter && ctx->verbose) {
        SF_display("civreghdfe: CUE reached max iterations\n");
    }

    /* Compute final residuals */
    cue_compute_residuals(y, ctx->X_all, beta, N, K_total, current_resid);
    if (resid) {
        memcpy(resid, current_resid, N * sizeof(ST_double));
    }

    /* Compute XZWZX_inv for VCE using the weighting matrix at the estimate
       (a failed allocation is an error: the caller must not post beta with
       another VCE) */
    if (XZWZX_inv_out) {
        /* Build S from final residuals and invert */
        ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        ST_double *ZOmegaZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));

        out_rc = (ZOmegaZ && ZOmegaZ_inv)
                     ? civ_omega_build(ctx->om, current_resid, 1, Z, K_iv, ZOmegaZ) : 920;
        if (!out_rc) {
            memcpy(ZOmegaZ_inv, ZOmegaZ, K_iv * K_iv * sizeof(ST_double));
            if (ctools_cholesky(ZOmegaZ_inv, K_iv) == 0) {
                out_rc = ctools_stata_rc(ctools_invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv));

                /* Compute (X'Z * W * Z'X)^{-1} where W = ZOmegaZ_inv */
                ST_double *XZW = (ST_double *)calloc((size_t)K_total * K_iv, sizeof(ST_double));
                ST_double *XZWZX = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
                if (!out_rc && (!XZW || !XZWZX)) out_rc = 920;

                if (!out_rc) {
                    for (i = 0; i < K_total; i++) {
                        for (j = 0; j < K_iv; j++) {
                            ST_double sum = 0.0;
                            for (k = 0; k < K_iv; k++) {
                                sum += ctx->ZtX[i * K_iv + k] * ZOmegaZ_inv[j * K_iv + k];
                            }
                            XZW[j * K_total + i] = sum;
                        }
                    }
                    for (i = 0; i < K_total; i++) {
                        for (j = 0; j < K_total; j++) {
                            ST_double sum = 0.0;
                            for (k = 0; k < K_iv; k++) {
                                sum += XZW[k * K_total + i] * ctx->ZtX[j * K_iv + k];
                            }
                            XZWZX[j * K_total + i] = sum;
                        }
                    }

                    if (ctools_cholesky(XZWZX, K_total) == 0) {
                        out_rc = ctools_stata_rc(ctools_invert_from_cholesky(XZWZX, K_total, XZWZX_inv_out));
                    }
                }
                if (XZW) free(XZW);
                if (XZWZX) free(XZWZX);
            }
        }
        if (ZOmegaZ) free(ZOmegaZ);
        if (ZOmegaZ_inv) free(ZOmegaZ_inv);
    }

cue_done:
    for (int t = 0; t < cue_nthreads; t++) {
        free(thread_beta[t]);
        free(thread_resid[t]);
        free(thread_grad[t]);
    }
    free(thread_beta);
    free(thread_resid);
    free(thread_grad);
    free(current_resid);
    free(gradient);
    free(hessian);
    free(direction);
    free(unit);

    return out_rc;
}

/*
    Wald statistic of the excluded instruments in one first-stage regression
    under the model's non-iid VCE (ivreg2's s_ffirst "robustified" Wald):
    Rb' (R Vr R')^-1 Rb with Vr = (Z'Z)^-1 S (Z'Z)^-1, S the raw moment
    covariance (civ_omega_build) of the first-stage residuals. Z holds the
    K_exog included instruments first; fs_coef are the first-stage
    coefficients and fitted holds Z*fs_coef on entry (overwritten with the
    residuals). Returns 0 or 920; *wald is NaN when R Vr R' is unusable.
*/
static ST_retcode civ_first_stage_robust_wald(
    const civ_omega *om, const ST_double *Z, const ST_double *X_e,
    ST_double *fitted, const ST_double *fs_coef, const ST_double *ZtZ_inv,
    ST_int N, ST_int K_iv, ST_int K_exog, ST_double *wald)
{
    ST_int L1 = K_iv - K_exog;
    *wald = NAN;
    if (L1 <= 0) return 0;
    for (ST_int i = 0; i < N; i++) fitted[i] = X_e[i] - fitted[i];

    ST_double *S = (ST_double *)calloc((size_t)K_iv * K_iv, sizeof(ST_double));
    ST_double *T = (ST_double *)calloc((size_t)K_iv * L1, sizeof(ST_double));
    ST_double *RVR = (ST_double *)calloc((size_t)L1 * L1, sizeof(ST_double));
    ST_double *sol = (ST_double *)calloc((size_t)L1, sizeof(ST_double));
    ST_retcode rc = (S && T && RVR && sol) ? 0 : 920;
    if (!rc && civ_omega_build(om, fitted, 1, Z, K_iv, S) != 0) rc = 920;
    if (!rc) {
        /* T = S (Z'Z)^-1 R', then R Vr R' = R (Z'Z)^-1 T */
        for (ST_int b = 0; b < L1; b++) {
            const ST_double *zinv_b = ZtZ_inv + (size_t)(K_exog + b) * K_iv;
            for (ST_int k = 0; k < K_iv; k++) {
                ST_double s = 0.0;
                for (ST_int j = 0; j < K_iv; j++) s += S[(size_t)j * K_iv + k] * zinv_b[j];
                T[(size_t)b * K_iv + k] = s;
            }
        }
        for (ST_int b = 0; b < L1; b++) {
            for (ST_int a = 0; a < L1; a++) {
                const ST_double *zinv_a = ZtZ_inv + (size_t)(K_exog + a) * K_iv;
                ST_double s = 0.0;
                for (ST_int k = 0; k < K_iv; k++) s += zinv_a[k] * T[(size_t)b * K_iv + k];
                RVR[(size_t)b * L1 + a] = s;
            }
        }
        if (civ_sym_solve(RVR, L1, fs_coef + K_exog, 1, sol) != 0) {
            rc = 920;
        } else {
            ST_double q = 0.0;
            for (ST_int a = 0; a < L1; a++) q += fs_coef[K_exog + a] * sol[a];
            if (isfinite(q)) *wald = q;
        }
    }
    free(S); free(T); free(RVR); free(sol);
    return rc;
}

/*
    Compute k-class IV estimation (includes 2SLS, LIML, Fuller, etc.)

    The k-class estimator is:
    beta = ((1-k)*X'X + k*X'P_Z X)^-1 ((1-k)*X'y + k*X'P_Z y)
    where P_Z = Z(Z'Z)^-1 Z' is the projection onto the instrument space

    When k=1, this reduces to 2SLS: beta = (X'P_Z X)^-1 X'P_Z y
    When k=lambda (min eigenvalue), this is LIML
    When k=lambda - alpha/(N-K_iv), this is Fuller LIML

    Parameters:
    - kclass: the k value (1.0 for 2SLS, lambda for LIML, etc.)
    - est_method: 0=2SLS, 1=LIML, 2=Fuller, 3=kclass, 4=GMM2S, 5=CUE
*/
/*
 * The coefficient phase of ivest_compute_2sls, without the Stata SPI (plan
 * 9.2): the k of the k-class family (LIML's lambda), the k-class solve, the
 * wmatrix() estimator, the moment covariance S of the chosen VCE, and the
 * GMM2S and CUE re-estimation. The full fit and the bootstrap replicates
 * both call it. Messages only with io->report; the arrays it allocates are
 * handed back in io, also on failure, for the caller to free.
 */
ST_retcode ivest_fit_core(civ_fit_io *io)
{
    const int report = io->report;
    const ST_double *y = io->y, *X_exog = io->X_exog, *X_endog = io->X_endog, *Z = io->Z;
    const ST_double *weights = io->weights;
    ST_int weight_type = io->weight_type, N = io->N;
    ST_double N_eff = io->N_eff;
    ST_int K_exog = io->K_exog, K_endog = io->K_endog, K_iv = io->K_iv, K_total = io->K_total;
    ST_double *beta = io->beta, *lambda_out = io->lambda_out;
    ST_int verbose = io->verbose, est_method = io->est_method;
    ST_double kclass_user = io->kclass_user, fuller_alpha = io->fuller_alpha;
    ST_int s_kind = io->s_kind, s_kernel = io->s_kernel, s_noniid = io->s_noniid;
    ST_double s_bw = io->s_bw;
    ST_int center = io->center;
    const ST_int *cluster_ids = io->cluster_ids, *cluster2_ids = io->cluster2_ids;
    ST_int num_clusters = io->num_clusters, num_clusters2 = io->num_clusters2;
    const civreghdfe_ts *ts = io->ts;
    const civreghdfe_extra *extra = io->extra;
    const ST_double *S0 = io->S0, *W0 = io->W0, *b0 = io->b0;
    ST_int dofminus = io->dofminus;
    ST_retcode est_rc = STATA_OK;
    ST_double *ZtZ = NULL;
    ST_double *ZtZ_inv = NULL;
    ST_double *ZtX = NULL;
    ST_double *Zty = NULL;
    ST_double *XtPzX = NULL;
    ST_double *XtPzy = NULL;
    ST_double *temp1 = NULL;
    ST_double *X_all = NULL;
    ST_double *resid = NULL;
    ST_double *X_hat = NULL;
    ST_double *ZtZ_inv_Zty = NULL;
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_double *XkX = NULL;
    ST_double *Xky = NULL;
    ST_double *XkX_copy = NULL;
    ST_double *beta_temp = NULL;
    ST_double *S_main = NULL;
    ST_double *gmm_hessian_inv = NULL;
    ST_double *W_bread = NULL;
    ST_double *W_zx = NULL;
    civ_omega om;
    memset(&om, 0, sizeof(om));
    ST_int i, j, k;
    ST_int effic1s = 0, gmmw_ineff = 0, s_fullrank = 1;
    (void)j;

    /* Determine k value based on estimation method */
    ST_double kclass = 1.0;  /* Default: 2SLS */
    ST_double lambda = 1.0;

    if (est_method == 1 || est_method == 2 || (est_method == 5 && !s_noniid)) {
        /* LIML, Fuller, or homoskedastic CUE: compute lambda
           Homoskedastic CUE = LIML (both minimize e'Pze / e'e) */
        lambda = compute_liml_lambda(y, X_endog, X_exog, Z, weights, weight_type, N, K_exog, K_endog, K_iv);
        if (isnan(lambda)) {
            /* never 2SLS under the LIML label */
            if (report) SF_error("civreghdfe: Memory allocation failed for the LIML eigenvalue\n");
            return 920;
        }
        if (lambda < 0.0 && K_iv != K_total) {
            /* a singular cross-product matrix: no LIML root, never 2SLS
               under the LIML label */
            if (report) SF_error("civreghdfe: the LIML eigenvalue cannot be computed (singular cross-product matrix)\n");
            return 498;
        }

        /* For exactly identified models, lambda should be 1 */
        if (K_iv == K_total) {
            lambda = 1.0;
        }

        if (est_method == 1 || (est_method == 5 && !s_noniid)) {
            /* LIML or homoskedastic CUE */
            kclass = lambda;
        } else {
            /* Fuller: k = lambda - alpha/(N_eff - K_iv)
               Use N_eff (sum of fweights for fweights, N otherwise) to match ivreg2 */
            if (fuller_alpha > (N_eff - K_iv)) {
                if (report) SF_error("civreghdfe: Invalid Fuller parameter\n");
                return 198;
            }
            kclass = lambda - fuller_alpha / (ST_double)(N_eff - K_iv);
        }

        if (lambda_out) *lambda_out = lambda;
    }
    else if (est_method == 3) {
        /* User-specified k-class */
        kclass = kclass_user;
    }
    /* est_method 0: k=1 (2SLS)
       est_method 4: GMM2S - will be computed after initial 2SLS
       est_method 5 with a non-iid S: CUE by minimizing the GMM objective */

    /* the k of the k-class estimator, ivreg2's e(kclass) */
    io->kclass = kclass;
    io->kclass_set = 1;

    const char *method_name = "2SLS";
    if (est_method == 1) method_name = "LIML";
    else if (est_method == 2) method_name = "Fuller LIML";
    else if (est_method == 3) method_name = "k-class";
    else if (est_method == 4) method_name = "GMM2S";
    else if (est_method == 5) method_name = "CUE";

    if (verbose) {
        char buf[256];
        snprintf(buf, sizeof(buf), "civreghdfe: Computing %s estimation\n", method_name);
        if (report) SF_display(buf);
        snprintf(buf, sizeof(buf), "  N=%d, K_exog=%d, K_endog=%d, K_iv=%d\n",
                 (int)N, (int)K_exog, (int)K_endog, (int)K_iv);
        if (report) SF_display(buf);
        if (est_method == 1 || est_method == 2) {
            snprintf(buf, sizeof(buf), "  lambda=%.6f, k=%.6f\n", lambda, kclass);
            if (report) SF_display(buf);
        } else if (est_method == 3) {
            snprintf(buf, sizeof(buf), "  k=%.6f\n", kclass);
            if (report) SF_display(buf);
        }
    }

    /* Check identification */
    if (K_iv < K_total) {
        if (report) SF_error("civreghdfe: Model is underidentified\n");
        return 198;
    }

    /* Allocate work arrays (with overflow checks) */
    ZtZ = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
    ZtZ_inv = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
    ZtX = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_total, sizeof(ST_double));
    Zty = (ST_double *)ctools_safe_calloc2((size_t)K_iv, sizeof(ST_double));
    XtPzX = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
    XtPzy = (ST_double *)ctools_safe_calloc2((size_t)K_total, sizeof(ST_double));
    temp1 = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_total, sizeof(ST_double));
    X_all = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_total, sizeof(ST_double));
    resid = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));

    if (!ZtZ || !ZtZ_inv || !ZtX || !Zty || !XtPzX || !XtPzy ||
        !temp1 || !X_all || !resid) {
        if (report) SF_error("civreghdfe: Memory allocation failed\n");
        est_rc = 920;
        goto est_cleanup;
    }

    /* Combine X_exog and X_endog into X_all */
    /* Layout: [X_exog (N x K_exog) | X_endog (N x K_endog)] */
    if (K_exog > 0 && X_exog) {
        size_t exog_size;
        if (ctools_safe_alloc_size((size_t)N, (size_t)K_exog, sizeof(ST_double), &exog_size) != 0) {
            if (report) SF_error("civreghdfe: Size overflow in memcpy\n");
            est_rc = 920;
            goto est_cleanup;
        }
        memcpy(X_all, X_exog, exog_size);
    }
    if (K_endog > 0 && X_endog) {
        size_t endog_size;
        if (ctools_safe_alloc_size((size_t)N, (size_t)K_endog, sizeof(ST_double), &endog_size) != 0) {
            if (report) SF_error("civreghdfe: Size overflow in memcpy\n");
            est_rc = 920;
            goto est_cleanup;
        }
        memcpy(X_all + (size_t)N * (size_t)K_exog, X_endog, endog_size);
    }

    /* Step 1: Compute Z'Z (weighted if needed) */
    if (weights && weight_type != 0) {
        matmul_atdb(Z, Z, weights, N, K_iv, K_iv, ZtZ);
    } else {
        ctools_matmul_ata(Z, N, K_iv, ZtZ);
    }

    /* Step 2: Invert Z'Z */
    {
        size_t ztz_size;
        ctools_safe_alloc_size((size_t)K_iv, (size_t)K_iv, sizeof(ST_double), &ztz_size);
        memcpy(ZtZ_inv, ZtZ, ztz_size);
    }
    if (ctools_cholesky(ZtZ_inv, K_iv) != 0) {
        if (report) SF_error("civreghdfe: Z'Z is singular (instruments may be collinear)\n");
        est_rc = 198;
        goto est_cleanup;
    }
    stata_retcode inverse_status_1 = ctools_invert_from_cholesky(ZtZ_inv, K_iv, ZtZ_inv);
    if (inverse_status_1 != STATA_OK) {
        if (report) SF_error("civreghdfe: Failed to invert Z'Z\n");
        est_rc = ctools_stata_rc(inverse_status_1);
        goto est_cleanup;
    }

    /* Step 3: Compute Z'X and Z'y */
    if (weights && weight_type != 0) {
        matmul_atdb(Z, X_all, weights, N, K_iv, K_total, ZtX);
        /* Z'y */
        for (i = 0; i < K_iv; i++) {
            ST_double sum = 0.0;
            const ST_double *z_col = Z + (size_t)i * N;
            for (k = 0; k < N; k++) {
                sum += z_col[k] * weights[k] * y[k];
            }
            Zty[i] = sum;
        }
    } else {
        matmul_atb(Z, X_all, N, K_iv, K_total, ZtX);
        /* Z'y */
        for (i = 0; i < K_iv; i++) {
            ST_double sum = 0.0;
            const ST_double *z_col = Z + (size_t)i * N;
            for (k = 0; k < N; k++) {
                sum += z_col[k] * y[k];
            }
            Zty[i] = sum;
        }
    }

    /* Step 4: Compute (Z'Z)^-1 * Z'X  -> temp1 (K_iv x K_total) */
    matmul_ab(ZtZ_inv, ZtX, K_iv, K_iv, K_total, temp1);

    /* Step 5: Compute X'P_Z X and X'P_Z y via fitted values (X_hat approach).
     * X_hat = Z * (Z'Z)^{-1} * Z'X = Z * temp1 (N x K_total)
     * X'PzX = X_hat' * X_hat using Kahan-compensated N-length dot products.
     *
     * This matches ivreg2's quadcross(X_hat, X_hat) precision. The compact
     * triple-product approach (Z'X)' * (Z'Z)^{-1} * Z'X only does K_iv-length
     * dot products, which don't benefit from compensated summation. With
     * condition number ~10^8 for X'PzX, the ~15 digit X'PzX precision after
     * compact triple-product leaves only ~7 sigfigs after inversion. The
     * X_hat approach gives ~19 digit X'PzX precision (Kahan N-length dots),
     * leaving ~11 sigfigs after inversion. */
    X_hat = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_total, sizeof(ST_double));
    if (!X_hat) {
        if (report) SF_error("civreghdfe: Memory allocation failed for X_hat\n");
        est_rc = 920;
        goto est_cleanup;
    }

    /* X_hat[i,j] = sum_{k=0}^{K_iv-1} Z[i,k] * temp1[k,j], accumulated in
     * that order for every element (row blocks in parallel).
     * For weighted: X_hat = sqrt(W) * Z * temp1 so X_hat'X_hat = X'W*PzX */
    ctools_matmul_xb(Z, N, K_iv, temp1, K_total, X_hat);
    if (weights && weight_type != 0) {
        #pragma omp parallel for schedule(static) if((size_t)N * K_total > 200000)
        for (j = 0; j < K_total; j++) {
            ST_double *xh = X_hat + (size_t)j * N;
            for (ST_int row = 0; row < N; row++) xh[row] = sqrt(weights[row]) * xh[row];
        }
    }

    /* X'PzX = X_hat' * X_hat using Kahan-compensated dot products (symmetric);
     * every element is an independent sum, so the columns run in parallel */
    #pragma omp parallel for schedule(dynamic) if(K_total > 1 && (size_t)N * K_total > 200000)
    for (j = 0; j < K_total; j++) {
        for (ST_int a = 0; a <= j; a++) {
            ST_double val = kahan_dot(X_hat + (size_t)a * N, X_hat + (size_t)j * N, N);
            XtPzX[j * K_total + a] = val;
            XtPzX[a * K_total + j] = val;
        }
    }

    /* X'Pzy = X_hat' * y_hat where y_hat = Z * (Z'Z)^{-1} * Z'y */
    /* Step 6: Compute (Z'Z)^-1 * Z'y */
    ZtZ_inv_Zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
    if (!ZtZ_inv_Zty) {
        if (report) SF_error("civreghdfe: Memory allocation failed for ZtZ_inv_Zty\n");
        est_rc = 920;
        goto est_cleanup;
    }
    for (i = 0; i < K_iv; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += ZtZ_inv[k * K_iv + i] * Zty[k];
        }
        ZtZ_inv_Zty[i] = sum;
    }

    /* Compute y_hat = Z * ZtZ_inv_Zty, then X'Pzy = X_hat' * y_hat via Kahan */
    {
        ST_double *y_hat = (ST_double *)malloc(N * sizeof(ST_double));
        if (!y_hat) {
            est_rc = 920;
            goto est_cleanup;
        }
        ctools_matmul_xb(Z, N, K_iv, ZtZ_inv_Zty, 1, y_hat);
        if (weights && weight_type != 0) {
            for (i = 0; i < N; i++) y_hat[i] = sqrt(weights[i]) * y_hat[i];
        }
        /* X'Pzy[a] = X_hat[:,a] . y_hat */
        #pragma omp parallel for schedule(static) if(K_total > 1 && (size_t)N * K_total > 200000)
        for (i = 0; i < K_total; i++) {
            XtPzy[i] = kahan_dot(X_hat + (size_t)i * N, y_hat, N);
        }
        free(y_hat);
    }
    free(X_hat);
    X_hat = NULL;

    /* Step 6b: For k-class estimation (k != 1), compute X'X and X'y */
    XkX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    Xky = (ST_double *)calloc(K_total, sizeof(ST_double));

    if (!XkX || !Xky) {
        if (report) SF_error("civreghdfe: Memory allocation failed for XkX/Xky\n");
        est_rc = 920;
        goto est_cleanup;
    }

    if (kclass != 1.0) {
        XtX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        Xty = (ST_double *)calloc(K_total, sizeof(ST_double));
        if (!XtX || !Xty) {
            if (report) SF_error("civreghdfe: Memory allocation failed for X'X/X'y\n");
            est_rc = 920;
            goto est_cleanup;
        }

        /* Compute X'X */
        if (weights && weight_type != 0) {
            matmul_atdb(X_all, X_all, weights, N, K_total, K_total, XtX);
        } else {
            matmul_atb(X_all, X_all, N, K_total, K_total, XtX);
        }

        /* Compute X'y */
        if (weights && weight_type != 0) {
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = X_all + (size_t)i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * weights[k] * y[k];
                }
                Xty[i] = sum;
            }
        } else {
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = X_all + (size_t)i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * y[k];
                }
                Xty[i] = sum;
            }
        }

        /* Compute k-class weighted matrices:
           XkX = (1-k)*X'X + k*X'P_Z X
           Xky = (1-k)*X'y + k*X'P_Z y */
        ST_double one_minus_k = 1.0 - kclass;
        for (i = 0; i < K_total * K_total; i++) {
            XkX[i] = one_minus_k * XtX[i] + kclass * XtPzX[i];
        }
        for (i = 0; i < K_total; i++) {
            Xky[i] = one_minus_k * Xty[i] + kclass * XtPzy[i];
        }
    } else {
        /* k = 1 (2SLS): just use X'P_Z X and X'P_Z y */
        memcpy(XkX, XtPzX, K_total * K_total * sizeof(ST_double));
        memcpy(Xky, XtPzy, K_total * sizeof(ST_double));
    }

    /* Step 7: Solve XkX * beta = Xky */
    XkX_copy = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
    if (!XkX_copy) {
        if (report) SF_error("civreghdfe: Memory allocation failed for XkX_copy\n");
        est_rc = 920;
        goto est_cleanup;
    }
    memcpy(XkX_copy, XkX, K_total * K_total * sizeof(ST_double));

    /* Use Cholesky solve */
    if (ctools_cholesky(XkX_copy, K_total) != 0) {
        if (report) SF_error("civreghdfe: XkX matrix is singular\n");
        est_rc = 198;
        goto est_cleanup;
    }

    /* Forward substitution */
    beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));
    if (!beta_temp) {
        if (report) SF_error("civreghdfe: Memory allocation failed for beta_temp\n");
        est_rc = 920;
        goto est_cleanup;
    }
    for (i = 0; i < K_total; i++) {
        ST_double sum = Xky[i];
        for (j = 0; j < i; j++) {
            sum -= XkX_copy[i * K_total + j] * beta_temp[j];
        }
        beta_temp[i] = sum / XkX_copy[i * K_total + i];
    }

    /* Backward substitution */
    for (i = K_total - 1; i >= 0; i--) {
        ST_double sum = beta_temp[i];
        for (j = i + 1; j < K_total; j++) {
            sum -= XkX_copy[j * K_total + i] * beta[j];
        }
        beta[i] = sum / XkX_copy[i * K_total + i];
    }

    if (verbose) {
        char buf[256];
        for (i = 0; i < K_total; i++) {
            snprintf(buf, sizeof(buf), "  beta[%d] = %g\n", (int)i, beta[i]);
            if (report) SF_display(buf);
        }
    }

    /* b0(): ivreg2 evaluates the CUE objective at the given coefficients
       (s_gmmcue with b0): they replace the estimate, and S, J and the
       efficient GMM VCE (X'Z S^-1 Z'X)^-1 are taken at their residuals */
    if (b0) memcpy(beta, b0, (size_t)K_total * sizeof(ST_double));

    /* ivreg2's first-step GMM: efficient (IV, or efficient GMM with a
       given S) unless the S is not iid or gmm2s is requested, and with a
       given S but no W. Otherwise a given W gives the inefficient
       first-step estimate (Z'X)'W(Z'X) b = (Z'X)'W Z'y, whose residuals
       give S (ivreg2's s_gmm1s); without W it is 2SLS. */
    effic1s = (!(est_method == 4) && !s_noniid) || (S0 && !W0);
    gmmw_ineff = W0 && !effic1s;
    if (gmmw_ineff) {
        /* W_zx = W Z'X (K_iv x K_total), W_bread = (X'Z W Z'X)^-1 */
        W_zx = (ST_double *)calloc((size_t)K_iv * K_total, sizeof(ST_double));
        W_bread = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
        ST_double *H = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
        ST_double *g = (ST_double *)calloc((size_t)K_total, sizeof(ST_double));
        ST_retcode wrc = (W_zx && W_bread && H && g) ? 0 : 920;
        if (!wrc) {
            ctools_matmul_ab(W0, ZtX, K_iv, K_iv, K_total, W_zx);
            for (ST_int c = 0; c < K_total; c++)
                for (ST_int r = 0; r < K_total; r++) {
                    ST_double sum = 0.0;
                    for (k = 0; k < K_iv; k++) sum += ZtX[(size_t)r * K_iv + k] * W_zx[(size_t)c * K_iv + k];
                    H[(size_t)c * K_total + r] = sum;
                }
            for (ST_int c = 0; c < K_total; c++)
                for (ST_int r = c + 1; r < K_total; r++) {
                    ST_double m = 0.5 * (H[(size_t)c * K_total + r] + H[(size_t)r * K_total + c]);
                    H[(size_t)c * K_total + r] = H[(size_t)r * K_total + c] = m;
                }
            for (ST_int r = 0; r < K_total; r++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_iv; k++) sum += W_zx[(size_t)r * K_iv + k] * Zty[k];
                g[r] = sum;
            }
            if (civ_sym_inverse(H, K_total, W_bread) != 0) wrc = 920;
        }
        if (!wrc) {
            for (ST_int r = 0; r < K_total; r++) {
                ST_double sum = 0.0;
                for (ST_int c = 0; c < K_total; c++) sum += W_bread[(size_t)c * K_total + r] * g[c];
                beta[r] = sum;
            }
        }
        free(H); free(g);
        if (wrc) {
            if (report) SF_error("civreghdfe: Memory allocation failed for the wmatrix() estimator\n");
            est_rc = wrc;
            goto est_cleanup;
        }
    }

    /* Step 8: Compute residuals using original X (not projected) */
    /* resid = y - X * beta */
    /* OPTIMIZED: Parallelized with OpenMP */
    #pragma omp parallel for schedule(static) if(N > 10000)
    for (i = 0; i < N; i++) {
        ST_double pred = 0.0;
        #pragma omp simd reduction(+:pred)
        for (ST_int kk = 0; kk < K_total; kk++) {
            pred += X_all[(size_t)kk * N + i] * beta[kk];
        }
        resid[i] = y[i] - pred;
    }

    /* Moment covariance S of the chosen VCE (ivreg2's m_omega). S_main is
       the S behind the reported J and C statistics: from the 2SLS residuals
       for GMM2S (its first step), at the estimate for CUE, from the model's
       residuals otherwise, and sigma^2 Z'WZ when iid (set after the RSS). */
    /* bw(auto): Newey and West's (1994) bandwidth from the first-step
       residuals, as ivreg2's abw selects it */
    if (s_kernel > 0 && s_bw < 0 && extra && extra->abw_z && ts) {
        s_bw = civ_auto_bandwidth(resid, extra->abw_z, N, ts, extra->abw_tobs, s_kernel);
        if (!(s_bw > 0)) {
            if (report) SF_error("civreghdfe: automatic bandwidth selection failed\n");
            est_rc = 498;
            goto est_cleanup;
        }
        io->bw_auto = s_bw;
        io->bw_auto_set = 1;
    }

    S_main = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
    ST_retcode om_rc = civ_omega_init(&om, s_kind, s_kernel, s_bw, N, (ST_double)N_eff,
                                      weights, weight_type, cluster_ids, num_clusters,
                                      cluster2_ids, num_clusters2, ts);
    om.center = center;
    om.psd = extra ? extra->psd : 0;   /* psd0 / psda */
    om.mult = io->mult;
    om.mult_groups = io->mult_groups;
    if (!om_rc && !S_main) om_rc = 920;
    if (!om_rc && S0) {
        /* the given S, in the raw scale (N times ivreg2's) */
        for (i = 0; i < K_iv * K_iv; i++) S_main[i] = S0[i] * (ST_double)N_eff;
    } else if (!om_rc && s_noniid) {
        om_rc = civ_omega_build(&om, resid, 1, Z, K_iv, S_main);
    } else if (!om_rc && W0 && est_method == 4) {
        /* gmm2s after a given W with iid errors: sigma^2 Z'WZ from the
           first-step residuals */
        ST_double ssr = 0.0;
        for (i = 0; i < N; i++) ssr += (weights && weight_type != 0 ? weights[i] : 1.0) * resid[i] * resid[i];
        for (i = 0; i < K_iv * K_iv; i++) S_main[i] = ssr / (ST_double)(N_eff - dofminus) * ZtZ[i];
    }
    if (om_rc) if (report) SF_error("civreghdfe: Cannot compute the moment covariance matrix\n");

    /* As ivreg2, a non-iid S must have full rank (Stata's invsym rank, too
       few clusters being the usual cause): efficient GMM (GMM2S, CUE) stops
       with r(506); the other estimators do not report the J, C and
       endogeneity statistics */
    s_fullrank = 1;
    if (!om_rc && (s_noniid || S0)) {
        ST_int s_rank = civ_sym_rank(S_main, K_iv, 1.0 / (ST_double)N_eff);
        if (s_rank < 0) {
            if (report) SF_error("civreghdfe: Memory allocation failed for the moment covariance rank\n");
            om_rc = 920;
        }
        s_fullrank = (s_rank == K_iv);
        if (!om_rc && !s_fullrank && !(S0 && effic1s) && (est_method == 4 || est_method == 5 || b0)) {
            if (report) SF_error("civreghdfe: estimated covariance matrix of moment conditions not of full rank,\n");
            if (report) SF_error("            and optimal GMM weighting matrix not unique\n");
            om_rc = 506;
        }
    }
    if (om_rc) {
        est_rc = om_rc;
        goto est_cleanup;
    }

    /* Allocate GMM/CUE Hessian inverse for VCE computation */

    /* Step 8b: For GMM2S, re-estimate with optimal weighting matrix (also
       the efficient GMM estimator of a given S, and two-step GMM after a
       given W) */
    if ((est_method == 4 && (s_noniid || W0 || S0)) || (S0 && effic1s)) {
        /*
            Two-step efficient GMM (ivreg2's s_egmm)
            Step 1: 2SLS residuals give S_main (built above with the chosen
                    VCE: robust, cluster, two-way, HAC, AC, DK or Kiefer)
            Step 2: Re-estimate β = (X'ZWZ'X)^-1 X'ZWZ'y with W = S_main^-1
        */

        /* Allocate storage for GMM Hessian inverse */
        gmm_hessian_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (!gmm_hessian_inv) {
            if (report) SF_error("civreghdfe: Memory allocation failed for the GMM Hessian\n");
            est_rc = 920;
            goto est_cleanup;
        }

        /* Create estimation context for GMM2S computation */
        IVEstContext gmm_ctx;
        gmm_ctx.N = N;
        gmm_ctx.K_exog = K_exog;
        gmm_ctx.K_endog = K_endog;
        gmm_ctx.K_iv = K_iv;
        gmm_ctx.K_total = K_total;
        gmm_ctx.weights = weights;
        gmm_ctx.weight_type = weight_type;
        gmm_ctx.Z = Z;
        gmm_ctx.y = y;
        gmm_ctx.ZtZ = ZtZ;
        gmm_ctx.ZtZ_inv = ZtZ_inv;
        gmm_ctx.ZtX = ZtX;
        gmm_ctx.Zty = Zty;
        gmm_ctx.XtPzX = XtPzX;
        gmm_ctx.XtPzy = XtPzy;
        gmm_ctx.temp_kiv_ktotal = temp1;
        gmm_ctx.X_all = X_all;
        gmm_ctx.om = &om;
        gmm_ctx.verbose = report ? verbose : 0;
        gmm_ctx.report = report;
        gmm_ctx.alloc_failed = 0;

        ST_retcode gmm_rc = ivest_compute_gmm2s(
            &gmm_ctx, y, S_main, beta, resid, gmm_hessian_inv
        );

        if (gmm_rc == 920) {
            /* beta may already be the GMM estimate: never post it with the
               2SLS VCE */
            if (report) SF_error("civreghdfe: Memory allocation failed in GMM2S\n");
            est_rc = gmm_rc;
            goto est_cleanup;
        } else if (gmm_rc != STATA_OK && !report) {
            /* A bootstrap replicate: its GMM2S estimate does not exist, so
               the replicate fails rather than counting 2SLS as GMM2S */
            est_rc = gmm_rc;
            goto est_cleanup;
        } else if (gmm_rc != STATA_OK) {
            if (verbose) SF_display("civreghdfe: GMM2S re-estimation failed, using 2SLS\n");
            if (gmm_hessian_inv) { free(gmm_hessian_inv); gmm_hessian_inv = NULL; }
        } else {
            if (verbose) {
                char buf[256];
                if (report) SF_display("civreghdfe: GMM2S estimates:\n");
                for (i = 0; i < K_total; i++) {
                    snprintf(buf, sizeof(buf), "  beta[%d] = %g\n", (int)i, beta[i]);
                    if (report) SF_display(buf);
                }
            }
        }

        /* Note: gmm_ctx uses pointers to existing arrays, no need to free context members */
    }

    /* Step 8c: CUE (Continuously Updated Estimator) */
    if (est_method == 5 && s_noniid) {
        /*
            Non-iid CUE: minimize g(β)'S(β)^-1 g(β) with S(β) from the chosen
            VCE (robust, cluster, two-way, HAC, AC, DK, Kiefer).
            For homoskedastic CUE, already handled above via LIML.
        */

        /* Allocate storage for CUE Hessian inverse */
        gmm_hessian_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (!gmm_hessian_inv) {
            if (report) SF_error("civreghdfe: Memory allocation failed for the CUE Hessian\n");
            est_rc = 920;
            goto est_cleanup;
        }

        /* Create estimation context for CUE computation */
        IVEstContext cue_ctx;
        cue_ctx.N = N;
        cue_ctx.K_exog = K_exog;
        cue_ctx.K_endog = K_endog;
        cue_ctx.K_iv = K_iv;
        cue_ctx.K_total = K_total;
        cue_ctx.weights = weights;
        cue_ctx.weight_type = weight_type;
        cue_ctx.Z = Z;
        cue_ctx.y = y;
        cue_ctx.ZtZ = ZtZ;
        cue_ctx.ZtZ_inv = ZtZ_inv;
        cue_ctx.ZtX = ZtX;
        cue_ctx.Zty = Zty;
        cue_ctx.XtPzX = XtPzX;
        cue_ctx.XtPzy = XtPzy;
        cue_ctx.temp_kiv_ktotal = temp1;
        cue_ctx.X_all = X_all;
        cue_ctx.om = &om;
        cue_ctx.verbose = report ? verbose : 0;
        cue_ctx.report = report;
        cue_ctx.alloc_failed = 0;

        const ST_int max_cue_iter = 100;
        const ST_double cue_tol = 1e-10;

        ST_retcode cue_rc = ivest_compute_cue(
            &cue_ctx, y, beta,
            beta, resid, max_cue_iter, cue_tol, gmm_hessian_inv
        );

        if (cue_rc != STATA_OK) {
            /* beta may already hold a Newton iterate: never post it with the
               2SLS VCE (ivest_compute_cue fails only when memory runs out) */
            if (report) SF_error("civreghdfe: Memory allocation failed in CUE estimation\n");
            est_rc = cue_rc;
            goto est_cleanup;
        } else {
            if (verbose) {
                char buf[256];
                if (report) SF_display("civreghdfe: CUE final estimates:\n");
                for (i = 0; i < K_total; i++) {
                    snprintf(buf, sizeof(buf), "  beta[%d] = %g\n", (int)i, beta[i]);
                    if (report) SF_display(buf);
                }
            }
            /* S at the CUE estimate: the J statistic is the CUE objective
               there (ivreg2's s_gmmcue), and the C statistic uses this S */
            om_rc = civ_omega_build(&om, resid, 1, Z, K_iv, S_main);
            if (!om_rc) {
                ST_int s_rank = civ_sym_rank(S_main, K_iv, 1.0 / (ST_double)N_eff);
                if (s_rank < 0) om_rc = 920;
                s_fullrank = (s_rank == K_iv);
            }
        }

        /* Note: cue_ctx uses pointers to existing arrays, no need to free context members */
    }
    if (om_rc) {
        if (report) SF_error("civreghdfe: Cannot compute the moment covariance matrix\n");
        est_rc = om_rc;
        goto est_cleanup;
    }

est_cleanup:
    io->ZtZ = ZtZ;
    io->ZtZ_inv = ZtZ_inv;
    io->ZtX = ZtX;
    io->Zty = Zty;
    io->XtPzX = XtPzX;
    io->XtPzy = XtPzy;
    io->temp1 = temp1;
    io->X_all = X_all;
    io->resid = resid;
    io->X_hat = X_hat;
    io->ZtZ_inv_Zty = ZtZ_inv_Zty;
    io->XtX = XtX;
    io->Xty = Xty;
    io->XkX = XkX;
    io->Xky = Xky;
    io->XkX_copy = XkX_copy;
    io->beta_temp = beta_temp;
    io->S_main = S_main;
    io->gmm_hessian_inv = gmm_hessian_inv;
    io->W_bread = W_bread;
    io->W_zx = W_zx;
    io->om = om;
    io->effic1s = effic1s;
    io->gmmw_ineff = gmmw_ineff;
    io->s_fullrank = s_fullrank;
    io->s_bw = s_bw;
    return est_rc;
}

/*
    ivreg2's e(S), e(W), e(condzz) and e(condxx) of the main fit, into
    __civreghdfe_S and __civreghdfe_W (top-left K_iv x K_iv; created by the
    caller over the instruments before collinearity) and two scalars.
    S = S_main / N. W: none for LIML and k-class (est_method 1-3); the given
    W0 for GMM with wmatrix(); for 2SLS (iv_w) (1/sigma^2) (Z'WZ/N)^-1 with
    sigma^2 = e'We/(N - dofminus), the first-step weighting matrix ivreg2
    posts also with a non-iid S; the given W0 with wmatrix() (also for
    two-step GMM, as ivreg2 posts it); otherwise (GMM2S, CUE, b0(),
    smatrix()) the efficient W = (S/N)^-1. cond() is the ratio of the largest to the
    smallest eigenvalue (missing when singular), as Mata's cond().
*/
static ST_double civ_cond_sym(const ST_double *A, ST_int n)
{
    if (n <= 0) return SV_missval;
    ST_double *ev = (ST_double *)malloc((size_t)n * sizeof(ST_double));
    ST_double *Q = (ST_double *)malloc((size_t)n * n * sizeof(ST_double));
    if (!ev || !Q || civ_sym_eigvec(A, n, ev, Q) != 0) { free(ev); free(Q); return SV_missval; }
    free(Q);
    ST_double lo = fabs(ev[0]), hi = fabs(ev[0]);
    for (ST_int i = 1; i < n; i++) {
        if (fabs(ev[i]) < lo) lo = fabs(ev[i]);
        if (fabs(ev[i]) > hi) hi = fabs(ev[i]);
    }
    free(ev);
    return (lo > 0.0 && isfinite(hi / lo)) ? hi / lo : SV_missval;
}

static ST_retcode civ_store_sw(const ST_double *S_main, const ST_double *ZtZ, const ST_double *ZtZ_inv,
                               const ST_double *X_all, const ST_double *weights, ST_int weight_type,
                               ST_int N, ST_double N_eff, ST_int K_iv, ST_int K_total,
                               ST_int est_method, const ST_double *W0, ST_int iv_w,
                               ST_double rss, ST_int dofminus)
{
    ST_retcode rc = 0;
    if (!S_main || K_iv <= 0) return 0;
    if (SF_row("__civreghdfe_S") < K_iv || SF_col("__civreghdfe_S") < K_iv) return 0;
    for (ST_int j = 0; j < K_iv && !rc; j++)
        for (ST_int i = 0; i < K_iv && !rc; i++)
            rc = ctools_mat_store("__civreghdfe_S", i + 1, j + 1, S_main[(size_t)j * K_iv + i] / N_eff);
    if (!rc && !(est_method >= 1 && est_method <= 3) &&
        SF_row("__civreghdfe_W") >= K_iv && SF_col("__civreghdfe_W") >= K_iv) {
        ST_double *W = (ST_double *)malloc((size_t)K_iv * K_iv * sizeof(ST_double));
        if (!W) return 920;
        if (W0) memcpy(W, W0, (size_t)K_iv * K_iv * sizeof(ST_double));
        else if (iv_w && ZtZ_inv) {
            ST_double f = (rss > 0.0) ? N_eff * (N_eff - (ST_double)dofminus) / rss : SV_missval;
            for (size_t k = 0; k < (size_t)K_iv * K_iv; k++) W[k] = f * ZtZ_inv[k];
        } else {
            if (civ_sym_inverse(S_main, K_iv, W) != 0) { free(W); return 920; }
            for (size_t k = 0; k < (size_t)K_iv * K_iv; k++) W[k] *= N_eff;
        }
        for (ST_int j = 0; j < K_iv && !rc; j++)
            for (ST_int i = 0; i < K_iv && !rc; i++)
                rc = ctools_mat_store("__civreghdfe_W", i + 1, j + 1, W[(size_t)j * K_iv + i]);
        free(W);
    }
    if (!rc && ZtZ) rc = ctools_scal_save("__civreghdfe_condzz", civ_cond_sym(ZtZ, K_iv));
    if (!rc && X_all && K_total > 0) {
        ST_double *XX = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
        if (!XX) return 920;
        if (weights && weight_type != 0) ctools_matmul_atdb(X_all, X_all, weights, N, K_total, K_total, XX);
        else ctools_matmul_ata(X_all, N, K_total, XX);
        rc = ctools_scal_save("__civreghdfe_condxx", civ_cond_sym(XX, K_total));
        free(XX);
    }
    return rc;
}

ST_retcode ivest_compute_2sls(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights (may exceed INT_MAX), N otherwise */
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_double *beta,
    ST_double *V,
    ST_double *first_stage_F,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int df_a,
    ST_int nested_adj,
    ST_int verbose,
    ST_int est_method,
    ST_double kclass_user,
    ST_double fuller_alpha,
    ST_double *lambda_out,
    ST_int kernel_type,
    ST_double bw,
    ST_int kiefer,
    const civreghdfe_ts *ts,
    ST_int sdofminus,
    ST_int center,
    const civreghdfe_test_cols *test_cols,
    ST_int n_absorbed_x,
    ST_int n_absorbed_iv,
    ST_int n_absorbed_iv_all,
    ST_int dofminus, ST_int small,
    const ST_double *b0, ST_int n_excl_dropped,
    const civreghdfe_extra *extra
)
{
    /* Heap state released at est_cleanup; every error path after the first
       allocation leaves through it */
    ST_retcode est_rc = STATA_OK;
    ST_double *ZtZ = NULL, *ZtZ_inv = NULL, *ZtX = NULL, *Zty = NULL;
    ST_double *XtPzX = NULL, *XtPzy = NULL, *temp1 = NULL, *X_all = NULL, *resid = NULL;
    ST_double *X_hat = NULL, *ZtZ_inv_Zty = NULL;
    ST_double *XtX = NULL, *Xty = NULL, *XkX = NULL, *Xky = NULL;
    ST_double *XkX_copy = NULL, *beta_temp = NULL, *S_main = NULL;
    ST_double *gmm_hessian_inv = NULL, *XkX_inv = NULL, *coviv_bread = NULL;
    civ_omega om;
    memset(&om, 0, sizeof(om));

    ST_int K_total = K_exog + K_endog;  /* Total regressors */
    /* Regressor and instrument counts for the degrees of freedom: they
       include the fixed-effect-absorbed columns that ivreghdfe keeps in its
       rank counts (see civreghdfe_impl.c), which are not estimated here */
    ST_int K_dof = K_total + n_absorbed_x;
    ST_int Kiv_dof = K_iv + n_absorbed_iv;
    ST_int i, j, k;

    /* smatrix()/wmatrix() (ivreg2's S0 and W0): they apply to IV and
       two-step GMM; LIML, k-class and CUE ignore them, as ivreg2 does */
    const ST_double *S0 = (extra && (est_method == 0 || est_method == 4) && !b0) ? extra->S0 : NULL;
    const ST_double *W0 = (extra && (est_method == 0 || est_method == 4) && !b0) ? extra->W0 : NULL;
    ST_double *W_bread = NULL, *W_zx = NULL;

    /* sdofminus is the FULL df_a (including all absorbed FE), used for test statistics DOF.
       This differs from the df_a parameter which is df_a_for_vce (excluding nested FE).
       For VCE computation, use df_a. For test statistics, use sdofminus. */
    ST_int df_a_full = (sdofminus > 0) ? sdofminus : df_a;

    /* Structure of the moment covariance S (ivreg2's m_omega), shared by the
       VCE, GMM2S/CUE, the J and C statistics: bw()/kernel() alone is AC,
       with robust HAC; dkraay() clusters on time with a kernel; kiefer is AC
       with the truncated kernel over every within-panel lag. */
    ST_int s_kind = CIV_OMEGA_IID, s_kernel = 0;
    ST_double s_bw = bw;
    if (kiefer) {
        /* heteroskedastic too with pweights (ivreg2 sets robust after its
           kiefer checks) */
        if (vce_type == CIVREGHDFE_VCE_ROBUST) s_kind = CIV_OMEGA_ROBUST;
        s_kernel = CIVREGHDFE_KERNEL_TRUNCATED;
        s_bw = ts ? ts->tmax_idx : 0;
    } else if (vce_type == CIVREGHDFE_VCE_DKRAAY) {
        s_kind = CIV_OMEGA_DKRAAY;
        s_kernel = (kernel_type > 0) ? kernel_type : CIVREGHDFE_KERNEL_BARTLETT;
    } else if (vce_type == CIVREGHDFE_VCE_CLUSTER2) {
        /* with a kernel: clusters on panel and time, Driscoll-Kraay over time */
        s_kind = CIV_OMEGA_TWOWAY;
        s_kernel = (bw != 0) ? kernel_type : 0;
    } else if (vce_type == CIVREGHDFE_VCE_CLUSTER) {
        s_kind = CIV_OMEGA_CLUSTER;
    } else {
        s_kind = (vce_type == CIVREGHDFE_VCE_ROBUST) ? CIV_OMEGA_ROBUST : CIV_OMEGA_IID;
        s_kernel = (bw != 0) ? kernel_type : 0;
    }
    /* Non-iid S: 2SLS/k-class are then inefficient and GMM2S/CUE reweight */
    ST_int s_noniid = (s_kind != CIV_OMEGA_IID) || (s_kernel > 0);
    ST_int s_cluster = (s_kind == CIV_OMEGA_CLUSTER || s_kind == CIV_OMEGA_TWOWAY ||
                        s_kind == CIV_OMEGA_DKRAAY);
    /* Clusters for the small-sample factor G/(G-1): time periods for
       Driscoll-Kraay, min(G1, G2) for two-way clustering */
    ST_int s_nclust = num_clusters;
    if (s_kind == CIV_OMEGA_TWOWAY && num_clusters2 < num_clusters) s_nclust = num_clusters2;
    if ((s_kernel > 0 || s_kind == CIV_OMEGA_DKRAAY) && !ts) {
        SF_error("civreghdfe: kernel-based VCE requires the tsset time variable\n");
        return 198;
    }

    if (verbose) {
        char buf[256];
        snprintf(buf, sizeof(buf), "civreghdfe: df_a=%d, sdofminus=%d, df_a_full=%d\n",
                 (int)df_a, (int)sdofminus, (int)df_a_full);
        SF_display(buf);
    }

    /* The coefficient phase (ivest_fit_core): k, the k-class solve, the
       wmatrix() estimator, S, GMM2S and CUE */
    civ_fit_io fio;
    memset(&fio, 0, sizeof(fio));
    fio.report = 1;
    fio.y = y; fio.X_exog = X_exog; fio.X_endog = X_endog; fio.Z = Z;
    fio.weights = weights; fio.weight_type = weight_type; fio.N = N; fio.N_eff = N_eff;
    fio.K_exog = K_exog; fio.K_endog = K_endog; fio.K_iv = K_iv; fio.K_total = K_total;
    fio.beta = beta; fio.lambda_out = lambda_out; fio.verbose = verbose;
    fio.est_method = est_method; fio.kclass_user = kclass_user; fio.fuller_alpha = fuller_alpha;
    fio.s_kind = s_kind; fio.s_kernel = s_kernel; fio.s_noniid = s_noniid; fio.s_bw = s_bw;
    fio.center = center;
    fio.cluster_ids = cluster_ids; fio.num_clusters = num_clusters;
    fio.cluster2_ids = cluster2_ids; fio.num_clusters2 = num_clusters2;
    fio.ts = ts; fio.extra = extra; fio.S0 = S0; fio.W0 = W0; fio.b0 = b0;
    fio.dofminus = dofminus;
    est_rc = ivest_fit_core(&fio);
    ZtZ = fio.ZtZ;
    ZtZ_inv = fio.ZtZ_inv;
    ZtX = fio.ZtX;
    Zty = fio.Zty;
    XtPzX = fio.XtPzX;
    XtPzy = fio.XtPzy;
    temp1 = fio.temp1;
    X_all = fio.X_all;
    resid = fio.resid;
    X_hat = fio.X_hat;
    ZtZ_inv_Zty = fio.ZtZ_inv_Zty;
    XtX = fio.XtX;
    Xty = fio.Xty;
    XkX = fio.XkX;
    Xky = fio.Xky;
    XkX_copy = fio.XkX_copy;
    beta_temp = fio.beta_temp;
    S_main = fio.S_main;
    gmm_hessian_inv = fio.gmm_hessian_inv;
    W_bread = fio.W_bread;
    W_zx = fio.W_zx;
    om = fio.om;
    if (fio.kclass_set) ctools_scal_save("__civreghdfe_k", fio.kclass);
    if (fio.bw_auto_set) ctools_scal_save("__civreghdfe_bw_auto", fio.bw_auto);
    if (est_rc) goto est_cleanup;
    s_bw = fio.s_bw;
    ST_int effic1s = fio.effic1s, gmmw_ineff = fio.gmmw_ineff, s_fullrank = fio.s_fullrank;

    /* Step 9: Compute VCE */
    /* For k-class estimators, the VCE is: sigma^2 * (XkX)^-1 */
    /* where sigma^2 = RSS / (N - K) using actual residuals */

    /* First, invert XkX (the k-class weighted matrix) */
    memcpy(XkX_copy, XkX, K_total * K_total * sizeof(ST_double));
    if (ctools_cholesky(XkX_copy, K_total) != 0) {
        SF_error("civreghdfe: Cannot compute VCE (XkX singular)\n");
        est_rc = 198;
        goto est_cleanup;
    }

    XkX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_int inverse_rc = XkX_inv ? ctools_invert_from_cholesky(XkX_copy, K_total, XkX_inv) : -1;
    if (inverse_rc != 0) {
        SF_error("civreghdfe: Memory allocation failed for VCE inverse\n");
        est_rc = 920;
        goto est_cleanup;
    }

    /* Compute RSS (needed for Sargan test; VCE helper computes its own internally) */
    ST_double rss = 0.0;
    if (weights && weight_type != 0) {
        for (i = 0; i < N; i++) rss += weights[i] * resid[i] * resid[i];
    } else {
        for (i = 0; i < N; i++) rss += resid[i] * resid[i];
    }

    /* iid S = sigma^2 Z'WZ with sigma^2 = e'We/(N-dofminus) from the model's
       residuals */
    if (!s_noniid && !S0 && !(W0 && est_method == 4)) {
        ST_double sigma2_s = rss / (ST_double)(N_eff - dofminus);
        for (i = 0; i < K_iv * K_iv; i++) S_main[i] = sigma2_s * ZtZ[i];
    }

    /* b0(): the efficient GMM VCE at S(b0), as s_gmmcue */
    if (b0) {
        ST_double *S_inv = (ST_double *)calloc((size_t)K_iv * K_iv, sizeof(ST_double));
        ST_double *A = (ST_double *)calloc((size_t)K_iv * K_total, sizeof(ST_double));
        ST_double *H = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
        gmm_hessian_inv = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
        ST_retcode b0_rc = (S_inv && A && H && gmm_hessian_inv) ? 0 : 920;
        if (!b0_rc && civ_sym_inverse(S_main, K_iv, S_inv) != 0) b0_rc = 920;
        if (!b0_rc) {
            /* A = S^-1 Z'X, H = X'Z S^-1 Z'X */
            ctools_matmul_ab(S_inv, ZtX, K_iv, K_iv, K_total, A);
            for (ST_int c = 0; c < K_total; c++)
                for (ST_int r = 0; r < K_total; r++) {
                    ST_double sum = 0.0;
                    for (ST_int k = 0; k < K_iv; k++) sum += ZtX[(size_t)r * K_iv + k] * A[(size_t)c * K_iv + k];
                    H[(size_t)c * K_total + r] = sum;
                }
            for (ST_int c = 0; c < K_total; c++)
                for (ST_int r = c + 1; r < K_total; r++) {
                    ST_double m = 0.5 * (H[(size_t)c * K_total + r] + H[(size_t)r * K_total + c]);
                    H[(size_t)c * K_total + r] = H[(size_t)r * K_total + c] = m;
                }
            if (civ_sym_inverse(H, K_total, gmm_hessian_inv) != 0) b0_rc = 920;
        }
        free(S_inv); free(A); free(H);
        if (b0_rc) {
            SF_error("civreghdfe: Memory allocation failed for the b0() VCE\n");
            est_rc = b0_rc;
            goto est_cleanup;
        }
    }

    /* ivreg2's small-sample factor for V: N/(N-K-df_a), or with clusters
       (incl. Driscoll-Kraay periods and two-way min(G1,G2))
       (N-1)/(N-K-df_a) * G/(G-1), where df_a = 1 when the fixed effects are
       nested in the clusters (ivreghdfe's absorb_ct) */
    ST_double dof_adj_s;
    if (s_cluster) {
        ST_int eff_df_a = (df_a == 0 && nested_adj == 1) ? 1 : df_a;
        ST_double denom = (ST_double)(N_eff - K_dof - eff_df_a);
        if (denom <= 0) denom = 1.0;
        ST_double G = (ST_double)s_nclust;
        dof_adj_s = ((ST_double)(N_eff - 1) / denom) * (G > 1.0 ? G / (G - 1.0) : 1.0);
    } else {
        /* (N-dofminus)/(N-K-dofminus-sdofminus); the iid S above carries
           the 1/(N-dofminus) */
        ST_double df_r_s = N_eff - K_dof - df_a - dofminus;
        if (df_r_s <= 0) df_r_s = 1;
        dof_adj_s = (ST_double)(N_eff - dofminus) / (ST_double)df_r_s;
    }

    if (!small) dof_adj_s = 1.0;

    /* coviv: LIML, Fuller and k-class take the 2SLS bread (X'PzX)^-1 in
       place of the k-class (X'(I - kMz)X)^-1, in the iid VCE and in every
       sandwich (ivreg2's s_liml); other estimators ignore it, as ivreg2 */
    ST_double *vce_bread_k = XkX_inv;
    {
        ST_double coviv_flag = 0;
        SF_scal_use("__civreghdfe_coviv", &coviv_flag);
        if (coviv_flag != 0 && est_method >= 1 && est_method <= 3) {
            coviv_bread = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
            if (!coviv_bread || civ_sym_inverse(XtPzX, K_total, coviv_bread) != 0) {
                SF_error("civreghdfe: Memory allocation failed for the coviv VCE\n");
                est_rc = 920;
                goto est_cleanup;
            }
            vce_bread_k = coviv_bread;
        }
    }

    if (gmm_hessian_inv == NULL && gmmw_ineff) {
        /* inefficient GMM with a given W (ivreg2's s_iegmm):
           V = B (X'Z) W S W (Z'X) B, B = (X'Z W Z'X)^-1 */
        ST_retcode vce_rc = ivvce_compute_sandwich(W_zx, S_main, W_bread, K_total, K_iv, dof_adj_s, V);
        if (vce_rc) {
            est_rc = vce_rc;
            goto est_cleanup;
        }
    } else if (gmm_hessian_inv == NULL && center && s_noniid && est_method >= 1 && est_method <= 3) {
        /* center: the k-class sandwich with the centered S (the 2SLS VCE
           does not change: its moments are orthogonal to the mean) */
        ST_retcode vce_rc = ivvce_compute_sandwich(temp1, S_main, vce_bread_k, K_total, K_iv, dof_adj_s, V);
        if (vce_rc) {
            est_rc = vce_rc;
            goto est_cleanup;
        }
    } else if (gmm_hessian_inv != NULL) {
        /*
           Efficient GMM2S / CUE VCE (ivreg2's s_egmm, s_gmmcue): with the
           optimal weights W = S^-1 (S of the first step for GMM2S, at the
           estimate for CUE) the VCE is (X'ZWZ'X)^-1, times the small-sample
           factor of the chosen VCE.
        */
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = gmm_hessian_inv[i] * dof_adj_s;
        }
    } else if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && cluster2_ids != NULL && num_clusters2 > 0 &&
               s_kernel <= 0) {
        /*
           Two-way clustered VCE using Cameron-Gelbach-Miller (2011) formula:
           V = V1 + V2 - V_intersection
           Apply nested_adj correction: when FE is nested in a cluster variable,
           df_a_for_vce is 0 but we still need to account for the partialled-out
           constant (effective_df_a = 1), matching ivreg2's sdofminus convention.
        */
        ST_int effective_df_a = df_a;
        if (df_a == 0 && nested_adj == 1) effective_df_a = 1;
        effective_df_a += n_absorbed_x;
        ST_retcode vce_rc = ivvce_compute_twoway(
            Z, resid, temp1, vce_bread_k,
            weights, weight_type,
            N, N_eff, K_total, K_iv,
            cluster_ids, num_clusters,
            cluster2_ids, num_clusters2,
            effective_df_a, small,
            V
        );
        if (vce_rc) {
            SF_error("civreghdfe: Memory allocation failed for the two-way cluster VCE\n");
            est_rc = vce_rc;
            goto est_cleanup;
        }
    } else if (s_kernel > 0 || s_kind == CIV_OMEGA_DKRAAY) {
        /*
           Kernel VCEs: HAC (robust + bw), AC (bw without robust),
           Driscoll-Kraay and Kiefer (AC with the truncated kernel over every
           within-panel lag), with lags by time within the tsset panels.
           Sandwich with the kernel S of the model's residuals (ivreg2's
           s_iegmm / s_liml): V = XkX^-1 A'SA XkX^-1, A = (Z'Z)^-1 Z'X.
           df_a equals the full absorbed DOF for Kiefer (civreghdfe_impl.c
           undoes the FE-nested-in-cluster subtraction).
        */
        if (verbose) {
            char buf[256];
            snprintf(buf, sizeof(buf), "civreghdfe: Kernel VCE (kernel=%d, lags=%d)\n",
                     (int)s_kernel, (int)om.max_lag);
            SF_display(buf);
        }
        ST_retcode vce_rc = ivvce_compute_sandwich(temp1, S_main, vce_bread_k, K_total, K_iv, dof_adj_s, V);
        if (vce_rc) {
            est_rc = vce_rc;
            goto est_cleanup;
        }
    } else {
        /*
           Standard VCE: Use refactored helper from civreghdfe_vce.c
           Handles unadjusted, robust (HC), and clustered VCE types.
        */
        /* For homoskedastic CUE, use 2SLS projection (X'PzX)^{-1} not (XkX)^{-1}
           because CUE VCE = σ²(X'Z(Z'Z)^{-1}Z'X)^{-1} even when coefficients = LIML */
        ST_double *vce_bread = vce_bread_k;
        ST_double *XtPzX_inv_cue = NULL;
        if (est_method == 5 && !s_noniid) {
            XtPzX_inv_cue = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
            /* a failed allocation must not leave the LIML bread in place */
            ST_int cue_inv_rc = XtPzX_inv_cue ? 0 : 920;
            if (XtPzX_inv_cue) {
                memcpy(XtPzX_inv_cue, XtPzX, K_total * K_total * sizeof(ST_double));
                if (ctools_cholesky(XtPzX_inv_cue, K_total) == 0) {
                    cue_inv_rc = ctools_stata_rc(ctools_invert_from_cholesky(XtPzX_inv_cue, K_total, XtPzX_inv_cue));
                    vce_bread = XtPzX_inv_cue;
                } else {
                    free(XtPzX_inv_cue);
                    XtPzX_inv_cue = NULL;
                }
            }
            if (cue_inv_rc) {
                SF_error("civreghdfe: Memory allocation failed for the CUE VCE\n");
                free(XtPzX_inv_cue);
                est_rc = cue_inv_rc;
                goto est_cleanup;
            }
        }

        ST_retcode vce_rc = ivvce_compute_full(
            Z, resid, temp1, vce_bread,
            weights, weight_type,
            N, N_eff, K_total, K_iv,
            vce_type, cluster_ids, num_clusters,
            /* the absorbed columns count like absorbed dof, after the
               nested-fixed-effect adjustment */
            ((vce_type == CIVREGHDFE_VCE_CLUSTER && df_a == 0 && nested_adj == 1) ? 1 : df_a) + n_absorbed_x,
            0, dofminus, small,
            V
        );

        if (XtPzX_inv_cue) free(XtPzX_inv_cue);
        if (vce_rc) {
            est_rc = vce_rc;
            goto est_cleanup;
        }
    }

    /* dofminus rescales the nonclustered moment covariance by N/(N-d).
     * This is independent of small and applies to HC, HAC and GMM alike. */
    if (s_noniid && !s_cluster && dofminus != 0 && !S0) {
        ST_double adjustment = N_eff / (N_eff - dofminus);
        for (i=0; i<K_total*K_total; i++) V[i] *= adjustment;
        for (i=0; i<K_iv*K_iv; i++) S_main[i] *= adjustment;
    }

    /* Step 10: Compute first-stage F statistics */
    /* For each endogenous variable, compute F-stat from first stage regression */
    /* The F tests whether excluded instruments are jointly significant,
       controlling for exogenous regressors. Uses partial R² approach:
       F = ((R²_full - R²_reduced) / L) / ((1 - R²_full) / df_resid) */
    if (first_stage_F) {
        for (ST_int e = 0; e < K_endog; e++) {
            const ST_double *X_e = X_endog + (size_t)e * N;

            /* Compute X_e'X_e */
            ST_double xx = 0.0;
            for (i = 0; i < N; i++) {
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                xx += w * X_e[i] * X_e[i];
            }
            if (xx <= 0) xx = 1.0;

            /* Compute R²_full: R² from projecting X_e onto all instruments Z.
               The first-stage coefficients (Z'Z)^-1 Z'X_e do not depend on
               the observation, so they are formed once; P_Z X_e then
               accumulates over the instruments in order. */
            ST_double xpx_full = 0.0;
            ST_double wald_robust = NAN;
            {
                ST_double *fs_coef = (ST_double *)calloc((size_t)(K_iv > 0 ? K_iv : 1), sizeof(ST_double));
                ST_double *pz_xe = (ST_double *)ctools_safe_malloc2((size_t)(N > 0 ? N : 1), sizeof(ST_double));
                if (!fs_coef || !pz_xe) {
                    free(fs_coef); free(pz_xe);
                    est_rc = 920;
                    goto est_cleanup;
                }
                for (k = 0; k < K_iv; k++) {
                    ST_double ziz_inv_zx = 0.0;
                    for (j = 0; j < K_iv; j++) {
                        ziz_inv_zx += ZtZ_inv[k * K_iv + j] * ZtX[(K_exog + e) * K_iv + j];
                    }
                    fs_coef[k] = ziz_inv_zx;
                }
                ctools_matmul_xb(Z, N, K_iv, fs_coef, 1, pz_xe);
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    xpx_full += w * X_e[i] * pz_xe[i];
                }
                /* ivreg2 robustifies the first-stage F under robust, cluster
                   and kernel VCEs (s_ffirst) */
                ST_retcode fs_rc = s_noniid
                    ? civ_first_stage_robust_wald(&om, Z, X_e, pz_xe, fs_coef, ZtZ_inv,
                                                  N, K_iv, K_exog, &wald_robust)
                    : 0;
                free(fs_coef); free(pz_xe);
                if (fs_rc) {
                    SF_error("civreghdfe: Memory allocation failed for the first-stage F\n");
                    est_rc = fs_rc;
                    goto est_cleanup;
                }
            }
            ST_double r2_full = xpx_full / xx;
            if (r2_full > 1.0) r2_full = 1.0;
            if (r2_full < 0.0) r2_full = 0.0;

            /* Compute R²_reduced: R² from projecting X_e onto only exogenous regressors */
            ST_double r2_reduced = 0.0;
            if (K_exog > 0) {
                /* Build Z_exog'Z_exog (K_exog x K_exog) from ZtZ */
                ST_double *ZeZe = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                ST_double *ZeZe_inv = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                ST_double *ZeXe = (ST_double *)calloc(K_exog, sizeof(ST_double));
                if (!ZeZe || !ZeZe_inv || !ZeXe) {
                    /* R2_reduced = 0 would overstate the first-stage F */
                    free(ZeZe); free(ZeZe_inv); free(ZeXe);
                    SF_error("civreghdfe: Memory allocation failed for the first-stage F\n");
                    est_rc = 920;
                    goto est_cleanup;
                }

                if (ZeZe && ZeZe_inv && ZeXe) {
                    /* Extract Z_exog'Z_exog from ZtZ */
                    for (ST_int l1 = 0; l1 < K_exog; l1++) {
                        for (ST_int l2 = 0; l2 < K_exog; l2++) {
                            ZeZe[l2 * K_exog + l1] = ZtZ[l2 * K_iv + l1];
                        }
                    }

                    /* Compute Z_exog'X_e */
                    for (ST_int l = 0; l < K_exog; l++) {
                        ST_double sum = 0.0;
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            sum += w * Z[(size_t)l * N + i] * X_e[i];
                        }
                        ZeXe[l] = sum;
                    }

                    /* Invert Z_exog'Z_exog */
                    memcpy(ZeZe_inv, ZeZe, K_exog * K_exog * sizeof(ST_double));
                    ST_int ze_ok = (ctools_cholesky(ZeZe_inv, K_exog) == 0);
                    stata_retcode ze_status = ze_ok ? ctools_invert_from_cholesky(ZeZe_inv, K_exog, ZeZe_inv) : STATA_OK;
                    if (ze_status != STATA_OK) {
                        free(ZeZe); free(ZeZe_inv); free(ZeXe);
                        SF_error("civreghdfe: Memory allocation failed for the first-stage F\n");
                        est_rc = ctools_stata_rc(ze_status);
                        goto est_cleanup;
                    }
                    if (ze_ok) {

                        /* Compute X_e'P_exog X_e = ZeXe' * ZeZe_inv * ZeXe */
                        ST_double xpx_reduced = 0.0;
                        for (ST_int l1 = 0; l1 < K_exog; l1++) {
                            ST_double temp = 0.0;
                            for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                temp += ZeZe_inv[l2 * K_exog + l1] * ZeXe[l2];
                            }
                            xpx_reduced += ZeXe[l1] * temp;
                        }
                        r2_reduced = xpx_reduced / xx;
                        if (r2_reduced > 1.0) r2_reduced = 1.0;
                        if (r2_reduced < 0.0) r2_reduced = 0.0;
                    }
                }

                free(ZeZe);
                free(ZeZe_inv);
                free(ZeXe);
            }

            /* Partial R² = R²_full - R²_reduced */
            ST_double partial_r2 = r2_full - r2_reduced;
            if (partial_r2 < 0.0) partial_r2 = 0.0;

            /* F = (partial_R² / L) / ((1 - R²_full) / df_resid) */
            ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
            if (L <= 0) L = 1;
            ST_double denom_df = N_eff - Kiv_dof - df_a - dofminus;
            if (denom_df <= 0) denom_df = 1;

            ST_double denom = (1.0 - r2_full) / (ST_double)denom_df;
            if (denom <= 0.0) denom = 1e-10;

            first_stage_F[e] = (partial_r2 / (ST_double)L) / denom;
            if (first_stage_F[e] < 0) first_stage_F[e] = 0;

            /* Non-iid VCE: ivreg2's robust Wald, scaled as its e(first)["F",]:
               Wald/(N-dofminus)*df_r/L, or with clusters
               Wald/(N-1)*(N-L_all-sdofminus)*(G-1)/G/L. The moment covariance
               carries the N/(N-dofminus) factor of m_omega without clusters. */
            if (s_noniid) {
                ST_double fs_df_a = (vce_type == CIVREGHDFE_VCE_CLUSTER && df_a == 0 && nested_adj == 1) ? 1 : df_a;
                ST_double fs_resid_df = N_eff - Kiv_dof - fs_df_a;
                ST_double wald = wald_robust;
                if (!s_cluster) {
                    if (dofminus != 0) wald *= (N_eff - dofminus) / N_eff;
                    first_stage_F[e] = wald / (N_eff - dofminus) * (fs_resid_df - dofminus) / (ST_double)L;
                } else {
                    first_stage_F[e] = wald / (N_eff - 1.0) * fs_resid_df
                                       * (ST_double)(s_nclust - 1) / (ST_double)s_nclust / (ST_double)L;
                }
                if (!isfinite(first_stage_F[e])) first_stage_F[e] = SV_missval;
            }
            /* denominator degrees of freedom of the reported F */
            char fdf_name[64];
            snprintf(fdf_name, sizeof(fdf_name), "__civreghdfe_F1df_r_%d", (int)(e + 1));
            ctools_scal_save(fdf_name, s_cluster ? (ST_double)(s_nclust - 1)
                                       : (ST_double)(N_eff - Kiv_dof - df_a - dofminus));

            /* Save the partial R² for ffirst display: ivreg2's pr2 from the
               classical Wald statistic, (R²_full - R²_reduced)/(1 - R²_reduced) */
            char r2_name[64];
            snprintf(r2_name, sizeof(r2_name), "__civreghdfe_partial_r2_%d", (int)(e + 1));
            ctools_scal_save(r2_name, (r2_reduced < 1.0) ? partial_r2 / (1.0 - r2_reduced) : 0.0);

            if (verbose) {
                char buf[256];
                snprintf(buf, sizeof(buf), "  First-stage F[%d] = %g (partial R² = %g, R²_full = %g)\n",
                         (int)e, first_stage_F[e], partial_r2, r2_full);
                SF_display(buf);
            }
        }
    }

    /* ivreg2's first-stage statistics (e(first), Anderson-Rubin, Stock-Wright
       and the system of first-stage and reduced-form equations) */
    if (extra && extra->ffirst && K_endog > 0 && K_iv > K_exog) {
        civ_ffirst_dof fd;
        fd.N = N_eff;
        fd.dofminus = dofminus;
        fd.sdofminus = (s_cluster && df_a == 0 && nested_adj == 1) ? 1 : df_a;
        fd.L_dof = Kiv_dof;
        fd.cluster = s_cluster;
        fd.nclust = s_nclust;
        fd.robust = s_noniid;
        est_rc = civ_ffirst(y, X_endog, K_endog, Z, K_iv, K_exog,
                            (weights && weight_type != 0) ? weights : NULL, N, &om, &fd);
        if (est_rc) {
            SF_error("civreghdfe: Memory allocation failed for the first-stage statistics\n");
            goto est_cleanup;
        }
    }

    /* ranktest's S, for the identification and redundancy statistics and the
       endogtest re-estimation: ivreghdfe passes bw() but never kernel() to
       ranktest and to the re-estimated model, so a kernel becomes Bartlett;
       kiefer passes neither bw() nor robust, so iid */
    civ_omega rk_om = civ_omega_with_kernel(&om, (!kiefer && s_kernel > 0) ? CIVREGHDFE_KERNEL_BARTLETT : 0);
    rk_om.center = 0;   /* ivreg2 does not pass center to ranktest */
    rk_om.psd = 0;      /* nor psd0 / psda */
    ST_int rk_cluster = s_cluster && !kiefer;

    /* Step 10b: Compute underidentification test and weak instrument stats */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double underid_stat = 0.0;
    ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
    ST_int underid_df = L;
    ST_double cd_f = 0.0;       /* Cragg-Donald Wald F (homoskedastic) */
    ST_double kp_f = 0.0;       /* Kleibergen-Paap rk Wald F (robust) */

    /* The tests return 920 when memory runs out; a statistic they cannot
       compute otherwise is posted as missing */
    est_rc = civreghdfe_compute_underid_test(
        X_endog, Z,
        weights, weight_type, N, N_eff, K_exog, K_endog, K_iv, df_a_full + n_absorbed_iv_all,
        dofminus, &rk_om, rk_cluster, s_nclust,
        &underid_stat, &underid_df, &cd_f, &kp_f
    );
    if (est_rc) {
        SF_error("civreghdfe: Memory allocation failed for the diagnostic tests\n");
        goto est_cleanup;
    }

    /* Excluded instruments dropped as collinear (with nocollin) or absorbed
       by the fixed effects still count in ivreg2's excluded-instrument
       count L2, which divides the Cragg-Donald and KP Wald F, enters their
       residual dof and the identification test dof (ivreg2 takes them
       from its variable lists) */
    if (n_excl_dropped > 0 && L > 0) {
        ST_double dfa_t = (ST_double)(df_a_full + n_absorbed_iv_all);
        ST_double dr = N_eff - K_iv - dfa_t - dofminus, dc = N_eff - K_iv - dfa_t;
        ST_double ls = (ST_double)L / (ST_double)(L + n_excl_dropped);
        if (dr > n_excl_dropped && dc > n_excl_dropped) {
            if (cd_f < SV_missval) cd_f *= (dr - n_excl_dropped) / dr * ls;
            if (kp_f < SV_missval) kp_f *= (rk_cluster ? (dc - n_excl_dropped) / dc : (dr - n_excl_dropped) / dr) * ls;
        }
        underid_df = L + n_excl_dropped - K_endog + 1;
        if (underid_df < 1) underid_df = 1;
    }

    /* Store underidentification test result */
    ctools_scal_save("__civreghdfe_underid", underid_stat);
    ctools_scal_save("__civreghdfe_underid_df", (ST_double)underid_df);

    /* Step 11: Compute Sargan/Hansen J overidentification test with S_main:
       the model's residuals when it is efficient for S (iid, GMM2S, CUE,
       whose J is its objective at the estimate), otherwise the efficient
       GMM residuals for S (ivreg2's s_iegmm / s_liml) */
    ST_double sargan_stat = 0.0;
    ST_int overid_df = K_iv - K_total;

    /* ivreg2's rankS: the rank of S as invsym() counts it (K_iv when iid) */
    {
        ST_int rank_s = (s_noniid || S0) ? civ_sym_rank(S_main, K_iv, 1.0 / (ST_double)N_eff) : K_iv;
        if (rank_s < 0) {
            SF_error("civreghdfe: Memory allocation failed for the moment covariance rank\n");
            est_rc = 920;
            goto est_cleanup;
        }
        ctools_scal_save("__civreghdfe_rankS", (ST_double)rank_s);
    }

    /* b0(): the CUE objective at b0, also when the model is exactly
       identified (ivreg2 reports it then too) */
    est_rc = civreghdfe_compute_hansen_j(
        y, X_all, b0 ? 0 : K_total, Z, K_iv,
        weights, weight_type, N, S_main,
        ((s_noniid || gmmw_ineff) && est_method != 4 && est_method != 5 && !b0 && !(S0 && effic1s)) ? NULL : resid,
        &sargan_stat, &overid_df
    );
    if (b0) overid_df = K_iv - K_total;
    if (est_rc) {
        SF_error("civreghdfe: Memory allocation failed for the diagnostic tests\n");
        goto est_cleanup;
    }

    /* Store diagnostic statistics as Stata scalars (J not reported when S
       is rank deficient) */
    ctools_scal_save("__civreghdfe_sargan", s_fullrank ? sargan_stat : SV_missval);
    overid_df += n_absorbed_iv - n_absorbed_x;   /* ivreghdfe's rankzz - rankxx */
    ctools_scal_save("__civreghdfe_sargan_df", (ST_double)overid_df);
    ctools_scal_save("__civreghdfe_cd_f", cd_f);
    ctools_scal_save("__civreghdfe_kp_f", kp_f);

    /* Step 12: Compute Durbin-Wu-Hausman endogeneity test */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double endog_chi2 = 0.0;
    ST_double endog_f = 0.0;
    ST_int endog_df = K_endog;

    est_rc = civreghdfe_compute_dwh_test(
        y, X_exog, X_endog, Z, temp1, weights, weight_type,
        N, N_eff, K_exog, K_endog, K_iv, df_a_full,
        &endog_chi2, &endog_f, &endog_df
    );
    if (est_rc) {
        SF_error("civreghdfe: Memory allocation failed for the diagnostic tests\n");
        goto est_cleanup;
    }

    ctools_scal_save("__civreghdfe_endog_chi2", endog_chi2);
    ctools_scal_save("__civreghdfe_endog_f", endog_f);
    ctools_scal_save("__civreghdfe_endog_df", (ST_double)endog_df);

    /* Step 13: Compute optional diagnostic tests (orthog, endogtest, redundant)
       on the columns selected in test_cols (already mapped to the final Z/X) */
    if (test_cols && test_cols->n_orthog + test_cols->n_orthog_x > 0) {
        ST_double cstat = 0.0;
        ST_int cstat_df = test_cols->n_orthog + test_cols->n_orthog_x;

        /* C = J - J_r with the S behind J (ivreg2's smatrix()) */
        est_rc = civreghdfe_compute_cstat(
            y, X_all, K_total, Z, K_iv, K_exog,
            weights, weight_type, N, S_main, sargan_stat,
            test_cols->orthog, test_cols->n_orthog,
            test_cols->orthog_x, test_cols->n_orthog_x,
            &cstat, &cstat_df
        );
        if (est_rc) {
            SF_error("civreghdfe: Memory allocation failed for the diagnostic tests\n");
            goto est_cleanup;
        }

        ctools_scal_save("__civreghdfe_cstat", s_fullrank ? cstat : SV_missval);
        ctools_scal_save("__civreghdfe_cstat_df", (ST_double)cstat_df);
    }

    if (test_cols && test_cols->n_endogtest > 0) {
        ST_double endogtest_stat = 0.0;
        ST_int endogtest_df_out = test_cols->n_endogtest;

        /* ivreg2 re-estimates with liml (LIML, Fuller) or gmm2s; 2SLS
           otherwise (cue and kclass are not passed) */
        ST_int endog_est = (est_method == 1 || est_method == 2) ? 1 : (est_method == 4 ? 4 : 0);
        est_rc = civreghdfe_compute_endogtest_subset(
            y, X_exog, X_endog, Z, N, weights, weight_type,
            K_exog, K_endog, K_iv, endog_est, &rk_om,
            test_cols->endogtest, test_cols->n_endogtest,
            &endogtest_stat, &endogtest_df_out
        );
        if (est_rc) {
            SF_error("civreghdfe: Memory allocation failed for the diagnostic tests\n");
            goto est_cleanup;
        }

        ctools_scal_save("__civreghdfe_endogtest_stat", s_fullrank ? endogtest_stat : SV_missval);
        ctools_scal_save("__civreghdfe_endogtest_df", (ST_double)endogtest_df_out);
    }

    if (test_cols && test_cols->n_redundant > 0) {
        ST_double redund_stat = 0.0;
        ST_int redund_df = K_endog * test_cols->n_redundant;

        est_rc = civreghdfe_compute_redundant(
            X_endog, Z, N, weights, weight_type,
            K_exog, K_endog, K_iv, &rk_om,
            test_cols->redundant, test_cols->n_redundant,
            &redund_stat, &redund_df
        );
        if (est_rc) {
            SF_error("civreghdfe: Memory allocation failed for the diagnostic tests\n");
            goto est_cleanup;
        }

        ctools_scal_save("__civreghdfe_redund_stat", redund_stat);
        ctools_scal_save("__civreghdfe_redund_df", (ST_double)redund_df);
    }

    /* ivreg2's e(S) and e(W) (its scale: divided by N) and the condition
       numbers e(condzz), e(condxx) of Z'WZ and X'WX, for the main fit */
    if (extra && extra->store_sw) {
        est_rc = civ_store_sw(S_main, ZtZ, ZtZ_inv, X_all, weights, weight_type, N, N_eff,
                              K_iv, K_total, est_method, (W0 && !b0) ? W0 : NULL,
                              (est_method == 0 && !S0 && !W0 && !b0), rss, dofminus);
        if (est_rc) goto est_cleanup;
    }

est_cleanup:
    /* Single exit: releases whatever the stage reached had allocated */
    civ_omega_free(&om);
    free(S_main);
    free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
    free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
    free(X_hat); free(ZtZ_inv_Zty);
    free(XkX_copy); free(beta_temp); free(XkX_inv); free(gmm_hessian_inv); free(coviv_bread);
    free(XtX); free(Xty); free(XkX); free(Xky);
    free(W_bread); free(W_zx);

    return est_rc;
}
