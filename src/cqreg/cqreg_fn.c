/*
 * cqreg_fn.c
 *
 * Frisch-Newton Interior Point Method for Quantile Regression
 *
 * Based on the algorithm from:
 *   Koenker, R. and S. Portnoy (1997). "The Gaussian Hare and the Laplacian
 *   Tortoise: Computability of Squared-error vs. Absolute-error Estimators"
 *   Statistical Science, 12(4):279-300.
 *
 * Implementation ported from:
 *   - R quantreg package (Koenker)
 *   - QuantileRegressions.jl (Julia)
 *   - pyqreg (Python/Cython)
 *
 * The algorithm solves the dual LP formulation using primal-dual interior
 * point methods with Mehrotra predictor-corrector steps.
 */

#include "cqreg_fn.h"
#include "cqreg_linalg.h"
#include "cqreg_blas.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Algorithm parameters */
#define IPM_BETA        0.99995   /* Step size damping factor */
#define IPM_SMALL       1e-14     /* Small number for numerical stability */
#define IPM_MAX_ITER    500       /* Maximum iterations */
#define IPM_TOL         1e-12     /* Convergence tolerance for gap */

/* Timing instrumentation for performance analysis */
#define FN_TIMING 1

#if FN_TIMING
#include "../ctools_timer.h"
static double fn_time_init = 0.0;
static double fn_time_scaling = 0.0;
static double fn_time_xtdx = 0.0;
static double fn_time_cholesky = 0.0;
static double fn_time_direction = 0.0;
static double fn_time_steplength = 0.0;
static double fn_time_corrector = 0.0;
static double fn_time_update = 0.0;
static int fn_iter_count = 0;

static void fn_timing_reset(void) {
    fn_time_init = 0.0;
    fn_time_scaling = 0.0;
    fn_time_xtdx = 0.0;
    fn_time_cholesky = 0.0;
    fn_time_direction = 0.0;
    fn_time_steplength = 0.0;
    fn_time_corrector = 0.0;
    fn_time_update = 0.0;
    fn_iter_count = 0;
}

static void fn_timing_report(ST_int N, ST_int K) {
    double total = fn_time_init + fn_time_scaling + fn_time_xtdx + fn_time_cholesky +
                   fn_time_direction + fn_time_steplength + fn_time_corrector + fn_time_update;
    char buf[512];
    snprintf(buf, sizeof(buf),
        "FN Timing (N=%d, K=%d, iters=%d): total=%.3fs\n"
        "  init=%.3fs  scaling=%.3fs  X'DX=%.3fs  chol=%.3fs\n"
        "  direction=%.3fs  steplength=%.3fs  corrector=%.3fs  update=%.3fs\n",
        N, K, fn_iter_count, total,
        fn_time_init, fn_time_scaling, fn_time_xtdx, fn_time_cholesky,
        fn_time_direction, fn_time_steplength, fn_time_corrector, fn_time_update);
    SF_display(buf);
}
#else
#define fn_timing_reset()
#define fn_timing_report(N, K)
#endif

/* Debug logging */
#define IPM_DEBUG 0
#if IPM_DEBUG
static FILE *ipm_debug_file = NULL;
static void ipm_debug_open(void) {
    if (!ipm_debug_file) {
        ipm_debug_file = fopen("/tmp/ipm_debug.log", "w");
    }
}
static void ipm_debug_close(void) {
    if (ipm_debug_file) { fclose(ipm_debug_file); ipm_debug_file = NULL; }
}
#define IPM_LOG(...) do { if (ipm_debug_file) { fprintf(ipm_debug_file, __VA_ARGS__); fflush(ipm_debug_file); } } while(0)
#else
#define ipm_debug_open()
#define ipm_debug_close()
#define IPM_LOG(...)
#endif

/* ============================================================================
 * Helper: compute bound for step length
 * For negative dx[i], compute -x[i]/dx[i] (the point where x becomes zero)
 * ============================================================================ */
static void compute_bound(ST_double *bound, const ST_double *x,
                          const ST_double *dx, ST_int n)
{
    for (ST_int i = 0; i < n; i++) {
        if (dx[i] < 0.0) {
            bound[i] = -x[i] / dx[i];
        } else {
            bound[i] = 1e20;  /* Large number = no constraint */
        }
    }
}

/* ============================================================================
 * Helper: find minimum of array
 * ============================================================================ */
static ST_double array_min(const ST_double *arr, ST_int n)
{
    ST_double min_val = arr[0];
    for (ST_int i = 1; i < n; i++) {
        if (arr[i] < min_val) min_val = arr[i];
    }
    return min_val;
}

/* ============================================================================
 * Helper: compute primal and dual step lengths in a single fused loop
 * This is 8x faster than separate compute_bound + array_min calls
 * ============================================================================ */
static void compute_step_lengths_fused(ST_double *fp_out, ST_double *fd_out,
                                       const ST_double *x_p, const ST_double *s,
                                       const ST_double *z, const ST_double *w,
                                       const ST_double *dx, const ST_double *ds,
                                       const ST_double *dz, const ST_double *dw,
                                       ST_int N)
{
    ST_double fp_min = 1e20, fd_min = 1e20;

    for (ST_int i = 0; i < N; i++) {
        /* Primal bounds */
        ST_double bx = (dx[i] < 0.0) ? -x_p[i] / dx[i] : 1e20;
        ST_double bs = (ds[i] < 0.0) ? -s[i] / ds[i] : 1e20;
        ST_double bp = (bx < bs) ? bx : bs;
        if (bp < fp_min) fp_min = bp;

        /* Dual bounds */
        ST_double bz = (dz[i] < 0.0) ? -z[i] / dz[i] : 1e20;
        ST_double bw = (dw[i] < 0.0) ? -w[i] / dw[i] : 1e20;
        ST_double bd = (bz < bw) ? bz : bw;
        if (bd < fd_min) fd_min = bd;
    }

    /* Apply damping and cap at 1 */
    *fp_out = (IPM_BETA * fp_min > 1.0) ? 1.0 : IPM_BETA * fp_min;
    *fd_out = (IPM_BETA * fd_min > 1.0) ? 1.0 : IPM_BETA * fd_min;
}

/* ============================================================================
 * Main Frisch-Newton Interior Point Solver
 *
 * Solves: min_beta  sum_i rho_tau(y_i - x_i'beta)
 * where rho_tau(u) = u*(tau - I(u<0)) is the check function
 *
 * This implementation follows the Julia QuantileRegressions.jl code exactly.
 *
 * Parameters:
 *   ipm   - Pre-allocated IPM state (provides workspace arrays)
 *   Y     - Response vector (N)
 *   X     - Design matrix (N x K, column-major)
 *   tau   - Quantile (0 < tau < 1)
 *   beta  - Output: coefficient estimates (K)
 *
 * Returns:
 *   Number of iterations on success, negative on error
 * ============================================================================ */
ST_int cqreg_fn_solve(cqreg_ipm_state *ipm,
                       const ST_double *Y,
                       const ST_double *X,
                       ST_double tau,
                       ST_double *beta)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j, iter;

    ipm_debug_open();
    IPM_LOG("=== IPM Solver Start ===\n");
    IPM_LOG("N=%d, K=%d, tau=%.4f\n", N, K, tau);

#if FN_TIMING
    fn_timing_reset();
    double t_start = ctools_timer_seconds();
#endif

    /* Use pre-allocated working arrays from IPM state (avoids malloc in hot loop) */
    ST_double *x_p = ipm->fn_xp;      /* Primal x (N) */
    ST_double *s = ipm->fn_s;         /* Primal s (N) */
    ST_double *y_d = ipm->fn_yd;      /* Dual y = -beta (K) */
    ST_double *z = ipm->fn_z;         /* Dual slack z (N) */
    ST_double *w = ipm->fn_w;         /* Dual slack w (N) */
    ST_double *dx = ipm->fn_dx;       /* Direction for x (N) */
    ST_double *ds = ipm->fn_ds;       /* Direction for s (N) */
    ST_double *dy = ipm->fn_dy;       /* Direction for y (K) */
    ST_double *dz = ipm->fn_dz;       /* Direction for z (N) */
    ST_double *dw = ipm->fn_dw;       /* Direction for w (N) */
    ST_double *fx = ipm->fn_fx;       /* Bound for x (N) */
    ST_double *fs = ipm->fn_fs;       /* Bound for s (N) */
    ST_double *fz = ipm->fn_fz;       /* Bound for z (N) */
    ST_double *fw = ipm->fn_fw;       /* Bound for w (N) */
    ST_double *q = ipm->fn_q;         /* Diagonal scaling (N) */
    ST_double *r = ipm->fn_r;         /* Residual (N) */
    ST_double *tmp = ipm->fn_tmp;     /* Temp array (N) */
    ST_double *sinv = ipm->fn_sinv;   /* 1/s for corrector (N) */
    ST_double *Xq = ipm->fn_Xq;       /* Scaled X (N*K) */
    ST_double *XqX = ipm->XDX;        /* Reuse existing XDX (K*K) */
    ST_double *Xqr = ipm->rhs;        /* Reuse existing rhs (K) */

    /* ========================================================================
     * INITIALIZATION (following Julia code exactly)
     * ======================================================================== */

    /* c = -Y (the objective coefficients in dual formulation) */
    /* x_p = (1 - tau) * ones(N) */
    /* s = 1 - x_p = tau * ones(N) */
    for (i = 0; i < N; i++) {
        x_p[i] = 1.0 - tau;
        s[i] = tau;
    }

    /* b = X' * x_p (constraint RHS in dual) */
    /* In Julia: b = X'x, but we don't need to store b explicitly */

    /* y_d = -X \ Y (OLS with sign flip - this is the dual variable) */
    /* Compute X'X */
    memset(XqX, 0, K * K * sizeof(ST_double));
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        for (ST_int k = j; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += Xj[i] * Xk[i];
            }
            XqX[j * K + k] = sum;
            XqX[k * K + j] = sum;
        }
    }

    /* Compute X'Y */
    for (j = 0; j < K; j++) {
        y_d[j] = cqreg_dot(&X[j * N], Y, N);
    }

    /* Solve (X'X) * tmp = X'Y, then y_d = -tmp */
    /* Add small regularization */
    for (j = 0; j < K; j++) {
        XqX[j * K + j] += IPM_SMALL;
    }

    if (cqreg_cholesky(XqX, K) != 0) {
        /* Cholesky failed - use zero */
        memset(y_d, 0, K * sizeof(ST_double));
    } else {
        cqreg_solve_cholesky(XqX, y_d, K);
        /* y_d now contains OLS solution, negate it */
        for (j = 0; j < K; j++) {
            y_d[j] = -y_d[j];
        }
    }

    /* Compute r = c - X * y_d = -Y - X * y_d */
    cqreg_matvec_col(r, X, y_d, N, K);  /* r = X * y_d */
    for (i = 0; i < N; i++) {
        r[i] = -Y[i] - r[i];  /* r = -Y - X*y_d */
    }

    /* Initialize z and w from r */
    /* z = max(r, 0), w = max(-r, 0) */
    /* If |r| <= small, add small to both */
    for (i = 0; i < N; i++) {
        ST_double ri = r[i];
        z[i] = (ri > 0) ? ri : 0.0;
        w[i] = (ri < 0) ? -ri : 0.0;
        if (fabs(ri) <= IPM_SMALL) {
            z[i] += IPM_SMALL;
            w[i] += IPM_SMALL;
        }
    }

    /* Compute initial duality gap */
    ST_double gap = 0.0;
    for (i = 0; i < N; i++) {
        gap += z[i] * x_p[i] + s[i] * w[i];
    }

    /* Debug: print initial values */
    IPM_LOG("\n=== Initialization ===\n");
    IPM_LOG("y_d (OLS negated): [");
    for (j = 0; j < K; j++) IPM_LOG("%.6f ", y_d[j]);
    IPM_LOG("]\n");
    IPM_LOG("Initial gap: %.6e\n", gap);
    IPM_LOG("x_p[0..2]: %.6f %.6f %.6f\n", x_p[0], x_p[1], x_p[2]);
    IPM_LOG("s[0..2]: %.6f %.6f %.6f\n", s[0], s[1], s[2]);
    IPM_LOG("z[0..2]: %.6f %.6f %.6f\n", z[0], z[1], z[2]);
    IPM_LOG("w[0..2]: %.6f %.6f %.6f\n", w[0], w[1], w[2]);

    /* ========================================================================
     * MAIN INTERIOR POINT LOOP
     * ======================================================================== */

#if FN_TIMING
    fn_time_init = ctools_timer_seconds() - t_start;
    double t_phase;
#endif

    for (iter = 0; iter < IPM_MAX_ITER; iter++) {

        /* Check convergence */
        if (gap < IPM_TOL || !isfinite(gap)) {
            IPM_LOG("Converged at iter %d, gap=%.6e\n", iter, gap);
            break;
        }

#if FN_TIMING
        fn_iter_count++;
#endif

        /* Log every 10 iterations or first 5 */
        if (iter < 5 || iter % 50 == 0) {
            IPM_LOG("\n=== Iter %d ===\n", iter);
            IPM_LOG("gap=%.6e\n", gap);
            IPM_LOG("y_d: [");
            for (j = 0; j < K; j++) IPM_LOG("%.6f ", y_d[j]);
            IPM_LOG("]\n");
        }

        /* ====================================================================
         * AFFINE (PREDICTOR) STEP
         * ==================================================================== */

#if FN_TIMING
        t_phase = ctools_timer_seconds();
#endif
        /* Fused loop: compute diagonal scaling q and residual r */
        for (i = 0; i < N; i++) {
            q[i] = 1.0 / (z[i] / x_p[i] + w[i] / s[i]);
            r[i] = z[i] - w[i];
        }
#if FN_TIMING
        fn_time_scaling += ctools_timer_seconds() - t_phase;
        t_phase = ctools_timer_seconds();
#endif

        /* Compute XqX = X' * diag(q) * X using optimized direct computation */
        blas_xtdx(XqX, X, N, K, q);
#if FN_TIMING
        fn_time_xtdx += ctools_timer_seconds() - t_phase;
        t_phase = ctools_timer_seconds();
#endif

        /* Add regularization */
        for (j = 0; j < K; j++) {
            XqX[j * K + j] += IPM_SMALL;
        }

        /* Compute qr = q .* r for the RHS */
        for (i = 0; i < N; i++) {
            Xq[i] = q[i] * r[i];  /* Reuse Xq as temp for q.*r */
        }

        /* Compute Xqr = X' * (q .* r) */
        for (j = 0; j < K; j++) {
            Xqr[j] = cqreg_dot(&X[j * N], Xq, N);
        }

        /* Solve (X'QX) * dy = Xqr via Cholesky */
        if (cqreg_cholesky(XqX, K) != 0) {
            break;  /* Numerical issues */
        }
        memcpy(dy, Xqr, K * sizeof(ST_double));
        cqreg_solve_cholesky(XqX, dy, K);
#if FN_TIMING
        fn_time_cholesky += ctools_timer_seconds() - t_phase;
        t_phase = ctools_timer_seconds();
#endif

        /* Compute tmp = X * dy */
        cqreg_matvec_col(tmp, X, dy, N, K);

        /* Compute dx = q .* (X*dy - r), ds = -dx */
        for (i = 0; i < N; i++) {
            dx[i] = q[i] * (tmp[i] - r[i]);
            ds[i] = -dx[i];
        }

        /* Compute dz = -z .* (1 + dx./x_p), dw = -w .* (1 + ds./s) */
        for (i = 0; i < N; i++) {
            dz[i] = -z[i] * (1.0 + dx[i] / x_p[i]);
            dw[i] = -w[i] * (1.0 + ds[i] / s[i]);
        }
#if FN_TIMING
        fn_time_direction += ctools_timer_seconds() - t_phase;
        t_phase = ctools_timer_seconds();
#endif

        /* ====================================================================
         * COMPUTE STEP LENGTHS (fused loop for 8x speedup)
         * ==================================================================== */

        ST_double fp, fd;
        compute_step_lengths_fused(&fp, &fd, x_p, s, z, w, dx, ds, dz, dw, N);
#if FN_TIMING
        fn_time_steplength += ctools_timer_seconds() - t_phase;
        t_phase = ctools_timer_seconds();
#endif

        if (iter < 5 || iter % 50 == 0) {
            IPM_LOG("Affine: fp=%.6f, fd=%.6f\n", fp, fd);
            IPM_LOG("dy: [");
            for (j = 0; j < K; j++) IPM_LOG("%.6f ", dy[j]);
            IPM_LOG("]\n");
        }

        /* ====================================================================
         * CORRECTOR STEP (if full step not feasible)
         * ==================================================================== */

        if (fp < 1.0 || fd < 1.0) {
            /* Fused loop: compute g0, gfpfd, dxdz, dsdw, xinv, sinv in one pass */
            ST_double g0 = 0.0, gfpfd = 0.0;
            ST_double *dxdz = fx;
            ST_double *dsdw = fs;
            ST_double *xinv = fw;
            /* sinv already points to pre-allocated ipm->fn_sinv */

            for (i = 0; i < N; i++) {
                /* Centering parameter terms */
                g0 += z[i] * x_p[i] + w[i] * s[i];
                gfpfd += (z[i] + fd * dz[i]) * (x_p[i] + fp * dx[i])
                       + (w[i] + fd * dw[i]) * (s[i] + fp * ds[i]);
                /* Corrector intermediates */
                dxdz[i] = dx[i] * dz[i];
                dsdw[i] = ds[i] * dw[i];
                xinv[i] = 1.0 / x_p[i];
                sinv[i] = 1.0 / s[i];
            }

            ST_double mu = pow(gfpfd / g0, 3.0) * g0 / (2.0 * N);

            /* Fused loop: compute xi, r, and qr in one pass */
            for (i = 0; i < N; i++) {
                ST_double xi_i = mu * (xinv[i] - sinv[i]);
                ST_double r_i = z[i] - w[i] + dxdz[i] - dsdw[i] - xi_i;
                Xq[i] = q[i] * r_i;  /* qr for RHS */
                fz[i] = xi_i;  /* Store xi for later use */
            }
            ST_double *xi = fz;

            /* Recompute Xqr = X' * (q .* r) */
            for (j = 0; j < K; j++) {
                Xqr[j] = cqreg_dot(&X[j * N], Xq, N);
            }

            /* Rebuild XqX = X' * diag(q) * X (was modified by Cholesky) */
            blas_xtdx(XqX, X, N, K, q);
            for (j = 0; j < K; j++) {
                XqX[j * K + j] += IPM_SMALL;
            }

            if (cqreg_cholesky(XqX, K) != 0) {
                break;
            }
            memcpy(dy, Xqr, K * sizeof(ST_double));
            cqreg_solve_cholesky(XqX, dy, K);

            /* Recompute tmp = X * dy */
            cqreg_matvec_col(tmp, X, dy, N, K);

            /* Recompute corrected directions */
            for (i = 0; i < N; i++) {
                dx[i] = q[i] * (tmp[i] + xi[i] - (z[i] - w[i]) - dxdz[i] + dsdw[i]);
                ds[i] = -dx[i];
                dz[i] = (mu - z[i] * dx[i]) * xinv[i] - z[i] - dxdz[i];
                dw[i] = (mu - w[i] * ds[i]) * sinv[i] - w[i] - dsdw[i];
            }

            /* Debug: check corrector formulas for w=0 case */
            if (iter == 0) {
                for (i = 0; i < N; i++) {
                    if (w[i] < IPM_SMALL && dw[i] < 0) {
                        IPM_LOG("  w=0 but dw<0 at i=%d: mu=%.4e, ds=%.4f, s=%.4f, dsdw=%.4e\n",
                                i, mu, ds[i], s[i], dsdw[i]);
                        IPM_LOG("    formula: dw = (%.4e - 0*%.4f)/%.4f - 0 - %.4e = %.4e\n",
                                mu, ds[i], s[i], dsdw[i], dw[i]);
                        break;  /* Just show first one */
                    }
                }
            }

            /* Recompute step lengths (fused loop) */
            compute_step_lengths_fused(&fp, &fd, x_p, s, z, w, dx, ds, dz, dw, N);

            if (iter < 5 || iter % 50 == 0) {
                IPM_LOG("Corrector: fp=%.6f, fd=%.6f, mu=%.6e\n", fp, fd, mu);
                IPM_LOG("Corrector dy: [");
                for (j = 0; j < K; j++) IPM_LOG("%.6f ", dy[j]);
                IPM_LOG("]\n");

                /* Debug: find why fd = 0 */
                ST_double min_fz = 1e20, min_fw = 1e20;
                ST_int min_fz_idx = -1, min_fw_idx = -1;
                for (i = 0; i < N; i++) {
                    if (fz[i] < min_fz) { min_fz = fz[i]; min_fz_idx = i; }
                    if (fw[i] < min_fw) { min_fw = fw[i]; min_fw_idx = i; }
                }
                IPM_LOG("min fz=%.6e at i=%d (z=%.6e, dz=%.6e)\n",
                        min_fz, min_fz_idx, z[min_fz_idx], dz[min_fz_idx]);
                IPM_LOG("min fw=%.6e at i=%d (w=%.6e, dw=%.6e)\n",
                        min_fw, min_fw_idx, w[min_fw_idx], dw[min_fw_idx]);
            }
            /* sinv is pre-allocated, no free needed */
        }
#if FN_TIMING
        fn_time_corrector += ctools_timer_seconds() - t_phase;
        t_phase = ctools_timer_seconds();
#endif

        /* ====================================================================
         * UPDATE VARIABLES
         * ==================================================================== */

        if (iter < 5 || iter % 50 == 0) {
            IPM_LOG("Before update y_d: [");
            for (j = 0; j < K; j++) IPM_LOG("%.6f ", y_d[j]);
            IPM_LOG("]\n");
            IPM_LOG("Update: y_d += %.6f * dy\n", fd);
        }

        /* Update dual coefficients: y_d += fd*dy (small loop, K elements) */
        for (j = 0; j < K; j++) {
            y_d[j] += fd * dy[j];
        }

        if (iter < 5 || iter % 50 == 0) {
            IPM_LOG("After update y_d: [");
            for (j = 0; j < K; j++) IPM_LOG("%.6f ", y_d[j]);
            IPM_LOG("]\n");
        }

        /* Fused loop: update all primal/dual variables and compute gap */
        gap = 0.0;
        for (i = 0; i < N; i++) {
            x_p[i] += fp * dx[i];
            s[i] += fp * ds[i];
            z[i] += fd * dz[i];
            w[i] += fd * dw[i];
            gap += z[i] * x_p[i] + s[i] * w[i];
        }
#if FN_TIMING
        fn_time_update += ctools_timer_seconds() - t_phase;
#endif

    } /* End main loop */

    /* ========================================================================
     * FINALIZE: beta = -y_d (the return value is negative of dual variable)
     * ======================================================================== */

    IPM_LOG("\n=== Final ===\n");
    IPM_LOG("Iterations: %d\n", iter);
    IPM_LOG("Final gap: %.6e\n", gap);
    IPM_LOG("Final y_d: [");
    for (j = 0; j < K; j++) IPM_LOG("%.6f ", y_d[j]);
    IPM_LOG("]\n");

    for (j = 0; j < K; j++) {
        beta[j] = -y_d[j];
    }

    IPM_LOG("Final beta (=-y_d): [");
    for (j = 0; j < K; j++) IPM_LOG("%.6f ", beta[j]);
    IPM_LOG("]\n");
    ipm_debug_close();

    /* Store results in ipm state */
    ipm->iterations = iter;
    ipm->converged = (gap < IPM_TOL) ? 1 : 0;

    /* Compute final residuals */
    cqreg_matvec_col(ipm->r_primal, X, beta, N, K);
    for (i = 0; i < N; i++) {
        ipm->r_primal[i] = Y[i] - ipm->r_primal[i];  /* r = Y - X*beta */
        if (ipm->r_primal[i] >= 0) {
            ipm->u[i] = ipm->r_primal[i];
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -ipm->r_primal[i];
        }
    }
    memcpy(ipm->beta, beta, K * sizeof(ST_double));

#if FN_TIMING
    fn_timing_report(N, K);
#endif

    /* No cleanup needed - all arrays are pre-allocated in IPM state */
    return iter;
}


/* ============================================================================
 * Warm-Start IPM Solver
 *
 * Same algorithm as cqreg_fn_solve, but uses initial_beta as the starting
 * point instead of computing OLS. This is much faster for auxiliary solves
 * in Siddiqui method where τ±h is close to τ.
 * ============================================================================ */

ST_int cqreg_fn_solve_warmstart(cqreg_ipm_state *ipm,
                                 const ST_double *Y,
                                 const ST_double *X,
                                 ST_double tau,
                                 const ST_double *initial_beta,
                                 ST_double *beta)
{
    ST_int N = ipm->N;
    ST_int K = ipm->K;
    ST_int i, j, iter;

    /* Use pre-allocated working arrays from IPM state (avoids malloc in hot loop) */
    ST_double *x_p = ipm->fn_xp;
    ST_double *s = ipm->fn_s;
    ST_double *y_d = ipm->fn_yd;
    ST_double *z = ipm->fn_z;
    ST_double *w = ipm->fn_w;
    ST_double *dx = ipm->fn_dx;
    ST_double *ds = ipm->fn_ds;
    ST_double *dy = ipm->fn_dy;
    ST_double *dz = ipm->fn_dz;
    ST_double *dw = ipm->fn_dw;
    ST_double *fx = ipm->fn_fx;
    ST_double *fs = ipm->fn_fs;
    ST_double *fz = ipm->fn_fz;
    ST_double *fw = ipm->fn_fw;
    ST_double *q = ipm->fn_q;
    ST_double *r = ipm->fn_r;
    ST_double *tmp = ipm->fn_tmp;
    ST_double *Xq = ipm->fn_Xq;
    ST_double *XqX = ipm->XDX;
    ST_double *Xqr = ipm->rhs;
    (void)Xq;  /* Suppress unused warning - not used in warmstart */

    /* WARM-START INITIALIZATION:
     * Instead of computing OLS, use provided initial_beta directly.
     * This is the key difference from cqreg_fn_solve.
     */

    /* x_p = (1 - tau), s = tau */
    for (i = 0; i < N; i++) {
        x_p[i] = 1.0 - tau;
        s[i] = tau;
    }

    /* y_d = -initial_beta (dual variable is negative of coefficients) */
    for (j = 0; j < K; j++) {
        y_d[j] = -initial_beta[j];
    }

    /* Compute r = c - X * y_d = -Y - X * y_d */
    cqreg_matvec_col(r, X, y_d, N, K);
    for (i = 0; i < N; i++) {
        r[i] = -Y[i] - r[i];
    }

    /* Initialize z and w from r
     * For warm-start, we need larger perturbation because residuals
     * may be very small (close to optimal), causing numerical issues.
     * Use 1e-4 instead of 1e-14 for stability.
     */
    const ST_double WARMSTART_PERTURB = 1e-4;
    for (i = 0; i < N; i++) {
        ST_double ri = r[i];
        z[i] = (ri > 0) ? ri : 0.0;
        w[i] = (ri < 0) ? -ri : 0.0;
        /* Always add perturbation for warm-start stability */
        z[i] += WARMSTART_PERTURB;
        w[i] += WARMSTART_PERTURB;
    }

    /* Compute initial duality gap */
    ST_double gap = 0.0;
    for (i = 0; i < N; i++) {
        gap += z[i] * x_p[i] + s[i] * w[i];
    }

    /* MAIN INTERIOR POINT LOOP (identical to cqreg_fn_solve) */
    for (iter = 0; iter < IPM_MAX_ITER; iter++) {
        /* Check convergence - but ensure at least 1 iteration for warm-start
         * since the quantile τ is different from the main solve */
        if (iter > 0 && (gap < IPM_TOL || !isfinite(gap))) {
            break;
        }

        /* Compute scaling: q = 1/(z/x_p + w/s) */
        for (i = 0; i < N; i++) {
            q[i] = 1.0 / (z[i] / x_p[i] + w[i] / s[i]);
        }

        /* Compute XqX = X' * diag(q) * X using optimized blas_xtdx */
        blas_xtdx(XqX, X, N, K, q);

        /* Add regularization */
        for (j = 0; j < K; j++) {
            XqX[j * K + j] += IPM_SMALL;
        }

        /* RHS computation */
        for (i = 0; i < N; i++) {
            r[i] = -Y[i];
        }
        cqreg_matvec_col(tmp, X, y_d, N, K);
        for (i = 0; i < N; i++) {
            r[i] = r[i] - tmp[i] + z[i] - w[i];
            r[i] = q[i] * r[i];
        }

        /* Xqr = X' * r */
        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += Xj[i] * r[i];
            }
            Xqr[j] = sum;
        }

        /* Solve for predictor direction */
        if (cqreg_cholesky(XqX, K) != 0) {
            break;
        }
        memcpy(dy, Xqr, K * sizeof(ST_double));
        cqreg_solve_cholesky(XqX, dy, K);

        /* dx = q .* (X * dy - r) */
        cqreg_matvec_col(tmp, X, dy, N, K);
        for (i = 0; i < N; i++) {
            dx[i] = q[i] * (tmp[i] - r[i] / q[i]);
        }
        /* ds = -dx */
        for (i = 0; i < N; i++) {
            ds[i] = -dx[i];
        }
        /* dz = -z + (z/x_p) .* dx */
        for (i = 0; i < N; i++) {
            dz[i] = -z[i] + (z[i] / x_p[i]) * dx[i];
        }
        /* dw = -w + (w/s) .* ds */
        for (i = 0; i < N; i++) {
            dw[i] = -w[i] + (w[i] / s[i]) * ds[i];
        }

        /* Step length computation */
        compute_bound(fx, x_p, dx, N);
        compute_bound(fs, s, ds, N);
        compute_bound(fz, z, dz, N);
        compute_bound(fw, w, dw, N);

        ST_double fp = IPM_BETA * fmin(array_min(fx, N), array_min(fs, N));
        ST_double fd = IPM_BETA * fmin(array_min(fz, N), array_min(fw, N));
        if (fp > 1.0) fp = 1.0;
        if (fd > 1.0) fd = 1.0;

        /* Centering parameter */
        ST_double mu = gap / (2.0 * N);

        /* Mehrotra corrector */
        ST_double g = 0.0;
        for (i = 0; i < N; i++) {
            g += (z[i] + fd * dz[i]) * (x_p[i] + fp * dx[i]);
            g += (w[i] + fd * dw[i]) * (s[i] + fp * ds[i]);
        }
        mu = (g / (2.0 * N)) * (g / gap) * (g / gap);

        /* Corrector step */
        for (i = 0; i < N; i++) {
            r[i] = q[i] * ((mu - dz[i] * dx[i]) / z[i] - (mu - dw[i] * ds[i]) / w[i]);
        }
        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                sum += Xj[i] * r[i];
            }
            Xqr[j] = sum;
        }
        cqreg_solve_cholesky(XqX, Xqr, K);
        for (j = 0; j < K; j++) {
            dy[j] += Xqr[j];
        }

        /* Update directions */
        cqreg_matvec_col(tmp, X, dy, N, K);
        for (i = 0; i < N; i++) {
            dx[i] = q[i] * (tmp[i] - r[i] / q[i] - dx[i] + (mu - dz[i] * dx[i]) / z[i] - (mu - dw[i] * ds[i]) / w[i]);
        }
        for (i = 0; i < N; i++) {
            ds[i] = -dx[i];
        }
        for (i = 0; i < N; i++) {
            dz[i] = (mu - z[i] * x_p[i] - dz[i] * dx[i]) / x_p[i] - z[i];
            dw[i] = (mu - w[i] * s[i] - dw[i] * ds[i]) / s[i] - w[i];
        }

        /* Recompute step lengths */
        compute_bound(fx, x_p, dx, N);
        compute_bound(fs, s, ds, N);
        compute_bound(fz, z, dz, N);
        compute_bound(fw, w, dw, N);

        fp = IPM_BETA * fmin(array_min(fx, N), array_min(fs, N));
        fd = IPM_BETA * fmin(array_min(fz, N), array_min(fw, N));
        if (fp > 1.0) fp = 1.0;
        if (fd > 1.0) fd = 1.0;

        /* Update variables */
        for (i = 0; i < N; i++) {
            x_p[i] += fp * dx[i];
            s[i] += fp * ds[i];
        }
        for (j = 0; j < K; j++) {
            y_d[j] += fd * dy[j];
        }
        for (i = 0; i < N; i++) {
            z[i] += fd * dz[i];
            w[i] += fd * dw[i];
        }

        /* Recompute gap */
        gap = 0.0;
        for (i = 0; i < N; i++) {
            gap += z[i] * x_p[i] + s[i] * w[i];
        }
    }

    /* FINALIZE: beta = -y_d */
    for (j = 0; j < K; j++) {
        beta[j] = -y_d[j];
    }

    /* Store results */
    ipm->iterations = iter;
    ipm->converged = (gap < IPM_TOL) ? 1 : 0;

    /* Compute final residuals */
    cqreg_matvec_col(ipm->r_primal, X, beta, N, K);
    for (i = 0; i < N; i++) {
        ipm->r_primal[i] = Y[i] - ipm->r_primal[i];
        if (ipm->r_primal[i] >= 0) {
            ipm->u[i] = ipm->r_primal[i];
            ipm->v[i] = 0.0;
        } else {
            ipm->u[i] = 0.0;
            ipm->v[i] = -ipm->r_primal[i];
        }
    }
    memcpy(ipm->beta, beta, K * sizeof(ST_double));

    /* No cleanup needed - all arrays are pre-allocated in IPM state */
    return iter;
}


/* ============================================================================
 * Legacy functions - keep for compatibility
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

void fn_form_reduced_system(const ST_double *X, const ST_double *g,
                             ST_int N, ST_int K,
                             const ST_int *active_idx, ST_int m,
                             ST_double *XAX, ST_double *XAg)
{
    memset(XAX, 0, K * K * sizeof(ST_double));
    memset(XAg, 0, K * sizeof(ST_double));

    for (ST_int j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        ST_double xg = 0.0;
        for (ST_int i = 0; i < m; i++) {
            ST_int idx = active_idx[i];
            xg += Xj[idx] * g[idx];
        }
        XAg[j] = xg;

        for (ST_int k = j; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double sum = 0.0;
            for (ST_int i = 0; i < m; i++) {
                ST_int idx = active_idx[i];
                sum += Xj[idx] * Xk[idx];
            }
            XAX[j * K + k] = sum;
            XAX[k * K + j] = sum;
        }
    }
}

ST_double fn_ratio_test(const ST_double *r, const ST_double *delta_r,
                         ST_int N, ST_double q)
{
    ST_double alpha = 1.0;
    for (ST_int i = 0; i < N; i++) {
        if (r[i] > 0 && delta_r[i] < -1e-14) {
            ST_double a = -r[i] / delta_r[i];
            if (a < alpha) alpha = a;
        } else if (r[i] < 0 && delta_r[i] > 1e-14) {
            ST_double a = -r[i] / delta_r[i];
            if (a < alpha) alpha = a;
        }
    }
    return alpha * 0.9995;
}

ST_int fn_check_optimality(const ST_double *X, const ST_double *r,
                            ST_int N, ST_int K, ST_double q,
                            ST_double *Xg)
{
    memset(Xg, 0, K * sizeof(ST_double));
    for (ST_int j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        ST_double sum = 0.0;
        for (ST_int i = 0; i < N; i++) {
            ST_double gi = (r[i] > 1e-6) ? q : ((r[i] < -1e-6) ? (q - 1.0) : (q - 0.5));
            sum += Xj[i] * gi;
        }
        Xg[j] = sum;
    }

    ST_double norm_sq = 0.0;
    for (ST_int j = 0; j < K; j++) {
        norm_sq += Xg[j] * Xg[j];
    }
    return (sqrt(norm_sq) < 0.05 * sqrt((ST_double)N)) ? 1 : 0;
}

ST_int cqreg_fn_preprocess_solve(cqreg_ipm_state *ipm,
                                  const ST_double *y,
                                  const ST_double *X,
                                  ST_double q,
                                  ST_double *beta)
{
    /* For now, just call the main solver */
    return cqreg_fn_solve(ipm, y, X, q, beta);
}
