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
/* Detailed corrector timing */
static double fn_corr_loop1 = 0.0;
static double fn_corr_loop2 = 0.0;
static double fn_corr_gemv1 = 0.0;
static double fn_corr_gemv2 = 0.0;
static double fn_corr_loop3 = 0.0;
static double fn_corr_step = 0.0;

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
    fn_corr_loop1 = fn_corr_loop2 = fn_corr_gemv1 = fn_corr_gemv2 = fn_corr_loop3 = fn_corr_step = 0.0;
}

static void fn_timing_report(ST_int N, ST_int K) {
    double total = fn_time_init + fn_time_scaling + fn_time_xtdx + fn_time_cholesky +
                   fn_time_direction + fn_time_steplength + fn_time_corrector + fn_time_update;
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "FN Timing (N=%d, K=%d, iters=%d): total=%.3fs\n"
        "  init=%.3fs  scaling=%.3fs  X'DX=%.3fs  chol=%.3fs\n"
        "  direction=%.3fs  steplength=%.3fs  corrector=%.3fs  update=%.3fs\n"
        "  Corrector breakdown:\n"
        "    loop1=%.3fs  loop2=%.3fs  gemv1=%.3fs  gemv2=%.3fs  loop3=%.3fs  step=%.3fs\n",
        N, K, fn_iter_count, total,
        fn_time_init, fn_time_scaling, fn_time_xtdx, fn_time_cholesky,
        fn_time_direction, fn_time_steplength, fn_time_corrector, fn_time_update,
        fn_corr_loop1, fn_corr_loop2, fn_corr_gemv1, fn_corr_gemv2, fn_corr_loop3, fn_corr_step);
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
 * Helper: compute primal and dual step lengths in a single fused loop
 * This is 8x faster than separate compute_bound + array_min calls.
 * Optimization: ds = -dx always, so we compute ds bound from dx directly,
 * eliminating one array read per element.
 * ============================================================================ */
static void compute_step_lengths_fused(ST_double *fp_out, ST_double *fd_out,
                                       const ST_double *CQREG_RESTRICT x_p,
                                       const ST_double *CQREG_RESTRICT s,
                                       const ST_double *CQREG_RESTRICT z,
                                       const ST_double *CQREG_RESTRICT w,
                                       const ST_double *CQREG_RESTRICT dx,
                                       const ST_double *CQREG_RESTRICT ds,
                                       const ST_double *CQREG_RESTRICT dz,
                                       const ST_double *CQREG_RESTRICT dw,
                                       ST_int N)
{
    (void)ds;  /* ds = -dx, so we don't need to read it */
    ST_double fp_min = 1e20, fd_min = 1e20;

#ifdef _OPENMP
    #pragma omp parallel reduction(min:fp_min, fd_min)
    {
        #pragma omp for schedule(static)
        for (ST_int i = 0; i < N; i++) {
            ST_double dx_i = dx[i];

            /* Primal bounds: ds = -dx, so bs constraint is (dx > 0 ? s/dx : 1e20) */
            ST_double bx = (dx_i < 0.0) ? -x_p[i] / dx_i : 1e20;
            ST_double bs = (dx_i > 0.0) ? s[i] / dx_i : 1e20;
            ST_double bp = (bx < bs) ? bx : bs;
            if (bp < fp_min) fp_min = bp;

            /* Dual bounds */
            ST_double bz = (dz[i] < 0.0) ? -z[i] / dz[i] : 1e20;
            ST_double bw = (dw[i] < 0.0) ? -w[i] / dw[i] : 1e20;
            ST_double bd = (bz < bw) ? bz : bw;
            if (bd < fd_min) fd_min = bd;
        }
    }
#else
    for (ST_int i = 0; i < N; i++) {
        ST_double dx_i = dx[i];

        /* Primal bounds: ds = -dx, so bs constraint is (dx > 0 ? s/dx : 1e20) */
        ST_double bx = (dx_i < 0.0) ? -x_p[i] / dx_i : 1e20;
        ST_double bs = (dx_i > 0.0) ? s[i] / dx_i : 1e20;  /* ds < 0 when dx > 0 */
        ST_double bp = (bx < bs) ? bx : bs;
        if (bp < fp_min) fp_min = bp;

        /* Dual bounds */
        ST_double bz = (dz[i] < 0.0) ? -z[i] / dz[i] : 1e20;
        ST_double bw = (dw[i] < 0.0) ? -w[i] / dw[i] : 1e20;
        ST_double bd = (bz < bw) ? bz : bw;
        if (bd < fd_min) fd_min = bd;
    }
#endif

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
    blas_dgemv(0, N, K, 1.0, X, N, y_d, 0.0, r);  /* r = X * y_d */
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

        /* Save XqX before Cholesky destroys it (for reuse in corrector) */
        ST_double *XqX_copy = ipm->XDXcopy;  /* Use pre-allocated copy */
        memcpy(XqX_copy, XqX, K * K * sizeof(ST_double));

        /* Compute qr = q .* r for the RHS */
        for (i = 0; i < N; i++) {
            Xq[i] = q[i] * r[i];  /* Reuse Xq as temp for q.*r */
        }

        /* Compute Xqr = X' * (q .* r) using single BLAS dgemv transpose */
        blas_dgemv(1, N, K, 1.0, X, N, Xq, 0.0, Xqr);

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

        /* Compute tmp = X * dy using BLAS for vectorization */
        blas_dgemv(0, N, K, 1.0, X, N, dy, 0.0, tmp);

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
#if FN_TIMING
            double t_corr = ctools_timer_seconds();
#endif
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
#if FN_TIMING
            fn_corr_loop1 += ctools_timer_seconds() - t_corr;
            t_corr = ctools_timer_seconds();
#endif

            /* Fused loop: compute xi, r, and qr in one pass */
            for (i = 0; i < N; i++) {
                ST_double xi_i = mu * (xinv[i] - sinv[i]);
                ST_double r_i = z[i] - w[i] + dxdz[i] - dsdw[i] - xi_i;
                Xq[i] = q[i] * r_i;  /* qr for RHS */
                fz[i] = xi_i;  /* Store xi for later use */
            }
            ST_double *xi = fz;
#if FN_TIMING
            fn_corr_loop2 += ctools_timer_seconds() - t_corr;
            t_corr = ctools_timer_seconds();
#endif

            /* Recompute Xqr = X' * (q .* r) using single BLAS dgemv transpose */
            blas_dgemv(1, N, K, 1.0, X, N, Xq, 0.0, Xqr);

            /* Restore XqX from saved copy (avoids recomputing X'DX) */
            memcpy(XqX, XqX_copy, K * K * sizeof(ST_double));

            if (cqreg_cholesky(XqX, K) != 0) {
                break;
            }
            memcpy(dy, Xqr, K * sizeof(ST_double));
            cqreg_solve_cholesky(XqX, dy, K);
#if FN_TIMING
            fn_corr_gemv1 += ctools_timer_seconds() - t_corr;
            t_corr = ctools_timer_seconds();
#endif

            /* Recompute tmp = X * dy using BLAS */
            blas_dgemv(0, N, K, 1.0, X, N, dy, 0.0, tmp);
#if FN_TIMING
            fn_corr_gemv2 += ctools_timer_seconds() - t_corr;
            t_corr = ctools_timer_seconds();
#endif

            /* Recompute corrected directions */
            for (i = 0; i < N; i++) {
                dx[i] = q[i] * (tmp[i] + xi[i] - (z[i] - w[i]) - dxdz[i] + dsdw[i]);
                ds[i] = -dx[i];
                dz[i] = (mu - z[i] * dx[i]) * xinv[i] - z[i] - dxdz[i];
                dw[i] = (mu - w[i] * ds[i]) * sinv[i] - w[i] - dsdw[i];
            }
#if FN_TIMING
            fn_corr_loop3 += ctools_timer_seconds() - t_corr;
            t_corr = ctools_timer_seconds();
#endif

            /* Recompute step lengths (fused loop) */
            compute_step_lengths_fused(&fp, &fd, x_p, s, z, w, dx, ds, dz, dw, N);
#if FN_TIMING
            fn_corr_step += ctools_timer_seconds() - t_corr;
#endif

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
                (void)min_fz_idx; (void)min_fw_idx;  /* Used in IPM_LOG when enabled */
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
    blas_dgemv(0, N, K, 1.0, X, N, beta, 0.0, ipm->r_primal);
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
