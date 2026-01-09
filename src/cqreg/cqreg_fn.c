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
#define IPM_TOL_DEFAULT 1e-12     /* Default convergence tolerance for gap */

/* SIMD support for ARM NEON and x86 SSE/AVX */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define CQREG_USE_NEON 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define CQREG_USE_SSE2 1
#endif

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

static void fn_timing_report(ST_int N, ST_int K, ST_int verbose) {
    if (!verbose) return;  /* Only report timing in verbose mode */
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
#define fn_timing_report(N, K, verbose)
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
 * SIMD-optimized corrector loop helpers
 * These optimize the hottest loops in the corrector step (loop1, loop3)
 * ============================================================================ */

#if CQREG_USE_NEON
/*
 * NEON-optimized loop1: compute g0, gfpfd, and corrector intermediates
 * Processes 2 doubles per iteration using NEON float64x2 vectors
 */
static void corrector_loop1_neon(
    ST_double *g0_out, ST_double *gfpfd_out,
    ST_double *CQREG_RESTRICT dxdz,
    ST_double *CQREG_RESTRICT dsdw,
    ST_double *CQREG_RESTRICT xinv,
    ST_double *CQREG_RESTRICT sinv,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT x_p,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT s,
    const ST_double *CQREG_RESTRICT dz,
    const ST_double *CQREG_RESTRICT dx,
    const ST_double *CQREG_RESTRICT dw,
    const ST_double *CQREG_RESTRICT ds,
    ST_double fp, ST_double fd, ST_int N)
{
    float64x2_t g0_vec = vdupq_n_f64(0.0);
    float64x2_t gfpfd_vec = vdupq_n_f64(0.0);
    float64x2_t fp_vec = vdupq_n_f64(fp);
    float64x2_t fd_vec = vdupq_n_f64(fd);
    float64x2_t one_vec = vdupq_n_f64(1.0);

    ST_int N2 = N - (N % 2);
    ST_int i;

    for (i = 0; i < N2; i += 2) {
        /* Load vectors */
        float64x2_t z_v = vld1q_f64(&z[i]);
        float64x2_t xp_v = vld1q_f64(&x_p[i]);
        float64x2_t w_v = vld1q_f64(&w[i]);
        float64x2_t s_v = vld1q_f64(&s[i]);
        float64x2_t dz_v = vld1q_f64(&dz[i]);
        float64x2_t dx_v = vld1q_f64(&dx[i]);
        float64x2_t dw_v = vld1q_f64(&dw[i]);
        float64x2_t ds_v = vld1q_f64(&ds[i]);

        /* g0 += z[i] * x_p[i] + w[i] * s[i] */
        g0_vec = vfmaq_f64(g0_vec, z_v, xp_v);
        g0_vec = vfmaq_f64(g0_vec, w_v, s_v);

        /* gfpfd terms */
        float64x2_t z_fd_dz = vfmaq_f64(z_v, fd_vec, dz_v);  /* z + fd*dz */
        float64x2_t xp_fp_dx = vfmaq_f64(xp_v, fp_vec, dx_v); /* xp + fp*dx */
        float64x2_t w_fd_dw = vfmaq_f64(w_v, fd_vec, dw_v);   /* w + fd*dw */
        float64x2_t s_fp_ds = vfmaq_f64(s_v, fp_vec, ds_v);   /* s + fp*ds */

        gfpfd_vec = vfmaq_f64(gfpfd_vec, z_fd_dz, xp_fp_dx);
        gfpfd_vec = vfmaq_f64(gfpfd_vec, w_fd_dw, s_fp_ds);

        /* Corrector intermediates */
        vst1q_f64(&dxdz[i], vmulq_f64(dx_v, dz_v));
        vst1q_f64(&dsdw[i], vmulq_f64(ds_v, dw_v));
        vst1q_f64(&xinv[i], vdivq_f64(one_vec, xp_v));
        vst1q_f64(&sinv[i], vdivq_f64(one_vec, s_v));
    }

    /* Horizontal sum */
    ST_double g0 = vgetq_lane_f64(g0_vec, 0) + vgetq_lane_f64(g0_vec, 1);
    ST_double gfpfd = vgetq_lane_f64(gfpfd_vec, 0) + vgetq_lane_f64(gfpfd_vec, 1);

    /* Handle remainder */
    for (; i < N; i++) {
        g0 += z[i] * x_p[i] + w[i] * s[i];
        gfpfd += (z[i] + fd * dz[i]) * (x_p[i] + fp * dx[i])
               + (w[i] + fd * dw[i]) * (s[i] + fp * ds[i]);
        dxdz[i] = dx[i] * dz[i];
        dsdw[i] = ds[i] * dw[i];
        xinv[i] = 1.0 / x_p[i];
        sinv[i] = 1.0 / s[i];
    }

    *g0_out = g0;
    *gfpfd_out = gfpfd;
}

/*
 * NEON-optimized loop3: recompute corrected directions
 * This is the most expensive loop in the corrector step
 */
static void corrector_loop3_neon(
    ST_double *CQREG_RESTRICT dx,
    ST_double *CQREG_RESTRICT ds,
    ST_double *CQREG_RESTRICT dz,
    ST_double *CQREG_RESTRICT dw,
    const ST_double *CQREG_RESTRICT q,
    const ST_double *CQREG_RESTRICT tmp,
    const ST_double *CQREG_RESTRICT xi,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT dxdz,
    const ST_double *CQREG_RESTRICT dsdw,
    const ST_double *CQREG_RESTRICT xinv,
    const ST_double *CQREG_RESTRICT sinv,
    ST_double mu, ST_int N)
{
    float64x2_t mu_vec = vdupq_n_f64(mu);
    ST_int N2 = N - (N % 2);
    ST_int i;

    for (i = 0; i < N2; i += 2) {
        /* Load vectors */
        float64x2_t q_v = vld1q_f64(&q[i]);
        float64x2_t tmp_v = vld1q_f64(&tmp[i]);
        float64x2_t xi_v = vld1q_f64(&xi[i]);
        float64x2_t z_v = vld1q_f64(&z[i]);
        float64x2_t w_v = vld1q_f64(&w[i]);
        float64x2_t dxdz_v = vld1q_f64(&dxdz[i]);
        float64x2_t dsdw_v = vld1q_f64(&dsdw[i]);
        float64x2_t xinv_v = vld1q_f64(&xinv[i]);
        float64x2_t sinv_v = vld1q_f64(&sinv[i]);

        /* dx[i] = q[i] * (tmp[i] + xi[i] - (z[i] - w[i]) - dxdz[i] + dsdw[i]) */
        float64x2_t z_minus_w = vsubq_f64(z_v, w_v);
        float64x2_t inner = vsubq_f64(vaddq_f64(tmp_v, xi_v), z_minus_w);
        inner = vsubq_f64(inner, dxdz_v);
        inner = vaddq_f64(inner, dsdw_v);
        float64x2_t dx_v = vmulq_f64(q_v, inner);

        /* ds[i] = -dx[i] */
        float64x2_t ds_v = vnegq_f64(dx_v);

        /* dz[i] = (mu - z[i] * dx[i]) * xinv[i] - z[i] - dxdz[i] */
        float64x2_t z_dx = vmulq_f64(z_v, dx_v);
        float64x2_t mu_minus_zdx = vsubq_f64(mu_vec, z_dx);
        float64x2_t dz_v = vfmsq_f64(vsubq_f64(vmulq_f64(mu_minus_zdx, xinv_v), z_v), vdupq_n_f64(1.0), dxdz_v);

        /* dw[i] = (mu - w[i] * ds[i]) * sinv[i] - w[i] - dsdw[i] */
        float64x2_t w_ds = vmulq_f64(w_v, ds_v);
        float64x2_t mu_minus_wds = vsubq_f64(mu_vec, w_ds);
        float64x2_t dw_v = vfmsq_f64(vsubq_f64(vmulq_f64(mu_minus_wds, sinv_v), w_v), vdupq_n_f64(1.0), dsdw_v);

        /* Store results */
        vst1q_f64(&dx[i], dx_v);
        vst1q_f64(&ds[i], ds_v);
        vst1q_f64(&dz[i], dz_v);
        vst1q_f64(&dw[i], dw_v);
    }

    /* Handle remainder */
    for (; i < N; i++) {
        dx[i] = q[i] * (tmp[i] + xi[i] - (z[i] - w[i]) - dxdz[i] + dsdw[i]);
        ds[i] = -dx[i];
        dz[i] = (mu - z[i] * dx[i]) * xinv[i] - z[i] - dxdz[i];
        dw[i] = (mu - w[i] * ds[i]) * sinv[i] - w[i] - dsdw[i];
    }
}
#endif /* CQREG_USE_NEON */

/* Scalar fallback versions (also used when N is small) */
static void corrector_loop1_scalar(
    ST_double *g0_out, ST_double *gfpfd_out,
    ST_double *CQREG_RESTRICT dxdz,
    ST_double *CQREG_RESTRICT dsdw,
    ST_double *CQREG_RESTRICT xinv,
    ST_double *CQREG_RESTRICT sinv,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT x_p,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT s,
    const ST_double *CQREG_RESTRICT dz,
    const ST_double *CQREG_RESTRICT dx,
    const ST_double *CQREG_RESTRICT dw,
    const ST_double *CQREG_RESTRICT ds,
    ST_double fp, ST_double fd, ST_int N)
{
    ST_double g0 = 0.0, gfpfd = 0.0;

    for (ST_int i = 0; i < N; i++) {
        g0 += z[i] * x_p[i] + w[i] * s[i];
        gfpfd += (z[i] + fd * dz[i]) * (x_p[i] + fp * dx[i])
               + (w[i] + fd * dw[i]) * (s[i] + fp * ds[i]);
        dxdz[i] = dx[i] * dz[i];
        dsdw[i] = ds[i] * dw[i];
        xinv[i] = 1.0 / x_p[i];
        sinv[i] = 1.0 / s[i];
    }

    *g0_out = g0;
    *gfpfd_out = gfpfd;
}

static void corrector_loop3_scalar(
    ST_double *CQREG_RESTRICT dx,
    ST_double *CQREG_RESTRICT ds,
    ST_double *CQREG_RESTRICT dz,
    ST_double *CQREG_RESTRICT dw,
    const ST_double *CQREG_RESTRICT q,
    const ST_double *CQREG_RESTRICT tmp,
    const ST_double *CQREG_RESTRICT xi,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT dxdz,
    const ST_double *CQREG_RESTRICT dsdw,
    const ST_double *CQREG_RESTRICT xinv,
    const ST_double *CQREG_RESTRICT sinv,
    ST_double mu, ST_int N)
{
    for (ST_int i = 0; i < N; i++) {
        dx[i] = q[i] * (tmp[i] + xi[i] - (z[i] - w[i]) - dxdz[i] + dsdw[i]);
        ds[i] = -dx[i];
        dz[i] = (mu - z[i] * dx[i]) * xinv[i] - z[i] - dxdz[i];
        dw[i] = (mu - w[i] * ds[i]) * sinv[i] - w[i] - dsdw[i];
    }
}

/* Dispatcher functions that select SIMD or scalar based on availability */
static inline void corrector_loop1(
    ST_double *g0_out, ST_double *gfpfd_out,
    ST_double *dxdz, ST_double *dsdw, ST_double *xinv, ST_double *sinv,
    const ST_double *z, const ST_double *x_p, const ST_double *w, const ST_double *s,
    const ST_double *dz, const ST_double *dx, const ST_double *dw, const ST_double *ds,
    ST_double fp, ST_double fd, ST_int N)
{
#if CQREG_USE_NEON
    if (N >= 64) {
        corrector_loop1_neon(g0_out, gfpfd_out, dxdz, dsdw, xinv, sinv,
                             z, x_p, w, s, dz, dx, dw, ds, fp, fd, N);
        return;
    }
#endif
    corrector_loop1_scalar(g0_out, gfpfd_out, dxdz, dsdw, xinv, sinv,
                           z, x_p, w, s, dz, dx, dw, ds, fp, fd, N);
}

static inline void corrector_loop3(
    ST_double *dx, ST_double *ds, ST_double *dz, ST_double *dw,
    const ST_double *q, const ST_double *tmp, const ST_double *xi,
    const ST_double *z, const ST_double *w,
    const ST_double *dxdz, const ST_double *dsdw,
    const ST_double *xinv, const ST_double *sinv,
    ST_double mu, ST_int N)
{
#if CQREG_USE_NEON
    if (N >= 64) {
        corrector_loop3_neon(dx, ds, dz, dw, q, tmp, xi, z, w,
                             dxdz, dsdw, xinv, sinv, mu, N);
        return;
    }
#endif
    corrector_loop3_scalar(dx, ds, dz, dw, q, tmp, xi, z, w,
                           dxdz, dsdw, xinv, sinv, mu, N);
}

/* ============================================================================
 * SIMD-optimized scaling loop: compute q and r
 * q[i] = 1.0 / (z[i]/x_p[i] + w[i]/s[i])
 * r[i] = z[i] - w[i]
 * ============================================================================ */
#if CQREG_USE_NEON
static void compute_scaling_neon(
    ST_double *CQREG_RESTRICT q,
    ST_double *CQREG_RESTRICT r,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT x_p,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT s,
    ST_int N)
{
    ST_int N2 = N - (N % 2);
    ST_int i;

    for (i = 0; i < N2; i += 2) {
        float64x2_t z_v = vld1q_f64(&z[i]);
        float64x2_t xp_v = vld1q_f64(&x_p[i]);
        float64x2_t w_v = vld1q_f64(&w[i]);
        float64x2_t s_v = vld1q_f64(&s[i]);

        /* q = 1 / (z/x_p + w/s) */
        float64x2_t zx = vdivq_f64(z_v, xp_v);
        float64x2_t ws = vdivq_f64(w_v, s_v);
        float64x2_t sum = vaddq_f64(zx, ws);
        float64x2_t q_v = vdivq_f64(vdupq_n_f64(1.0), sum);

        /* r = z - w */
        float64x2_t r_v = vsubq_f64(z_v, w_v);

        vst1q_f64(&q[i], q_v);
        vst1q_f64(&r[i], r_v);
    }

    for (; i < N; i++) {
        q[i] = 1.0 / (z[i] / x_p[i] + w[i] / s[i]);
        r[i] = z[i] - w[i];
    }
}
#endif

static void compute_scaling_scalar(
    ST_double *CQREG_RESTRICT q,
    ST_double *CQREG_RESTRICT r,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT x_p,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT s,
    ST_int N)
{
    for (ST_int i = 0; i < N; i++) {
        q[i] = 1.0 / (z[i] / x_p[i] + w[i] / s[i]);
        r[i] = z[i] - w[i];
    }
}

static inline void compute_scaling(
    ST_double *q, ST_double *r,
    const ST_double *z, const ST_double *x_p,
    const ST_double *w, const ST_double *s, ST_int N)
{
#if CQREG_USE_NEON
    if (N >= 64) {
        compute_scaling_neon(q, r, z, x_p, w, s, N);
        return;
    }
#endif
    compute_scaling_scalar(q, r, z, x_p, w, s, N);
}

/* ============================================================================
 * SIMD-optimized direction computation (fused loops)
 * Computes dx, ds, dz, dw in a single pass
 * ============================================================================ */
#if CQREG_USE_NEON
static void compute_directions_neon(
    ST_double *CQREG_RESTRICT dx,
    ST_double *CQREG_RESTRICT ds,
    ST_double *CQREG_RESTRICT dz,
    ST_double *CQREG_RESTRICT dw,
    const ST_double *CQREG_RESTRICT q,
    const ST_double *CQREG_RESTRICT tmp,
    const ST_double *CQREG_RESTRICT r,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT x_p,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT s,
    ST_int N)
{
    float64x2_t one = vdupq_n_f64(1.0);
    float64x2_t neg_one = vdupq_n_f64(-1.0);
    ST_int N2 = N - (N % 2);
    ST_int i;

    for (i = 0; i < N2; i += 2) {
        float64x2_t q_v = vld1q_f64(&q[i]);
        float64x2_t tmp_v = vld1q_f64(&tmp[i]);
        float64x2_t r_v = vld1q_f64(&r[i]);
        float64x2_t z_v = vld1q_f64(&z[i]);
        float64x2_t xp_v = vld1q_f64(&x_p[i]);
        float64x2_t w_v = vld1q_f64(&w[i]);
        float64x2_t s_v = vld1q_f64(&s[i]);

        /* dx = q * (tmp - r) */
        float64x2_t dx_v = vmulq_f64(q_v, vsubq_f64(tmp_v, r_v));
        /* ds = -dx */
        float64x2_t ds_v = vmulq_f64(neg_one, dx_v);

        /* dz = -z * (1 + dx/x_p) */
        float64x2_t dz_v = vmulq_f64(vnegq_f64(z_v),
                                      vaddq_f64(one, vdivq_f64(dx_v, xp_v)));
        /* dw = -w * (1 + ds/s) */
        float64x2_t dw_v = vmulq_f64(vnegq_f64(w_v),
                                      vaddq_f64(one, vdivq_f64(ds_v, s_v)));

        vst1q_f64(&dx[i], dx_v);
        vst1q_f64(&ds[i], ds_v);
        vst1q_f64(&dz[i], dz_v);
        vst1q_f64(&dw[i], dw_v);
    }

    for (; i < N; i++) {
        dx[i] = q[i] * (tmp[i] - r[i]);
        ds[i] = -dx[i];
        dz[i] = -z[i] * (1.0 + dx[i] / x_p[i]);
        dw[i] = -w[i] * (1.0 + ds[i] / s[i]);
    }
}
#endif

static void compute_directions_scalar(
    ST_double *CQREG_RESTRICT dx,
    ST_double *CQREG_RESTRICT ds,
    ST_double *CQREG_RESTRICT dz,
    ST_double *CQREG_RESTRICT dw,
    const ST_double *CQREG_RESTRICT q,
    const ST_double *CQREG_RESTRICT tmp,
    const ST_double *CQREG_RESTRICT r,
    const ST_double *CQREG_RESTRICT z,
    const ST_double *CQREG_RESTRICT x_p,
    const ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT s,
    ST_int N)
{
    for (ST_int i = 0; i < N; i++) {
        ST_double dx_i = q[i] * (tmp[i] - r[i]);
        ST_double ds_i = -dx_i;
        dx[i] = dx_i;
        ds[i] = ds_i;
        dz[i] = -z[i] * (1.0 + dx_i / x_p[i]);
        dw[i] = -w[i] * (1.0 + ds_i / s[i]);
    }
}

static inline void compute_directions(
    ST_double *dx, ST_double *ds, ST_double *dz, ST_double *dw,
    const ST_double *q, const ST_double *tmp, const ST_double *r,
    const ST_double *z, const ST_double *x_p,
    const ST_double *w, const ST_double *s, ST_int N)
{
#if CQREG_USE_NEON
    if (N >= 64) {
        compute_directions_neon(dx, ds, dz, dw, q, tmp, r, z, x_p, w, s, N);
        return;
    }
#endif
    compute_directions_scalar(dx, ds, dz, dw, q, tmp, r, z, x_p, w, s, N);
}

/* ============================================================================
 * SIMD-optimized update loop with parallel gap reduction
 * Updates x_p, s, z, w and computes gap = sum(z*x_p + s*w)
 * ============================================================================ */
#if CQREG_USE_NEON
static ST_double update_variables_neon(
    ST_double *CQREG_RESTRICT x_p,
    ST_double *CQREG_RESTRICT s,
    ST_double *CQREG_RESTRICT z,
    ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT dx,
    const ST_double *CQREG_RESTRICT ds,
    const ST_double *CQREG_RESTRICT dz,
    const ST_double *CQREG_RESTRICT dw,
    ST_double fp, ST_double fd, ST_int N)
{
    float64x2_t fp_v = vdupq_n_f64(fp);
    float64x2_t fd_v = vdupq_n_f64(fd);
    float64x2_t gap_v = vdupq_n_f64(0.0);
    ST_int N2 = N - (N % 2);
    ST_int i;

    for (i = 0; i < N2; i += 2) {
        float64x2_t xp_v = vld1q_f64(&x_p[i]);
        float64x2_t s_v = vld1q_f64(&s[i]);
        float64x2_t z_v = vld1q_f64(&z[i]);
        float64x2_t w_v = vld1q_f64(&w[i]);

        float64x2_t dx_v = vld1q_f64(&dx[i]);
        float64x2_t ds_v = vld1q_f64(&ds[i]);
        float64x2_t dz_v = vld1q_f64(&dz[i]);
        float64x2_t dw_v = vld1q_f64(&dw[i]);

        /* Update */
        xp_v = vfmaq_f64(xp_v, fp_v, dx_v);
        s_v = vfmaq_f64(s_v, fp_v, ds_v);
        z_v = vfmaq_f64(z_v, fd_v, dz_v);
        w_v = vfmaq_f64(w_v, fd_v, dw_v);

        vst1q_f64(&x_p[i], xp_v);
        vst1q_f64(&s[i], s_v);
        vst1q_f64(&z[i], z_v);
        vst1q_f64(&w[i], w_v);

        /* Accumulate gap: z*x_p + s*w */
        gap_v = vfmaq_f64(gap_v, z_v, xp_v);
        gap_v = vfmaq_f64(gap_v, s_v, w_v);
    }

    ST_double gap = vgetq_lane_f64(gap_v, 0) + vgetq_lane_f64(gap_v, 1);

    for (; i < N; i++) {
        x_p[i] += fp * dx[i];
        s[i] += fp * ds[i];
        z[i] += fd * dz[i];
        w[i] += fd * dw[i];
        gap += z[i] * x_p[i] + s[i] * w[i];
    }

    return gap;
}
#endif

static ST_double update_variables_scalar(
    ST_double *CQREG_RESTRICT x_p,
    ST_double *CQREG_RESTRICT s,
    ST_double *CQREG_RESTRICT z,
    ST_double *CQREG_RESTRICT w,
    const ST_double *CQREG_RESTRICT dx,
    const ST_double *CQREG_RESTRICT ds,
    const ST_double *CQREG_RESTRICT dz,
    const ST_double *CQREG_RESTRICT dw,
    ST_double fp, ST_double fd, ST_int N)
{
    ST_double gap = 0.0;
    for (ST_int i = 0; i < N; i++) {
        x_p[i] += fp * dx[i];
        s[i] += fp * ds[i];
        z[i] += fd * dz[i];
        w[i] += fd * dw[i];
        gap += z[i] * x_p[i] + s[i] * w[i];
    }
    return gap;
}

static inline ST_double update_variables(
    ST_double *x_p, ST_double *s, ST_double *z, ST_double *w,
    const ST_double *dx, const ST_double *ds,
    const ST_double *dz, const ST_double *dw,
    ST_double fp, ST_double fd, ST_int N)
{
#if CQREG_USE_NEON
    if (N >= 64) {
        return update_variables_neon(x_p, s, z, w, dx, ds, dz, dw, fp, fd, N);
    }
#endif
    return update_variables_scalar(x_p, s, z, w, dx, ds, dz, dw, fp, fd, N);
}

/* ============================================================================
 * Helper: compute primal and dual step lengths using chunked parallel reduction
 *
 * Optimization strategy:
 * 1. OpenMP parallel for with reduction(min:...) handles thread-private mins
 * 2. 4-way unrolling for better instruction pipelining
 * 3. Tree reduction within each 4-element group reduces comparisons
 *
 * ds = -dx always, so we compute ds bound from dx directly.
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
    ST_int N4 = N - (N % 4);

#ifdef _OPENMP
    /*
     * OpenMP parallel for with reduction clause.
     * Key improvement: 4-way unrolling inside the loop for better pipelining.
     * schedule(static) gives each thread a contiguous chunk for cache locality.
     */
    #pragma omp parallel for schedule(static) reduction(min:fp_min, fd_min)
    for (ST_int i = 0; i < N4; i += 4) {
        /* Load and compute bounds for 4 elements */
        ST_double dx_0 = dx[i];
        ST_double bx_0 = (dx_0 < 0.0) ? -x_p[i] / dx_0 : 1e20;
        ST_double bs_0 = (dx_0 > 0.0) ? s[i] / dx_0 : 1e20;
        ST_double bp_0 = (bx_0 < bs_0) ? bx_0 : bs_0;
        ST_double bz_0 = (dz[i] < 0.0) ? -z[i] / dz[i] : 1e20;
        ST_double bw_0 = (dw[i] < 0.0) ? -w[i] / dw[i] : 1e20;
        ST_double bd_0 = (bz_0 < bw_0) ? bz_0 : bw_0;

        ST_double dx_1 = dx[i+1];
        ST_double bx_1 = (dx_1 < 0.0) ? -x_p[i+1] / dx_1 : 1e20;
        ST_double bs_1 = (dx_1 > 0.0) ? s[i+1] / dx_1 : 1e20;
        ST_double bp_1 = (bx_1 < bs_1) ? bx_1 : bs_1;
        ST_double bz_1 = (dz[i+1] < 0.0) ? -z[i+1] / dz[i+1] : 1e20;
        ST_double bw_1 = (dw[i+1] < 0.0) ? -w[i+1] / dw[i+1] : 1e20;
        ST_double bd_1 = (bz_1 < bw_1) ? bz_1 : bw_1;

        ST_double dx_2 = dx[i+2];
        ST_double bx_2 = (dx_2 < 0.0) ? -x_p[i+2] / dx_2 : 1e20;
        ST_double bs_2 = (dx_2 > 0.0) ? s[i+2] / dx_2 : 1e20;
        ST_double bp_2 = (bx_2 < bs_2) ? bx_2 : bs_2;
        ST_double bz_2 = (dz[i+2] < 0.0) ? -z[i+2] / dz[i+2] : 1e20;
        ST_double bw_2 = (dw[i+2] < 0.0) ? -w[i+2] / dw[i+2] : 1e20;
        ST_double bd_2 = (bz_2 < bw_2) ? bz_2 : bw_2;

        ST_double dx_3 = dx[i+3];
        ST_double bx_3 = (dx_3 < 0.0) ? -x_p[i+3] / dx_3 : 1e20;
        ST_double bs_3 = (dx_3 > 0.0) ? s[i+3] / dx_3 : 1e20;
        ST_double bp_3 = (bx_3 < bs_3) ? bx_3 : bs_3;
        ST_double bz_3 = (dz[i+3] < 0.0) ? -z[i+3] / dz[i+3] : 1e20;
        ST_double bw_3 = (dw[i+3] < 0.0) ? -w[i+3] / dw[i+3] : 1e20;
        ST_double bd_3 = (bz_3 < bw_3) ? bz_3 : bw_3;

        /* Tree reduction: find min of 4 in 3 comparisons instead of 4 */
        ST_double bp_01 = (bp_0 < bp_1) ? bp_0 : bp_1;
        ST_double bp_23 = (bp_2 < bp_3) ? bp_2 : bp_3;
        ST_double bp_min4 = (bp_01 < bp_23) ? bp_01 : bp_23;
        if (bp_min4 < fp_min) fp_min = bp_min4;

        ST_double bd_01 = (bd_0 < bd_1) ? bd_0 : bd_1;
        ST_double bd_23 = (bd_2 < bd_3) ? bd_2 : bd_3;
        ST_double bd_min4 = (bd_01 < bd_23) ? bd_01 : bd_23;
        if (bd_min4 < fd_min) fd_min = bd_min4;
    }
#else
    /* Sequential version with 4-way unrolling */
    for (ST_int i = 0; i < N4; i += 4) {
        ST_double dx_0 = dx[i];
        ST_double bx_0 = (dx_0 < 0.0) ? -x_p[i] / dx_0 : 1e20;
        ST_double bs_0 = (dx_0 > 0.0) ? s[i] / dx_0 : 1e20;
        ST_double bp_0 = (bx_0 < bs_0) ? bx_0 : bs_0;
        ST_double bz_0 = (dz[i] < 0.0) ? -z[i] / dz[i] : 1e20;
        ST_double bw_0 = (dw[i] < 0.0) ? -w[i] / dw[i] : 1e20;
        ST_double bd_0 = (bz_0 < bw_0) ? bz_0 : bw_0;

        ST_double dx_1 = dx[i+1];
        ST_double bx_1 = (dx_1 < 0.0) ? -x_p[i+1] / dx_1 : 1e20;
        ST_double bs_1 = (dx_1 > 0.0) ? s[i+1] / dx_1 : 1e20;
        ST_double bp_1 = (bx_1 < bs_1) ? bx_1 : bs_1;
        ST_double bz_1 = (dz[i+1] < 0.0) ? -z[i+1] / dz[i+1] : 1e20;
        ST_double bw_1 = (dw[i+1] < 0.0) ? -w[i+1] / dw[i+1] : 1e20;
        ST_double bd_1 = (bz_1 < bw_1) ? bz_1 : bw_1;

        ST_double dx_2 = dx[i+2];
        ST_double bx_2 = (dx_2 < 0.0) ? -x_p[i+2] / dx_2 : 1e20;
        ST_double bs_2 = (dx_2 > 0.0) ? s[i+2] / dx_2 : 1e20;
        ST_double bp_2 = (bx_2 < bs_2) ? bx_2 : bs_2;
        ST_double bz_2 = (dz[i+2] < 0.0) ? -z[i+2] / dz[i+2] : 1e20;
        ST_double bw_2 = (dw[i+2] < 0.0) ? -w[i+2] / dw[i+2] : 1e20;
        ST_double bd_2 = (bz_2 < bw_2) ? bz_2 : bw_2;

        ST_double dx_3 = dx[i+3];
        ST_double bx_3 = (dx_3 < 0.0) ? -x_p[i+3] / dx_3 : 1e20;
        ST_double bs_3 = (dx_3 > 0.0) ? s[i+3] / dx_3 : 1e20;
        ST_double bp_3 = (bx_3 < bs_3) ? bx_3 : bs_3;
        ST_double bz_3 = (dz[i+3] < 0.0) ? -z[i+3] / dz[i+3] : 1e20;
        ST_double bw_3 = (dw[i+3] < 0.0) ? -w[i+3] / dw[i+3] : 1e20;
        ST_double bd_3 = (bz_3 < bw_3) ? bz_3 : bw_3;

        ST_double bp_01 = (bp_0 < bp_1) ? bp_0 : bp_1;
        ST_double bp_23 = (bp_2 < bp_3) ? bp_2 : bp_3;
        ST_double bp_min4 = (bp_01 < bp_23) ? bp_01 : bp_23;
        if (bp_min4 < fp_min) fp_min = bp_min4;

        ST_double bd_01 = (bd_0 < bd_1) ? bd_0 : bd_1;
        ST_double bd_23 = (bd_2 < bd_3) ? bd_2 : bd_3;
        ST_double bd_min4 = (bd_01 < bd_23) ? bd_01 : bd_23;
        if (bd_min4 < fd_min) fd_min = bd_min4;
    }
#endif

    /* Handle remainder (up to 3 elements) */
    for (ST_int i = N4; i < N; i++) {
        ST_double dx_i = dx[i];
        ST_double bx = (dx_i < 0.0) ? -x_p[i] / dx_i : 1e20;
        ST_double bs = (dx_i > 0.0) ? s[i] / dx_i : 1e20;
        ST_double bp = (bx < bs) ? bx : bs;
        if (bp < fp_min) fp_min = bp;

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

    /* Get tolerance from config (use default if not set) */
    ST_double ipm_tol = ipm->config.tol_gap;
    if (ipm_tol <= 0.0 || ipm_tol > 1.0) {
        ipm_tol = IPM_TOL_DEFAULT;
    }

    for (iter = 0; iter < IPM_MAX_ITER; iter++) {

        /* Check convergence using configurable tolerance */
        if (gap < ipm_tol || !isfinite(gap)) {
            IPM_LOG("Converged at iter %d, gap=%.6e (tol=%.6e)\n", iter, gap, ipm_tol);
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
        /* SIMD-optimized scaling: compute diagonal scaling q and residual r */
        compute_scaling(q, r, z, x_p, w, s, N);
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

        /* SIMD-optimized fused direction computation:
         * dx = q .* (X*dy - r), ds = -dx
         * dz = -z .* (1 + dx./x_p), dw = -w .* (1 + ds./s)
         */
        compute_directions(dx, ds, dz, dw, q, tmp, r, z, x_p, w, s, N);
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
            /* SIMD-optimized loop1: compute g0, gfpfd, and corrector intermediates */
            ST_double g0 = 0.0, gfpfd = 0.0;
            ST_double *dxdz = fx;
            ST_double *dsdw = fs;
            ST_double *xinv = fw;
            /* sinv already points to pre-allocated ipm->fn_sinv */

            corrector_loop1(&g0, &gfpfd, dxdz, dsdw, xinv, sinv,
                           z, x_p, w, s, dz, dx, dw, ds, fp, fd, N);

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

            /* SIMD-optimized loop3: recompute corrected directions */
            corrector_loop3(dx, ds, dz, dw, q, tmp, xi, z, w,
                           dxdz, dsdw, xinv, sinv, mu, N);
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

        /* SIMD-optimized update: update all primal/dual variables and compute gap */
        gap = update_variables(x_p, s, z, w, dx, ds, dz, dw, fp, fd, N);
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
    ipm->converged = (gap < ipm_tol) ? 1 : 0;

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
    fn_timing_report(N, K, ipm->config.verbose);
#endif

    /* No cleanup needed - all arrays are pre-allocated in IPM state */
    return iter;
}
