/*
 * creghdfe_ols.c
 *
 * OLS computations: X'X, X'y, Cholesky decomposition, collinearity detection
 * Highly optimized with OpenMP parallelization
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_ols.h"
#include "../ctools_unroll.h"

/* ========================================================================
 * Compensated dot product using FMA-based error-free transformations.
 * Achieves near-quad-precision accumulation for OLS numerical stability.
 * ======================================================================== */

#include <math.h>

ST_double kahan_dot(const ST_double * RESTRICT x,
                    const ST_double * RESTRICT y,
                    ST_int N)
{
    ST_double s = 0.0;
    ST_double c = 0.0;
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double p = x[i] * y[i];
        ST_double pe = fma(x[i], y[i], -p);
        ST_double t = s + p;
        ST_double v = t - s;
        ST_double te = (s - (t - v)) + (p - v);
        c += te + pe;
        s = t;
    }
    return s + c;
}

/* ========================================================================
 * Double-double (quad-precision) accumulation functions
 * These match Stata's quadcross() precision for RSS/TSS calculations
 * ======================================================================== */

ST_double dd_sum_sq(const ST_double *x, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        dd_real prod = two_prod(x[i], x[i]);
        /* Add both hi and lo parts to accumulator */
        sum = dd_add_d(sum, prod.hi);
        sum = dd_add_d(sum, prod.lo);
    }
    return sum.hi + sum.lo;
}

ST_double dd_sum_sq_weighted(const ST_double *x, const ST_double *w, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        /* Compute x[i]^2 with error tracking */
        dd_real sq = two_prod(x[i], x[i]);
        /* Multiply by weight: w * sq.hi and w * sq.lo */
        dd_real wsq_hi = two_prod(w[i], sq.hi);
        dd_real wsq_lo = two_prod(w[i], sq.lo);
        /* Accumulate all parts */
        sum = dd_add_d(sum, wsq_hi.hi);
        sum = dd_add_d(sum, wsq_hi.lo);
        sum = dd_add_d(sum, wsq_lo.hi);
        sum = dd_add_d(sum, wsq_lo.lo);
    }
    return sum.hi + sum.lo;
}

ST_double dd_sum_sq_dev(const ST_double *x, ST_double mean, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double d = x[i] - mean;
        dd_real prod = two_prod(d, d);
        sum = dd_add_d(sum, prod.hi);
        sum = dd_add_d(sum, prod.lo);
    }
    return sum.hi + sum.lo;
}

ST_double dd_sum_sq_dev_weighted(const ST_double *x, ST_double mean,
                                  const ST_double *w, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double d = x[i] - mean;
        dd_real sq = two_prod(d, d);
        dd_real wsq_hi = two_prod(w[i], sq.hi);
        dd_real wsq_lo = two_prod(w[i], sq.lo);
        sum = dd_add_d(sum, wsq_hi.hi);
        sum = dd_add_d(sum, wsq_hi.lo);
        sum = dd_add_d(sum, wsq_lo.hi);
        sum = dd_add_d(sum, wsq_lo.lo);
    }
    return sum.hi + sum.lo;
}

/* ========================================================================
 * BLAS-like optimized dot product - uses ctools K-way unrolling abstraction
 * ======================================================================== */

ST_double fast_dot(const ST_double * RESTRICT x,
                   const ST_double * RESTRICT y,
                   ST_int N)
{
    return ctools_dot_unrolled(x, y, N);
}

/* ========================================================================
 * Compute X'X (K x K) and X'y (K x 1) - BLAS-like implementation
 * ======================================================================== */

void compute_xtx_xty(
    const ST_double *data,  /* N x K matrix in column-major order */
    ST_int N,
    ST_int K,               /* K includes y as first column */
    ST_double *xtx,         /* Output: (K-1) x (K-1) */
    ST_double *xty          /* Output: (K-1) x 1 */
)
{
    ST_int j, k;
    ST_int K_x = K - 1;  /* Number of X columns (excluding y) */
    const ST_double *y = &data[0];  /* Column 0 is y */

    /* Initialize output */
    memset(xtx, 0, K_x * K_x * sizeof(ST_double));
    memset(xty, 0, K_x * sizeof(ST_double));

    /* Always use Kahan-compensated summation for OLS dot products.
     * reghdfe uses Mata's quadcross() (quad-precision), so we need
     * compensated summation to match its numerical precision. */
    int use_kahan = 1;

#ifdef _OPENMP
    /* Parallel computation of X'y and diagonal of X'X */
    #pragma omp parallel for schedule(static) if(K_x > 2)
#endif
    for (j = 0; j < K_x; j++) {
        const ST_double *xj = &data[(j + 1) * N];

        /* X'y: j-th element */
        if (use_kahan) {
            xty[j] = kahan_dot(xj, y, N);
        } else {
            xty[j] = fast_dot(xj, y, N);
        }

        /* Diagonal of X'X */
        if (use_kahan) {
            xtx[j * K_x + j] = kahan_dot(xj, xj, N);
        } else {
            xtx[j * K_x + j] = fast_dot(xj, xj, N);
        }
    }

    /* Off-diagonal elements (symmetric) */
    for (j = 0; j < K_x; j++) {
        const ST_double *xj = &data[(j + 1) * N];
        for (k = j + 1; k < K_x; k++) {
            const ST_double *xk = &data[(k + 1) * N];
            ST_double val;
            if (use_kahan) {
                val = kahan_dot(xj, xk, N);
            } else {
                val = fast_dot(xj, xk, N);
            }
            xtx[j * K_x + k] = val;
            xtx[k * K_x + j] = val;  /* Symmetric */
        }
    }
}

/* ========================================================================
 * Compute weighted X'WX and X'Wy - BLAS-like implementation
 * W = diag(weights), so X'WX = sum_i (w_i * x_i * x_i')
 * For aweight/pweight: normalizes weights to sum to N (matches reghdfe.mata line 3598, 3609)
 * ======================================================================== */

void compute_xtx_xty_weighted(
    const ST_double *data,     /* N x K matrix in column-major order */
    const ST_double *weights,  /* N x 1 weight vector */
    ST_int weight_type,        /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_int N,
    ST_int K,                  /* K includes y as first column */
    ST_double *xtx,            /* Output: (K-1) x (K-1) */
    ST_double *xty             /* Output: (K-1) x 1 */
)
{
    ST_int i, j, k;
    ST_int K_x = K - 1;  /* Number of X columns (excluding y) */
    const ST_double *y = &data[0];  /* Column 0 is y */

    /* For aweight/pweight: normalize weights to sum to N (reghdfe.mata line 3598) */
    ST_double weight_scale = 1.0;
    if (weight_type == 1 || weight_type == 3) {
        ST_double sum_w = 0.0;
        for (i = 0; i < N; i++) sum_w += weights[i];
        weight_scale = (ST_double)N / sum_w;
    }

    /* Initialize output */
    memset(xtx, 0, K_x * K_x * sizeof(ST_double));
    memset(xty, 0, K_x * sizeof(ST_double));

    /* Compute X'WX and X'Wy where W = diag(weights) */
    for (i = 0; i < N; i++) {
        ST_double w = weights[i];
        /* Normalize for aweight/pweight */
        if (weight_type == 1 || weight_type == 3) {
            w = w * weight_scale;
        }
        ST_double yi = y[i];

        for (j = 0; j < K_x; j++) {
            ST_double xj = data[(j + 1) * N + i];
            ST_double wxj = w * xj;

            /* X'Wy: accumulate w * x_j * y */
            xty[j] += wxj * yi;

            /* Diagonal of X'WX */
            xtx[j * K_x + j] += wxj * xj;

            /* Off-diagonal (upper triangle, then copy to lower) */
            for (k = j + 1; k < K_x; k++) {
                ST_double xk = data[(k + 1) * N + i];
                ST_double val = wxj * xk;
                xtx[j * K_x + k] += val;
                xtx[k * K_x + j] += val;  /* Symmetric */
            }
        }
    }
}

/* Cholesky decomposition and matrix inversion are now in ctools_ols.c
 * The header file defines macros to redirect cholesky and invert_from_cholesky
 * to ctools_cholesky and ctools_invert_from_cholesky */

/* ========================================================================
 * Detect collinear variables using Cholesky decomposition
 * Returns number of collinear variables detected (-1 on error)
 * Algorithm:
 *  1. Attempt Cholesky decomposition
 *  2. Variables with near-zero diagonal are collinear
 * ======================================================================== */

ST_int detect_collinearity(const ST_double *xx, ST_int K, ST_int *is_collinear, ST_int verbose)
{
    ST_int i, j, k, num_collinear;
    ST_double *L, sum, diag_elem;
    const ST_double tol = 1e-14;

    (void)verbose;  /* Unused for now */

    L = (ST_double *)malloc(K * K * sizeof(ST_double));
    if (!L) return -1;

    memcpy(L, xx, K * K * sizeof(ST_double));

    for (i = 0; i < K; i++) {
        is_collinear[i] = 0;
    }

    /* Modified Cholesky that detects collinearity */
    for (k = 0; k < K; k++) {
        /* Compute L[k,k] */
        sum = L[k * K + k];
        for (j = 0; j < k; j++) {
            sum -= L[k * K + j] * L[k * K + j];
        }

        /* Check if diagonal is effectively zero (collinear) */
        if (sum < tol) {
            is_collinear[k] = 1;
            L[k * K + k] = 0.0;
            for (i = k + 1; i < K; i++) {
                L[i * K + k] = 0.0;
            }
            continue;
        }

        diag_elem = sqrt(sum);
        L[k * K + k] = diag_elem;

        for (i = k + 1; i < K; i++) {
            sum = L[i * K + k];
            for (j = 0; j < k; j++) {
                sum -= L[i * K + j] * L[k * K + j];
            }
            L[i * K + k] = sum / diag_elem;
        }
    }

    num_collinear = 0;
    for (i = 0; i < K; i++) {
        if (is_collinear[i]) num_collinear++;
    }

    free(L);
    return num_collinear;
}
