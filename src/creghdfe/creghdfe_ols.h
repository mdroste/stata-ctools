/*
 * creghdfe_ols.h
 *
 * OLS computations: X'X, X'y, Cholesky decomposition, collinearity detection
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_OLS_H
#define CREGHDFE_OLS_H

#include "creghdfe_types.h"
#include "../ctools_ols.h"

/* Backward compatibility macros - use shared implementations */
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/* ========================================================================
 * Double-double arithmetic for quad-precision accumulation
 * Represents a number as hi + lo where |lo| <= 0.5 * ulp(hi)
 * ======================================================================== */

typedef struct {
    ST_double hi;
    ST_double lo;
} dd_real;

/* Error-free transformations */
static inline dd_real two_sum(ST_double a, ST_double b) {
    ST_double s = a + b;
    ST_double v = s - a;
    ST_double e = (a - (s - v)) + (b - v);
    return (dd_real){s, e};
}

static inline dd_real two_prod(ST_double a, ST_double b) {
    ST_double p = a * b;
    ST_double e = fma(a, b, -p);
    return (dd_real){p, e};
}

/* Double-double addition: (a.hi + a.lo) + b */
static inline dd_real dd_add_d(dd_real a, ST_double b) {
    dd_real s = two_sum(a.hi, b);
    s.lo += a.lo;
    /* Renormalize */
    ST_double t = s.hi + s.lo;
    s.lo = s.lo - (t - s.hi);
    s.hi = t;
    return s;
}

/* Double-double sum of squares: sum(x[i]^2) with quad precision */
ST_double dd_sum_sq(const ST_double *x, ST_int N);

/* Double-double weighted sum of squares: sum(w[i] * x[i]^2) */
ST_double dd_sum_sq_weighted(const ST_double *x, const ST_double *w, ST_int N);

/* Double-double sum of squared deviations: sum((x[i] - mean)^2) */
ST_double dd_sum_sq_dev(const ST_double *x, ST_double mean, ST_int N);

/* Double-double weighted sum of squared deviations */
ST_double dd_sum_sq_dev_weighted(const ST_double *x, ST_double mean,
                                  const ST_double *w, ST_int N);

/* ========================================================================
 * Dot product variants for numerical precision
 * ======================================================================== */

/* Kahan-compensated dot product for numerical stability */
ST_double kahan_dot(const ST_double * RESTRICT x,
                    const ST_double * RESTRICT y,
                    ST_int N);

/* BLAS-like optimized dot product with 8-way unrolling */
ST_double fast_dot(const ST_double * RESTRICT x,
                   const ST_double * RESTRICT y,
                   ST_int N);

/* ========================================================================
 * Matrix computations
 * ======================================================================== */

/* Compute X'X (K x K) and X'y (K x 1) */
void compute_xtx_xty(
    const ST_double *data,  /* N x K matrix in column-major order */
    ST_int N,
    ST_int K,               /* K includes y as first column */
    ST_double *xtx,         /* Output: (K-1) x (K-1) */
    ST_double *xty          /* Output: (K-1) x 1 */
);

/* Compute weighted X'WX and X'Wy where W = diag(weights)
 * For aweight/pweight: normalizes weights to sum to N (matches reghdfe.mata line 3598, 3609) */
void compute_xtx_xty_weighted(
    const ST_double *data,     /* N x K matrix in column-major order */
    const ST_double *weights,  /* N x 1 weight vector */
    ST_int weight_type,        /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_int N,
    ST_int K,                  /* K includes y as first column */
    ST_double *xtx,            /* Output: (K-1) x (K-1) */
    ST_double *xty             /* Output: (K-1) x 1 */
);

/* Cholesky decomposition and inversion are provided by ctools_ols.h
 * via the macros defined above */

/* ========================================================================
 * Collinearity detection
 * ======================================================================== */

/* Detect collinear variables using Cholesky decomposition
 * Returns number of collinear variables detected (-1 on error) */
ST_int detect_collinearity(const ST_double *xx, ST_int K, ST_int *is_collinear, ST_int verbose);

#endif /* CREGHDFE_OLS_H */
