/*
 * creghdfe_ols.h
 *
 * OLS computations: X'X, X'y, Cholesky decomposition, collinearity detection
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_OLS_H
#define CREGHDFE_OLS_H

#include "creghdfe_types.h"

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

/* Compute weighted X'WX and X'Wy where W = diag(weights) */
void compute_xtx_xty_weighted(
    const ST_double *data,     /* N x K matrix in column-major order */
    const ST_double *weights,  /* N x 1 weight vector */
    ST_int N,
    ST_int K,                  /* K includes y as first column */
    ST_double *xtx,            /* Output: (K-1) x (K-1) */
    ST_double *xty             /* Output: (K-1) x 1 */
);

/* Cholesky decomposition: A = L * L' (in-place, returns L in lower triangle) */
ST_int cholesky(ST_double *A, ST_int n);

/* Matrix inversion via Cholesky (for positive definite matrices) */
ST_int invert_from_cholesky(const ST_double *L, ST_int n, ST_double *inv);

/* ========================================================================
 * Collinearity detection
 * ======================================================================== */

/* Detect collinear variables using Cholesky decomposition
 * Returns number of collinear variables detected (-1 on error) */
ST_int detect_collinearity(const ST_double *xx, ST_int K, ST_int *is_collinear, ST_int verbose);

#endif /* CREGHDFE_OLS_H */
