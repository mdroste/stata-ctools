/*
 * ctools_ols.h
 *
 * Shared OLS utilities: Cholesky decomposition, matrix inversion,
 * precision arithmetic, dot products, collinearity detection.
 * Used by both creghdfe and civreghdfe.
 * Part of the ctools Stata plugin suite.
 */

#ifndef CTOOLS_OLS_H
#define CTOOLS_OLS_H

#include "stplugin.h"
#include "ctools_config.h"
#include "ctools_types.h"
#include <math.h>

/* Compiler hints for restrict keyword */
#ifndef RESTRICT
#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif
#endif

/* ========================================================================
 * Cholesky decomposition: A = L * L' (in-place, returns L in lower triangle)
 * Dense row-major n x n, leading dimension n; upper triangle is zeroed.
 * Ordinary double accumulation in index order. No scratch allocation.
 * n==0 is a no-op; on failure A can contain a partial factor.
 * Returns STATA_OK or STATA_ERR_INVALID_INPUT / STATA_ERR_SINGULAR
 * ======================================================================== */
ST_int ctools_cholesky(ST_double *A, ST_int n);

/* ========================================================================
 * Matrix inversion via Cholesky (for positive definite matrices)
 * Input: L (lower triangular Cholesky factor)
 * Output: inv (inverse of L*L')
 * Both row-major n x n; L is read-only and must not alias inv. Allocates n*n
 * doubles for the inverse factor. Requires a nonsingular factor from a
 * successful factorization; no new factorization is performed here.
 * Returns STATA_OK or an internal stata_retcode (memory/input/singularity)
 * ======================================================================== */
ST_int ctools_invert_from_cholesky(const ST_double *L, ST_int n, ST_double *inv);

/* ========================================================================
 * Solve linear system Ax = b using Cholesky decomposition
 * A must be symmetric positive definite
 * Row-major n x n; A and b are read-only, x receives n doubles. Allocates a
 * factor copy and n-vector scratch, then factors the copy (unlike the cqreg
 * solve API, which consumes an existing lower factor and overwrites b).
 * Returns STATA_OK or an internal stata_retcode; caller converts at the boundary
 * ======================================================================== */
ST_int ctools_solve_cholesky(const ST_double *A, const ST_double *b, ST_int n, ST_double *x);

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

/* ========================================================================
 * Quad-precision accumulation functions
 * Match Stata's quadcross() precision for RSS/TSS calculations
 * ======================================================================== */

/* Double-double sum of squares: sum(x[i]^2) with quad precision */
ST_double dd_sum_sq(const ST_double *x, ST_int N);

/* Double-double weighted sum of squares: sum(w[i] * x[i]^2) */
ST_double dd_sum_sq_weighted(const ST_double *x, const ST_double *w, ST_int N);

/* Double-double sum of squared deviations: sum((x[i] - mean)^2) */
ST_double dd_sum_sq_dev(const ST_double *x, ST_double mean, ST_int N);

/* Double-double weighted sum of squared deviations */
ST_double dd_sum_sq_dev_weighted(const ST_double *x, ST_double mean,
                                  const ST_double *w, ST_int N);

/* Double-double weighted dot product: sum((w[i] * x[i]) * y[i]), with w == NULL
 * meaning unit weights. The product w[i] * x[i] is rounded to double first;
 * accumulation is sequential in index order. Returned unnormalized so a
 * caller can keep accumulating; the value is hi + lo. */
dd_real dd_dot_weighted(const ST_double *x, const ST_double *y, const ST_double *w, ST_int N);

/* The same sums for n right-hand vectors at once: out[t] is exactly
 * dd_dot_weighted(x, ys[t], w, N). The accumulations are independent, so up
 * to four are interleaved per pass over the data, which several times
 * shortens the wait on each sum's dependency chain. */
void dd_dot_weighted_multi(const ST_double *x, const ST_double *const *ys, ST_int n,
                           const ST_double *w, ST_int N, dd_real *out);

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

/* ========================================================================
 * Collinearity detection
 * ======================================================================== */

/* Detect collinear variables using Cholesky decomposition
 * Returns number of collinear variables detected (-1 on error) */
ST_int detect_collinearity(const ST_double *xx, ST_int K, ST_int *is_collinear, ST_int verbose);

/* Scale-invariant collinearity detection on a symmetric cross-product matrix.
 * xx (K x K) is modified in place: rows and columns of the variables flagged
 * in is_collinear on entry are zeroed, the matrix is equilibrated to a unit
 * diagonal, and detect_collinearity() runs on the result. Variables flagged
 * on entry stay flagged. Returns what detect_collinearity() returns for the
 * equilibrated matrix (-1 on allocation failure, is_collinear then unchanged). */
ST_int detect_collinearity_scaled(ST_double *xx, ST_int K, ST_int *is_collinear, ST_int verbose);

/* Inverse cross-product matrix of [X, 1] from that of mean-centred X.
 * inv is (X~'W X~)^-1 for the K centred regressors (K x K, symmetric), means
 * their K (weighted) means and sum_w the total weight 1'W1 (N if unweighted).
 * out is (K+1) x (K+1) with the constant last, by the block-partition
 * formula; it must not overlap inv. K == 0 gives the 1 x 1 matrix 1/sum_w. */
void ctools_inverse_add_constant(const ST_double *inv, ST_int K, const ST_double *means,
                                 ST_double sum_w, ST_double *out);

/* Sweep-based generalized inverse for symmetric matrices (matches Mata invsym):
 * pivots on the largest remaining diagonal and zeros a column whose remaining
 * diagonal is not positive, below 1e-9 of its original diagonal, or below 1e-19.
 * Returns 0 on success, -1 on error. */
ST_int ctools_invsym(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out);

#endif /* CTOOLS_OLS_H */
