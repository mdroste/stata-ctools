#include "ctools_linalg.h"
/*
 * ctools_ols.c
 *
 * Shared OLS utilities: Cholesky decomposition, matrix inversion,
 * precision arithmetic, dot products, collinearity detection.
 * Used by both creghdfe and civreghdfe.
 * Part of the ctools Stata plugin suite.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ctools_ols.h"
#include "ctools_config.h"
#include "ctools_unroll.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* ========================================================================
 * Cholesky decomposition: A = L * L' (in-place, returns L in lower triangle)
 * ======================================================================== */

ST_int ctools_cholesky(ST_double *A, ST_int n)
{
    if (n == 0) return STATA_OK;
    if (!A || n < 0) return STATA_ERR_INVALID_INPUT;
    ST_int i, j, k;

    for (j = 0; j < n; j++) {
        ST_double sum = A[(size_t)j * n + j];
        for (k = 0; k < j; k++) {
            sum -= A[(size_t)j * n + k] * A[(size_t)j * n + k];
        }
        if (!isfinite(sum) || sum <= 0.0) return STATA_ERR_SINGULAR;  /* Not positive definite */
        A[(size_t)j * n + j] = sqrt(sum);

        for (i = j + 1; i < n; i++) {
            sum = A[(size_t)i * n + j];
            for (k = 0; k < j; k++) {
                sum -= A[(size_t)i * n + k] * A[(size_t)j * n + k];
            }
            A[(size_t)i * n + j] = sum / A[(size_t)j * n + j];
        }
    }

    /* Zero upper triangle */
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            A[(size_t)i * n + j] = 0.0;
        }
    }

    return 0;
}

/* ========================================================================
 * Matrix inversion via Cholesky (for positive definite matrices)
 * Input: L (lower triangular Cholesky factor)
 * Output: inv (inverse of L*L')
 * ======================================================================== */

ST_int ctools_invert_from_cholesky(const ST_double *L, ST_int n, ST_double *inv)
{
    if (n == 0) return STATA_OK;
    if (!L || !inv || n < 0) return STATA_ERR_INVALID_INPUT;
    ST_int i, j, k;
    ST_double *L_inv = (ST_double *)ctools_safe_calloc3((size_t)n, (size_t)n, sizeof(ST_double));

    if (!L_inv) return STATA_ERR_MEMORY;

    /* Invert L (lower triangular) */
    for (i = 0; i < n; i++) {
        L_inv[(size_t)i * n + i] = 1.0 / L[(size_t)i * n + i];
        for (j = i + 1; j < n; j++) {
            ST_double sum = 0.0;
            for (k = i; k < j; k++) {
                sum -= L[(size_t)j * n + k] * L_inv[(size_t)k * n + i];
            }
            L_inv[(size_t)j * n + i] = sum / L[(size_t)j * n + j];
        }
    }

    /* Compute inv = L_inv' * L_inv */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            inv[(size_t)i * n + j] = 0.0;
            for (k = (i > j ? i : j); k < n; k++) {
                inv[(size_t)i * n + j] += L_inv[(size_t)k * n + i] * L_inv[(size_t)k * n + j];
            }
        }
    }

    free(L_inv);
    return 0;
}

/* ========================================================================
 * Solve linear system Ax = b using Cholesky decomposition
 * ======================================================================== */

ST_int ctools_solve_cholesky(const ST_double *A, const ST_double *b, ST_int n, ST_double *x)
{
    if (n == 0) return STATA_OK;
    if (!A || !b || !x || n < 0) return STATA_ERR_INVALID_INPUT;
    ST_double *L = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    ST_double *y = (ST_double *)ctools_safe_malloc2((size_t)n, sizeof(ST_double));

    if (!L || !y) {
        free(L);
        free(y);
        return STATA_ERR_MEMORY;
    }

    memcpy(L, A, (size_t)n * (size_t)n * sizeof(ST_double));

    stata_retcode status = ctools_cholesky(L, n);
    if (status != STATA_OK) {
        free(L);
        free(y);
        return status;
    }

    ctools_solve_lower_factor(L, b, n, y);
    ctools_solve_lower_transpose(L, y, n, x);

    free(L);
    free(y);
    return 0;
}

/* ========================================================================
 * Compensated dot product using FMA-based error-free transformations.
 * Achieves near-quad-precision accumulation for OLS numerical stability.
 * ======================================================================== */

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

dd_real dd_dot_weighted(const ST_double *x, const ST_double *y, const ST_double *w, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double wx = w ? w[i] * x[i] : x[i];
        dd_real prod = two_prod(wx, y[i]);
        sum = dd_add_d(sum, prod.hi);
        sum = dd_add_d(sum, prod.lo);
    }
    return sum;
}

void dd_dot_weighted_multi(const ST_double *x, const ST_double *const *ys, ST_int n,
                           const ST_double *w, ST_int N, dd_real *out)
{
    ST_int t = 0;

    for (; t + 3 < n; t += 4) {
        const ST_double *y0 = ys[t], *y1 = ys[t + 1], *y2 = ys[t + 2], *y3 = ys[t + 3];
        dd_real s0 = {0.0, 0.0}, s1 = {0.0, 0.0}, s2 = {0.0, 0.0}, s3 = {0.0, 0.0};
        for (ST_int i = 0; i < N; i++) {
            ST_double wx = w ? w[i] * x[i] : x[i];
            dd_real p0 = two_prod(wx, y0[i]);
            dd_real p1 = two_prod(wx, y1[i]);
            dd_real p2 = two_prod(wx, y2[i]);
            dd_real p3 = two_prod(wx, y3[i]);
            s0 = dd_add_d(dd_add_d(s0, p0.hi), p0.lo);
            s1 = dd_add_d(dd_add_d(s1, p1.hi), p1.lo);
            s2 = dd_add_d(dd_add_d(s2, p2.hi), p2.lo);
            s3 = dd_add_d(dd_add_d(s3, p3.hi), p3.lo);
        }
        out[t] = s0;
        out[t + 1] = s1;
        out[t + 2] = s2;
        out[t + 3] = s3;
    }
    for (; t < n; t++) {
        out[t] = dd_dot_weighted(x, ys[t], w, N);
    }
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
    memset(xtx, 0, (size_t)K_x * K_x * sizeof(ST_double));
    memset(xty, 0, K_x * sizeof(ST_double));

    /* Kahan-compensated dot products for numerical stability.
     * reghdfe uses Mata's quadcross() (quad-precision), so we need
     * compensated summation to match after FE partialling out.
     * BLAS dsyrk would be faster but lacks the precision needed here. */

#ifdef _OPENMP
    /* Parallel computation of X'y and diagonal of X'X */
    #pragma omp parallel for schedule(static) if(K_x > 2)
#endif
    for (j = 0; j < K_x; j++) {
        const ST_double *xj = &data[(size_t)(j + 1) * N];
        xty[j] = kahan_dot(xj, y, N);
        xtx[(size_t)j * K_x + j] = kahan_dot(xj, xj, N);
    }

    /* Off-diagonal elements (symmetric, trivially parallelizable) */
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(K_x > 4) private(k)
#endif
    for (j = 0; j < K_x; j++) {
        const ST_double *xj = &data[(size_t)(j + 1) * N];
        for (k = j + 1; k < K_x; k++) {
            const ST_double *xk = &data[(size_t)(k + 1) * N];
            ST_double val = kahan_dot(xj, xk, N);
            xtx[(size_t)j * K_x + k] = val;
            xtx[(size_t)k * K_x + j] = val;  /* Symmetric */
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
    memset(xtx, 0, (size_t)K_x * K_x * sizeof(ST_double));
    memset(xty, 0, K_x * sizeof(ST_double));

    /* Pre-normalize weights if needed to avoid branch per observation */
    ST_double *w_scaled = NULL;
    if (weight_type == 1 || weight_type == 3) {
        w_scaled = (ST_double *)malloc(N * sizeof(ST_double));
        if (w_scaled) {
            for (i = 0; i < N; i++) w_scaled[i] = weights[i] * weight_scale;
        }
    }
    const ST_double *w_eff = w_scaled ? w_scaled : weights;

    /* xty[j] = sum_i (w_i x_ij) y_i and xtx[j,k] = sum_i (w_i x_ij) x_ik for
     * k >= j, every element accumulated over the observations in order with
     * plain double arithmetic. The elements are independent sums, so the rows
     * of the upper triangle run in parallel, four right-hand columns per
     * pass over the data. */
    #pragma omp parallel for schedule(dynamic) private(i, k) if(K_x > 1 && (size_t)N * (size_t)K_x > 100000)
    for (j = 0; j < K_x; j++) {
        const ST_double *xj = data + (size_t)(j + 1) * N;
        /* right-hand columns: y, then x_j .. x_{K_x-1} */
        const ST_int n_rhs = K_x - j + 1;
        for (ST_int t = 0; t < n_rhs; t += 4) {
            const ST_double *r[4];
            ST_double acc[4] = {0.0, 0.0, 0.0, 0.0};
            const ST_int n_now = (n_rhs - t < 4) ? n_rhs - t : 4;
            for (ST_int q = 0; q < 4; q++) {
                ST_int rhs = t + (q < n_now ? q : 0);
                r[q] = (rhs == 0) ? y : data + (size_t)(j + rhs) * N;
            }
            for (i = 0; i < N; i++) {
                ST_double wxj = w_eff[i] * xj[i];
                acc[0] += wxj * r[0][i];
                acc[1] += wxj * r[1][i];
                acc[2] += wxj * r[2][i];
                acc[3] += wxj * r[3][i];
            }
            for (ST_int q = 0; q < n_now; q++) {
                ST_int rhs = t + q;
                if (rhs == 0) xty[j] = acc[q];
                else xtx[(size_t)j * K_x + (j + rhs - 1)] = acc[q];
            }
        }
    }

    /* Mirror upper triangle to lower */
    for (j = 0; j < K_x; j++) {
        for (k = j + 1; k < K_x; k++) {
            xtx[(size_t)k * K_x + j] = xtx[(size_t)j * K_x + k];
        }
    }

    if (w_scaled) free(w_scaled);
}

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
    const ST_double abs_tol = 1e-14;
    const ST_double rel_tol = 1e-10;

    (void)verbose;  /* Unused for now */

    L = (ST_double *)ctools_safe_malloc3((size_t)K, (size_t)K, sizeof(ST_double));
    if (!L) return -1;

    memcpy(L, xx, (size_t)K * K * sizeof(ST_double));

    for (i = 0; i < K; i++) {
        is_collinear[i] = 0;
    }

    /* Modified Cholesky that detects collinearity */
    for (k = 0; k < K; k++) {
        /* Save original diagonal before subtraction for relative check */
        ST_double orig_diag = L[(size_t)k * K + k];

        /* Compute L[k,k] */
        sum = L[(size_t)k * K + k];
        for (j = 0; j < k; j++) {
            sum -= L[(size_t)k * K + j] * L[(size_t)k * K + j];
        }

        /* Check if diagonal is effectively zero (collinear)
         * Use both absolute and relative tolerance to handle varying scales.
         * Absolute: catches zero-variance columns
         * Relative: catches collinearity where floating-point cancellation
         *   leaves a residual of O(N * eps * diag), which can exceed abs_tol */
        if (sum < abs_tol || (orig_diag > 0 && sum / orig_diag < rel_tol)) {
            is_collinear[k] = 1;
            L[(size_t)k * K + k] = 0.0;
            for (i = k + 1; i < K; i++) {
                L[(size_t)i * K + k] = 0.0;
            }
            continue;
        }

        diag_elem = sqrt(sum);
        L[(size_t)k * K + k] = diag_elem;

        for (i = k + 1; i < K; i++) {
            sum = L[(size_t)i * K + k];
            for (j = 0; j < k; j++) {
                sum -= L[(size_t)i * K + j] * L[(size_t)k * K + j];
            }
            L[(size_t)i * K + k] = sum / diag_elem;
        }
    }

    num_collinear = 0;
    for (i = 0; i < K; i++) {
        if (is_collinear[i]) num_collinear++;
    }

    free(L);
    return num_collinear;
}

/* ========================================================================
 * Scale-invariant collinearity detection (see ctools_ols.h)
 * ======================================================================== */
ST_int detect_collinearity_scaled(ST_double *xx, ST_int K, ST_int *is_collinear, ST_int verbose)
{
    ST_int i, j, num_collinear;
    ST_int *excluded = (ST_int *)malloc((size_t)(K > 0 ? K : 1) * sizeof(ST_int));

    if (!excluded) return -1;
    memcpy(excluded, is_collinear, (size_t)K * sizeof(ST_int));

    /* Excluded variables get zero rows and columns, so the equilibration
     * leaves them at zero and the factorization flags them without involving
     * the other columns */
    for (i = 0; i < K; i++) {
        if (!excluded[i]) continue;
        for (j = 0; j < K; j++) {
            xx[(size_t)j * K + i] = 0.0;
            xx[(size_t)i * K + j] = 0.0;
        }
    }

    /* Equilibrate to a correlation matrix, which keeps mixed-scale
     * variables from producing false positives */
    for (i = 0; i < K; i++) {
        ST_double d = xx[(size_t)i * K + i];
        if (d > 0.0) {
            ST_double scale = 1.0 / sqrt(d);
            for (j = 0; j < K; j++) {
                xx[(size_t)j * K + i] *= scale;
                xx[(size_t)i * K + j] *= scale;
            }
        }
    }

    /* detect_collinearity() resets the flags; on failure they are restored */
    num_collinear = detect_collinearity(xx, K, is_collinear, verbose);
    if (num_collinear < 0) {
        memcpy(is_collinear, excluded, (size_t)K * sizeof(ST_int));
    } else {
        for (i = 0; i < K; i++) {
            if (excluded[i]) is_collinear[i] = 1;
        }
    }
    free(excluded);
    return num_collinear;
}

/* ========================================================================
 * Inverse cross-product matrix with a constant (see ctools_ols.h)
 * ======================================================================== */
void ctools_inverse_add_constant(const ST_double *inv, ST_int K, const ST_double *means,
                                 ST_double sum_w, ST_double *out)
{
    ST_int i, j;
    const ST_int Kc = K + 1;

    /* side = -inv * means, kept in the last column of out */
    for (j = 0; j < K; j++) {
        ST_double side = 0.0;
        for (i = 0; i < K; i++) {
            side -= means[i] * inv[(size_t)i * K + j];
        }
        out[(size_t)j * Kc + K] = side;
    }

    ST_double corner = 1.0 / sum_w;
    for (i = 0; i < K; i++) {
        corner -= means[i] * out[(size_t)i * Kc + K];
    }

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            out[(size_t)i * Kc + j] = inv[(size_t)i * K + j];
        }
        out[(size_t)K * Kc + i] = out[(size_t)i * Kc + K];
    }
    out[(size_t)K * Kc + K] = corner;
}

/* ========================================================================
 * Sweep-based generalized inverse for symmetric matrices.
 * Matches Mata's invsym() (checked in Stata): the next pivot is the column
 * with the largest remaining diagonal (lowest index on ties); a column whose
 * remaining diagonal is not positive, has fallen below 1e-9 of its own
 * original diagonal, or is below 1e-19 in absolute value is dependent and
 * zeroed, and the sweep continues with the other columns.
 * Returns 0 on success, -1 on memory error.
 * Output: inv (K x K), rank_out (number of non-singular pivots)
 * ======================================================================== */
ST_int ctools_invsym(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out)
{
    if (K == 0) { if (rank_out) *rank_out = 0; return STATA_OK; }
    if (!A || !inv || K < 0) return STATA_ERR_INVALID_INPUT;
    ST_int i, j;

    ST_double *W = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));
    /* 0 = not yet pivoted, 1 = swept, 2 = dependent (zeroed) */
    ST_int *swept = (ST_int *)calloc(K, sizeof(ST_int));
    if (!W || !swept) { free(W); free(swept); return STATA_ERR_MEMORY; }

    memcpy(W, A, (size_t)K * K * sizeof(ST_double));

    ST_int rank = 0;
    for (;;) {
        /* Largest remaining diagonal among the columns not yet pivoted */
        ST_int p = -1;
        for (i = 0; i < K; i++) {
            if (!swept[i] && (p < 0 || W[(size_t)i * K + i] > W[(size_t)p * K + p])) p = i;
        }
        if (p < 0) break;

        ST_double d = W[(size_t)p * K + p];
        ST_double d0 = A[(size_t)p * K + p];
        if (!(d > 0.0) || !(d0 > 0.0) || d < 1e-9 * d0 || d < 1e-19) {
            swept[p] = 2;
            continue;
        }

        /* Sweep on column/row p (Gauss-Jordan pivot for symmetric matrix).
         * L is small (number of excluded instruments, typically 2-10),
         * so parallelization would be pure overhead — do not add OMP here. */
        ST_double inv_d = 1.0 / W[(size_t)p * K + p];

        /* Update off-diagonal block: W[i,j] -= W[i,p] * W[p,j] / d
         * W[.,p] and W[p,.] are read-only during this loop (only W[i!=p, j!=p] written),
         * so we precompute scaled pivot column to avoid repeated division. */
        for (i = 0; i < K; i++) {
            if (i == p) continue;
            ST_double wip = W[(size_t)p * K + i] * inv_d;  /* W[i,p] / d, constant for inner loop */
            for (j = i; j < K; j++) {
                if (j == p) continue;
                W[(size_t)j * K + i] -= wip * W[(size_t)j * K + p];
                if (i != j) W[(size_t)i * K + j] = W[(size_t)j * K + i];
            }
        }

        /* Update pivot row/col: W[i,p] /= d, W[p,i] /= d */
        for (i = 0; i < K; i++) {
            if (i != p) {
                W[(size_t)p * K + i] *= inv_d;
                W[(size_t)i * K + p] *= inv_d;
            }
        }

        /* Pivot diagonal: the symmetric sweep keeps the swept block as
         * -inverse (the later updates use W[p,p] with that sign), so the
         * result is negated below. */
        W[(size_t)p * K + p] = -inv_d;
        swept[p] = 1;
        rank++;
    }

    /* Copy result: swept positions have -W (the inverse of the swept
     * submatrix, as invsym gives), zero out dependent columns */
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            if (swept[i] == 1 && swept[j] == 1) {
                inv[(size_t)j * K + i] = -W[(size_t)j * K + i];
            } else {
                inv[(size_t)j * K + i] = 0.0;
            }
        }
    }

    if (rank_out) *rank_out = rank;
    free(W);
    free(swept);
    return 0;
}

/* ========================================================================
 * Pseudo-inverse for symmetric positive semidefinite matrices
 * Uses Jacobi eigenvalue decomposition (exact for small matrices).
 * Eigenvalues below rel_tol * max_eigenvalue are treated as zero.
 * Returns 0 on success, -1 on error.
 * Output: inv (K x K), rank (effective rank of matrix)
 * ======================================================================== */
