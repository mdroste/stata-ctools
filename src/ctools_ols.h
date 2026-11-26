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
 * Both row-major n x n; L is read-only. L may be the same array as inv: L is
 * read completely into the inverse factor before inv is written (callers
 * invert in place). Allocates n*n doubles for the inverse factor. Requires a nonsingular factor from a
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

/* The exact rounding error a*b - p of p = fl(a*b). With hardware FMA this
 * is one fused instruction. The baseline x86-64 targets have no FMA, and
 * there fma() is a library call; MinGW-w64's software fma() is not
 * correctly rounded either, so the "error-free" product was inexact on
 * Windows. Dekker's split product is exact (and so bitwise equal to a
 * correctly rounded fma) whenever no split overflows and the error does not
 * underflow; the rare operands outside that range use fma(). */
#if defined(__FMA__) || defined(__aarch64__) || defined(_M_ARM64)
#define CTOOLS_PROD_ERR(a, b, p) __builtin_fma((a), (b), -(p))
#else
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
/* Operands outside the split's range: the FMA instruction where the CPU has
 * one (every x86-64 CPU since 2013), so that they too are exact. */
__attribute__((target("fma"), noinline, unused))
static ST_double ctools_fma_insn(ST_double a, ST_double b, ST_double c) {
    return __builtin_fma(a, b, c);
}
#define CTOOLS_FMA_RARE(a, b, c) \
    (__builtin_cpu_supports("fma") ? ctools_fma_insn((a), (b), (c)) : fma((a), (b), (c)))
#else
#define CTOOLS_FMA_RARE(a, b, c) fma((a), (b), (c))
#endif
static inline ST_double ctools_prod_err(ST_double a, ST_double b, ST_double p) {
    if (fabs(a) < 0x1p995 && fabs(b) < 0x1p995 && fabs(p) > 0x1p-968) {
        const ST_double split = 134217729.0; /* 2^27 + 1 */
        ST_double ca = split * a, cb = split * b;
        ST_double ah = ca - (ca - a), al = a - ah;
        ST_double bh = cb - (cb - b), bl = b - bh;
        return ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    }
    return (p == 0.0 && (a == 0.0 || b == 0.0)) ? 0.0 : CTOOLS_FMA_RARE(a, b, -p);
}
#define CTOOLS_PROD_ERR(a, b, p) ctools_prod_err((a), (b), (p))
#endif

static inline dd_real two_prod(ST_double a, ST_double b) {
    ST_double p = a * b;
    ST_double e = CTOOLS_PROD_ERR(a, b, p);
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
 * Kernels of error-free products, compiled twice on baseline x86-64
 *
 * The split product above costs seventeen operations per element, the FMA
 * instruction two. A kernel that takes "fused" as a constant first argument
 * and forms its products with two_prod_as() is compiled once with the split
 * and once for the FMA instruction, and the copy the CPU can run is chosen
 * at run time. Both copies compute the exact error of every product (the
 * split is exact in its range and uses the instruction outside it), so the
 * choice changes no result. This relies, like the rest of the code, on
 * -ffp-contract=off: no other multiply-add may fuse.
 * ======================================================================== */
#if !defined(__FMA__) && (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
#define CTOOLS_FMA_DISPATCH 1
#if defined(__clang__)
#define CTOOLS_FMA_TARGET __attribute__((target("fma")))
#else
/* 128-bit vectors at most: the FMA instruction brings AVX with it, and GCC
 * does not keep 256-bit spills aligned on the Windows x64 stack */
#define CTOOLS_FMA_TARGET __attribute__((target("fma,prefer-vector-width=128")))
#endif
#define CTOOLS_CPU_HAS_FMA() (__builtin_cpu_supports("fma") != 0)
#else
#define CTOOLS_FMA_DISPATCH 0
#define CTOOLS_FMA_TARGET
#define CTOOLS_CPU_HAS_FMA() 0
#endif

/* two_prod() inside such a kernel. fused may be nonzero only in a function
 * declared CTOOLS_FMA_TARGET. */
static CTOOLS_INLINE dd_real two_prod_as(const int fused, ST_double a, ST_double b) {
    ST_double p = a * b;
    ST_double e = fused ? __builtin_fma(a, b, -p) : CTOOLS_PROD_ERR(a, b, p);
    return (dd_real){p, e};
}

/* Define ret name params from static CTOOLS_INLINE ret name##_body(const int
 * fused, ...): the arguments after the parameter list are those passed on. */
#if CTOOLS_FMA_DISPATCH
#define CTOOLS_FMA_VARIANTS(linkage, ret, name, params, ...)                                       \
    static ret name##_split params { return name##_body(0, __VA_ARGS__); }                         \
    CTOOLS_FMA_TARGET static ret name##_fused params { return name##_body(1, __VA_ARGS__); }       \
    linkage ret name params                                                                        \
    {                                                                                              \
        return CTOOLS_CPU_HAS_FMA() ? name##_fused(__VA_ARGS__) : name##_split(__VA_ARGS__);       \
    }
#define CTOOLS_FMA_VARIANTS_VOID(linkage, name, params, ...)                                       \
    static void name##_split params { name##_body(0, __VA_ARGS__); }                               \
    CTOOLS_FMA_TARGET static void name##_fused params { name##_body(1, __VA_ARGS__); }             \
    linkage void name params                                                                       \
    {                                                                                              \
        if (CTOOLS_CPU_HAS_FMA()) name##_fused(__VA_ARGS__);                                       \
        else name##_split(__VA_ARGS__);                                                            \
    }
#else
#define CTOOLS_FMA_VARIANTS(linkage, ret, name, params, ...)                                       \
    linkage ret name params { return name##_body(0, __VA_ARGS__); }
#define CTOOLS_FMA_VARIANTS_VOID(linkage, name, params, ...)                                       \
    linkage void name params { name##_body(0, __VA_ARGS__); }
#endif

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

/* The same sums over fixed blocks of rows (each in row order, then the
 * blocks in order), computed in parallel: equal to dd_dot_weighted_multi to
 * double-double accuracy and independent of the thread count. Returns 0, or
 * 920 on allocation failure. */
int dd_dot_weighted_multi_blocked(const ST_double *x, const ST_double *const *ys, ST_int n,
                                  const ST_double *w, ST_int N, dd_real *out);

/* ========================================================================
 * Dot product variants for numerical precision
 * ======================================================================== */

/* Kahan-compensated dot product for numerical stability: the rows in order,
 * each product with its exact error */
ST_double kahan_dot(const ST_double * RESTRICT x,
                    const ST_double * RESTRICT y,
                    ST_int N);

/* Compensated cross products of the K columns of X (N x K, column-major):
 * xtx (K x K, both triangles) receives X'X and, when y is not NULL, xty (K)
 * receives X'y. Up to CTOOLS_KAHAN_BLOCK rows every element is exactly
 * kahan_dot() of its two columns. Beyond that each block of
 * CTOOLS_KAHAN_BLOCK rows is summed as kahan_dot() does and the blocks are
 * added in order in double-double, which equals kahan_dot() to its own
 * accuracy; the blocks are fixed by N, so an element is the same at every
 * thread count and depends only on its two columns. The blocks run in
 * parallel, all pairs of columns in one pass over the rows. No failure mode:
 * a scratch allocation that fails selects a slower arrangement of the same
 * sums. */
#define CTOOLS_KAHAN_BLOCK 32768
void ctools_kahan_cross(const ST_double *X, const ST_double *y, ST_int N, ST_int K,
                        ST_double *xtx, ST_double *xty);

/* The same sums for X' diag(w) X: element (j, k) adds (w[i] * X[i,j]) * X[i,k],
 * the weighted left factor rounded to double first and that product taken
 * with its exact error. Only the upper triangle upper[j * K + k], k >= j, is
 * written, each element as the long double sum of its two parts (the running
 * sum and its compensation for one block, the double-double total for more),
 * for callers that go on in extended precision. */
void ctools_kahan_cross_weighted(const ST_double *X, const ST_double *w, ST_int N, ST_int K,
                                 long double *upper);

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

/* ctools_invsym with the columns swept in the order 0..K-1, as Mata's
 * invsym(A, 1..K); a column that fails the tolerance when its turn comes is
 * zeroed. */
ST_int ctools_invsym_ordered(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out);

#endif /* CTOOLS_OLS_H */
