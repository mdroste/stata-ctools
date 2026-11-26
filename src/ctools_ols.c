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
 * Compensated (Kahan) products: every product enters with its exact error.
 * Near-quad-precision accumulation for OLS numerical stability.
 * ======================================================================== */

/* (s, c) += x * y: the running sum in s; in c the rounding errors of the
 * additions and of the products */
#define KAHAN_STEP(s, c, x, y, fused)                                                              \
    do {                                                                                           \
        const ST_double p_ = (x) * (y);                                                            \
        const ST_double pe_ =                                                                      \
            (fused) ? __builtin_fma((x), (y), -p_) : CTOOLS_PROD_ERR((x), (y), p_);                \
        const ST_double t_ = (s) + p_;                                                             \
        const ST_double v_ = t_ - (s);                                                             \
        const ST_double te_ = ((s) - (t_ - v_)) + (p_ - v_);                                       \
        (c) += te_ + pe_;                                                                          \
        (s) = t_;                                                                                  \
    } while (0)

/* n rows of x against nr (1 to 4) right-hand columns. The sums advance
 * together, which keeps up to four dependency chains in flight instead of
 * one; each is the sequential sum over the rows. sc holds (s, c) per
 * right-hand column and is continued, not reset. */
static CTOOLS_INLINE void kahan_rows_body(const int fused, const int nr, const ST_double *x,
                                          const ST_double *const *r, ST_int n, ST_double *sc)
{
    const ST_double *r0 = r[0], *r1 = r[nr > 1 ? 1 : 0];
    const ST_double *r2 = r[nr > 2 ? 2 : 0], *r3 = r[nr > 3 ? 3 : 0];
    ST_double s0 = sc[0], c0 = sc[1];
    ST_double s1 = nr > 1 ? sc[2] : 0.0, c1 = nr > 1 ? sc[3] : 0.0;
    ST_double s2 = nr > 2 ? sc[4] : 0.0, c2 = nr > 2 ? sc[5] : 0.0;
    ST_double s3 = nr > 3 ? sc[6] : 0.0, c3 = nr > 3 ? sc[7] : 0.0;

    for (ST_int i = 0; i < n; i++) {
        const ST_double xi = x[i];
        KAHAN_STEP(s0, c0, xi, r0[i], fused);
        if (nr > 1) KAHAN_STEP(s1, c1, xi, r1[i], fused);
        if (nr > 2) KAHAN_STEP(s2, c2, xi, r2[i], fused);
        if (nr > 3) KAHAN_STEP(s3, c3, xi, r3[i], fused);
    }
    sc[0] = s0; sc[1] = c0;
    if (nr > 1) { sc[2] = s1; sc[3] = c1; }
    if (nr > 2) { sc[4] = s2; sc[5] = c2; }
    if (nr > 3) { sc[6] = s3; sc[7] = c3; }
}

#define KAHAN_ROWS_PARAMS (const ST_double *x, const ST_double *const *r, ST_int n, ST_double *sc)
#define KAHAN_ROWS(nr)                                                                             \
    static CTOOLS_INLINE void kahan_rows##nr##_body(const int fused, const ST_double *x,           \
                                                    const ST_double *const *r, ST_int n,           \
                                                    ST_double *sc)                                 \
    {                                                                                              \
        kahan_rows_body(fused, nr, x, r, n, sc);                                                   \
    }                                                                                              \
    CTOOLS_FMA_VARIANTS_VOID(static, kahan_rows##nr, KAHAN_ROWS_PARAMS, x, r, n, sc)
KAHAN_ROWS(1)
KAHAN_ROWS(2)
KAHAN_ROWS(3)
KAHAN_ROWS(4)

static void kahan_rows(int nr, const ST_double *x, const ST_double *const *r, ST_int n,
                       ST_double *sc)
{
    switch (nr) {
    case 1: kahan_rows1(x, r, n, sc); break;
    case 2: kahan_rows2(x, r, n, sc); break;
    case 3: kahan_rows3(x, r, n, sc); break;
    default: kahan_rows4(x, r, n, sc); break;
    }
}

ST_double kahan_dot(const ST_double * RESTRICT x,
                    const ST_double * RESTRICT y,
                    ST_int N)
{
    ST_double sc[2] = {0.0, 0.0};
    const ST_double *r[1] = {y};

    kahan_rows1(x, r, N, sc);
    return sc[0] + sc[1];
}

/* ------------------------------------------------------------------------
 * ctools_kahan_cross: X'X and X'y, compensated, blocked and in parallel
 *
 * The pairs of left column j are (x_j, y) when y is given, then (x_j, x_k)
 * for k = j .. K-1, in that order: the right-hand columns of x_j, taken
 * four at a time.
 * ------------------------------------------------------------------------ */

/* Rows per pass inside a block of the block arrangement: the columns of
 * that many rows stay in the cache while every pair is taken. The sums are
 * continued from pass to pass, so this changes no result. */
#define CTOOLS_KAHAN_ROWS 2048
/* Blocks summed at once by the block arrangement, and the (s, c) pairs it
 * may keep per block: beyond that the column arrangement runs */
#define CTOOLS_KAHAN_WAVE 64
#define CTOOLS_KAHAN_PAIRS ((size_t)1 << 16)
/* Element products below which the work stays on the calling thread */
#define CTOOLS_KAHAN_PARALLEL 400000.0

static inline const ST_double *kahan_cross_rhs(const ST_double *X, const ST_double *y, ST_int N,
                                               ST_int j, ST_int rhs)
{
    if (y) {
        if (rhs == 0) return y;
        rhs--;
    }
    return X + (size_t)(j + rhs) * N;
}

/* Where the sums go: doubles (X'X in both triangles, and X'y), or the upper
 * triangle of X'X in long double */
typedef struct {
    ST_double *xtx, *xty;
    long double *wide;
    int has_y;
    ST_int K;
} kahan_cross_out;

/* The pair's value is hi + lo */
static inline void kahan_cross_store(const kahan_cross_out *o, ST_int j, ST_int rhs, ST_double hi,
                                     ST_double lo)
{
    if (o->has_y) {
        if (rhs == 0) {
            o->xty[j] = hi + lo;
            return;
        }
        rhs--;
    }
    const ST_int k = j + rhs;
    if (o->wide) {
        o->wide[(size_t)j * o->K + k] = (long double)hi + (long double)lo;
        return;
    }
    const ST_double value = hi + lo;
    o->xtx[(size_t)j * o->K + k] = value;
    o->xtx[(size_t)k * o->K + j] = value;
}

/* Rows lo .. lo+len-1 of left column j against its right-hand columns
 * t .. t+n_now-1, the (s, c) pairs in sc continued. With w the left column
 * is w[i] * x_j[i], rounded before it is multiplied. */
static inline void kahan_cross_rows(const ST_double *X, const ST_double *y, const ST_double *w,
                                    ST_int N, ST_int j, ST_int t, int n_now, size_t lo, ST_int len,
                                    ST_double *sc)
{
    ST_double wx[CTOOLS_KAHAN_ROWS];
    const ST_double *col[4], *r[4];

    for (int q = 0; q < n_now; q++) col[q] = kahan_cross_rhs(X, y, N, j, t + q);
    for (ST_int r0 = 0; r0 < len; r0 += CTOOLS_KAHAN_ROWS) {
        const ST_int rows = (len - r0 < CTOOLS_KAHAN_ROWS) ? len - r0 : CTOOLS_KAHAN_ROWS;
        const size_t at = lo + (size_t)r0;
        const ST_double *x = X + (size_t)j * N + at;
        if (w) {
            const ST_double *wi = w + at;
            for (ST_int i = 0; i < rows; i++) wx[i] = wi[i] * x[i];
            x = wx;
        }
        for (int q = 0; q < n_now; q++) r[q] = col[q] + at;
        kahan_rows(n_now, x, r, rows, sc);
    }
}

/* Rows lo .. lo+len-1, every pair: sc receives (s, c) per pair, the pairs in
 * the order of the left columns and their right-hand columns */
static void kahan_cross_block(const ST_double *X, const ST_double *y, const ST_double *w, ST_int N,
                              ST_int K, ST_int lo, ST_int len, ST_double *sc, size_t npairs)
{
    memset(sc, 0, 2 * npairs * sizeof(ST_double));
    for (ST_int r0 = 0; r0 < len; r0 += CTOOLS_KAHAN_ROWS) {
        const ST_int rows = (len - r0 < CTOOLS_KAHAN_ROWS) ? len - r0 : CTOOLS_KAHAN_ROWS;
        ST_double *p = sc;
        for (ST_int j = 0; j < K; j++) {
            const ST_int n_rhs = K - j + (y != NULL);
            for (ST_int t = 0; t < n_rhs; t += 4) {
                const int n_now = (n_rhs - t < 4) ? (int)(n_rhs - t) : 4;
                kahan_cross_rows(X, y, w, N, j, t, n_now, (size_t)lo + (size_t)r0, rows, p);
                p += 2 * n_now;
            }
        }
    }
}

static void kahan_cross_core(const ST_double *X, const ST_double *y, const ST_double *w, ST_int N,
                             ST_int K, const kahan_cross_out *out)
{
    const int has_y = out->has_y;
    /* N + CTOOLS_KAHAN_BLOCK - 1 would overflow int for N near INT_MAX */
    const ST_int nblocks = N / CTOOLS_KAHAN_BLOCK + (N % CTOOLS_KAHAN_BLOCK != 0);
    const size_t npairs = (size_t)K * ((size_t)K + 1) / 2 + (has_y ? (size_t)K : 0);
    const int threaded = (double)N * (double)npairs >= CTOOLS_KAHAN_PARALLEL;
    const ST_int wave = (nblocks < CTOOLS_KAHAN_WAVE) ? nblocks : CTOOLS_KAHAN_WAVE;
    ST_double *part = NULL;
    dd_real *total = NULL;
    (void)threaded; /* read only by the OpenMP clauses */

    /* Block arrangement: the blocks of a wave in parallel, every pair in one
     * pass over a block's rows, then the wave's blocks added in order */
    if (nblocks > 1 && npairs <= CTOOLS_KAHAN_PAIRS) {
        part = (ST_double *)malloc((size_t)wave * 2 * npairs * sizeof(ST_double));
        total = (dd_real *)malloc(npairs * sizeof(dd_real));
        if (!part || !total) {
            free(part);
            free(total);
            part = NULL;
            total = NULL;
        }
    }
    if (part) {
        for (size_t p = 0; p < npairs; p++) total[p] = (dd_real){0.0, 0.0};
        for (ST_int b0 = 0; b0 < nblocks; b0 += wave) {
            const ST_int nb = (nblocks - b0 < wave) ? nblocks - b0 : wave;
            ST_int b;
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 1) if(threaded)
#endif
            for (b = 0; b < nb; b++) {
                const ST_int lo = (b0 + b) * CTOOLS_KAHAN_BLOCK;
                const ST_int len = (N - lo < CTOOLS_KAHAN_BLOCK) ? N - lo : CTOOLS_KAHAN_BLOCK;
                kahan_cross_block(X, y, w, N, K, lo, len, part + (size_t)b * 2 * npairs, npairs);
            }
            for (b = 0; b < nb; b++) {
                const ST_double *sc = part + (size_t)b * 2 * npairs;
                for (size_t p = 0; p < npairs; p++) {
                    total[p] = dd_add_d(dd_add_d(total[p], sc[2 * p]), sc[2 * p + 1]);
                }
            }
        }
        size_t p = 0;
        for (ST_int j = 0; j < K; j++) {
            const ST_int n_rhs = K - j + has_y;
            for (ST_int rhs = 0; rhs < n_rhs; rhs++, p++) {
                kahan_cross_store(out, j, rhs, total[p].hi, total[p].lo);
            }
        }
        free(part);
        free(total);
        return;
    }

    /* Column arrangement (one block, very many pairs, or no memory for the
     * block sums): the left columns in parallel, each taking its blocks in
     * order. The same sums, with no scratch beyond the stack. */
    ST_int j;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) if(threaded && K > 1)
#endif
    for (j = 0; j < K; j++) {
        const ST_int n_rhs = K - j + has_y;
        for (ST_int t = 0; t < n_rhs; t += 4) {
            const int n_now = (n_rhs - t < 4) ? (int)(n_rhs - t) : 4;
            dd_real sum[4] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
            ST_double sc[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            for (ST_int b = 0; b < nblocks; b++) {
                const ST_int lo = b * CTOOLS_KAHAN_BLOCK;
                const ST_int len = (N - lo < CTOOLS_KAHAN_BLOCK) ? N - lo : CTOOLS_KAHAN_BLOCK;
                for (int q = 0; q < 8; q++) sc[q] = 0.0;
                kahan_cross_rows(X, y, w, N, j, t, n_now, (size_t)lo, len, sc);
                if (nblocks > 1) {
                    for (int q = 0; q < n_now; q++) {
                        sum[q] = dd_add_d(dd_add_d(sum[q], sc[2 * q]), sc[2 * q + 1]);
                    }
                }
            }
            for (int q = 0; q < n_now; q++) {
                /* one block: the sequential sum's own (s, c) */
                if (nblocks > 1) kahan_cross_store(out, j, t + q, sum[q].hi, sum[q].lo);
                else kahan_cross_store(out, j, t + q, sc[2 * q], sc[2 * q + 1]);
            }
        }
    }
}

void ctools_kahan_cross(const ST_double *X, const ST_double *y, ST_int N, ST_int K,
                        ST_double *xtx, ST_double *xty)
{
    if (K <= 0) return;
    if (N < 0) N = 0;
    if (!xty) y = NULL;
    const kahan_cross_out out = {xtx, xty, NULL, y != NULL, K};
    kahan_cross_core(X, y, NULL, N, K, &out);
}

void ctools_kahan_cross_weighted(const ST_double *X, const ST_double *w, ST_int N, ST_int K,
                                 long double *upper)
{
    if (K <= 0) return;
    if (N < 0) N = 0;
    const kahan_cross_out out = {NULL, NULL, upper, 0, K};
    kahan_cross_core(X, NULL, w, N, K, &out);
}

/* ========================================================================
 * Double-double (quad-precision) accumulation functions
 * These match Stata's quadcross() precision for RSS/TSS calculations
 * ======================================================================== */

static CTOOLS_INLINE ST_double dd_sum_sq_body(const int fused, const ST_double *x, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        dd_real prod = two_prod_as(fused, x[i], x[i]);
        /* Add both hi and lo parts to accumulator */
        sum = dd_add_d(sum, prod.hi);
        sum = dd_add_d(sum, prod.lo);
    }
    return sum.hi + sum.lo;
}
CTOOLS_FMA_VARIANTS(, ST_double, dd_sum_sq, (const ST_double *x, ST_int N), x, N)

static CTOOLS_INLINE ST_double dd_sum_sq_weighted_body(const int fused, const ST_double *x,
                                                       const ST_double *w, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        /* Compute x[i]^2 with error tracking */
        dd_real sq = two_prod_as(fused, x[i], x[i]);
        /* Multiply by weight: w * sq.hi and w * sq.lo */
        dd_real wsq_hi = two_prod_as(fused, w[i], sq.hi);
        dd_real wsq_lo = two_prod_as(fused, w[i], sq.lo);
        /* Accumulate all parts */
        sum = dd_add_d(sum, wsq_hi.hi);
        sum = dd_add_d(sum, wsq_hi.lo);
        sum = dd_add_d(sum, wsq_lo.hi);
        sum = dd_add_d(sum, wsq_lo.lo);
    }
    return sum.hi + sum.lo;
}
CTOOLS_FMA_VARIANTS(, ST_double, dd_sum_sq_weighted,
                    (const ST_double *x, const ST_double *w, ST_int N), x, w, N)

static CTOOLS_INLINE ST_double dd_sum_sq_dev_body(const int fused, const ST_double *x,
                                                  ST_double mean, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double d = x[i] - mean;
        dd_real prod = two_prod_as(fused, d, d);
        sum = dd_add_d(sum, prod.hi);
        sum = dd_add_d(sum, prod.lo);
    }
    return sum.hi + sum.lo;
}
CTOOLS_FMA_VARIANTS(, ST_double, dd_sum_sq_dev, (const ST_double *x, ST_double mean, ST_int N), x,
                    mean, N)

static CTOOLS_INLINE ST_double dd_sum_sq_dev_weighted_body(const int fused, const ST_double *x,
                                                           ST_double mean, const ST_double *w,
                                                           ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double d = x[i] - mean;
        dd_real sq = two_prod_as(fused, d, d);
        dd_real wsq_hi = two_prod_as(fused, w[i], sq.hi);
        dd_real wsq_lo = two_prod_as(fused, w[i], sq.lo);
        sum = dd_add_d(sum, wsq_hi.hi);
        sum = dd_add_d(sum, wsq_hi.lo);
        sum = dd_add_d(sum, wsq_lo.hi);
        sum = dd_add_d(sum, wsq_lo.lo);
    }
    return sum.hi + sum.lo;
}
CTOOLS_FMA_VARIANTS(, ST_double, dd_sum_sq_dev_weighted,
                    (const ST_double *x, ST_double mean, const ST_double *w, ST_int N), x, mean, w,
                    N)

static CTOOLS_INLINE dd_real dd_dot_weighted_body(const int fused, const ST_double *x,
                                                  const ST_double *y, const ST_double *w, ST_int N)
{
    dd_real sum = {0.0, 0.0};
    ST_int i;

    for (i = 0; i < N; i++) {
        ST_double wx = w ? w[i] * x[i] : x[i];
        dd_real prod = two_prod_as(fused, wx, y[i]);
        sum = dd_add_d(sum, prod.hi);
        sum = dd_add_d(sum, prod.lo);
    }
    return sum;
}
CTOOLS_FMA_VARIANTS(, dd_real, dd_dot_weighted,
                    (const ST_double *x, const ST_double *y, const ST_double *w, ST_int N), x, y,
                    w, N)

static CTOOLS_INLINE void dd_dot_weighted_multi_body(const int fused, const ST_double *x,
                                                     const ST_double *const *ys, ST_int n,
                                                     const ST_double *w, ST_int N, dd_real *out)
{
    ST_int t = 0;

    for (; t + 3 < n; t += 4) {
        const ST_double *y0 = ys[t], *y1 = ys[t + 1], *y2 = ys[t + 2], *y3 = ys[t + 3];
        dd_real s0 = {0.0, 0.0}, s1 = {0.0, 0.0}, s2 = {0.0, 0.0}, s3 = {0.0, 0.0};
        for (ST_int i = 0; i < N; i++) {
            ST_double wx = w ? w[i] * x[i] : x[i];
            dd_real p0 = two_prod_as(fused, wx, y0[i]);
            dd_real p1 = two_prod_as(fused, wx, y1[i]);
            dd_real p2 = two_prod_as(fused, wx, y2[i]);
            dd_real p3 = two_prod_as(fused, wx, y3[i]);
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
    /* The one to three sums left over, in one pass as well: a double-double
     * sum is one long chain of dependent additions, so sums taken together
     * cost little more than one taken alone (each is still its own sum over
     * the rows in order: the same values as one pass per sum). */
    if (n - t == 3) {
        const ST_double *y0 = ys[t], *y1 = ys[t + 1], *y2 = ys[t + 2];
        dd_real s0 = {0.0, 0.0}, s1 = {0.0, 0.0}, s2 = {0.0, 0.0};
        for (ST_int i = 0; i < N; i++) {
            ST_double wx = w ? w[i] * x[i] : x[i];
            dd_real p0 = two_prod_as(fused, wx, y0[i]);
            dd_real p1 = two_prod_as(fused, wx, y1[i]);
            dd_real p2 = two_prod_as(fused, wx, y2[i]);
            s0 = dd_add_d(dd_add_d(s0, p0.hi), p0.lo);
            s1 = dd_add_d(dd_add_d(s1, p1.hi), p1.lo);
            s2 = dd_add_d(dd_add_d(s2, p2.hi), p2.lo);
        }
        out[t] = s0;
        out[t + 1] = s1;
        out[t + 2] = s2;
    } else if (n - t == 2) {
        const ST_double *y0 = ys[t], *y1 = ys[t + 1];
        dd_real s0 = {0.0, 0.0}, s1 = {0.0, 0.0};
        for (ST_int i = 0; i < N; i++) {
            ST_double wx = w ? w[i] * x[i] : x[i];
            dd_real p0 = two_prod_as(fused, wx, y0[i]);
            dd_real p1 = two_prod_as(fused, wx, y1[i]);
            s0 = dd_add_d(dd_add_d(s0, p0.hi), p0.lo);
            s1 = dd_add_d(dd_add_d(s1, p1.hi), p1.lo);
        }
        out[t] = s0;
        out[t + 1] = s1;
    } else if (n - t == 1) {
        out[t] = dd_dot_weighted_body(fused, x, ys[t], w, N);
    }
}
CTOOLS_FMA_VARIANTS_VOID(, dd_dot_weighted_multi,
                         (const ST_double *x, const ST_double *const *ys, ST_int n,
                          const ST_double *w, ST_int N, dd_real *out),
                         x, ys, n, w, N, out)

/* dd_dot_weighted_multi over fixed blocks of CTOOLS_DD_BLOCK rows: each
 * block's sums run in row order and the blocks are added in block order, so
 * the result is independent of the thread count (and equal to the
 * sequential sums to double-double accuracy). Blocks run in parallel unless
 * the caller is already in a parallel region. Returns 0 or 920. */
#define CTOOLS_DD_BLOCK 32768
int dd_dot_weighted_multi_blocked(const ST_double *x, const ST_double *const *ys, ST_int n,
                                  const ST_double *w, ST_int N, dd_real *out)
{
    ST_int nblocks = (N + CTOOLS_DD_BLOCK - 1) / CTOOLS_DD_BLOCK;
    if (n < 1) return 0;
    if (nblocks <= 1) {
        dd_dot_weighted_multi(x, ys, n, w, N, out);
        return 0;
    }
    dd_real *part = (dd_real *)malloc((size_t)nblocks * (size_t)n * sizeof(dd_real));
    if (!part) return 920;
    int failed = 0;
#ifdef _OPENMP
    #pragma omp parallel if (nblocks > 8 && !omp_in_parallel())
#endif
    {
        const ST_double **ysb = (const ST_double **)malloc((size_t)n * sizeof(*ysb));
        if (!ysb) {
#ifdef _OPENMP
            #pragma omp atomic write
#endif
            failed = 1;
        }
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (ST_int b = 0; b < nblocks; b++) {
            if (!ysb) continue;
            ST_int lo = b * CTOOLS_DD_BLOCK;
            ST_int len = N - lo < CTOOLS_DD_BLOCK ? N - lo : CTOOLS_DD_BLOCK;
            for (ST_int t = 0; t < n; t++) ysb[t] = ys[t] + lo;
            dd_dot_weighted_multi(x + lo, ysb, n, w ? w + lo : NULL, len, part + (size_t)b * n);
        }
        free(ysb);
    }
    if (failed) { free(part); return 920; }
    for (ST_int t = 0; t < n; t++) {
        dd_real s = {0.0, 0.0};
        for (ST_int b = 0; b < nblocks; b++) {
            s = dd_add_d(s, part[(size_t)b * n + t].hi);
            s = dd_add_d(s, part[(size_t)b * n + t].lo);
        }
        out[t] = s;
    }
    free(part);
    return 0;
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
    /* Kahan-compensated dot products for numerical stability.
     * reghdfe uses Mata's quadcross() (quad-precision), so we need
     * compensated summation to match after FE partialling out.
     * BLAS dsyrk would be faster but lacks the precision needed here. */
    ctools_kahan_cross(data + (size_t)N, data, N, K - 1, xtx, xty);
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
    const int scaled = (weight_type == 1 || weight_type == 3);
    ST_double weight_scale = 1.0;
    if (scaled) {
        ST_double sum_w = 0.0;
        for (i = 0; i < N; i++) sum_w += weights[i];
        weight_scale = (ST_double)N / sum_w;
    }

    /* Initialize output */
    memset(xtx, 0, (size_t)K_x * K_x * sizeof(ST_double));
    memset(xty, 0, K_x * sizeof(ST_double));

    /* xty[j] = sum_i (w_i x_ij) y_i and xtx[j,k] = sum_i (w_i x_ij) x_ik for
     * k >= j, every element accumulated over the observations in order with
     * plain double arithmetic, w_i = weights[i] * weight_scale rounded to a
     * double (formed in the loop: a stored copy of the scaled weights would
     * cost a pass over the rows and N doubles). The elements are independent
     * sums, so the rows of the upper triangle run in parallel, four
     * right-hand columns per pass over the data. */
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
            if (scaled) {
                for (i = 0; i < N; i++) {
                    const ST_double w_i = weights[i] * weight_scale;
                    ST_double wxj = w_i * xj[i];
                    acc[0] += wxj * r[0][i];
                    acc[1] += wxj * r[1][i];
                    acc[2] += wxj * r[2][i];
                    acc[3] += wxj * r[3][i];
                }
            } else {
                for (i = 0; i < N; i++) {
                    ST_double wxj = weights[i] * xj[i];
                    acc[0] += wxj * r[0][i];
                    acc[1] += wxj * r[1][i];
                    acc[2] += wxj * r[2][i];
                    acc[3] += wxj * r[3][i];
                }
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
static ST_int invsym_sweep(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out,
                           int in_order);

ST_int ctools_invsym(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out)
{
    return invsym_sweep(A, K, inv, rank_out, 0);
}

ST_int ctools_invsym_ordered(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out)
{
    return invsym_sweep(A, K, inv, rank_out, 1);
}

/* in_order: sweep the columns 0..K-1 in turn (Mata's invsym(A, 1..K));
 * otherwise pivot on the largest remaining diagonal */
static ST_int invsym_sweep(const ST_double *A, ST_int K, ST_double *inv, ST_int *rank_out,
                           int in_order)
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
            if (in_order && p >= 0) break;
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
