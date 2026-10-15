#include "ctools_linalg.h"
/*
 * cqreg_linalg.c
 *
 * Linear algebra operations for the quantile regression solver.
 * Part of the ctools suite.
 */

#include "cqreg_linalg.h"
#include "../ctools_unroll.h"
#include <math.h>
#include <string.h>

/* ============================================================================
 * Dot Products
 * ============================================================================ */

ST_double cqreg_dot(const ST_double * CQREG_RESTRICT x,
                    const ST_double * CQREG_RESTRICT y,
                    ST_int N)
{
    return ctools_dot_unrolled(x, y, N);
}

ST_double cqreg_dot_self(const ST_double * CQREG_RESTRICT x, ST_int N)
{
    return ctools_dot_self_unrolled(x, N);
}

/* ============================================================================
 * Cholesky Decomposition and Solve
 * ============================================================================ */

ST_int cqreg_cholesky(ST_double * CQREG_RESTRICT A, ST_int K)
{
    ST_int j, k, i;

    for (j = 0; j < K; j++) {
        ST_double sum = A[j * K + j];

        /* Subtract L[j,k]^2 for k < j */
        for (k = 0; k < j; k++) {
            ST_double Ljk = A[j * K + k];
            sum -= Ljk * Ljk;
        }

        if (sum <= 0.0) {
            /* Matrix is not positive definite */
            return -1;
        }

        A[j * K + j] = sqrt(sum);

        /* Update column j below diagonal */
        for (i = j + 1; i < K; i++) {
            sum = A[i * K + j];

            for (k = 0; k < j; k++) {
                sum -= A[i * K + k] * A[j * K + k];
            }

            A[i * K + j] = sum / A[j * K + j];
        }
    }

    return 0;
}

static void cqreg_solve_lower(const ST_double * CQREG_RESTRICT L,
                              ST_double * CQREG_RESTRICT b,
                              ST_int K)
{
    ctools_solve_lower_factor(L, b, K, b);
}

static void cqreg_solve_lower_t(const ST_double * CQREG_RESTRICT L,
                                ST_double * CQREG_RESTRICT b,
                                ST_int K)
{
    ctools_solve_lower_transpose(L, b, K, b);
}

void cqreg_solve_cholesky(const ST_double * CQREG_RESTRICT L,
                          ST_double * CQREG_RESTRICT b,
                          ST_int K)
{
    /* Solve A * x = b where A = L * L' */
    /* Step 1: Solve L * y = b */
    cqreg_solve_lower(L, b, K);
    /* Step 2: Solve L' * x = y */
    cqreg_solve_lower_t(L, b, K);
}

void cqreg_invert_cholesky(ST_double * CQREG_RESTRICT Ainv,
                           const ST_double * CQREG_RESTRICT L,
                           ST_int K)
{
    ST_int j;

    /* Compute A^{-1} by solving A * Ainv[:,j] = e_j for each column */
    memset(Ainv, 0, (size_t)K * K * sizeof(ST_double));

    for (j = 0; j < K; j++) {
        /* Set column j to e_j */
        ST_double *col = &Ainv[(size_t)j * K];
        col[j] = 1.0;

        /* Solve L * L' * col = e_j */
        cqreg_solve_cholesky(L, col, K);
    }

    /* Note: result is stored column-major but should be symmetric */
}

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

void cqreg_add_regularization(ST_double * CQREG_RESTRICT A, ST_int K, ST_double lambda)
{
    ST_int i;

    for (i = 0; i < K; i++) {
        A[i * K + i] += lambda;
    }
}

/* ============================================================================
 * Quantile and Statistics Functions
 * ============================================================================ */

/* Comparison function for qsort */
static int compare_double(const void *a, const void *b)
{
    ST_double da = *(const ST_double *)a;
    ST_double db = *(const ST_double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

ST_double cqreg_select(ST_double *values, ST_int N, ST_int k)
{
    ST_int lo = 0, hi = N;
    size_t budget = (size_t)N;
    budget = budget <= SIZE_MAX / 8 ? 8 * budget : SIZE_MAX;

    /* Three-way partitioning removes all copies of the pivot in one pass.
     * Two-way quickselect can take quadratic time on tied residuals. Bound
     * the total partition work, falling back to sorting only the remainder
     * if median-of-three pivots repeatedly fail to shrink the search. */
    while (hi - lo > 1) {
        size_t count = (size_t)(hi - lo);
        if (count > budget) {
            qsort(values + lo, count, sizeof(*values), compare_double);
            return values[k];
        }
        budget -= count;
        ST_double a = values[lo], b = values[lo + (hi - lo) / 2];
        ST_double c = values[hi - 1], tmp;
        if (a > b) { tmp = a; a = b; b = tmp; }
        if (b > c) { b = c; }
        ST_double pivot = a > b ? a : b;
        ST_int lt = lo, i = lo, gt = hi;
        while (i < gt) {
            if (values[i] < pivot) {
                tmp = values[lt]; values[lt++] = values[i]; values[i++] = tmp;
            } else if (values[i] > pivot) {
                tmp = values[--gt]; values[gt] = values[i]; values[i] = tmp;
            } else {
                i++;
            }
        }
        if (k < lt) hi = lt;
        else if (k >= gt) lo = gt;
        else return values[k];
    }
    return values[k];
}

ST_double cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau)
{
    if (N <= 0) return 0.0;
    /* Select from a copy so the response stays in observation order. */
    ST_double *y_work = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    if (y_work == NULL) {
        return NAN;  /* the caller reports r(920); 0 would pass as a quantile */
    }

    memcpy(y_work, y, (size_t)N * sizeof(ST_double));

    /*
     * Compute quantile using Stata's qreg method:
     * k = floor((N + 1) * tau)
     *
     * Examples:
     * - N=74, tau=0.25: floor(75*0.25) = floor(18.75) = 18, return y[18]
     * - N=74, tau=0.50: floor(75*0.50) = floor(37.5) = 37, return y[37]
     * - N=74, tau=0.75: floor(75*0.75) = floor(56.25) = 56, return y[56]
     * - N=15, tau=0.50: floor(16*0.50) = floor(8.0) = 8, return y[8]
     */
    ST_double pos = ((ST_double)N + 1.0) * tau;
    ST_int k = (ST_int)floor(pos);

    /* Handle boundary cases */
    if (k < 1) k = 1;
    if (k > N) k = N;

    ST_double q = cqreg_select(y_work, N, k - 1);  /* Convert to 0-indexed */

    free(y_work);
    return q;
}

ST_double cqreg_sum_raw_deviations(const ST_double *y, ST_int N, ST_double q, ST_double tau)
{
    ST_double sum = 0.0;
    ST_int i;

    /* Check function: rho_tau(u) = u * (tau - I(u < 0))
     * = tau * u if u >= 0
     * = (tau - 1) * u if u < 0
     */
    for (i = 0; i < N; i++) {
        ST_double u = y[i] - q;
        if (u >= 0.0) {
            sum += tau * u;
        } else {
            sum += (tau - 1.0) * u;
        }
    }

    return sum;
}
