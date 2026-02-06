/*
 * cqreg_linalg.c
 *
 * Optimized linear algebra operations for quantile regression IPM solver.
 * Part of the ctools suite.
 */

#include "cqreg_linalg.h"
#include "../ctools_unroll.h"
#include "../ctools_select.h"
#include "../ctools_simd.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * Dot Products
 * ============================================================================ */

ST_double cqreg_dot(const ST_double * CQREG_RESTRICT x,
                    const ST_double * CQREG_RESTRICT y,
                    ST_int N)
{
    return ctools_dot_unrolled(x, y, N);
}

ST_double cqreg_dot_kahan(const ST_double * CQREG_RESTRICT x,
                          const ST_double * CQREG_RESTRICT y,
                          ST_int N)
{
    ST_int i;
    ST_double sum = 0.0;
    ST_double c = 0.0;  /* Compensation for lost low-order bits */

    for (i = 0; i < N; i++) {
        ST_double prod = x[i] * y[i] - c;
        ST_double t = sum + prod;
        c = (t - sum) - prod;
        sum = t;
    }

    return sum;
}

ST_double cqreg_dot_weighted(const ST_double * CQREG_RESTRICT x,
                             const ST_double * CQREG_RESTRICT y,
                             const ST_double * CQREG_RESTRICT w,
                             ST_int N)
{
    /* Use SIMD-accelerated weighted dot product (AVX2/NEON with FMA) */
    return ctools_simd_dot_weighted(w, x, y, (size_t)N);
}

ST_double cqreg_dot_self(const ST_double * CQREG_RESTRICT x, ST_int N)
{
    return ctools_dot_self_unrolled(x, N);
}

/* ============================================================================
 * Vector Operations
 * ============================================================================ */

void cqreg_vcopy(ST_double * CQREG_RESTRICT dst,
                 const ST_double * CQREG_RESTRICT src,
                 ST_int N)
{
    memcpy(dst, src, N * sizeof(ST_double));
}

void cqreg_vscale(ST_double * CQREG_RESTRICT x, ST_double alpha, ST_int N)
{
    ST_int i;
    ST_int N4 = N - (N % 4);

    #pragma omp simd
    for (i = 0; i < N4; i += 4) {
        x[i]     *= alpha;
        x[i + 1] *= alpha;
        x[i + 2] *= alpha;
        x[i + 3] *= alpha;
    }

    for (; i < N; i++) {
        x[i] *= alpha;
    }
}

void cqreg_vaxpy(ST_double *y,
                 ST_double alpha,
                 const ST_double *x,
                 ST_int N)
{
    /* Use SIMD-accelerated AXPY (AVX2/NEON with FMA) */
    ctools_simd_axpy(y, x, alpha, (size_t)N);
}

void cqreg_vmul(ST_double * CQREG_RESTRICT z,
                const ST_double * CQREG_RESTRICT x,
                const ST_double * CQREG_RESTRICT y,
                ST_int N)
{
    ST_int i;

    #pragma omp simd
    for (i = 0; i < N; i++) {
        z[i] = x[i] * y[i];
    }
}

void cqreg_vdiv(ST_double * CQREG_RESTRICT z,
                const ST_double * CQREG_RESTRICT x,
                const ST_double * CQREG_RESTRICT y,
                ST_int N)
{
    ST_int i;

    #pragma omp simd
    for (i = 0; i < N; i++) {
        z[i] = x[i] / y[i];
    }
}

void cqreg_vset(ST_double * CQREG_RESTRICT x, ST_double val, ST_int N)
{
    ST_int i;

    #pragma omp simd
    for (i = 0; i < N; i++) {
        x[i] = val;
    }
}

ST_double cqreg_vsum(const ST_double * CQREG_RESTRICT x, ST_int N)
{
    ST_int i;
    ST_double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    ST_int N4 = N - (N % 4);

    for (i = 0; i < N4; i += 4) {
        sum0 += x[i];
        sum1 += x[i + 1];
        sum2 += x[i + 2];
        sum3 += x[i + 3];
    }

    for (; i < N; i++) {
        sum0 += x[i];
    }

    return (sum0 + sum1) + (sum2 + sum3);
}

ST_double cqreg_vmax_abs(const ST_double * CQREG_RESTRICT x, ST_int N)
{
    ST_int i;
    ST_int N4 = N - (N & 3);
    ST_double max0 = 0.0, max1 = 0.0, max2 = 0.0, max3 = 0.0;

    /* 4-way unrolled loop for better ILP */
    for (i = 0; i < N4; i += 4) {
        ST_double abs0 = fabs(x[i]);
        ST_double abs1 = fabs(x[i + 1]);
        ST_double abs2 = fabs(x[i + 2]);
        ST_double abs3 = fabs(x[i + 3]);
        if (abs0 > max0) max0 = abs0;
        if (abs1 > max1) max1 = abs1;
        if (abs2 > max2) max2 = abs2;
        if (abs3 > max3) max3 = abs3;
    }

    /* Remainder */
    for (; i < N; i++) {
        ST_double abs_val = fabs(x[i]);
        if (abs_val > max0) max0 = abs_val;
    }

    /* Reduce */
    if (max1 > max0) max0 = max1;
    if (max2 > max0) max0 = max2;
    if (max3 > max0) max0 = max3;

    return max0;
}

/* ============================================================================
 * Matrix-Vector Operations
 * ============================================================================ */

void cqreg_matvec(ST_double * CQREG_RESTRICT y,
                  const ST_double * CQREG_RESTRICT A,
                  const ST_double * CQREG_RESTRICT x,
                  ST_int M, ST_int N)
{
    ST_int i, j;

    for (i = 0; i < M; i++) {
        ST_double sum = 0.0;
        const ST_double *row = &A[i * N];
        for (j = 0; j < N; j++) {
            sum += row[j] * x[j];
        }
        y[i] = sum;
    }
}

void cqreg_matvec_t(ST_double * CQREG_RESTRICT y,
                    const ST_double * CQREG_RESTRICT A,
                    const ST_double * CQREG_RESTRICT x,
                    ST_int M, ST_int N)
{
    ST_int i, j;

    memset(y, 0, N * sizeof(ST_double));

    for (i = 0; i < M; i++) {
        const ST_double *row = &A[i * N];
        ST_double xi = x[i];
        for (j = 0; j < N; j++) {
            y[j] += row[j] * xi;
        }
    }
}

void cqreg_matvec_col(ST_double * CQREG_RESTRICT y,
                      const ST_double * CQREG_RESTRICT A,
                      const ST_double * CQREG_RESTRICT x,
                      ST_int M, ST_int N)
{
    ST_int i, j;
    ST_int M8 = M - (M & 7);  /* M rounded down to multiple of 8 */

    /* Initialize y to zero */
    memset(y, 0, M * sizeof(ST_double));

    /* y = sum_j x[j] * A[:,j]
     * OPTIMIZED: 8x loop unrolling for better vectorization */
    for (j = 0; j < N; j++) {
        const ST_double *col = &A[j * M];
        ST_double xj = x[j];

        /* 8x unrolled main loop */
        for (i = 0; i < M8; i += 8) {
            y[i]     += xj * col[i];
            y[i + 1] += xj * col[i + 1];
            y[i + 2] += xj * col[i + 2];
            y[i + 3] += xj * col[i + 3];
            y[i + 4] += xj * col[i + 4];
            y[i + 5] += xj * col[i + 5];
            y[i + 6] += xj * col[i + 6];
            y[i + 7] += xj * col[i + 7];
        }
        /* Remainder */
        for (; i < M; i++) {
            y[i] += xj * col[i];
        }
    }
}

void cqreg_matvec_t_col(ST_double * CQREG_RESTRICT y,
                        const ST_double * CQREG_RESTRICT A,
                        const ST_double * CQREG_RESTRICT x,
                        ST_int M, ST_int N)
{
    ST_int j;

    /* y[j] = A[:,j]' * x = dot(A[:,j], x) */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(N >= 2 && M > 500)
    #endif
    for (j = 0; j < N; j++) {
        y[j] = cqreg_dot(&A[j * M], x, M);
    }
}

void cqreg_xtdx(ST_double * CQREG_RESTRICT XDX,
                const ST_double * CQREG_RESTRICT X,
                const ST_double * CQREG_RESTRICT D,
                ST_int N, ST_int K)
{
    ST_int j, k;

    /* Initialize to zero */
    memset(XDX, 0, K * K * sizeof(ST_double));

    /* Compute X' * D * X using SIMD-accelerated weighted dot products
     * XDX[j,k] = sum_i X[i,j] * D[i] * X[i,k]
     *          = sum_i D[i] * X[j*N + i] * X[k*N + i]
     * = ctools_simd_dot_weighted(D, Xj, Xk, N)
     */

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(K >= 2 && N > 500)
    #endif
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];

        /* Diagonal element: XDX[j,j] = X[:,j]' * D * X[:,j] */
        XDX[j * K + j] = ctools_simd_dot_weighted(D, Xj, Xj, (size_t)N);

        /* Off-diagonal elements (upper triangle) */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double total = ctools_simd_dot_weighted(D, Xj, Xk, (size_t)N);
            XDX[j * K + k] = total;
            XDX[k * K + j] = total;  /* Symmetric */
        }
    }
}

void cqreg_xtv(ST_double * CQREG_RESTRICT result,
               const ST_double * CQREG_RESTRICT X,
               const ST_double * CQREG_RESTRICT v,
               ST_int N, ST_int K)
{
    ST_int j;

    /* result[j] = X[:,j]' * v */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(K >= 2 && N > 500)
    #endif
    for (j = 0; j < K; j++) {
        result[j] = cqreg_dot(&X[j * N], v, N);
    }
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

void cqreg_solve_lower(const ST_double * CQREG_RESTRICT L,
                       ST_double * CQREG_RESTRICT b,
                       ST_int K)
{
    ST_int i, j;

    /* Forward substitution: L * x = b */
    for (i = 0; i < K; i++) {
        ST_double sum = b[i];
        for (j = 0; j < i; j++) {
            sum -= L[i * K + j] * b[j];
        }
        b[i] = sum / L[i * K + i];
    }
}

void cqreg_solve_lower_t(const ST_double * CQREG_RESTRICT L,
                         ST_double * CQREG_RESTRICT b,
                         ST_int K)
{
    ST_int i, j;

    /* Back substitution: L' * x = b */
    for (i = K - 1; i >= 0; i--) {
        ST_double sum = b[i];
        for (j = i + 1; j < K; j++) {
            sum -= L[j * K + i] * b[j];
        }
        b[i] = sum / L[i * K + i];
    }
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
    memset(Ainv, 0, K * K * sizeof(ST_double));

    for (j = 0; j < K; j++) {
        /* Set column j to e_j */
        ST_double *col = &Ainv[j * K];
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

ST_int cqreg_is_symmetric(const ST_double * CQREG_RESTRICT A, ST_int K, ST_double tol)
{
    ST_int i, j;

    for (i = 0; i < K; i++) {
        for (j = i + 1; j < K; j++) {
            if (fabs(A[i * K + j] - A[j * K + i]) > tol) {
                return 0;
            }
        }
    }

    return 1;
}

void cqreg_symmetrize_lower(ST_double * CQREG_RESTRICT A, ST_int K)
{
    ST_int i, j;

    for (i = 0; i < K; i++) {
        for (j = i + 1; j < K; j++) {
            A[i * K + j] = A[j * K + i];
        }
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

ST_double cqreg_compute_quantile(const ST_double *y, ST_int N, ST_double tau)
{
    /* Make a copy of y for sorting */
    ST_double *y_sorted = (ST_double *)malloc(N * sizeof(ST_double));
    if (y_sorted == NULL) {
        return 0.0;  /* Error: return 0 */
    }

    memcpy(y_sorted, y, N * sizeof(ST_double));

    /* Sort using qsort (O(N log N)) */
    qsort(y_sorted, N, sizeof(ST_double), compare_double);

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

    ST_double q = y_sorted[k - 1];  /* Convert to 0-indexed */

    free(y_sorted);
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
