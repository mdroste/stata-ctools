/*
 * cqreg_blas.c
 *
 * BLAS/LAPACK abstraction layer implementation.
 * Uses Apple Accelerate on macOS, OpenBLAS on Linux/Windows,
 * with pure C fallbacks for all operations.
 *
 * Part of the ctools suite.
 */

#include "cqreg_blas.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif


/* ============================================================================
 * Utility
 * ============================================================================ */

const char* blas_get_impl_name(void)
{
    return BLAS_IMPL_NAME;
}


/* ============================================================================
 * BLAS Level 1 Operations
 * ============================================================================ */

ST_double blas_ddot(ST_int N, const ST_double *x, const ST_double *y)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    return cblas_ddot(N, x, 1, y, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    return cblas_ddot(N, x, 1, y, 1);
#else
    /* Pure C fallback with 8-way unrolling */
    ST_int i;
    ST_double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    ST_double sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
    ST_int N8 = N - (N % 8);

    for (i = 0; i < N8; i += 8) {
        sum0 += x[i]     * y[i];
        sum1 += x[i + 1] * y[i + 1];
        sum2 += x[i + 2] * y[i + 2];
        sum3 += x[i + 3] * y[i + 3];
        sum4 += x[i + 4] * y[i + 4];
        sum5 += x[i + 5] * y[i + 5];
        sum6 += x[i + 6] * y[i + 6];
        sum7 += x[i + 7] * y[i + 7];
    }

    for (; i < N; i++) {
        sum0 += x[i] * y[i];
    }

    return ((sum0 + sum4) + (sum1 + sum5)) + ((sum2 + sum6) + (sum3 + sum7));
#endif
}


void blas_daxpy(ST_int N, ST_double alpha, const ST_double *x, ST_double *y)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_daxpy(N, alpha, x, 1, y, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_daxpy(N, alpha, x, 1, y, 1);
#else
    ST_int i;
    for (i = 0; i < N; i++) {
        y[i] += alpha * x[i];
    }
#endif
}


void blas_dscal(ST_int N, ST_double alpha, ST_double *x)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_dscal(N, alpha, x, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_dscal(N, alpha, x, 1);
#else
    ST_int i;
    ST_int N4 = N - (N % 4);
    for (i = 0; i < N4; i += 4) {
        x[i]     *= alpha;
        x[i + 1] *= alpha;
        x[i + 2] *= alpha;
        x[i + 3] *= alpha;
    }
    for (; i < N; i++) {
        x[i] *= alpha;
    }
#endif
}


void blas_dcopy(ST_int N, const ST_double *x, ST_double *y)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_dcopy(N, x, 1, y, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_dcopy(N, x, 1, y, 1);
#else
    memcpy(y, x, N * sizeof(ST_double));
#endif
}


ST_double blas_dasum(ST_int N, const ST_double *x)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    return cblas_dasum(N, x, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    return cblas_dasum(N, x, 1);
#else
    ST_int i;
    ST_double sum = 0.0;
    for (i = 0; i < N; i++) {
        sum += fabs(x[i]);
    }
    return sum;
#endif
}


ST_double blas_dnrm2(ST_int N, const ST_double *x)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    return cblas_dnrm2(N, x, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    return cblas_dnrm2(N, x, 1);
#else
    return sqrt(blas_ddot(N, x, x));
#endif
}


/* ============================================================================
 * BLAS Level 2 Operations
 * ============================================================================ */

void blas_dgemv(int trans, ST_int M, ST_int N,
                ST_double alpha, const ST_double *A, ST_int lda,
                const ST_double *x,
                ST_double beta, ST_double *y)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_dgemv(CblasColMajor,
                trans ? CblasTrans : CblasNoTrans,
                M, N, alpha, A, lda, x, 1, beta, y, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_dgemv(CblasColMajor,
                trans ? CblasTrans : CblasNoTrans,
                M, N, alpha, A, lda, x, 1, beta, y, 1);
#else
    ST_int i, j;

    if (trans == 0) {
        /* y = alpha * A * x + beta * y, A is M x N column-major */
        /* y[i] = alpha * sum_j A[i + j*lda] * x[j] + beta * y[i] */
        for (i = 0; i < M; i++) {
            y[i] *= beta;
        }
        for (j = 0; j < N; j++) {
            ST_double axj = alpha * x[j];
            const ST_double *Aj = &A[j * lda];
            for (i = 0; i < M; i++) {
                y[i] += axj * Aj[i];
            }
        }
    } else {
        /* y = alpha * A' * x + beta * y, A is M x N column-major */
        /* y[j] = alpha * sum_i A[i + j*lda] * x[i] + beta * y[j] */
        for (j = 0; j < N; j++) {
            const ST_double *Aj = &A[j * lda];
            ST_double sum = 0.0;
            for (i = 0; i < M; i++) {
                sum += Aj[i] * x[i];
            }
            y[j] = alpha * sum + beta * y[j];
        }
    }
#endif
}


void blas_dsymv(ST_int N,
                ST_double alpha, const ST_double *A, ST_int lda,
                const ST_double *x,
                ST_double beta, ST_double *y)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_dsymv(CblasColMajor, CblasLower, N, alpha, A, lda, x, 1, beta, y, 1);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_dsymv(CblasColMajor, CblasLower, N, alpha, A, lda, x, 1, beta, y, 1);
#else
    /* Fallback: treat as general matrix-vector multiply */
    blas_dgemv(0, N, N, alpha, A, lda, x, beta, y);
#endif
}


/* ============================================================================
 * BLAS Level 3 Operations
 * ============================================================================ */

void blas_dgemm(int transA, int transB,
                ST_int M, ST_int N, ST_int K,
                ST_double alpha, const ST_double *A, ST_int lda,
                const ST_double *B, ST_int ldb,
                ST_double beta, ST_double *C, ST_int ldc)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_dgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                transB ? CblasTrans : CblasNoTrans,
                M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_dgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                transB ? CblasTrans : CblasNoTrans,
                M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#else
    /* Pure C fallback - naive triple loop */
    ST_int i, j, k;

    /* Initialize C with beta scaling */
    for (j = 0; j < N; j++) {
        for (i = 0; i < M; i++) {
            C[i + j * ldc] *= beta;
        }
    }

    /* C += alpha * op(A) * op(B) */
    if (transA == 0 && transB == 0) {
        /* C = alpha * A * B + beta * C */
        /* A is M x K, B is K x N, C is M x N */
        for (j = 0; j < N; j++) {
            for (k = 0; k < K; k++) {
                ST_double alphaBkj = alpha * B[k + j * ldb];
                for (i = 0; i < M; i++) {
                    C[i + j * ldc] += alphaBkj * A[i + k * lda];
                }
            }
        }
    } else if (transA == 1 && transB == 0) {
        /* C = alpha * A' * B + beta * C */
        /* A is K x M, B is K x N, C is M x N */
        for (j = 0; j < N; j++) {
            for (i = 0; i < M; i++) {
                ST_double sum = 0.0;
                for (k = 0; k < K; k++) {
                    sum += A[k + i * lda] * B[k + j * ldb];
                }
                C[i + j * ldc] += alpha * sum;
            }
        }
    } else if (transA == 0 && transB == 1) {
        /* C = alpha * A * B' + beta * C */
        for (j = 0; j < N; j++) {
            for (k = 0; k < K; k++) {
                ST_double alphaBjk = alpha * B[j + k * ldb];
                for (i = 0; i < M; i++) {
                    C[i + j * ldc] += alphaBjk * A[i + k * lda];
                }
            }
        }
    } else {
        /* C = alpha * A' * B' + beta * C */
        for (j = 0; j < N; j++) {
            for (i = 0; i < M; i++) {
                ST_double sum = 0.0;
                for (k = 0; k < K; k++) {
                    sum += A[k + i * lda] * B[j + k * ldb];
                }
                C[i + j * ldc] += alpha * sum;
            }
        }
    }
#endif
}


void blas_dsyrk(int trans, ST_int N, ST_int K,
                ST_double alpha, const ST_double *A, ST_int lda,
                ST_double beta, ST_double *C, ST_int ldc)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    cblas_dsyrk(CblasColMajor, CblasLower,
                trans ? CblasTrans : CblasNoTrans,
                N, K, alpha, A, lda, beta, C, ldc);
    /* Symmetrize - copy lower to upper */
    for (ST_int j = 0; j < N; j++) {
        for (ST_int i = j + 1; i < N; i++) {
            C[j + i * ldc] = C[i + j * ldc];
        }
    }
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    cblas_dsyrk(CblasColMajor, CblasLower,
                trans ? CblasTrans : CblasNoTrans,
                N, K, alpha, A, lda, beta, C, ldc);
    /* Symmetrize */
    for (ST_int j = 0; j < N; j++) {
        for (ST_int i = j + 1; i < N; i++) {
            C[j + i * ldc] = C[i + j * ldc];
        }
    }
#else
    /* Use dgemm fallback */
    if (trans == 0) {
        /* C = alpha * A * A' + beta * C */
        blas_dgemm(0, 1, N, N, K, alpha, A, lda, A, lda, beta, C, ldc);
    } else {
        /* C = alpha * A' * A + beta * C */
        blas_dgemm(1, 0, N, N, K, alpha, A, lda, A, lda, beta, C, ldc);
    }
#endif
}


/* ============================================================================
 * LAPACK Operations
 * ============================================================================ */

ST_int lapack_dpotrf(ST_int N, ST_double *A, ST_int lda)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    /* Apple Accelerate uses CLAPACK interface */
    char uplo = 'L';
    int info = 0;
    int n = (int)N;
    int ld = (int)lda;
    dpotrf_(&uplo, &n, A, &ld, &info);
    return info;
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    return LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', N, A, lda);
#else
    /* Pure C Cholesky factorization */
    ST_int i, j, k;

    for (j = 0; j < N; j++) {
        /* Compute diagonal element */
        ST_double sum = A[j + j * lda];
        for (k = 0; k < j; k++) {
            ST_double Ljk = A[j + k * lda];
            sum -= Ljk * Ljk;
        }

        if (sum <= 0.0) {
            return j + 1;  /* Not positive definite */
        }

        A[j + j * lda] = sqrt(sum);

        /* Compute column j below diagonal */
        for (i = j + 1; i < N; i++) {
            sum = A[i + j * lda];
            for (k = 0; k < j; k++) {
                sum -= A[i + k * lda] * A[j + k * lda];
            }
            A[i + j * lda] = sum / A[j + j * lda];
        }
    }

    return 0;
#endif
}


ST_int lapack_dpotrs(ST_int N, ST_int NRHS,
                     const ST_double *A, ST_int lda,
                     ST_double *B, ST_int ldb)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    char uplo = 'L';
    int info = 0;
    int n = (int)N;
    int nrhs = (int)NRHS;
    int ld_a = (int)lda;
    int ld_b = (int)ldb;
    dpotrs_(&uplo, &n, &nrhs, (ST_double *)A, &ld_a, B, &ld_b, &info);
    return info;
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    return LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', N, NRHS, A, lda, B, ldb);
#else
    /* Pure C triangular solve */
    ST_int i, j, k;

    for (k = 0; k < NRHS; k++) {
        ST_double *b = &B[k * ldb];

        /* Forward substitution: L * y = b */
        for (i = 0; i < N; i++) {
            ST_double sum = b[i];
            for (j = 0; j < i; j++) {
                sum -= A[i + j * lda] * b[j];
            }
            b[i] = sum / A[i + i * lda];
        }

        /* Back substitution: L' * x = y */
        for (i = N - 1; i >= 0; i--) {
            ST_double sum = b[i];
            for (j = i + 1; j < N; j++) {
                sum -= A[j + i * lda] * b[j];
            }
            b[i] = sum / A[i + i * lda];
        }
    }

    return 0;
#endif
}


ST_int lapack_dpotri(ST_int N, ST_double *A, ST_int lda)
{
#if USE_BLAS && defined(HAVE_ACCELERATE)
    char uplo = 'L';
    int info = 0;
    int n = (int)N;
    int ld = (int)lda;
    dpotri_(&uplo, &n, A, &ld, &info);
    /* Symmetrize */
    for (ST_int j = 0; j < N; j++) {
        for (ST_int i = j + 1; i < N; i++) {
            A[j + i * lda] = A[i + j * lda];
        }
    }
    return info;
#elif USE_BLAS && defined(HAVE_OPENBLAS)
    int info = LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', N, A, lda);
    /* Symmetrize */
    for (ST_int j = 0; j < N; j++) {
        for (ST_int i = j + 1; i < N; i++) {
            A[j + i * lda] = A[i + j * lda];
        }
    }
    return info;
#else
    /* Pure C inverse via column-by-column solve */
    ST_int j, k;
    ST_double *col = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *Ainv = (ST_double *)malloc(N * N * sizeof(ST_double));

    if (!col || !Ainv) {
        free(col);
        free(Ainv);
        return -1;
    }

    for (j = 0; j < N; j++) {
        /* Set column to e_j */
        memset(col, 0, N * sizeof(ST_double));
        col[j] = 1.0;

        /* Solve L * L' * col = e_j */
        /* Forward substitution */
        for (k = 0; k < N; k++) {
            ST_double sum = col[k];
            for (ST_int m = 0; m < k; m++) {
                sum -= A[k + m * lda] * col[m];
            }
            col[k] = sum / A[k + k * lda];
        }

        /* Back substitution */
        for (k = N - 1; k >= 0; k--) {
            ST_double sum = col[k];
            for (ST_int m = k + 1; m < N; m++) {
                sum -= A[m + k * lda] * col[m];
            }
            col[k] = sum / A[k + k * lda];
        }

        /* Store column j of inverse */
        for (k = 0; k < N; k++) {
            Ainv[k + j * N] = col[k];
        }
    }

    /* Copy result back to A */
    memcpy(A, Ainv, N * N * sizeof(ST_double));

    free(col);
    free(Ainv);

    return 0;
#endif
}


/* ============================================================================
 * High-Level Operations
 * ============================================================================ */

void blas_xtdx(ST_double *XDX,
               const ST_double *X, ST_int N, ST_int K,
               const ST_double *D)
{
    /*
     * Compute X' * diag(D) * X = X' * (D .* X)
     *
     * Strategy: Form W = D .* X (scaled X), then compute X' * W
     * This is O(N*K) for scaling + O(N*K^2) for matrix multiply
     */

#if USE_BLAS && (defined(HAVE_ACCELERATE) || defined(HAVE_OPENBLAS))
    /*
     * For large N, use BLAS dgemm: XDX = X' * W
     * where W[:,j] = D .* X[:,j]
     */
    ST_int j, i;

    /* Allocate scaled matrix W */
    ST_double *W = (ST_double *)malloc(N * K * sizeof(ST_double));
    if (!W) {
        /* Fallback to direct computation */
        goto fallback;
    }

    /* Form W = diag(D) * X, i.e., W[:,j] = D .* X[:,j] */
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(K > 2)
    #endif
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];
        ST_double *Wj = &W[j * N];
        for (i = 0; i < N; i++) {
            Wj[i] = D[i] * Xj[i];
        }
    }

    /* XDX = X' * W using BLAS dgemm */
    /* X is N x K (col-major), W is N x K, XDX is K x K */
    /* XDX = 1.0 * X' * W + 0.0 * XDX */
    blas_dgemm(1, 0, K, K, N, 1.0, X, N, W, N, 0.0, XDX, K);

    free(W);
    return;

fallback:
#endif
    {
        /* Pure C implementation with OpenMP parallelization */
        ST_int j, k, i;

        memset(XDX, 0, K * K * sizeof(ST_double));

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) if(K > 2)
        #endif
        for (j = 0; j < K; j++) {
            const ST_double *Xj = &X[j * N];

            /* Diagonal element */
            ST_double diag = 0.0;
            for (i = 0; i < N; i++) {
                diag += D[i] * Xj[i] * Xj[i];
            }
            XDX[j + j * K] = diag;

            /* Off-diagonal (upper triangle, then symmetrize) */
            for (k = j + 1; k < K; k++) {
                const ST_double *Xk = &X[k * N];
                ST_double sum = 0.0;
                for (i = 0; i < N; i++) {
                    sum += D[i] * Xj[i] * Xk[i];
                }
                XDX[j + k * K] = sum;
                XDX[k + j * K] = sum;
            }
        }
    }
}


void blas_xtv(ST_double *result,
              const ST_double *X, ST_int N, ST_int K,
              const ST_double *v)
{
    /*
     * Compute result = X' * v
     * X is N x K (column-major), v is N, result is K
     */

#if USE_BLAS && (defined(HAVE_ACCELERATE) || defined(HAVE_OPENBLAS))
    /* Use BLAS dgemv: result = X' * v */
    blas_dgemv(1, N, K, 1.0, X, N, v, 0.0, result);
#else
    ST_int j;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(K > 2)
    #endif
    for (j = 0; j < K; j++) {
        result[j] = blas_ddot(N, &X[j * N], v);
    }
#endif
}


ST_int blas_solve_weighted_ls(const ST_double *X, ST_int N, ST_int K,
                               const ST_double *D,
                               const ST_double *y,
                               ST_double *beta)
{
    /*
     * Solve (X'DX) * beta = X'Dy via Cholesky
     */

    /* Allocate workspace */
    ST_double *XDX = (ST_double *)malloc(K * K * sizeof(ST_double));
    ST_double *Dy = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *rhs = (ST_double *)malloc(K * sizeof(ST_double));

    if (!XDX || !Dy || !rhs) {
        free(XDX);
        free(Dy);
        free(rhs);
        return -1;
    }

    /* Compute X'DX */
    blas_xtdx(XDX, X, N, K, D);

    /* Add small regularization for numerical stability */
    for (ST_int j = 0; j < K; j++) {
        XDX[j + j * K] += 1e-14;
    }

    /* Compute Dy = D .* y */
    for (ST_int i = 0; i < N; i++) {
        Dy[i] = D[i] * y[i];
    }

    /* Compute rhs = X' * Dy */
    blas_xtv(rhs, X, N, K, Dy);

    /* Cholesky factorization */
    ST_int info = lapack_dpotrf(K, XDX, K);
    if (info != 0) {
        free(XDX);
        free(Dy);
        free(rhs);
        return info;
    }

    /* Solve */
    memcpy(beta, rhs, K * sizeof(ST_double));
    info = lapack_dpotrs(K, 1, XDX, K, beta, K);

    free(XDX);
    free(Dy);
    free(rhs);

    return info;
}
