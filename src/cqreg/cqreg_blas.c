/*
 * cqreg_blas.c
 *
 * BLAS abstraction layer implementation (dgemv and the X'DX kernel).
 * Uses Apple Accelerate on macOS, OpenBLAS on Linux/Windows,
 * with pure C fallbacks for all operations.
 *
 * Part of the ctools suite.
 */

#include "cqreg_blas.h"
#include "../ctools_unroll.h"
#include "../ctools_config.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif


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
            y[i] = beta == 0.0 ? 0.0 : beta * y[i];
        }
        for (j = 0; j < N; j++) {
            ST_double axj = alpha * x[j];
            const ST_double *Aj = &A[(size_t)j * lda];
            for (i = 0; i < M; i++) {
                y[i] += axj * Aj[i];
            }
        }
    } else {
        /* y = alpha * A' * x + beta * y, A is M x N column-major */
        /* y[j] = alpha * sum_i A[i + j*lda] * x[i] + beta * y[j] */
        for (j = 0; j < N; j++) {
            const ST_double *Aj = &A[(size_t)j * lda];
            ST_double sum = 0.0;
            for (i = 0; i < M; i++) {
                sum += Aj[i] * x[i];
            }
            y[j] = beta == 0.0 ? alpha * sum : alpha * sum + beta * y[j];
        }
    }
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
     * Compute X' * diag(D) * X
     *
     * For small K (typical in regression), direct computation is faster
     * than BLAS dgemm because it avoids allocating N*K temporary arrays.
     *
     * Parallelization strategy:
     * - Single parallel region over N to minimize thread creation overhead
     * - Each thread accumulates its own K*K partial sums
     * - Final reduction combines partial sums
     * - Uses 4-way unrolling for better pipelining
     */

    ST_int j, k;

    /* Initialize output to zero */
    memset(XDX, 0, (size_t)K * K * sizeof(ST_double));

/*
     * OpenMP parallelization strategy:
     * Parallelize over (j,k) matrix element pairs - each element is independent.
     * This avoids thread-local storage and works well with small K.
     * For typical K=2-10 in regression, this gives ~K*(K+1)/2 parallel tasks.
     */
#ifdef _OPENMP
    if (N > 100000 && K >= 2) {
        /* Each element XDX[j,k] = sum_i D[i] * X[i,j] * X[i,k] of the upper
         * triangle is an independent sum, accumulated four rows at a time
         * exactly as the sequential version below does. For a fixed j the
         * columns k are taken XTDX_GROUP at a time, so D * Xj is formed once
         * per row and the data are read in far fewer passes than one per
         * element; the (j, group) tasks run in parallel. */
        enum { XTDX_GROUP = 8 };
        const ST_int N4 = N - (N % 4);
        const ST_int groups_per_row = (K + XTDX_GROUP - 1) / XTDX_GROUP;
        ST_int task;
        #pragma omp parallel for schedule(dynamic)
        for (task = 0; task < K * groups_per_row; task++) {
            const ST_int jj = task / groups_per_row;
            const ST_int k0 = jj + (task % groups_per_row) * XTDX_GROUP;
            if (k0 >= K) continue;
            const ST_int nk = (K - k0 < XTDX_GROUP) ? K - k0 : XTDX_GROUP;
            const ST_double *Xj = &X[(size_t)jj * N];
            const ST_double *Xk[XTDX_GROUP];
            ST_double sum[XTDX_GROUP];
            for (ST_int q = 0; q < XTDX_GROUP; q++) {
                /* unused slots repeat the first column */
                Xk[q] = &X[(size_t)(k0 + (q < nk ? q : 0)) * N];
                sum[q] = 0.0;
            }

            for (ST_int i = 0; i < N4; i += 4) {
                const ST_double dx0 = D[i]   * Xj[i];
                const ST_double dx1 = D[i+1] * Xj[i+1];
                const ST_double dx2 = D[i+2] * Xj[i+2];
                const ST_double dx3 = D[i+3] * Xj[i+3];
                for (ST_int q = 0; q < XTDX_GROUP; q++) {
                    const ST_double *xk = Xk[q];
                    sum[q] += dx0 * xk[i] + dx1 * xk[i+1] + dx2 * xk[i+2] + dx3 * xk[i+3];
                }
            }
            for (ST_int i = N4; i < N; i++) {
                const ST_double dx = D[i] * Xj[i];
                for (ST_int q = 0; q < XTDX_GROUP; q++) sum[q] += dx * Xk[q][i];
            }

            for (ST_int q = 0; q < nk; q++) {
                const ST_int kk = k0 + q;
                XDX[jj + kk * K] = sum[q];
                if (kk > jj) {
                    XDX[kk + jj * K] = sum[q];  /* Symmetrize */
                }
            }
        }
    } else
#endif
    {
        /* Sequential version with 4-way unrolling */
        const ST_int N4 = N - (N % 4);

        for (ST_int i = 0; i < N4; i += 4) {
            ST_double D0 = D[i], D1 = D[i+1], D2 = D[i+2], D3 = D[i+3];

            for (j = 0; j < K; j++) {
                const ST_double *Xj = &X[(size_t)j * N];
                ST_double X0j = Xj[i];
                ST_double X1j = Xj[i + 1];
                ST_double X2j = Xj[i + 2];
                ST_double X3j = Xj[i + 3];

                ST_double DX0j = D0 * X0j;
                ST_double DX1j = D1 * X1j;
                ST_double DX2j = D2 * X2j;
                ST_double DX3j = D3 * X3j;

                /* Diagonal */
                XDX[j + j * K] += DX0j * X0j + DX1j * X1j + DX2j * X2j + DX3j * X3j;

                /* Off-diagonal (upper triangle) */
                for (k = j + 1; k < K; k++) {
                    const ST_double *Xk = &X[(size_t)k * N];
                    ST_double X0k = Xk[i];
                    ST_double X1k = Xk[i + 1];
                    ST_double X2k = Xk[i + 2];
                    ST_double X3k = Xk[i + 3];
                    XDX[j + k * K] += DX0j * X0k + DX1j * X1k + DX2j * X2k + DX3j * X3k;
                }
            }
        }

        /* Handle remainder */
        for (ST_int i = N4; i < N; i++) {
            ST_double Di = D[i];
            for (j = 0; j < K; j++) {
                ST_double DXij = Di * X[(size_t)j * N + i];
                XDX[j + j * K] += DXij * X[(size_t)j * N + i];
                for (k = j + 1; k < K; k++) {
                    XDX[j + k * K] += DXij * X[(size_t)k * N + i];
                }
            }
        }

        /* Symmetrize */
        for (j = 0; j < K; j++) {
            for (k = j + 1; k < K; k++) {
                XDX[k + j * K] = XDX[j + k * K];
            }
        }
    }
}


