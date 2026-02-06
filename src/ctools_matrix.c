/*
 * ctools_matrix.c
 *
 * Shared matrix operations for creghdfe and civreghdfe.
 * Provides optimized matrix multiplication with OpenMP and SIMD.
 * Part of the ctools Stata plugin suite.
 */

#include <stdlib.h>
#include <string.h>

#include "ctools_matrix.h"
#include "ctools_config.h"
#include "ctools_unroll.h"
#include "ctools_simd.h"

/*
 * Matrix multiply C = A' * B
 * A is N x K1, B is N x K2, Result C is K1 x K2
 *
 * OPTIMIZED: Uses K-way unrolled dot products and OpenMP parallelization
 */
void ctools_matmul_atb(const ST_double * restrict A, const ST_double * restrict B,
                       ST_int N, ST_int K1, ST_int K2,
                       ST_double * restrict C)
{
    ST_int j;

    memset(C, 0, K1 * K2 * sizeof(ST_double));

    /* Parallelize over output columns - i must be private to avoid race condition */
    #pragma omp parallel for schedule(static) if(K1 * K2 > 4)
    for (j = 0; j < K2; j++) {
        ST_int i;  /* Declare inside parallel region to make it private */
        const ST_double *b_col = B + j * N;
        for (i = 0; i < K1; i++) {
            const ST_double *a_col = A + i * N;
            C[j * K1 + i] = ctools_dot_unrolled(a_col, b_col, N);
        }
    }
}

/*
 * Matrix multiply C = A * B
 * A is K1 x K2, B is K2 x K3, Result C is K1 x K3
 *
 * OPTIMIZED: Uses cache-friendly loop order and OpenMP parallelization
 */
void ctools_matmul_ab(const ST_double * restrict A, const ST_double * restrict B,
                      ST_int K1, ST_int K2, ST_int K3,
                      ST_double * restrict C)
{
    ST_int j;

    memset(C, 0, K1 * K3 * sizeof(ST_double));

    /* C[i,j] = sum_k A[i,k] * B[k,j] */
    /* A is K1 x K2, stored column-major: A[i,k] = A[k*K1 + i] */
    /* B is K2 x K3, stored column-major: B[k,j] = B[j*K2 + k] */
    /* C is K1 x K3, stored column-major: C[i,j] = C[j*K1 + i] */

    /* Reorder loops for better cache locality: k-j-i instead of j-i-k */
    /* This allows sequential access to A columns and B columns */
    #pragma omp parallel for schedule(static) if(K1 * K3 > 16)
    for (j = 0; j < K3; j++) {
        ST_int i, k;  /* Private to each thread */
        const ST_double *b_col = B + j * K2;
        ST_double *c_col = C + j * K1;
        for (k = 0; k < K2; k++) {
            ST_double b_kj = b_col[k];
            const ST_double *a_col = A + k * K1;
            #pragma omp simd
            for (i = 0; i < K1; i++) {
                c_col[i] += a_col[i] * b_kj;
            }
        }
    }
}

/*
 * Weighted matrix multiply C = A' * diag(w) * B
 *
 * OPTIMIZED: Uses SIMD-accelerated weighted dot products (AVX2/NEON with FMA)
 * and OpenMP parallelization
 */
void ctools_matmul_atdb(const ST_double * restrict A, const ST_double * restrict B,
                        const ST_double * restrict w, ST_int N, ST_int K1, ST_int K2,
                        ST_double * restrict C)
{
    ST_int j;

    memset(C, 0, K1 * K2 * sizeof(ST_double));

    if (w) {
        /* Weighted case: use SIMD-accelerated weighted dot product */
        #pragma omp parallel for schedule(static) if(K1 * K2 > 4)
        for (j = 0; j < K2; j++) {
            ST_int i;  /* Private to each thread */
            const ST_double *b_col = B + j * N;
            for (i = 0; i < K1; i++) {
                const ST_double *a_col = A + i * N;
                C[j * K1 + i] = ctools_simd_dot_weighted(w, a_col, b_col, (size_t)N);
            }
        }
    } else {
        /* Unweighted case: use SIMD-accelerated dot product */
        #pragma omp parallel for schedule(static) if(K1 * K2 > 4)
        for (j = 0; j < K2; j++) {
            ST_int i;  /* Private to each thread */
            const ST_double *b_col = B + j * N;
            for (i = 0; i < K1; i++) {
                const ST_double *a_col = A + i * N;
                C[j * K1 + i] = ctools_simd_dot(a_col, b_col, (size_t)N);
            }
        }
    }
}
