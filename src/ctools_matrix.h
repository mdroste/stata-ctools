/*
 * ctools_matrix.h
 *
 * Shared matrix operations for creghdfe and civreghdfe.
 * Provides optimized matrix multiplication with OpenMP and SIMD.
 * Part of the ctools Stata plugin suite.
 */

#ifndef CTOOLS_MATRIX_H
#define CTOOLS_MATRIX_H

#include "stplugin.h"

/*
 * Matrix multiply C = A' * B
 * A is N x K1, B is N x K2, Result C is K1 x K2
 * All matrices stored in column-major order.
 *
 * OPTIMIZED: Uses K-way unrolled dot products and OpenMP parallelization
 */
void ctools_matmul_atb(const ST_double * restrict A, const ST_double * restrict B,
                       ST_int N, ST_int K1, ST_int K2,
                       ST_double * restrict C);

/*
 * Matrix multiply C = A * B
 * A is K1 x K2, B is K2 x K3, Result C is K1 x K3
 * All matrices stored in column-major order.
 *
 * OPTIMIZED: Uses cache-friendly loop order and OpenMP parallelization
 */
void ctools_matmul_ab(const ST_double * restrict A, const ST_double * restrict B,
                      ST_int K1, ST_int K2, ST_int K3,
                      ST_double * restrict C);

/*
 * Weighted matrix multiply C = A' * diag(w) * B
 * A is N x K1, B is N x K2, w is N x 1 (weights, may be NULL)
 * Result C is K1 x K2
 *
 * OPTIMIZED: Uses SIMD-accelerated weighted dot products and OpenMP parallelization
 */
void ctools_matmul_atdb(const ST_double * restrict A, const ST_double * restrict B,
                        const ST_double * restrict w, ST_int N, ST_int K1, ST_int K2,
                        ST_double * restrict C);

#endif /* CTOOLS_MATRIX_H */
