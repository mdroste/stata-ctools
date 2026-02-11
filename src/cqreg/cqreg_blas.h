/*
 * cqreg_blas.h
 *
 * BLAS/LAPACK abstraction layer for cqreg.
 * Uses Apple Accelerate on macOS, OpenBLAS on Linux/Windows.
 *
 * Part of the ctools suite.
 */

#ifndef CQREG_BLAS_H
#define CQREG_BLAS_H

#include "cqreg_types.h"

/* ============================================================================
 * Platform Detection and BLAS Selection
 * ============================================================================
 *
 * The build system defines:
 *   SYSTEM=APPLEMAC  - macOS (use Accelerate)
 *   SYSTEM=STUNIX    - Linux (use OpenBLAS)
 *   SYSTEM=STWIN32   - Windows (use OpenBLAS)
 *
 * Compile flags:
 *   -DUSE_BLAS=1     - Enable BLAS (default on all platforms)
 *   -DUSE_BLAS=0     - Force pure C fallback
 */

/* Default: enable BLAS on all platforms */
#ifndef USE_BLAS
#define USE_BLAS 1
#endif

/* Platform-specific BLAS headers */
#if USE_BLAS

#if defined(SYSTEM) && SYSTEM == APPLEMAC
    /* Apple Accelerate framework - built into macOS, highly optimized */
    #include <Accelerate/Accelerate.h>
    #define BLAS_IMPL_NAME "Apple Accelerate"
    #define HAVE_ACCELERATE 1

#elif defined(SYSTEM) && (SYSTEM == STUNIX || SYSTEM == STWIN32)
    /* OpenBLAS - open source, available on Linux and Windows */
    #ifdef HAVE_OPENBLAS
        #include <cblas.h>
        #include <lapacke.h>
        #define BLAS_IMPL_NAME "OpenBLAS"
    #else
        /* Fallback to pure C if OpenBLAS not available */
        #undef USE_BLAS
        #define USE_BLAS 0
        #define BLAS_IMPL_NAME "Pure C (no BLAS)"
    #endif

#else
    /* Unknown platform - use pure C */
    #undef USE_BLAS
    #define USE_BLAS 0
    #define BLAS_IMPL_NAME "Pure C (unknown platform)"
#endif

#else
    #define BLAS_IMPL_NAME "Pure C (BLAS disabled)"
#endif /* USE_BLAS */


/* ============================================================================
 * BLAS Level 2 Operations (Matrix-Vector)
 * ============================================================================ */

/*
 * General matrix-vector multiply: y = alpha * A * x + beta * y
 * A is M x N (column-major)
 *
 * trans = 0: y = alpha * A * x + beta * y    (A is M x N, x is N, y is M)
 * trans = 1: y = alpha * A' * x + beta * y   (A is M x N, x is M, y is N)
 */
void blas_dgemv(int trans, ST_int M, ST_int N,
                ST_double alpha, const ST_double *A, ST_int lda,
                const ST_double *x,
                ST_double beta, ST_double *y);


/* ============================================================================
 * High-Level Operations (Optimized combinations)
 * ============================================================================ */

/*
 * Compute X' * D * X where D is diagonal
 * X is N x K (column-major), D is N, result is K x K (column-major)
 *
 * This is a hot path in the IPM solver, heavily optimized.
 */
void blas_xtdx(ST_double *XDX,
               const ST_double *X, ST_int N, ST_int K,
               const ST_double *D);


/* ============================================================================
 * Vectorized Element-wise Operations (SIMD-optimized)
 * Uses vDSP on macOS for maximum performance.
 * ============================================================================ */

/*
 * Element-wise divide: C = A / B
 */
void blas_vdiv(ST_int N, const ST_double *A, const ST_double *B, ST_double *C);

/*
 * Element-wise add: C = A + B
 */
void blas_vadd(ST_int N, const ST_double *A, const ST_double *B, ST_double *C);

/*
 * Element-wise subtract: C = A - B
 */
void blas_vsub(ST_int N, const ST_double *A, const ST_double *B, ST_double *C);

/*
 * Element-wise multiply: C = A * B
 */
void blas_vmul(ST_int N, const ST_double *A, const ST_double *B, ST_double *C);

/*
 * Scalar add: C = A + scalar
 */
void blas_vsadd(ST_int N, const ST_double *A, ST_double scalar, ST_double *C);

/*
 * Scalar multiply: C = A * scalar
 */
void blas_vsmul(ST_int N, const ST_double *A, ST_double scalar, ST_double *C);


#endif /* CQREG_BLAS_H */
