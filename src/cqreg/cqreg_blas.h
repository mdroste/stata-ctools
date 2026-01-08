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
 * BLAS Level 1 Operations (Vector-Vector)
 * ============================================================================ */

/*
 * Dot product: result = x' * y
 */
ST_double blas_ddot(ST_int N, const ST_double *x, const ST_double *y);

/*
 * Scaled vector addition: y = alpha * x + y
 */
void blas_daxpy(ST_int N, ST_double alpha, const ST_double *x, ST_double *y);

/*
 * Vector scale: x = alpha * x
 */
void blas_dscal(ST_int N, ST_double alpha, ST_double *x);

/*
 * Vector copy: y = x
 */
void blas_dcopy(ST_int N, const ST_double *x, ST_double *y);

/*
 * Sum of absolute values: result = |x[0]| + |x[1]| + ... + |x[N-1]|
 */
ST_double blas_dasum(ST_int N, const ST_double *x);

/*
 * Euclidean norm: result = sqrt(x' * x)
 */
ST_double blas_dnrm2(ST_int N, const ST_double *x);


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

/*
 * Symmetric matrix-vector multiply: y = alpha * A * x + beta * y
 * A is N x N symmetric (column-major, lower triangle stored)
 */
void blas_dsymv(ST_int N,
                ST_double alpha, const ST_double *A, ST_int lda,
                const ST_double *x,
                ST_double beta, ST_double *y);


/* ============================================================================
 * BLAS Level 3 Operations (Matrix-Matrix)
 * ============================================================================ */

/*
 * General matrix-matrix multiply: C = alpha * A * B + beta * C
 * All matrices in column-major format
 *
 * transA = 0: use A as-is
 * transA = 1: use A'
 * transB = 0: use B as-is
 * transB = 1: use B'
 *
 * If transA=0, transB=0: A is M x K, B is K x N, C is M x N
 */
void blas_dgemm(int transA, int transB,
                ST_int M, ST_int N, ST_int K,
                ST_double alpha, const ST_double *A, ST_int lda,
                const ST_double *B, ST_int ldb,
                ST_double beta, ST_double *C, ST_int ldc);

/*
 * Symmetric rank-k update: C = alpha * A * A' + beta * C
 * A is N x K, C is N x N symmetric (column-major, lower triangle)
 *
 * trans = 0: C = alpha * A * A' + beta * C  (A is N x K)
 * trans = 1: C = alpha * A' * A + beta * C  (A is K x N)
 */
void blas_dsyrk(int trans, ST_int N, ST_int K,
                ST_double alpha, const ST_double *A, ST_int lda,
                ST_double beta, ST_double *C, ST_int ldc);


/* ============================================================================
 * LAPACK Operations
 * ============================================================================ */

/*
 * Cholesky factorization: A = L * L'
 * A is N x N symmetric positive definite (column-major)
 * On output, L is stored in lower triangle of A
 *
 * Returns: 0 on success, >0 if matrix not positive definite
 */
ST_int lapack_dpotrf(ST_int N, ST_double *A, ST_int lda);

/*
 * Solve using Cholesky factorization: A * X = B
 * A is N x N with Cholesky factor from dpotrf (column-major)
 * B is N x NRHS (column-major), solution X stored in B
 *
 * Returns: 0 on success
 */
ST_int lapack_dpotrs(ST_int N, ST_int NRHS,
                     const ST_double *A, ST_int lda,
                     ST_double *B, ST_int ldb);

/*
 * Compute inverse from Cholesky factorization
 * A is N x N with Cholesky factor from dpotrf (column-major)
 * On output, A contains A^{-1} (full symmetric matrix)
 *
 * Returns: 0 on success
 */
ST_int lapack_dpotri(ST_int N, ST_double *A, ST_int lda);


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

/*
 * Compute X' * v for column-major X
 * X is N x K, v is N, result is K
 */
void blas_xtv(ST_double *result,
              const ST_double *X, ST_int N, ST_int K,
              const ST_double *v);

/*
 * Solve (X'DX) * beta = X'y via Cholesky
 * X is N x K, D is N (diagonal), y is N, beta is K (output)
 *
 * Returns: 0 on success, <0 on error
 */
ST_int blas_solve_weighted_ls(const ST_double *X, ST_int N, ST_int K,
                               const ST_double *D,
                               const ST_double *y,
                               ST_double *beta);


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


/* ============================================================================
 * Utility
 * ============================================================================ */

/*
 * Get the name of the BLAS implementation being used
 */
const char* blas_get_impl_name(void);


#endif /* CQREG_BLAS_H */
