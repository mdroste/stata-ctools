/*
 * ctools_ols.h
 *
 * Shared OLS utilities: Cholesky decomposition, matrix inversion
 * Used by both creghdfe and civreghdfe
 * Part of the ctools Stata plugin suite
 */

#ifndef CTOOLS_OLS_H
#define CTOOLS_OLS_H

#include "stplugin.h"
#include "ctools_config.h"

/* ========================================================================
 * Cholesky decomposition: A = L * L' (in-place, returns L in lower triangle)
 * Returns 0 on success, -1 if not positive definite
 * ======================================================================== */
ST_int ctools_cholesky(ST_double *A, ST_int n);

/* ========================================================================
 * Matrix inversion via Cholesky (for positive definite matrices)
 * Input: L (lower triangular Cholesky factor)
 * Output: inv (inverse of L*L')
 * Returns 0 on success, -1 on memory allocation failure
 * ======================================================================== */
ST_int ctools_invert_from_cholesky(const ST_double *L, ST_int n, ST_double *inv);

/* ========================================================================
 * Combined solve: compute inverse of SPD matrix A using Cholesky
 * A_inv = (A)^{-1}
 * Returns 0 on success, -1 on failure
 * ======================================================================== */
ST_int ctools_invert_spd(ST_double *A, ST_int n, ST_double *A_inv);

/* ========================================================================
 * Solve linear system Ax = b using Cholesky decomposition
 * A must be symmetric positive definite
 * Returns 0 on success, -1 on failure
 * ======================================================================== */
ST_int ctools_solve_cholesky(const ST_double *A, const ST_double *b, ST_int n, ST_double *x);

#endif /* CTOOLS_OLS_H */
