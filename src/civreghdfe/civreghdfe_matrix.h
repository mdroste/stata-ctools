/*
    civreghdfe_matrix.h
    Matrix operations for IV regression

    Helper functions for matrix multiplication, decomposition,
    eigenvalue computation, and linear algebra operations.
*/

#ifndef CIVREGHDFE_MATRIX_H
#define CIVREGHDFE_MATRIX_H

#include "../stplugin.h"
#include "../ctools_matrix.h"

/*
    HAC Kernel weight functions for Newey-West style standard errors.

    Kernel types:
    0 = None (not HAC)
    1 = Bartlett / Newey-West (triangular)
    2 = Parzen
    3 = Quadratic Spectral
    4 = Truncated
    5 = Tukey-Hanning

    Returns kernel weight for lag j given bandwidth bw.
*/
ST_double civreghdfe_kernel_weight(ST_int kernel_type, ST_int j, ST_int bw);

/*
    Jacobi eigenvalue algorithm for symmetric matrices.
    Computes all eigenvalues of a K x K symmetric matrix A.
    If eigenvalues is NULL, returns just the minimum eigenvalue.
    If eigenvalues is provided, stores all K eigenvalues in it (sorted ascending).

    A is stored in column-major order and is NOT modified.
    Returns the minimum eigenvalue.
*/
ST_double civreghdfe_jacobi_eigenvalues(const ST_double *A, ST_int K, ST_double *eigenvalues);

/*
    Compute matrix square root inverse: M^(-1/2)
    For a symmetric positive definite matrix M, compute M^(-1/2)
    using eigendecomposition.

    M is K x K (column-major), not modified.
    Result is stored in Minv_sqrt (K x K, column-major)
    Returns 0 on success, -1 if M is not positive definite.
*/
ST_int civreghdfe_matpowersym_neg_half(const ST_double *M, ST_int K, ST_double *Minv_sqrt);

/*
    Compute LIML lambda (minimum eigenvalue for k-class estimation)

    Y = [y, X_endog]  (N x (1 + K_endog))
    Z = all instruments (N x K_iv)
    X_exog = included instruments / exogenous regressors (N x K_exog)

    Computes:
    QWW = Y'(I - Z(Z'Z)^-1 Z')Y = Y'M_Z Y
    QWW1 = Y'(I - X_exog(X_exog'X_exog)^-1 X_exog')Y = Y'M_{X_exog} Y
    lambda = min eigenvalue of QWW^(-1/2) * QWW1 * QWW^(-1/2)

    Returns lambda (1.0 on error, equivalent to 2SLS).
*/
ST_double civreghdfe_compute_liml_lambda(
    const ST_double *y,
    const ST_double *X_endog,
    const ST_double *X_exog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv);

#endif /* CIVREGHDFE_MATRIX_H */
