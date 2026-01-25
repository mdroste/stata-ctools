/*
 * ctools_ols.c
 *
 * Shared OLS utilities: Cholesky decomposition, matrix inversion
 * Used by both creghdfe and civreghdfe
 * Part of the ctools Stata plugin suite
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ctools_ols.h"
#include "ctools_config.h"

/* ========================================================================
 * Cholesky decomposition: A = L * L' (in-place, returns L in lower triangle)
 * ======================================================================== */

ST_int ctools_cholesky(ST_double *A, ST_int n)
{
    ST_int i, j, k;

    for (j = 0; j < n; j++) {
        ST_double sum = A[j * n + j];
        for (k = 0; k < j; k++) {
            sum -= A[j * n + k] * A[j * n + k];
        }
        if (sum <= 0.0) return -1;  /* Not positive definite */
        A[j * n + j] = sqrt(sum);

        for (i = j + 1; i < n; i++) {
            sum = A[i * n + j];
            for (k = 0; k < j; k++) {
                sum -= A[i * n + k] * A[j * n + k];
            }
            A[i * n + j] = sum / A[j * n + j];
        }
    }

    /* Zero upper triangle */
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            A[i * n + j] = 0.0;
        }
    }

    return 0;
}

/* ========================================================================
 * Matrix inversion via Cholesky (for positive definite matrices)
 * Input: L (lower triangular Cholesky factor)
 * Output: inv (inverse of L*L')
 * ======================================================================== */

ST_int ctools_invert_from_cholesky(const ST_double *L, ST_int n, ST_double *inv)
{
    ST_int i, j, k;
    ST_double *L_inv = (ST_double *)ctools_safe_calloc3((size_t)n, (size_t)n, sizeof(ST_double));

    if (!L_inv) return -1;

    /* Invert L (lower triangular) */
    for (i = 0; i < n; i++) {
        L_inv[i * n + i] = 1.0 / L[i * n + i];
        for (j = i + 1; j < n; j++) {
            ST_double sum = 0.0;
            for (k = i; k < j; k++) {
                sum -= L[j * n + k] * L_inv[k * n + i];
            }
            L_inv[j * n + i] = sum / L[j * n + j];
        }
    }

    /* Compute inv = L_inv' * L_inv */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            inv[i * n + j] = 0.0;
            for (k = (i > j ? i : j); k < n; k++) {
                inv[i * n + j] += L_inv[k * n + i] * L_inv[k * n + j];
            }
        }
    }

    free(L_inv);
    return 0;
}

/* ========================================================================
 * Combined solve: compute inverse of SPD matrix A using Cholesky
 * ======================================================================== */

ST_int ctools_invert_spd(ST_double *A, ST_int n, ST_double *A_inv)
{
    ST_double *L = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    if (!L) return -1;

    memcpy(L, A, (size_t)n * (size_t)n * sizeof(ST_double));

    if (ctools_cholesky(L, n) != 0) {
        free(L);
        return -1;
    }

    ST_int rc = ctools_invert_from_cholesky(L, n, A_inv);
    free(L);
    return rc;
}

/* ========================================================================
 * Solve linear system Ax = b using Cholesky decomposition
 * ======================================================================== */

ST_int ctools_solve_cholesky(const ST_double *A, const ST_double *b, ST_int n, ST_double *x)
{
    ST_int i, j;
    ST_double *L = (ST_double *)ctools_safe_malloc3((size_t)n, (size_t)n, sizeof(ST_double));
    ST_double *y = (ST_double *)ctools_safe_malloc2((size_t)n, sizeof(ST_double));

    if (!L || !y) {
        free(L);
        free(y);
        return -1;
    }

    memcpy(L, A, (size_t)n * (size_t)n * sizeof(ST_double));

    if (ctools_cholesky(L, n) != 0) {
        free(L);
        free(y);
        return -1;
    }

    /* Forward substitution: solve Ly = b */
    for (i = 0; i < n; i++) {
        ST_double sum = b[i];
        for (j = 0; j < i; j++) {
            sum -= L[i * n + j] * y[j];
        }
        y[i] = sum / L[i * n + i];
    }

    /* Backward substitution: solve L'x = y */
    for (i = n - 1; i >= 0; i--) {
        ST_double sum = y[i];
        for (j = i + 1; j < n; j++) {
            sum -= L[j * n + i] * x[j];
        }
        x[i] = sum / L[i * n + i];
    }

    free(L);
    free(y);
    return 0;
}
