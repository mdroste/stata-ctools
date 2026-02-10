/*
    civreghdfe_matrix.c
    Matrix operations for IV regression

    Helper functions for matrix multiplication, decomposition,
    eigenvalue computation, and linear algebra operations.

    OPTIMIZED: Uses ctools_dot_unrolled for K-way unrolled dot products,
    OpenMP parallelization, and restrict pointers for better performance.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "civreghdfe_matrix.h"
#include "../ctools_config.h"
#include "../ctools_ols.h"

/* Use shared Cholesky functions from ctools_ols */
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/*
    HAC Kernel weight functions

    Kernel weights for HAC (Newey-West style) estimators.
    bw = bandwidth parameter
    j = lag index

    Kernels are classified as:
    - "lag" windows (Bartlett, Parzen, Truncated, Tukey-Hanning): truncated at bw
    - "spectral" windows (Quadratic Spectral): weights computed for all lags

    Following ivreg2/livreg2's m_calckw() function exactly.
*/
ST_double civreghdfe_kernel_weight(ST_int kernel_type, ST_int j, ST_int bw)
{
    if (j == 0) return 1.0;
    if (bw <= 0) return 0.0;

    /* karg = tau / bw, following ivreg2's m_calckw() */
    ST_double x = (ST_double)j / (ST_double)bw;

    switch (kernel_type) {
        case 1:  /* Bartlett / Newey-West - lag window, truncate at bw */
            /* kw = 1 - karg */
            if (j > bw) return 0.0;
            return 1.0 - x;

        case 2:  /* Parzen - lag window, truncate at bw */
            if (j > bw) return 0.0;
            if (x <= 0.5) {
                return 1.0 - 6.0 * x * x + 6.0 * x * x * x;
            } else {
                ST_double t = 1.0 - x;
                return 2.0 * t * t * t;
            }

        case 3:  /* Quadratic Spectral - spectral window, no truncation */
            /* kw = 25/(12*pi^2*karg^2) * (sin(6*pi*karg/5)/(6*pi*karg/5) - cos(6*pi*karg/5)) */
            {
                ST_double pi = 3.14159265358979323846;
                ST_double z = 6.0 * pi * x / 5.0;
                if (fabs(z) < 1e-10) return 1.0;
                return 25.0 / (12.0 * pi * pi * x * x) *
                       (sin(z) / z - cos(z));
            }

        case 4:  /* Truncated - lag window, truncate at bw */
            if (j > bw) return 0.0;
            return 1.0;

        case 5:  /* Tukey-Hanning - lag window, truncate at bw */
            if (j > bw) return 0.0;
            {
                ST_double pi = 3.14159265358979323846;
                return 0.5 + 0.5 * cos(pi * x);
            }

        default:
            return 0.0;
    }
}

/*
    Jacobi eigenvalue algorithm for symmetric matrices.
*/
ST_double civreghdfe_jacobi_eigenvalues(const ST_double *A, ST_int K, ST_double *eigenvalues)
{
    const ST_int max_iter = 100;
    const ST_double tol = 1e-12;
    ST_int iter, p, q, i;

    /* Work with a copy */
    ST_double *D = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));
    if (!D) {
        /* Return identity eigenvalue on allocation failure */
        if (eigenvalues) {
            for (i = 0; i < K; i++) eigenvalues[i] = 1.0;
        }
        return 1.0;
    }
    memcpy(D, A, (size_t)K * K * sizeof(ST_double));

    for (iter = 0; iter < max_iter; iter++) {
        /* Find largest off-diagonal element */
        ST_double max_off = 0.0;
        p = 0; q = 1;
        for (i = 0; i < K; i++) {
            for (ST_int j = i + 1; j < K; j++) {
                ST_double val = fabs(D[j * K + i]);
                if (val > max_off) {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }

        /* Check convergence */
        if (max_off < tol) break;

        /* Compute rotation */
        ST_double app = D[p * K + p];
        ST_double aqq = D[q * K + q];
        ST_double apq = D[q * K + p];

        ST_double theta = 0.5 * atan2(2.0 * apq, aqq - app);
        ST_double c = cos(theta);
        ST_double s = sin(theta);

        /* Apply rotation */
        for (i = 0; i < K; i++) {
            if (i != p && i != q) {
                ST_double dip = D[p * K + i];
                ST_double diq = D[q * K + i];
                D[p * K + i] = c * dip - s * diq;
                D[i * K + p] = D[p * K + i];
                D[q * K + i] = s * dip + c * diq;
                D[i * K + q] = D[q * K + i];
            }
        }

        D[p * K + p] = c * c * app - 2 * s * c * apq + s * s * aqq;
        D[q * K + q] = s * s * app + 2 * s * c * apq + c * c * aqq;
        D[q * K + p] = 0.0;
        D[p * K + q] = 0.0;
    }

    /* Extract eigenvalues (diagonal elements) */
    ST_double min_eval = D[0];
    if (eigenvalues) {
        for (i = 0; i < K; i++) {
            eigenvalues[i] = D[i * K + i];
        }
        /* Sort eigenvalues in ascending order (simple bubble sort for small K) */
        for (i = 0; i < K - 1; i++) {
            for (ST_int j = i + 1; j < K; j++) {
                if (eigenvalues[j] < eigenvalues[i]) {
                    ST_double tmp = eigenvalues[i];
                    eigenvalues[i] = eigenvalues[j];
                    eigenvalues[j] = tmp;
                }
            }
        }
        min_eval = eigenvalues[0];
    } else {
        for (i = 1; i < K; i++) {
            if (D[i * K + i] < min_eval) {
                min_eval = D[i * K + i];
            }
        }
    }

    free(D);
    return min_eval;
}

/*
    Returns just the minimum eigenvalue of a symmetric matrix.
*/
ST_double civreghdfe_jacobi_min_eigenvalue(const ST_double *A, ST_int K)
{
    return civreghdfe_jacobi_eigenvalues(A, K, NULL);
}

/*
    Compute matrix square root inverse: M^(-1/2)
*/
ST_int civreghdfe_matpowersym_neg_half(const ST_double *M, ST_int K, ST_double *Minv_sqrt)
{
    const ST_int max_iter = 100;
    const ST_double tol = 1e-12;
    ST_int iter, p, q, i, j;

    /* Work matrices */
    ST_double *D = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));
    ST_double *V = (ST_double *)malloc((size_t)K * K * sizeof(ST_double));

    if (!D || !V) {
        free(D); free(V);
        return -1;  /* Allocation failure */
    }

    memcpy(D, M, (size_t)K * K * sizeof(ST_double));

    /* Initialize V as identity */
    memset(V, 0, (size_t)K * K * sizeof(ST_double));
    for (i = 0; i < K; i++) V[i * K + i] = 1.0;

    /* Jacobi iteration to diagonalize D, accumulate eigenvectors in V */
    for (iter = 0; iter < max_iter; iter++) {
        ST_double max_off = 0.0;
        p = 0; q = 1;
        for (i = 0; i < K; i++) {
            for (j = i + 1; j < K; j++) {
                ST_double val = fabs(D[j * K + i]);
                if (val > max_off) {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }

        if (max_off < tol) break;

        ST_double app = D[p * K + p];
        ST_double aqq = D[q * K + q];
        ST_double apq = D[q * K + p];

        ST_double theta = 0.5 * atan2(2.0 * apq, aqq - app);
        ST_double c = cos(theta);
        ST_double s = sin(theta);

        /* Apply rotation to D */
        for (i = 0; i < K; i++) {
            if (i != p && i != q) {
                ST_double dip = D[p * K + i];
                ST_double diq = D[q * K + i];
                D[p * K + i] = c * dip - s * diq;
                D[i * K + p] = D[p * K + i];
                D[q * K + i] = s * dip + c * diq;
                D[i * K + q] = D[q * K + i];
            }
        }
        D[p * K + p] = c * c * app - 2 * s * c * apq + s * s * aqq;
        D[q * K + q] = s * s * app + 2 * s * c * apq + c * c * aqq;
        D[q * K + p] = 0.0;
        D[p * K + q] = 0.0;

        /* Apply rotation to V */
        for (i = 0; i < K; i++) {
            ST_double vip = V[p * K + i];
            ST_double viq = V[q * K + i];
            V[p * K + i] = c * vip - s * viq;
            V[q * K + i] = s * vip + c * viq;
        }
    }

    /* Compute V * diag(1/sqrt(eigenvalues)) * V' */
    ST_double *VD = (ST_double *)malloc(K * K * sizeof(ST_double));
    if (!VD) {
        free(D); free(V);
        return -1;  /* Allocation failure */
    }
    for (j = 0; j < K; j++) {
        ST_double eval = D[j * K + j];
        if (eval <= 0) {
            free(D); free(V); free(VD);
            return -1;  /* Not positive definite */
        }
        ST_double inv_sqrt_eval = 1.0 / sqrt(eval);
        for (i = 0; i < K; i++) {
            VD[j * K + i] = V[j * K + i] * inv_sqrt_eval;
        }
    }

    /* Compute VD * V' */
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (ST_int k = 0; k < K; k++) {
                sum += VD[k * K + i] * V[k * K + j];
            }
            Minv_sqrt[j * K + i] = sum;
        }
    }

    free(D);
    free(V);
    free(VD);
    return 0;
}

/*
    Compute LIML lambda (minimum eigenvalue for k-class)
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
    ST_int K_iv)
{
    ST_int K_Y = 1 + K_endog;
    ST_int i, j;
    int use_weights = (weights != NULL && weight_type != 0);

    /* Initialize pointers to NULL for safe cleanup */
    ST_double *Y = NULL, *ZtZ = NULL, *ZtZ_inv = NULL, *ZtZ_L = NULL;
    ST_double *ZtY = NULL, *YtY = NULL, *ZtZ_inv_ZtY = NULL;
    ST_double *QWW = NULL, *ZtY_t_ZtZ_inv_ZtY = NULL;

    /* Allocate Y = [y, X_endog] */
    Y = (ST_double *)malloc((size_t)N * K_Y * sizeof(ST_double));
    if (!Y) return 1.0;
    memcpy(Y, y, (size_t)N * sizeof(ST_double));
    if (K_endog > 0) {
        memcpy(Y + N, X_endog, (size_t)N * K_endog * sizeof(ST_double));
    }

    /* Compute Z'Z and its inverse */
    ZtZ = (ST_double *)calloc((size_t)K_iv * K_iv, sizeof(ST_double));
    ZtZ_inv = (ST_double *)calloc((size_t)K_iv * K_iv, sizeof(ST_double));
    if (!ZtZ || !ZtZ_inv) {
        free(Y); free(ZtZ); free(ZtZ_inv);
        return 1.0;
    }
    if (use_weights)
        ctools_matmul_atdb(Z, Z, weights, N, K_iv, K_iv, ZtZ);
    else
        ctools_matmul_atb(Z, Z, N, K_iv, K_iv, ZtZ);

    /* Invert Z'Z using Cholesky */
    ZtZ_L = (ST_double *)malloc((size_t)K_iv * K_iv * sizeof(ST_double));
    if (!ZtZ_L) {
        free(Y); free(ZtZ); free(ZtZ_inv);
        return 1.0;
    }
    memcpy(ZtZ_L, ZtZ, (size_t)K_iv * K_iv * sizeof(ST_double));
    if (cholesky(ZtZ_L, K_iv) != 0) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtZ_L);
        return 1.0;
    }
    invert_from_cholesky(ZtZ_L, K_iv, ZtZ_inv);
    free(ZtZ_L);
    ZtZ_L = NULL;

    /* Compute Z'Y */
    ZtY = (ST_double *)calloc((size_t)K_iv * K_Y, sizeof(ST_double));
    if (!ZtY) {
        free(Y); free(ZtZ); free(ZtZ_inv);
        return 1.0;
    }
    if (use_weights)
        ctools_matmul_atdb(Z, Y, weights, N, K_iv, K_Y, ZtY);
    else
        ctools_matmul_atb(Z, Y, N, K_iv, K_Y, ZtY);

    /* Compute Y'Y */
    YtY = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    if (!YtY) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY);
        return 1.0;
    }
    if (use_weights)
        ctools_matmul_atdb(Y, Y, weights, N, K_Y, K_Y, YtY);
    else
        ctools_matmul_atb(Y, Y, N, K_Y, K_Y, YtY);

    /* Compute (Z'Z)^-1 * Z'Y */
    ZtZ_inv_ZtY = (ST_double *)calloc((size_t)K_iv * K_Y, sizeof(ST_double));
    if (!ZtZ_inv_ZtY) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
        return 1.0;
    }
    ctools_matmul_ab(ZtZ_inv, ZtY, K_iv, K_iv, K_Y, ZtZ_inv_ZtY);

    /* Compute QWW = Y'Y - (Z'Y)'(Z'Z)^-1(Z'Y) = Y'M_Z Y */
    QWW = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    ZtY_t_ZtZ_inv_ZtY = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    if (!QWW || !ZtY_t_ZtZ_inv_ZtY) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
        free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY);
        return 1.0;
    }

    for (i = 0; i < K_Y; i++) {
        for (j = 0; j < K_Y; j++) {
            ST_double sum = 0.0;
            for (ST_int k = 0; k < K_iv; k++) {
                sum += ZtY[i * K_iv + k] * ZtZ_inv_ZtY[j * K_iv + k];
            }
            ZtY_t_ZtZ_inv_ZtY[j * K_Y + i] = sum;
        }
    }

    for (i = 0; i < K_Y * K_Y; i++) {
        QWW[i] = YtY[i] - ZtY_t_ZtZ_inv_ZtY[i];
    }

    /* Compute QWW1 = Y'M_{X_exog} Y */
    ST_double *QWW1 = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    if (!QWW1) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
        free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY);
        return 1.0;
    }

    if (K_exog > 0) {
        ST_double *Z2tZ2 = (ST_double *)calloc((size_t)K_exog * K_exog, sizeof(ST_double));
        ST_double *Z2tZ2_inv = (ST_double *)calloc((size_t)K_exog * K_exog, sizeof(ST_double));
        ST_double *Z2tZ2_L = (ST_double *)malloc((size_t)K_exog * K_exog * sizeof(ST_double));
        if (!Z2tZ2 || !Z2tZ2_inv || !Z2tZ2_L) {
            free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
            free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY); free(QWW1);
            free(Z2tZ2); free(Z2tZ2_inv); free(Z2tZ2_L);
            return 1.0;
        }
        if (use_weights)
            ctools_matmul_atdb(X_exog, X_exog, weights, N, K_exog, K_exog, Z2tZ2);
        else
            ctools_matmul_atb(X_exog, X_exog, N, K_exog, K_exog, Z2tZ2);

        memcpy(Z2tZ2_L, Z2tZ2, (size_t)K_exog * K_exog * sizeof(ST_double));
        if (cholesky(Z2tZ2_L, K_exog) != 0) {
            free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
            free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY); free(QWW1);
            free(Z2tZ2); free(Z2tZ2_inv); free(Z2tZ2_L);
            return 1.0;
        }
        invert_from_cholesky(Z2tZ2_L, K_exog, Z2tZ2_inv);
        free(Z2tZ2_L);

        ST_double *Z2tY = (ST_double *)calloc((size_t)K_exog * K_Y, sizeof(ST_double));
        ST_double *Z2tZ2_inv_Z2tY = (ST_double *)calloc((size_t)K_exog * K_Y, sizeof(ST_double));
        ST_double *Z2tY_t_Z2tZ2_inv_Z2tY = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
        if (!Z2tY || !Z2tZ2_inv_Z2tY || !Z2tY_t_Z2tZ2_inv_Z2tY) {
            free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
            free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY); free(QWW1);
            free(Z2tZ2); free(Z2tZ2_inv);
            free(Z2tY); free(Z2tZ2_inv_Z2tY); free(Z2tY_t_Z2tZ2_inv_Z2tY);
            return 1.0;
        }
        if (use_weights)
            ctools_matmul_atdb(X_exog, Y, weights, N, K_exog, K_Y, Z2tY);
        else
            ctools_matmul_atb(X_exog, Y, N, K_exog, K_Y, Z2tY);
        ctools_matmul_ab(Z2tZ2_inv, Z2tY, K_exog, K_exog, K_Y, Z2tZ2_inv_Z2tY);

        for (i = 0; i < K_Y; i++) {
            for (j = 0; j < K_Y; j++) {
                ST_double sum = 0.0;
                for (ST_int k = 0; k < K_exog; k++) {
                    sum += Z2tY[i * K_exog + k] * Z2tZ2_inv_Z2tY[j * K_exog + k];
                }
                Z2tY_t_Z2tZ2_inv_Z2tY[j * K_Y + i] = sum;
            }
        }

        for (i = 0; i < K_Y * K_Y; i++) {
            QWW1[i] = YtY[i] - Z2tY_t_Z2tZ2_inv_Z2tY[i];
        }

        free(Z2tZ2); free(Z2tZ2_inv); free(Z2tY);
        free(Z2tZ2_inv_Z2tY); free(Z2tY_t_Z2tZ2_inv_Z2tY);
    }
    else {
        memcpy(QWW1, YtY, (size_t)K_Y * K_Y * sizeof(ST_double));
    }

    /* Compute lambda = min eigenvalue of QWW^(-1/2) * QWW1 * QWW^(-1/2) */
    ST_double *QWW_inv_sqrt = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    ST_double *temp1 = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    ST_double *M_for_eigen = (ST_double *)calloc((size_t)K_Y * K_Y, sizeof(ST_double));
    if (!QWW_inv_sqrt || !temp1 || !M_for_eigen) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
        free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY);
        free(QWW1); free(QWW_inv_sqrt); free(temp1); free(M_for_eigen);
        return 1.0;
    }

    if (civreghdfe_matpowersym_neg_half(QWW, K_Y, QWW_inv_sqrt) != 0) {
        free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
        free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY);
        free(QWW1); free(QWW_inv_sqrt); free(temp1); free(M_for_eigen);
        return 1.0;
    }

    ctools_matmul_ab(QWW_inv_sqrt, QWW1, K_Y, K_Y, K_Y, temp1);
    ctools_matmul_ab(temp1, QWW_inv_sqrt, K_Y, K_Y, K_Y, M_for_eigen);

    ST_double lambda = civreghdfe_jacobi_min_eigenvalue(M_for_eigen, K_Y);

    /* Cleanup */
    free(Y); free(ZtZ); free(ZtZ_inv); free(ZtY); free(YtY);
    free(ZtZ_inv_ZtY); free(QWW); free(ZtY_t_ZtZ_inv_ZtY);
    free(QWW1); free(QWW_inv_sqrt); free(temp1); free(M_for_eigen);

    return lambda;
}
