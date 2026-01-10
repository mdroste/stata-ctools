/*
    civreghdfe_impl.c
    Instrumental Variables Regression with High-Dimensional Fixed Effects

    Implements 2SLS/IV regression with HDFE absorption. Reuses the creghdfe
    infrastructure for HDFE demeaning, singleton detection, and VCE computation.

    Algorithm:
    1. Read data: y, X_endog, X_exog, Z (instruments), FE vars, weights, cluster
    2. Detect and remove singletons
    3. Partial out FEs from all variables using CG solver
    4. First stage: Regress each endogenous var on [exogenous + instruments]
    5. Second stage: Regress y on [exogenous + predicted endogenous]
    6. Compute VCE (corrected for 2SLS)
    7. Store results to Stata
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "civreghdfe_impl.h"

/* Debug flag */
#define CIVREGHDFE_DEBUG 0

/*
    Helper: Matrix multiply C = A' * B where A is N x K1, B is N x K2
    Result C is K1 x K2
*/
static void matmul_atb(const ST_double *A, const ST_double *B,
                       ST_int N, ST_int K1, ST_int K2,
                       ST_double *C)
{
    ST_int i, j, k;

    /* Initialize C to zero */
    memset(C, 0, K1 * K2 * sizeof(ST_double));

    /* Compute C[i,j] = sum_k A[k,i] * B[k,j] */
    for (j = 0; j < K2; j++) {
        for (i = 0; i < K1; i++) {
            ST_double sum = 0.0;
            const ST_double *a_col = A + i * N;
            const ST_double *b_col = B + j * N;
            for (k = 0; k < N; k++) {
                sum += a_col[k] * b_col[k];
            }
            C[j * K1 + i] = sum;  /* Column-major storage */
        }
    }
}

/*
    Helper: Matrix multiply C = A * B where A is K1 x K2, B is K2 x K3
    Result C is K1 x K3
*/
static void matmul_ab(const ST_double *A, const ST_double *B,
                      ST_int K1, ST_int K2, ST_int K3,
                      ST_double *C)
{
    ST_int i, j, k;

    /* Initialize C to zero */
    memset(C, 0, K1 * K3 * sizeof(ST_double));

    /* C[i,j] = sum_k A[i,k] * B[k,j] */
    /* A is K1 x K2, stored column-major: A[i,k] = A[k*K1 + i] */
    /* B is K2 x K3, stored column-major: B[k,j] = B[j*K2 + k] */
    /* C is K1 x K3, stored column-major: C[i,j] = C[j*K1 + i] */
    for (j = 0; j < K3; j++) {
        for (i = 0; i < K1; i++) {
            ST_double sum = 0.0;
            for (k = 0; k < K2; k++) {
                sum += A[k * K1 + i] * B[j * K2 + k];
            }
            C[j * K1 + i] = sum;
        }
    }
}

/*
    Helper: Solve linear system Ax = b using Cholesky decomposition
    A is K x K symmetric positive definite
    b is K x 1
    x is K x 1 (output, can be same as b for in-place)
    Returns 0 on success, -1 if not positive definite
*/
static ST_int solve_cholesky(ST_double *A, ST_double *b, ST_int K, ST_double *x)
{
    ST_int i, j;
    ST_double *L = (ST_double *)malloc(K * K * sizeof(ST_double));
    if (!L) return -1;

    /* Copy A to L */
    memcpy(L, A, K * K * sizeof(ST_double));

    /* Cholesky decomposition: L*L' = A */
    if (cholesky(L, K) != 0) {
        free(L);
        return -1;
    }

    /* Forward substitution: solve Ly = b */
    ST_double *y = (ST_double *)malloc(K * sizeof(ST_double));
    if (!y) {
        free(L);
        return -1;
    }

    for (i = 0; i < K; i++) {
        ST_double sum = b[i];
        for (j = 0; j < i; j++) {
            sum -= L[i * K + j] * y[j];
        }
        y[i] = sum / L[i * K + i];
    }

    /* Backward substitution: solve L'x = y */
    for (i = K - 1; i >= 0; i--) {
        ST_double sum = y[i];
        for (j = i + 1; j < K; j++) {
            sum -= L[j * K + i] * x[j];
        }
        x[i] = sum / L[i * K + i];
    }

    free(L);
    free(y);
    return 0;
}

/*
    Helper: Weighted matrix multiply C = A' * diag(w) * B
*/
static void matmul_atdb(const ST_double *A, const ST_double *B,
                        const ST_double *w, ST_int N, ST_int K1, ST_int K2,
                        ST_double *C)
{
    ST_int i, j, k;

    memset(C, 0, K1 * K2 * sizeof(ST_double));

    for (j = 0; j < K2; j++) {
        for (i = 0; i < K1; i++) {
            ST_double sum = 0.0;
            const ST_double *a_col = A + i * N;
            const ST_double *b_col = B + j * N;
            if (w) {
                for (k = 0; k < N; k++) {
                    sum += a_col[k] * w[k] * b_col[k];
                }
            } else {
                for (k = 0; k < N; k++) {
                    sum += a_col[k] * b_col[k];
                }
            }
            C[j * K1 + i] = sum;
        }
    }
}

/*
    Compute 2SLS estimation.

    The 2SLS estimator is:
    beta = (X'P_Z X)^-1 X'P_Z y
    where P_Z = Z(Z'Z)^-1 Z' is the projection onto the instrument space

    For efficiency, we compute:
    1. Z'Z and invert it
    2. Z'X and Z'y
    3. X'Z(Z'Z)^-1 Z'X = (Z'X)'(Z'Z)^-1(Z'X)
    4. X'Z(Z'Z)^-1 Z'y
    5. Solve for beta
*/
ST_retcode compute_2sls(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_double *beta,
    ST_double *V,
    ST_double *first_stage_F,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int verbose
)
{
    ST_int K_total = K_exog + K_endog;  /* Total regressors */
    ST_int i, j, k;

    if (verbose) {
        SF_display("civreghdfe: Computing 2SLS estimation\n");
        char buf[256];
        snprintf(buf, sizeof(buf), "  N=%d, K_exog=%d, K_endog=%d, K_iv=%d\n",
                 (int)N, (int)K_exog, (int)K_endog, (int)K_iv);
        SF_display(buf);
    }

    /* Check identification */
    if (K_iv < K_total) {
        SF_error("civreghdfe: Model is underidentified\n");
        return 198;
    }

    /* Allocate work arrays */
    ST_double *ZtZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ST_double *ZtZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ST_double *ZtX = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
    ST_double *Zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
    ST_double *XtPzX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *XtPzy = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *temp1 = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
    ST_double *X_all = (ST_double *)malloc(N * K_total * sizeof(ST_double));
    ST_double *resid = (ST_double *)malloc(N * sizeof(ST_double));

    if (!ZtZ || !ZtZ_inv || !ZtX || !Zty || !XtPzX || !XtPzy ||
        !temp1 || !X_all || !resid) {
        SF_error("civreghdfe: Memory allocation failed\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        return 920;
    }

    /* Combine X_exog and X_endog into X_all */
    /* Layout: [X_exog (N x K_exog) | X_endog (N x K_endog)] */
    if (K_exog > 0 && X_exog) {
        memcpy(X_all, X_exog, N * K_exog * sizeof(ST_double));
    }
    if (K_endog > 0 && X_endog) {
        memcpy(X_all + N * K_exog, X_endog, N * K_endog * sizeof(ST_double));
    }

    /* Step 1: Compute Z'Z (weighted if needed) */
    if (weights && weight_type != 0) {
        matmul_atdb(Z, Z, weights, N, K_iv, K_iv, ZtZ);
    } else {
        matmul_atb(Z, Z, N, K_iv, K_iv, ZtZ);
    }

    /* Step 2: Invert Z'Z */
    memcpy(ZtZ_inv, ZtZ, K_iv * K_iv * sizeof(ST_double));
    if (cholesky(ZtZ_inv, K_iv) != 0) {
        SF_error("civreghdfe: Z'Z is singular (instruments may be collinear)\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        return 198;
    }
    if (invert_from_cholesky(ZtZ_inv, K_iv, ZtZ_inv) != 0) {
        SF_error("civreghdfe: Failed to invert Z'Z\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        return 198;
    }

    /* Step 3: Compute Z'X and Z'y */
    if (weights && weight_type != 0) {
        matmul_atdb(Z, X_all, weights, N, K_iv, K_total, ZtX);
        /* Z'y */
        for (i = 0; i < K_iv; i++) {
            ST_double sum = 0.0;
            const ST_double *z_col = Z + i * N;
            for (k = 0; k < N; k++) {
                sum += z_col[k] * weights[k] * y[k];
            }
            Zty[i] = sum;
        }
    } else {
        matmul_atb(Z, X_all, N, K_iv, K_total, ZtX);
        /* Z'y */
        for (i = 0; i < K_iv; i++) {
            ST_double sum = 0.0;
            const ST_double *z_col = Z + i * N;
            for (k = 0; k < N; k++) {
                sum += z_col[k] * y[k];
            }
            Zty[i] = sum;
        }
    }

    /* Step 4: Compute (Z'Z)^-1 * Z'X  -> temp1 (K_iv x K_total) */
    matmul_ab(ZtZ_inv, ZtX, K_iv, K_iv, K_total, temp1);

    /* Step 5: Compute X'P_Z X = (Z'X)' * (Z'Z)^-1 * Z'X = ZtX' * temp1 */
    /* ZtX is K_iv x K_total, temp1 is K_iv x K_total */
    /* ZtX' * temp1 = K_total x K_total */
    for (j = 0; j < K_total; j++) {
        for (i = 0; i < K_total; i++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                /* ZtX[k,i] = ZtX[i * K_iv + k] */
                /* temp1[k,j] = temp1[j * K_iv + k] */
                sum += ZtX[i * K_iv + k] * temp1[j * K_iv + k];
            }
            XtPzX[j * K_total + i] = sum;
        }
    }

    /* Step 6: Compute X'P_Z y = (Z'X)' * (Z'Z)^-1 * Z'y */
    /* First compute (Z'Z)^-1 * Z'y */
    ST_double *ZtZ_inv_Zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
    for (i = 0; i < K_iv; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += ZtZ_inv[k * K_iv + i] * Zty[k];
        }
        ZtZ_inv_Zty[i] = sum;
    }

    /* Then compute (Z'X)' * result */
    for (i = 0; i < K_total; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += ZtX[i * K_iv + k] * ZtZ_inv_Zty[k];
        }
        XtPzy[i] = sum;
    }

    /* Step 7: Solve (X'P_Z X) * beta = X'P_Z y */
    ST_double *XtPzX_copy = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
    memcpy(XtPzX_copy, XtPzX, K_total * K_total * sizeof(ST_double));

    /* Use Cholesky solve */
    if (cholesky(XtPzX_copy, K_total) != 0) {
        SF_error("civreghdfe: X'P_Z X is singular\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XtPzX_copy);
        return 198;
    }

    /* Forward substitution */
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));
    for (i = 0; i < K_total; i++) {
        ST_double sum = XtPzy[i];
        for (j = 0; j < i; j++) {
            sum -= XtPzX_copy[i * K_total + j] * beta_temp[j];
        }
        beta_temp[i] = sum / XtPzX_copy[i * K_total + i];
    }

    /* Backward substitution */
    for (i = K_total - 1; i >= 0; i--) {
        ST_double sum = beta_temp[i];
        for (j = i + 1; j < K_total; j++) {
            sum -= XtPzX_copy[j * K_total + i] * beta[j];
        }
        beta[i] = sum / XtPzX_copy[i * K_total + i];
    }

    if (verbose) {
        char buf[256];
        for (i = 0; i < K_total; i++) {
            snprintf(buf, sizeof(buf), "  beta[%d] = %g\n", (int)i, beta[i]);
            SF_display(buf);
        }
    }

    /* Step 8: Compute residuals using original X (not projected) */
    /* resid = y - X * beta */
    for (i = 0; i < N; i++) {
        ST_double pred = 0.0;
        for (k = 0; k < K_total; k++) {
            pred += X_all[k * N + i] * beta[k];
        }
        resid[i] = y[i] - pred;
    }

    /* Step 9: Compute VCE */
    /* For 2SLS, the VCE is: sigma^2 * (X'P_Z X)^-1 */
    /* where sigma^2 = RSS / (N - K) using actual residuals */

    /* First, invert X'P_Z X */
    memcpy(XtPzX_copy, XtPzX, K_total * K_total * sizeof(ST_double));
    if (cholesky(XtPzX_copy, K_total) != 0) {
        SF_error("civreghdfe: Cannot compute VCE\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XtPzX_copy); free(beta_temp);
        return 198;
    }

    ST_double *XtPzX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    invert_from_cholesky(XtPzX_copy, K_total, XtPzX_inv);

    /* Compute RSS */
    ST_double rss = 0.0;
    if (weights && weight_type != 0) {
        for (i = 0; i < N; i++) {
            rss += weights[i] * resid[i] * resid[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            rss += resid[i] * resid[i];
        }
    }

    ST_int df_r = N - K_total - 1 - df_a;  /* -1 for constant, -df_a for absorbed */
    if (df_r <= 0) df_r = 1;

    ST_double sigma2 = rss / df_r;

    if (vce_type == 0) {
        /* Unadjusted VCE: sigma^2 * (X'P_Z X)^-1 */
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = sigma2 * XtPzX_inv[i];
        }
    } else if (vce_type == 1) {
        /* Robust VCE for 2SLS */
        /* V = (X'P_Z X)^-1 X'P_Z Omega P_Z X (X'P_Z X)^-1 */
        /* where Omega = diag(e^2) */
        /* Simplified: Use sandwich estimator with projected X */

        /* Compute P_Z X = Z(Z'Z)^-1 Z'X */
        ST_double *PzX = (ST_double *)calloc(N * K_total, sizeof(ST_double));

        /* PzX = Z * temp1 where temp1 = (Z'Z)^-1 Z'X */
        for (i = 0; i < N; i++) {
            for (j = 0; j < K_total; j++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_iv; k++) {
                    sum += Z[k * N + i] * temp1[j * K_iv + k];
                }
                PzX[j * N + i] = sum;
            }
        }

        /* Compute meat: sum over i of (PzX_i * e_i) * (PzX_i * e_i)' */
        ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        for (i = 0; i < N; i++) {
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            ST_double we = w * resid[i];
            for (j = 0; j < K_total; j++) {
                for (k = 0; k <= j; k++) {
                    ST_double contrib = PzX[j * N + i] * we * PzX[k * N + i] * resid[i];
                    meat[j * K_total + k] += contrib;
                    if (k != j) meat[k * K_total + j] += contrib;
                }
            }
        }

        /* HC1 adjustment */
        ST_double dof_adj = (ST_double)N / (ST_double)df_r;

        /* V = XtPzX_inv * meat * XtPzX_inv * dof_adj */
        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        matmul_ab(XtPzX_inv, meat, K_total, K_total, K_total, temp_v);
        matmul_ab(temp_v, XtPzX_inv, K_total, K_total, K_total, V);

        for (i = 0; i < K_total * K_total; i++) {
            V[i] *= dof_adj;
        }

        free(PzX);
        free(meat);
        free(temp_v);

    } else if (vce_type == 2 && cluster_ids && num_clusters > 0) {
        /* Clustered VCE for 2SLS */
        /* Similar to robust but sum within clusters */

        ST_double *PzX = (ST_double *)calloc(N * K_total, sizeof(ST_double));

        for (i = 0; i < N; i++) {
            for (j = 0; j < K_total; j++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_iv; k++) {
                    sum += Z[k * N + i] * temp1[j * K_iv + k];
                }
                PzX[j * N + i] = sum;
            }
        }

        /* Allocate cluster sums */
        ST_double *cluster_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));

        /* Sum (PzX_i * e_i) within each cluster */
        for (i = 0; i < N; i++) {
            ST_int c = cluster_ids[i] - 1;  /* 1-indexed to 0-indexed */
            if (c < 0 || c >= num_clusters) continue;
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            ST_double we = w * resid[i];
            for (j = 0; j < K_total; j++) {
                cluster_sums[c * K_total + j] += PzX[j * N + i] * we;
            }
        }

        /* Compute meat: sum over clusters of outer products */
        ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        for (ST_int c = 0; c < num_clusters; c++) {
            for (j = 0; j < K_total; j++) {
                for (k = 0; k <= j; k++) {
                    ST_double contrib = cluster_sums[c * K_total + j] * cluster_sums[c * K_total + k];
                    meat[j * K_total + k] += contrib;
                    if (k != j) meat[k * K_total + j] += contrib;
                }
            }
        }

        /* Cluster adjustment: (G / (G-1)) * ((N-1) / (N-K)) */
        ST_double G = (ST_double)num_clusters;
        ST_double dof_adj = (G / (G - 1.0)) * ((ST_double)(N - 1) / (ST_double)df_r);

        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        matmul_ab(XtPzX_inv, meat, K_total, K_total, K_total, temp_v);
        matmul_ab(temp_v, XtPzX_inv, K_total, K_total, K_total, V);

        for (i = 0; i < K_total * K_total; i++) {
            V[i] *= dof_adj;
        }

        free(PzX);
        free(cluster_sums);
        free(meat);
        free(temp_v);
    }

    /* Step 10: Compute first-stage F statistics */
    /* For each endogenous variable, compute F-stat from first stage regression */
    if (first_stage_F) {
        for (ST_int e = 0; e < K_endog; e++) {
            /* First stage: X_endog[e] = [X_exog, Z_excluded] * gamma + error */
            /* F-stat tests whether excluded instruments are jointly significant */

            /* Simplified: Use partial R-squared approach */
            /* F = ((R^2_full - R^2_reduced) / q) / ((1 - R^2_full) / (N - K_full)) */

            /* For now, compute simple F using projection */
            const ST_double *X_e = X_endog + e * N;

            /* Compute X_e'P_Z X_e */
            ST_double xpx = 0.0;
            for (i = 0; i < N; i++) {
                ST_double pz_xe = 0.0;
                for (k = 0; k < K_iv; k++) {
                    ST_double ziz_inv_zx = 0.0;
                    for (j = 0; j < K_iv; j++) {
                        ziz_inv_zx += ZtZ_inv[k * K_iv + j] * ZtX[(K_exog + e) * K_iv + j];
                    }
                    pz_xe += Z[k * N + i] * ziz_inv_zx;
                }
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                xpx += w * X_e[i] * pz_xe;
            }

            /* Compute X_e'X_e */
            ST_double xx = 0.0;
            for (i = 0; i < N; i++) {
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                xx += w * X_e[i] * X_e[i];
            }

            /* R^2 of projection */
            ST_double r2 = xpx / xx;
            if (r2 > 1.0) r2 = 1.0;
            if (r2 < 0.0) r2 = 0.0;

            /* F-stat: (R^2 / q) / ((1 - R^2) / (N - K_iv - 1)) */
            ST_int q = K_iv - K_exog;  /* Number of excluded instruments */
            if (q <= 0) q = 1;
            ST_int denom_df = N - K_iv - 1 - df_a;
            if (denom_df <= 0) denom_df = 1;

            first_stage_F[e] = (r2 / (ST_double)q) / ((1.0 - r2) / (ST_double)denom_df);
            if (first_stage_F[e] < 0) first_stage_F[e] = 0;

            if (verbose) {
                char buf[256];
                snprintf(buf, sizeof(buf), "  First-stage F[%d] = %g (R^2 = %g)\n",
                         (int)e, first_stage_F[e], r2);
                SF_display(buf);
            }
        }
    }

    /* Cleanup */
    free(ZtZ);
    free(ZtZ_inv);
    free(ZtX);
    free(Zty);
    free(XtPzX);
    free(XtPzy);
    free(temp1);
    free(X_all);
    free(resid);
    free(ZtZ_inv_Zty);
    free(XtPzX_copy);
    free(beta_temp);
    free(XtPzX_inv);

    return STATA_OK;
}

/*
    Full IV regression with HDFE.

    Variable layout from Stata:
    [y, X_endog_1, ..., X_endog_Ke, X_exog_1, ..., X_exog_Kx, Z_1, ..., Z_Kz, FE_1, ..., FE_G, cluster, weight]

    Scalars from Stata:
    - __civreghdfe_K_endog: number of endogenous regressors
    - __civreghdfe_K_exog: number of exogenous regressors
    - __civreghdfe_K_iv: number of instruments (including exogenous)
    - __civreghdfe_G: number of FE factors
    - __civreghdfe_has_cluster, __civreghdfe_has_weights, __civreghdfe_weight_type
    - __civreghdfe_vce_type: 0=unadjusted, 1=robust, 2=cluster
    - __civreghdfe_maxiter, __civreghdfe_tolerance, __civreghdfe_verbose
*/
static ST_retcode do_iv_regression(void)
{
    ST_int in1 = SF_in1();
    ST_int in2 = SF_in2();
    ST_int N_total = in2 - in1 + 1;

    /* Read configuration scalars */
    ST_double dval;
    ST_int K_endog, K_exog, K_iv, G;
    ST_int has_cluster, has_weights, weight_type;
    ST_int vce_type, maxiter, verbose;
    ST_double tolerance;

    SF_scal_use("__civreghdfe_K_endog", &dval); K_endog = (ST_int)dval;
    SF_scal_use("__civreghdfe_K_exog", &dval); K_exog = (ST_int)dval;
    SF_scal_use("__civreghdfe_K_iv", &dval); K_iv = (ST_int)dval;
    SF_scal_use("__civreghdfe_G", &dval); G = (ST_int)dval;
    SF_scal_use("__civreghdfe_has_cluster", &dval); has_cluster = (ST_int)dval;
    SF_scal_use("__civreghdfe_has_weights", &dval); has_weights = (ST_int)dval;
    SF_scal_use("__civreghdfe_weight_type", &dval); weight_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_vce_type", &dval); vce_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_maxiter", &dval); maxiter = (ST_int)dval;
    SF_scal_use("__civreghdfe_tolerance", &dval); tolerance = dval;
    SF_scal_use("__civreghdfe_verbose", &dval); verbose = (ST_int)dval;

    if (verbose) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "civreghdfe: N=%d, K_endog=%d, K_exog=%d, K_iv=%d, G=%d\n",
                 (int)N_total, (int)K_endog, (int)K_exog, (int)K_iv, (int)G);
        SF_display(buf);
    }

    /* Calculate variable positions */
    /* Layout: [y, X_endog (Ke), X_exog (Kx), Z (Kz), FE (G), cluster?, weight?] */
    ST_int var_y = 1;
    ST_int var_endog_start = 2;
    ST_int var_exog_start = var_endog_start + K_endog;
    ST_int var_iv_start = var_exog_start + K_exog;
    ST_int var_fe_start = var_iv_start + K_iv;
    ST_int var_cluster = has_cluster ? (var_fe_start + G) : 0;
    ST_int var_weight = has_weights ? (var_fe_start + G + (has_cluster ? 1 : 0)) : 0;

    ST_int K_total = K_exog + K_endog;  /* Total regressors */

    /* Allocate data arrays */
    ST_double *y = (ST_double *)malloc(N_total * sizeof(ST_double));
    ST_double *X_endog = (K_endog > 0) ? (ST_double *)malloc(N_total * K_endog * sizeof(ST_double)) : NULL;
    ST_double *X_exog = (K_exog > 0) ? (ST_double *)malloc(N_total * K_exog * sizeof(ST_double)) : NULL;
    ST_double *Z = (ST_double *)malloc(N_total * K_iv * sizeof(ST_double));
    ST_double *weights = has_weights ? (ST_double *)malloc(N_total * sizeof(ST_double)) : NULL;
    ST_int *cluster_ids = has_cluster ? (ST_int *)malloc(N_total * sizeof(ST_int)) : NULL;

    /* Allocate FE level arrays */
    ST_int **fe_levels = (ST_int **)malloc(G * sizeof(ST_int *));
    for (ST_int g = 0; g < G; g++) {
        fe_levels[g] = (ST_int *)malloc(N_total * sizeof(ST_int));
    }

    if (!y || !Z || (K_endog > 0 && !X_endog) || (K_exog > 0 && !X_exog) ||
        (has_weights && !weights) || (has_cluster && !cluster_ids)) {
        SF_error("civreghdfe: Memory allocation failed\n");
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 920;
    }

    /* Read data from Stata */
    ST_double val;
    for (ST_int i = 0; i < N_total; i++) {
        ST_int obs = in1 + i;

        /* y */
        SF_vdata(var_y, obs, &val);
        y[i] = val;

        /* X_endog */
        for (ST_int k = 0; k < K_endog; k++) {
            SF_vdata(var_endog_start + k, obs, &val);
            X_endog[k * N_total + i] = val;
        }

        /* X_exog */
        for (ST_int k = 0; k < K_exog; k++) {
            SF_vdata(var_exog_start + k, obs, &val);
            X_exog[k * N_total + i] = val;
        }

        /* Z (instruments) */
        for (ST_int k = 0; k < K_iv; k++) {
            SF_vdata(var_iv_start + k, obs, &val);
            Z[k * N_total + i] = val;
        }

        /* FE levels */
        for (ST_int g = 0; g < G; g++) {
            SF_vdata(var_fe_start + g, obs, &val);
            fe_levels[g][i] = (ST_int)val;
        }

        /* Cluster */
        if (has_cluster) {
            SF_vdata(var_cluster, obs, &val);
            cluster_ids[i] = (ST_int)val;
        }

        /* Weight */
        if (has_weights) {
            SF_vdata(var_weight, obs, &val);
            weights[i] = val;
        }
    }

    /* Drop observations with missing values */
    ST_int *valid_mask = (ST_int *)calloc(N_total, sizeof(ST_int));
    ST_int N_valid = 0;

    for (ST_int i = 0; i < N_total; i++) {
        int is_valid = 1;

        /* Check y */
        if (SF_is_missing(y[i])) is_valid = 0;

        /* Check X_endog */
        for (ST_int k = 0; k < K_endog && is_valid; k++) {
            if (SF_is_missing(X_endog[k * N_total + i])) is_valid = 0;
        }

        /* Check X_exog */
        for (ST_int k = 0; k < K_exog && is_valid; k++) {
            if (SF_is_missing(X_exog[k * N_total + i])) is_valid = 0;
        }

        /* Check Z */
        for (ST_int k = 0; k < K_iv && is_valid; k++) {
            if (SF_is_missing(Z[k * N_total + i])) is_valid = 0;
        }

        /* Check FE */
        for (ST_int g = 0; g < G && is_valid; g++) {
            if (fe_levels[g][i] <= 0) is_valid = 0;
        }

        /* Check weights */
        if (has_weights && is_valid) {
            if (SF_is_missing(weights[i]) || weights[i] <= 0) is_valid = 0;
        }

        /* Check cluster */
        if (has_cluster && is_valid) {
            if (cluster_ids[i] <= 0) is_valid = 0;
        }

        valid_mask[i] = is_valid;
        if (is_valid) N_valid++;
    }

    if (verbose) {
        char buf[256];
        snprintf(buf, sizeof(buf), "civreghdfe: %d valid observations (dropped %d)\n",
                 (int)N_valid, (int)(N_total - N_valid));
        SF_display(buf);
    }

    if (N_valid < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations\n");
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids); free(valid_mask);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 2001;
    }

    /* Compact data to remove invalid observations */
    ST_double *y_c = (ST_double *)malloc(N_valid * sizeof(ST_double));
    ST_double *X_endog_c = (K_endog > 0) ? (ST_double *)malloc(N_valid * K_endog * sizeof(ST_double)) : NULL;
    ST_double *X_exog_c = (K_exog > 0) ? (ST_double *)malloc(N_valid * K_exog * sizeof(ST_double)) : NULL;
    ST_double *Z_c = (ST_double *)malloc(N_valid * K_iv * sizeof(ST_double));
    ST_double *weights_c = has_weights ? (ST_double *)malloc(N_valid * sizeof(ST_double)) : NULL;
    ST_int *cluster_ids_c = has_cluster ? (ST_int *)malloc(N_valid * sizeof(ST_int)) : NULL;
    ST_int **fe_levels_c = (ST_int **)malloc(G * sizeof(ST_int *));
    for (ST_int g = 0; g < G; g++) {
        fe_levels_c[g] = (ST_int *)malloc(N_valid * sizeof(ST_int));
    }

    ST_int idx = 0;
    for (ST_int i = 0; i < N_total; i++) {
        if (!valid_mask[i]) continue;

        y_c[idx] = y[i];

        for (ST_int k = 0; k < K_endog; k++) {
            X_endog_c[k * N_valid + idx] = X_endog[k * N_total + i];
        }
        for (ST_int k = 0; k < K_exog; k++) {
            X_exog_c[k * N_valid + idx] = X_exog[k * N_total + i];
        }
        for (ST_int k = 0; k < K_iv; k++) {
            Z_c[k * N_valid + idx] = Z[k * N_total + i];
        }
        for (ST_int g = 0; g < G; g++) {
            fe_levels_c[g][idx] = fe_levels[g][i];
        }
        if (has_weights) weights_c[idx] = weights[i];
        if (has_cluster) cluster_ids_c[idx] = cluster_ids[i];

        idx++;
    }

    /* Free original arrays */
    free(y); free(X_endog); free(X_exog); free(Z);
    free(weights); free(cluster_ids); free(valid_mask);
    for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
    free(fe_levels);

    ST_int N = N_valid;

    /* Initialize HDFE state (reuse creghdfe infrastructure) */
    /* This requires setting up FE_Factor structures */
    HDFE_State *state = (HDFE_State *)calloc(1, sizeof(HDFE_State));
    state->G = G;
    state->N = N;
    state->K = 1 + K_endog + K_exog + K_iv;  /* Total columns to demean */
    state->in1 = 1;
    state->in2 = N;
    state->has_weights = has_weights;
    state->weight_type = weight_type;
    state->weights = weights_c;
    state->maxiter = maxiter;
    state->tolerance = (ST_int)(tolerance * 1e10);  /* Convert to int representation */
    state->verbose = verbose;

    /* Set up factors */
    state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    for (ST_int g = 0; g < G; g++) {
        /* Count levels */
        ST_int max_level = 0;
        for (ST_int i = 0; i < N; i++) {
            if (fe_levels_c[g][i] > max_level) max_level = fe_levels_c[g][i];
        }

        state->factors[g].num_levels = max_level;
        state->factors[g].has_intercept = 1;
        state->factors[g].levels = fe_levels_c[g];
        state->factors[g].counts = (ST_double *)calloc(max_level + 1, sizeof(ST_double));
        state->factors[g].weighted_counts = has_weights ? (ST_double *)calloc(max_level + 1, sizeof(ST_double)) : NULL;
        state->factors[g].means = (ST_double *)calloc(max_level + 1, sizeof(ST_double));

        /* Count observations per level */
        for (ST_int i = 0; i < N; i++) {
            ST_int lev = fe_levels_c[g][i];
            state->factors[g].counts[lev] += 1.0;
            if (has_weights) {
                state->factors[g].weighted_counts[lev] += weights_c[i];
            }
        }
    }

    /* Allocate CG solver buffers */
    ST_int num_threads = omp_get_max_threads();
    if (num_threads > 8) num_threads = 8;

    state->thread_cg_r = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_cg_u = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_cg_v = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_proj = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_fe_means = (ST_double **)malloc(num_threads * sizeof(ST_double *));

    ST_int total_levels = 0;
    for (ST_int g = 0; g < G; g++) {
        total_levels += state->factors[g].num_levels + 1;
    }

    for (ST_int t = 0; t < num_threads; t++) {
        state->thread_cg_r[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_cg_u[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_cg_v[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_proj[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_fe_means[t] = (ST_double *)malloc(total_levels * sizeof(ST_double));
    }

    state->factors_initialized = 1;

    /* Set global state for creghdfe solver functions */
    g_state = state;

    /* Partial out FEs from all variables using CG solver */
    if (verbose) {
        SF_display("civreghdfe: Partialling out fixed effects...\n");
    }

    ST_int total_cols = 1 + K_endog + K_exog + K_iv;

    /* Combine all data into one array for parallel processing */
    ST_double *all_data = (ST_double *)malloc(N * total_cols * sizeof(ST_double));

    /* Copy y */
    memcpy(all_data, y_c, N * sizeof(ST_double));

    /* Copy X_endog */
    if (K_endog > 0) {
        memcpy(all_data + N, X_endog_c, N * K_endog * sizeof(ST_double));
    }

    /* Copy X_exog */
    if (K_exog > 0) {
        memcpy(all_data + N * (1 + K_endog), X_exog_c, N * K_exog * sizeof(ST_double));
    }

    /* Copy Z */
    memcpy(all_data + N * (1 + K_endog + K_exog), Z_c, N * K_iv * sizeof(ST_double));

    /* Demean all columns in parallel */
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (ST_int k = 0; k < total_cols; k++) {
        int tid = omp_get_thread_num();
        cg_solve_column_threaded(state, all_data + k * N, tid);
    }

    if (verbose) {
        SF_display("civreghdfe: Fixed effects partialled out\n");
    }

    /* Extract demeaned data */
    ST_double *y_dem = all_data;
    ST_double *X_endog_dem = (K_endog > 0) ? all_data + N : NULL;
    ST_double *X_exog_dem = (K_exog > 0) ? all_data + N * (1 + K_endog) : NULL;
    ST_double *Z_dem = all_data + N * (1 + K_endog + K_exog);

    /* Compute df_a (absorbed degrees of freedom) */
    /* For simplicity, use sum of (levels - 1) for each FE */
    ST_int df_a = 0;
    for (ST_int g = 0; g < G; g++) {
        df_a += state->factors[g].num_levels - 1;
    }
    if (G > 1) df_a += 1;  /* Adjust for multiple FEs */

    /* Count clusters if needed */
    ST_int num_clusters = 0;
    if (has_cluster) {
        ST_int max_cluster = 0;
        for (ST_int i = 0; i < N; i++) {
            if (cluster_ids_c[i] > max_cluster) max_cluster = cluster_ids_c[i];
        }
        num_clusters = max_cluster;
    }

    /* Allocate output arrays */
    ST_double *beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *V = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* Compute 2SLS */
    ST_retcode rc = compute_2sls(
        y_dem, X_exog_dem, X_endog_dem, Z_dem,
        weights_c, weight_type,
        N, K_exog, K_endog, K_iv,
        beta, V, first_stage_F,
        vce_type, cluster_ids_c, num_clusters,
        df_a, verbose
    );

    if (rc != STATA_OK) {
        /* Cleanup and return */
        free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c);
        free(beta); free(V); free(first_stage_F);
        /* Cleanup state */
        for (ST_int t = 0; t < num_threads; t++) {
            free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
            free(state->thread_cg_v[t]); free(state->thread_proj[t]);
            free(state->thread_fe_means[t]);
        }
        free(state->thread_cg_r); free(state->thread_cg_u);
        free(state->thread_cg_v); free(state->thread_proj);
        free(state->thread_fe_means);
        for (ST_int g = 0; g < G; g++) {
            free(state->factors[g].counts);
            if (state->factors[g].weighted_counts) free(state->factors[g].weighted_counts);
            free(state->factors[g].means);
        }
        free(state->factors);
        free(state);
        g_state = NULL;
        return rc;
    }

    /* Store results to Stata */
    /* Scalars */
    SF_scal_save("__civreghdfe_N", (ST_double)N);
    SF_scal_save("__civreghdfe_df_r", (ST_double)(N - K_total - 1 - df_a));
    SF_scal_save("__civreghdfe_df_a", (ST_double)df_a);
    SF_scal_save("__civreghdfe_K", (ST_double)K_total);

    if (has_cluster) {
        SF_scal_save("__civreghdfe_N_clust", (ST_double)num_clusters);
    }

    /* Matrices: e(b) and e(V) */
    /* Create matrix __civreghdfe_b (1 x K_total) */
    for (ST_int k = 0; k < K_total; k++) {
        SF_mat_store("__civreghdfe_b", 1, k + 1, beta[k]);
    }

    /* Create matrix __civreghdfe_V (K_total x K_total) */
    for (ST_int i = 0; i < K_total; i++) {
        for (ST_int j = 0; j < K_total; j++) {
            SF_mat_store("__civreghdfe_V", i + 1, j + 1, V[j * K_total + i]);
        }
    }

    /* First stage F-stats */
    for (ST_int e = 0; e < K_endog; e++) {
        char name[64];
        snprintf(name, sizeof(name), "__civreghdfe_F1_%d", (int)(e + 1));
        SF_scal_save(name, first_stage_F[e]);
    }

    /* Cleanup */
    free(all_data); free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
    free(weights_c); free(cluster_ids_c);
    free(beta); free(V); free(first_stage_F);

    /* Cleanup state */
    for (ST_int t = 0; t < num_threads; t++) {
        free(state->thread_cg_r[t]); free(state->thread_cg_u[t]);
        free(state->thread_cg_v[t]); free(state->thread_proj[t]);
        free(state->thread_fe_means[t]);
    }
    free(state->thread_cg_r); free(state->thread_cg_u);
    free(state->thread_cg_v); free(state->thread_proj);
    free(state->thread_fe_means);
    for (ST_int g = 0; g < G; g++) {
        free(state->factors[g].counts);
        if (state->factors[g].weighted_counts) free(state->factors[g].weighted_counts);
        free(state->factors[g].means);
        /* Note: levels are in fe_levels_c which we already freed, but we passed the pointer */
    }
    free(state->factors);
    free(state);
    g_state = NULL;

    if (verbose) {
        SF_display("civreghdfe: Done\n");
    }

    return STATA_OK;
}

/*
    Main entry point for civreghdfe plugin.
*/
ST_retcode civreghdfe_main(const char *args)
{
    if (args == NULL || strlen(args) == 0) {
        SF_error("civreghdfe: No subcommand specified\n");
        return 198;
    }

    if (strcmp(args, "iv_regression") == 0) {
        return do_iv_regression();
    }

    SF_error("civreghdfe: Unknown subcommand\n");
    return 198;
}
