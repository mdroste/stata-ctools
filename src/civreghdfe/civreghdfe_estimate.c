/*
    civreghdfe_estimate.c
    Core IV Estimation: 2SLS, LIML, Fuller, GMM2S, CUE

    This module implements the k-class family of IV estimators.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "civreghdfe_estimate.h"
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"

/* Forward declarations from creghdfe_solver */
extern int cholesky(double *A, int K);
extern int invert_from_cholesky(const double *L, int K, double *A_inv);

/*
    Initialize IV estimation context and compute basic matrices.
*/
ST_retcode ivest_init_context(
    IVEstContext *ctx,
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
    ST_int verbose
)
{
    ST_int K_total = K_exog + K_endog;
    ST_int i, j, k;

    /* Store dimensions */
    ctx->N = N;
    ctx->K_exog = K_exog;
    ctx->K_endog = K_endog;
    ctx->K_iv = K_iv;
    ctx->K_total = K_total;
    ctx->weights = weights;
    ctx->weight_type = weight_type;
    ctx->verbose = verbose;

    /* Store pointers to original data (not owned) */
    ctx->Z = Z;
    ctx->y = y;

    /* Check identification */
    if (K_iv < K_total) {
        SF_error("civreghdfe: Model is underidentified\n");
        return 198;
    }

    /* Allocate matrices */
    ctx->ZtZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ctx->ZtZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ctx->ZtX = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
    ctx->Zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
    ctx->XtPzX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ctx->XtPzy = (ST_double *)calloc(K_total, sizeof(ST_double));
    ctx->temp_kiv_ktotal = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
    ctx->X_all = (ST_double *)malloc(N * K_total * sizeof(ST_double));

    if (!ctx->ZtZ || !ctx->ZtZ_inv || !ctx->ZtX || !ctx->Zty ||
        !ctx->XtPzX || !ctx->XtPzy || !ctx->temp_kiv_ktotal || !ctx->X_all) {
        SF_error("civreghdfe: Memory allocation failed\n");
        ivest_free_context(ctx);
        return 920;
    }

    /* Combine X_exog and X_endog into X_all */
    if (K_exog > 0 && X_exog) {
        memcpy(ctx->X_all, X_exog, N * K_exog * sizeof(ST_double));
    }
    if (K_endog > 0 && X_endog) {
        memcpy(ctx->X_all + N * K_exog, X_endog, N * K_endog * sizeof(ST_double));
    }

    /* Compute Z'Z (weighted if needed) */
    if (weights && weight_type != 0) {
        civreghdfe_matmul_atdb(Z, Z, weights, N, K_iv, K_iv, ctx->ZtZ);
    } else {
        civreghdfe_matmul_atb(Z, Z, N, K_iv, K_iv, ctx->ZtZ);
    }

    /* Invert Z'Z */
    memcpy(ctx->ZtZ_inv, ctx->ZtZ, K_iv * K_iv * sizeof(ST_double));
    if (cholesky(ctx->ZtZ_inv, K_iv) != 0) {
        SF_error("civreghdfe: Z'Z is singular (instruments may be collinear)\n");
        ivest_free_context(ctx);
        return 198;
    }
    if (invert_from_cholesky(ctx->ZtZ_inv, K_iv, ctx->ZtZ_inv) != 0) {
        SF_error("civreghdfe: Failed to invert Z'Z\n");
        ivest_free_context(ctx);
        return 198;
    }

    /* Compute Z'X and Z'y */
    if (weights && weight_type != 0) {
        civreghdfe_matmul_atdb(Z, ctx->X_all, weights, N, K_iv, K_total, ctx->ZtX);
        /* Z'y */
        for (i = 0; i < K_iv; i++) {
            ST_double sum = 0.0;
            const ST_double *z_col = Z + i * N;
            for (k = 0; k < N; k++) {
                sum += z_col[k] * weights[k] * y[k];
            }
            ctx->Zty[i] = sum;
        }
    } else {
        civreghdfe_matmul_atb(Z, ctx->X_all, N, K_iv, K_total, ctx->ZtX);
        /* Z'y */
        for (i = 0; i < K_iv; i++) {
            ST_double sum = 0.0;
            const ST_double *z_col = Z + i * N;
            for (k = 0; k < N; k++) {
                sum += z_col[k] * y[k];
            }
            ctx->Zty[i] = sum;
        }
    }

    /* Compute (Z'Z)^-1 * Z'X -> temp_kiv_ktotal */
    civreghdfe_matmul_ab(ctx->ZtZ_inv, ctx->ZtX, K_iv, K_iv, K_total, ctx->temp_kiv_ktotal);

    /* Compute X'P_Z X = (Z'X)' * (Z'Z)^-1 * Z'X */
    for (j = 0; j < K_total; j++) {
        for (i = 0; i < K_total; i++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += ctx->ZtX[i * K_iv + k] * ctx->temp_kiv_ktotal[j * K_iv + k];
            }
            ctx->XtPzX[j * K_total + i] = sum;
        }
    }

    /* Compute X'P_Z y = (Z'X)' * (Z'Z)^-1 * Z'y */
    ST_double *ZtZ_inv_Zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
    for (i = 0; i < K_iv; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += ctx->ZtZ_inv[k * K_iv + i] * ctx->Zty[k];
        }
        ZtZ_inv_Zty[i] = sum;
    }

    for (i = 0; i < K_total; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += ctx->ZtX[i * K_iv + k] * ZtZ_inv_Zty[k];
        }
        ctx->XtPzy[i] = sum;
    }
    free(ZtZ_inv_Zty);

    return STATA_OK;
}

/*
    Free resources allocated by ivest_init_context.
*/
void ivest_free_context(IVEstContext *ctx)
{
    if (ctx->ZtZ) { free(ctx->ZtZ); ctx->ZtZ = NULL; }
    if (ctx->ZtZ_inv) { free(ctx->ZtZ_inv); ctx->ZtZ_inv = NULL; }
    if (ctx->ZtX) { free(ctx->ZtX); ctx->ZtX = NULL; }
    if (ctx->Zty) { free(ctx->Zty); ctx->Zty = NULL; }
    if (ctx->XtPzX) { free(ctx->XtPzX); ctx->XtPzX = NULL; }
    if (ctx->XtPzy) { free(ctx->XtPzy); ctx->XtPzy = NULL; }
    if (ctx->temp_kiv_ktotal) { free(ctx->temp_kiv_ktotal); ctx->temp_kiv_ktotal = NULL; }
    if (ctx->X_all) { free(ctx->X_all); ctx->X_all = NULL; }
}

/*
    Compute k-class estimator coefficients.
*/
ST_retcode ivest_compute_kclass(
    IVEstContext *ctx,
    const ST_double *y,
    ST_double kclass,
    ST_double *beta,
    ST_double *resid
)
{
    ST_int N = ctx->N;
    ST_int K_total = ctx->K_total;
    ST_int i, j, k;

    /* Allocate work arrays */
    ST_double *XkX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *Xky = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *XkX_L = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));

    if (!XkX || !Xky || !XkX_L || !beta_temp) {
        SF_error("civreghdfe: Memory allocation failed\n");
        free(XkX); free(Xky); free(XkX_L); free(beta_temp);
        return 920;
    }

    if (kclass != 1.0) {
        /* Compute X'X and X'y for k-class */
        ST_double *XtX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        ST_double *Xty = (ST_double *)calloc(K_total, sizeof(ST_double));

        if (ctx->weights && ctx->weight_type != 0) {
            civreghdfe_matmul_atdb(ctx->X_all, ctx->X_all, ctx->weights, N, K_total, K_total, XtX);
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = ctx->X_all + i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * ctx->weights[k] * y[k];
                }
                Xty[i] = sum;
            }
        } else {
            civreghdfe_matmul_atb(ctx->X_all, ctx->X_all, N, K_total, K_total, XtX);
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = ctx->X_all + i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * y[k];
                }
                Xty[i] = sum;
            }
        }

        /* XkX = (1-k)*X'X + k*X'P_Z X */
        /* Xky = (1-k)*X'y + k*X'P_Z y */
        ST_double one_minus_k = 1.0 - kclass;
        for (i = 0; i < K_total * K_total; i++) {
            XkX[i] = one_minus_k * XtX[i] + kclass * ctx->XtPzX[i];
        }
        for (i = 0; i < K_total; i++) {
            Xky[i] = one_minus_k * Xty[i] + kclass * ctx->XtPzy[i];
        }

        free(XtX);
        free(Xty);
    } else {
        /* k = 1 (2SLS): just use X'P_Z X and X'P_Z y */
        memcpy(XkX, ctx->XtPzX, K_total * K_total * sizeof(ST_double));
        memcpy(Xky, ctx->XtPzy, K_total * sizeof(ST_double));
    }

    /* Solve XkX * beta = Xky using Cholesky */
    memcpy(XkX_L, XkX, K_total * K_total * sizeof(ST_double));

    if (cholesky(XkX_L, K_total) != 0) {
        SF_error("civreghdfe: XkX matrix is singular\n");
        free(XkX); free(Xky); free(XkX_L); free(beta_temp);
        return 198;
    }

    /* Forward substitution */
    for (i = 0; i < K_total; i++) {
        ST_double sum = Xky[i];
        for (j = 0; j < i; j++) {
            sum -= XkX_L[i * K_total + j] * beta_temp[j];
        }
        beta_temp[i] = sum / XkX_L[i * K_total + i];
    }

    /* Backward substitution */
    for (i = K_total - 1; i >= 0; i--) {
        ST_double sum = beta_temp[i];
        for (j = i + 1; j < K_total; j++) {
            sum -= XkX_L[j * K_total + i] * beta[j];
        }
        beta[i] = sum / XkX_L[i * K_total + i];
    }

    /* Compute residuals if requested */
    if (resid) {
        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_total; k++) {
                pred += ctx->X_all[k * N + i] * beta[k];
            }
            resid[i] = y[i] - pred;
        }
    }

    free(XkX);
    free(Xky);
    free(XkX_L);
    free(beta_temp);

    return STATA_OK;
}

/*
    Compute GMM2S (two-step efficient GMM) estimator.

    Two-step efficient GMM:
    Step 1: Initial residuals from 2SLS (provided as input)
    Step 2: Compute optimal weighting matrix W = (Z'ΩZ)^-1 where Ω = diag(e²)
    Step 3: Re-estimate β = (X'ZWZ'X)^-1 X'ZWZ'y
*/
ST_retcode ivest_compute_gmm2s(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *initial_resid,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *beta,
    ST_double *resid
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    ST_int K_total = ctx->K_total;
    const ST_double *Z = ctx->Z;
    ST_int i, j, k;

    /* Compute optimal weighting matrix Z'ΩZ */
    ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ST_double *ZOmegaZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));

    if (!ZOmegaZ || !ZOmegaZ_inv) {
        free(ZOmegaZ); free(ZOmegaZ_inv);
        return 920;
    }

    if (cluster_ids && num_clusters > 0) {
        /* Cluster-robust optimal weighting matrix */
        ST_double *cluster_ze = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
        if (!cluster_ze) {
            free(ZOmegaZ); free(ZOmegaZ_inv);
            return 920;
        }

        for (i = 0; i < N; i++) {
            ST_int c = cluster_ids[i] - 1;
            if (c < 0 || c >= num_clusters) continue;
            ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
            ST_double we = w * initial_resid[i];
            for (j = 0; j < K_iv; j++) {
                cluster_ze[c * K_iv + j] += Z[j * N + i] * we;
            }
        }

        for (ST_int c = 0; c < num_clusters; c++) {
            for (j = 0; j < K_iv; j++) {
                for (k = 0; k <= j; k++) {
                    ST_double contrib = cluster_ze[c * K_iv + j] * cluster_ze[c * K_iv + k];
                    ZOmegaZ[j * K_iv + k] += contrib;
                    if (k != j) ZOmegaZ[k * K_iv + j] += contrib;
                }
            }
        }
        free(cluster_ze);
    } else {
        /* Heteroskedastic optimal weighting matrix: Z'diag(e²)Z */
        for (i = 0; i < N; i++) {
            ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
            ST_double e2 = w * initial_resid[i] * initial_resid[i];
            for (j = 0; j < K_iv; j++) {
                ST_double z_j = Z[j * N + i];
                ST_double z_j_e2 = z_j * e2;
                for (k = 0; k <= j; k++) {
                    ST_double contrib = z_j_e2 * Z[k * N + i];
                    ZOmegaZ[j * K_iv + k] += contrib;
                    if (k != j) ZOmegaZ[k * K_iv + j] += contrib;
                }
            }
        }
    }

    /* Invert Z'ΩZ */
    memcpy(ZOmegaZ_inv, ZOmegaZ, K_iv * K_iv * sizeof(ST_double));
    if (cholesky(ZOmegaZ_inv, K_iv) != 0) {
        SF_error("civreghdfe: GMM2S optimal weighting matrix is singular\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        return 198;
    }
    if (invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv) != 0) {
        SF_error("civreghdfe: Failed to invert GMM2S weighting matrix\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        return 198;
    }

    /* Compute GMM estimator: β = (X'ZWZ'X)^-1 X'ZWZ'y */
    /* 1. WZtX = W * Z'X */
    ST_double *WZtX = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
    civreghdfe_matmul_ab(ZOmegaZ_inv, ctx->ZtX, K_iv, K_iv, K_total, WZtX);

    /* 2. XZW = (Z'X)' * W (K_total x K_iv) */
    ST_double *XZW = (ST_double *)calloc(K_total * K_iv, sizeof(ST_double));
    for (i = 0; i < K_total; i++) {
        for (j = 0; j < K_iv; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += ctx->ZtX[i * K_iv + k] * ZOmegaZ_inv[j * K_iv + k];
            }
            XZW[j * K_total + i] = sum;
        }
    }

    /* 3. XZWZX = XZW * Z'X (K_total x K_total) */
    ST_double *XZWZX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    for (i = 0; i < K_total; i++) {
        for (j = 0; j < K_total; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += XZW[k * K_total + i] * ctx->ZtX[j * K_iv + k];
            }
            XZWZX[j * K_total + i] = sum;
        }
    }

    /* 4. XZWZy = XZW * Z'y (K_total x 1) */
    ST_double *XZWZy = (ST_double *)calloc(K_total, sizeof(ST_double));
    for (i = 0; i < K_total; i++) {
        ST_double sum = 0.0;
        for (k = 0; k < K_iv; k++) {
            sum += XZW[k * K_total + i] * ctx->Zty[k];
        }
        XZWZy[i] = sum;
    }

    /* 5. Solve XZWZX * beta = XZWZy */
    ST_double *XZWZX_L = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));
    memcpy(XZWZX_L, XZWZX, K_total * K_total * sizeof(ST_double));

    if (cholesky(XZWZX_L, K_total) != 0) {
        SF_error("civreghdfe: GMM2S X'ZWZ'X matrix is singular\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        free(WZtX); free(XZW); free(XZWZX); free(XZWZy);
        free(XZWZX_L); free(beta_temp);
        return 198;
    }

    /* Forward substitution */
    for (i = 0; i < K_total; i++) {
        ST_double sum = XZWZy[i];
        for (j = 0; j < i; j++) {
            sum -= XZWZX_L[i * K_total + j] * beta_temp[j];
        }
        beta_temp[i] = sum / XZWZX_L[i * K_total + i];
    }

    /* Backward substitution */
    for (i = K_total - 1; i >= 0; i--) {
        ST_double sum = beta_temp[i];
        for (j = i + 1; j < K_total; j++) {
            sum -= XZWZX_L[j * K_total + i] * beta[j];
        }
        beta[i] = sum / XZWZX_L[i * K_total + i];
    }

    /* Compute residuals */
    if (resid) {
        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_total; k++) {
                pred += ctx->X_all[k * N + i] * beta[k];
            }
            resid[i] = y[i] - pred;
        }
    }

    free(ZOmegaZ); free(ZOmegaZ_inv);
    free(WZtX); free(XZW); free(XZWZX); free(XZWZy);
    free(XZWZX_L); free(beta_temp);

    return STATA_OK;
}

/*
    Compute CUE (Continuously Updated Estimator).

    CUE minimizes: Q(β) = g(β)' W(β)^-1 g(β)
    where g(β) = Z'(y - Xβ) and W(β) = Z'Ω(β)Z

    Uses iterative re-weighting:
    1. Start with initial beta (from 2SLS/GMM2S)
    2. Iterate: update W based on current residuals, re-estimate β
    3. Stop when β converges
*/
ST_retcode ivest_compute_cue(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *initial_beta,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *beta,
    ST_double *resid,
    ST_int max_iter,
    ST_double tol
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    ST_int K_total = ctx->K_total;
    const ST_double *Z = ctx->Z;
    ST_int i, j, k;

    /* Initialize beta from input */
    memcpy(beta, initial_beta, K_total * sizeof(ST_double));

    /* Allocate working arrays */
    ST_double *current_resid = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *beta_old = (ST_double *)malloc(K_total * sizeof(ST_double));
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));

    if (!current_resid || !beta_old || !beta_temp) {
        free(current_resid); free(beta_old); free(beta_temp);
        return 920;
    }

    if (ctx->verbose) {
        SF_display("civreghdfe: Computing CUE (iterative re-weighting)\n");
    }

    ST_int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* Compute current residuals */
        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_total; k++) {
                pred += ctx->X_all[k * N + i] * beta[k];
            }
            current_resid[i] = y[i] - pred;
        }

        /* Compute optimal weighting matrix based on current residuals */
        ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        ST_double *ZOmegaZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));

        if (!ZOmegaZ || !ZOmegaZ_inv) {
            free(ZOmegaZ); free(ZOmegaZ_inv);
            free(current_resid); free(beta_old); free(beta_temp);
            return 920;
        }

        if (cluster_ids && num_clusters > 0) {
            /* Cluster-robust */
            ST_double *cluster_ze = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
            if (!cluster_ze) {
                free(ZOmegaZ); free(ZOmegaZ_inv);
                free(current_resid); free(beta_old); free(beta_temp);
                return 920;
            }

            for (i = 0; i < N; i++) {
                ST_int c = cluster_ids[i] - 1;
                if (c < 0 || c >= num_clusters) continue;
                ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
                ST_double we = w * current_resid[i];
                for (j = 0; j < K_iv; j++) {
                    cluster_ze[c * K_iv + j] += Z[j * N + i] * we;
                }
            }

            for (ST_int c = 0; c < num_clusters; c++) {
                for (j = 0; j < K_iv; j++) {
                    for (k = 0; k <= j; k++) {
                        ST_double contrib = cluster_ze[c * K_iv + j] * cluster_ze[c * K_iv + k];
                        ZOmegaZ[j * K_iv + k] += contrib;
                        if (k != j) ZOmegaZ[k * K_iv + j] += contrib;
                    }
                }
            }
            free(cluster_ze);
        } else {
            /* Heteroskedastic */
            for (i = 0; i < N; i++) {
                ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
                ST_double e2 = w * current_resid[i] * current_resid[i];
                for (j = 0; j < K_iv; j++) {
                    ST_double z_j = Z[j * N + i];
                    ST_double z_j_e2 = z_j * e2;
                    for (k = 0; k <= j; k++) {
                        ST_double contrib = z_j_e2 * Z[k * N + i];
                        ZOmegaZ[j * K_iv + k] += contrib;
                        if (k != j) ZOmegaZ[k * K_iv + j] += contrib;
                    }
                }
            }
        }

        /* Invert Z'ΩZ */
        memcpy(ZOmegaZ_inv, ZOmegaZ, K_iv * K_iv * sizeof(ST_double));
        if (cholesky(ZOmegaZ_inv, K_iv) != 0) {
            if (ctx->verbose) SF_display("civreghdfe: CUE weighting matrix singular, stopping\n");
            free(ZOmegaZ); free(ZOmegaZ_inv);
            break;
        }
        invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv);

        /* Compute GMM update: β = (X'ZWZ'X)^-1 X'ZWZ'y */
        ST_double *XZW = (ST_double *)calloc(K_total * K_iv, sizeof(ST_double));
        for (i = 0; i < K_total; i++) {
            for (j = 0; j < K_iv; j++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_iv; k++) {
                    sum += ctx->ZtX[i * K_iv + k] * ZOmegaZ_inv[j * K_iv + k];
                }
                XZW[j * K_total + i] = sum;
            }
        }

        ST_double *XZWZX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        for (i = 0; i < K_total; i++) {
            for (j = 0; j < K_total; j++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_iv; k++) {
                    sum += XZW[k * K_total + i] * ctx->ZtX[j * K_iv + k];
                }
                XZWZX[j * K_total + i] = sum;
            }
        }

        ST_double *XZWZy = (ST_double *)calloc(K_total, sizeof(ST_double));
        for (i = 0; i < K_total; i++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += XZW[k * K_total + i] * ctx->Zty[k];
            }
            XZWZy[i] = sum;
        }

        /* Solve for new beta */
        ST_double *XZWZX_L = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
        memcpy(XZWZX_L, XZWZX, K_total * K_total * sizeof(ST_double));

        if (cholesky(XZWZX_L, K_total) != 0) {
            if (ctx->verbose) SF_display("civreghdfe: CUE X'ZWZ'X singular, stopping\n");
            free(ZOmegaZ); free(ZOmegaZ_inv);
            free(XZW); free(XZWZX); free(XZWZy); free(XZWZX_L);
            break;
        }

        /* Forward substitution */
        ST_double *beta_new = (ST_double *)calloc(K_total, sizeof(ST_double));
        for (i = 0; i < K_total; i++) {
            ST_double sum = XZWZy[i];
            for (j = 0; j < i; j++) {
                sum -= XZWZX_L[i * K_total + j] * beta_temp[j];
            }
            beta_temp[i] = sum / XZWZX_L[i * K_total + i];
        }

        /* Backward substitution */
        for (i = K_total - 1; i >= 0; i--) {
            ST_double sum = beta_temp[i];
            for (j = i + 1; j < K_total; j++) {
                sum -= XZWZX_L[j * K_total + i] * beta_new[j];
            }
            beta_new[i] = sum / XZWZX_L[i * K_total + i];
        }

        /* Check convergence */
        ST_double max_diff = 0.0;
        for (i = 0; i < K_total; i++) {
            ST_double diff = fabs(beta_new[i] - beta[i]);
            if (diff > max_diff) max_diff = diff;
        }

        /* Update beta */
        memcpy(beta_old, beta, K_total * sizeof(ST_double));
        memcpy(beta, beta_new, K_total * sizeof(ST_double));

        free(ZOmegaZ); free(ZOmegaZ_inv);
        free(XZW); free(XZWZX); free(XZWZy); free(XZWZX_L);
        free(beta_new);

        if (max_diff < tol) {
            if (ctx->verbose) {
                char buf[256];
                snprintf(buf, sizeof(buf), "civreghdfe: CUE converged in %d iterations (diff=%.2e)\n",
                         iter + 1, max_diff);
                SF_display(buf);
            }
            break;
        }
    }

    if (iter == max_iter && ctx->verbose) {
        SF_display("civreghdfe: CUE reached max iterations\n");
    }

    /* Compute final residuals */
    if (resid) {
        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_total; k++) {
                pred += ctx->X_all[k * N + i] * beta[k];
            }
            resid[i] = y[i] - pred;
        }
    }

    free(current_resid);
    free(beta_old);
    free(beta_temp);

    return STATA_OK;
}

/*
    Compute first-stage F-statistics.

    Uses the Z matrix stored in context to compute first-stage F statistics
    for each endogenous variable.
*/
ST_retcode ivest_first_stage_f(
    IVEstContext *ctx,
    const ST_double *X_endog,
    ST_int df_a,
    ST_double *first_stage_F
)
{
    ST_int N = ctx->N;
    ST_int K_exog = ctx->K_exog;
    ST_int K_endog = ctx->K_endog;
    ST_int K_iv = ctx->K_iv;
    const ST_double *Z = ctx->Z;
    ST_int i, j, k;

    for (ST_int e = 0; e < K_endog; e++) {
        const ST_double *X_e = X_endog + e * N;

        /* Compute X_e'P_Z X_e */
        ST_double xpx = 0.0;
        for (i = 0; i < N; i++) {
            ST_double pz_xe = 0.0;
            for (k = 0; k < K_iv; k++) {
                ST_double ziz_inv_zx = 0.0;
                for (j = 0; j < K_iv; j++) {
                    ziz_inv_zx += ctx->ZtZ_inv[k * K_iv + j] * ctx->ZtX[(K_exog + e) * K_iv + j];
                }
                pz_xe += Z[k * N + i] * ziz_inv_zx;
            }
            ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
            xpx += w * X_e[i] * pz_xe;
        }

        /* Compute X_e'X_e */
        ST_double xx = 0.0;
        for (i = 0; i < N; i++) {
            ST_double w = (ctx->weights && ctx->weight_type != 0) ? ctx->weights[i] : 1.0;
            xx += w * X_e[i] * X_e[i];
        }

        /* R^2 of projection */
        ST_double r2 = (xx > 0) ? xpx / xx : 0.0;
        if (r2 > 1.0) r2 = 1.0;
        if (r2 < 0.0) r2 = 0.0;

        /* F-stat: (R^2 / q) / ((1 - R^2) / (N - K_iv - df_a)) */
        ST_int q = K_iv - K_exog;  /* Number of excluded instruments */
        if (q <= 0) q = 1;
        ST_int denom_df = N - K_iv - df_a;
        if (denom_df <= 0) denom_df = 1;

        first_stage_F[e] = (r2 / (ST_double)q) / ((1.0 - r2) / (ST_double)denom_df);
        if (first_stage_F[e] < 0) first_stage_F[e] = 0;

        if (ctx->verbose) {
            char buf[256];
            snprintf(buf, sizeof(buf), "  First-stage F[%d] = %g (R^2 = %g)\n",
                     (int)e, first_stage_F[e], r2);
            SF_display(buf);
        }
    }

    return STATA_OK;
}
