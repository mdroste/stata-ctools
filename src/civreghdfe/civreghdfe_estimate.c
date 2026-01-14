/*
    civreghdfe_estimate.c
    Core IV Estimation: 2SLS, LIML, Fuller, GMM2S, CUE

    This module implements the k-class family of IV estimators.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    /* Mark parameters as intentionally unused - stub function */
    (void)ctx;
    (void)y;
    (void)initial_resid;
    (void)cluster_ids;
    (void)num_clusters;
    (void)beta;
    (void)resid;

    /* GMM2S requires Z matrix which is not stored in context.
       Use compute_2sls() directly for GMM2S estimation. */
    SF_error("civreghdfe: GMM2S requires Z matrix - use compute_2sls directly\n");
    return 198;
}

/*
    Compute CUE (Continuously Updated Estimator).
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
    /* Mark parameters as intentionally unused - stub function */
    (void)ctx;
    (void)y;
    (void)initial_beta;
    (void)cluster_ids;
    (void)num_clusters;
    (void)beta;
    (void)resid;
    (void)max_iter;
    (void)tol;

    /* CUE requires iterative optimization with Z matrix.
       Use compute_2sls() directly for CUE estimation. */
    SF_error("civreghdfe: CUE requires Z matrix - use compute_2sls directly\n");
    return 198;
}

/*
    Compute first-stage F-statistics.

    Note: This function requires Z matrix for full computation.
    Use civreghdfe_compute_first_stage_F() in civreghdfe_tests.c
    which takes Z as an explicit parameter.
*/
ST_retcode ivest_first_stage_f(
    IVEstContext *ctx,
    const ST_double *X_endog,
    ST_int df_a,
    ST_double *first_stage_F
)
{
    /* Mark parameters as intentionally unused - stub function */
    (void)ctx;
    (void)X_endog;
    (void)df_a;
    (void)first_stage_F;

    /* This function requires Z matrix which is not stored in context.
       Use civreghdfe_compute_first_stage_F() from civreghdfe_tests.c instead,
       which takes Z as an explicit parameter. */
    SF_error("civreghdfe: ivest_first_stage_f requires Z matrix - use civreghdfe_compute_first_stage_F\n");
    return 198;
}
