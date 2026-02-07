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
#include "civreghdfe_impl.h"  /* For g_civreghdfe_noreturn */
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"
#include "civreghdfe_tests.h"
#include "../ctools_config.h"
#include "../ctools_spi.h"  /* Error-checking SPI wrappers */

/* OpenMP for parallel residual computation */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Convenience aliases for the matrix functions */
#define matmul_atb  ctools_matmul_atb
#define matmul_ab   ctools_matmul_ab
#define matmul_atdb ctools_matmul_atdb
#define compute_liml_lambda civreghdfe_compute_liml_lambda

/* Shared OLS functions */
#include "../ctools_ols.h"
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

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

    /* Allocate matrices - with overflow checks */
    size_t kiv_kiv_size, kiv_ktotal_size, ktotal_ktotal_size, n_ktotal_size;
    size_t n_kexog_size, n_kendog_size;

    if (ctools_safe_alloc_size((size_t)K_iv, (size_t)K_iv, sizeof(ST_double), &kiv_kiv_size) != 0 ||
        ctools_safe_alloc_size((size_t)K_iv, (size_t)K_total, sizeof(ST_double), &kiv_ktotal_size) != 0 ||
        ctools_safe_alloc_size((size_t)K_total, (size_t)K_total, sizeof(ST_double), &ktotal_ktotal_size) != 0 ||
        ctools_safe_alloc_size((size_t)N, (size_t)K_total, sizeof(ST_double), &n_ktotal_size) != 0) {
        SF_error("civreghdfe: Matrix size overflow\n");
        return 920;
    }

    ctx->ZtZ = (ST_double *)calloc(1, kiv_kiv_size);
    ctx->ZtZ_inv = (ST_double *)calloc(1, kiv_kiv_size);
    ctx->ZtX = (ST_double *)calloc(1, kiv_ktotal_size);
    ctx->Zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
    ctx->XtPzX = (ST_double *)calloc(1, ktotal_ktotal_size);
    ctx->XtPzy = (ST_double *)calloc(K_total, sizeof(ST_double));
    ctx->temp_kiv_ktotal = (ST_double *)calloc(1, kiv_ktotal_size);
    ctx->X_all = (ST_double *)malloc(n_ktotal_size);

    if (!ctx->ZtZ || !ctx->ZtZ_inv || !ctx->ZtX || !ctx->Zty ||
        !ctx->XtPzX || !ctx->XtPzy || !ctx->temp_kiv_ktotal || !ctx->X_all) {
        SF_error("civreghdfe: Memory allocation failed\n");
        ivest_free_context(ctx);
        return 920;
    }

    /* Combine X_exog and X_endog into X_all - with overflow checks */
    if (K_exog > 0 && X_exog) {
        if (ctools_safe_alloc_size((size_t)N, (size_t)K_exog, sizeof(ST_double), &n_kexog_size) != 0) {
            SF_error("civreghdfe: Matrix size overflow\n");
            ivest_free_context(ctx);
            return 920;
        }
        memcpy(ctx->X_all, X_exog, n_kexog_size);
    }
    if (K_endog > 0 && X_endog) {
        if (ctools_safe_alloc_size((size_t)N, (size_t)K_endog, sizeof(ST_double), &n_kendog_size) != 0) {
            SF_error("civreghdfe: Matrix size overflow\n");
            ivest_free_context(ctx);
            return 920;
        }
        memcpy(ctx->X_all + (size_t)N * K_exog, X_endog, n_kendog_size);
    }

    /* Compute Z'Z (weighted if needed) */
    if (weights && weight_type != 0) {
        ctools_matmul_atdb(Z, Z, weights, N, K_iv, K_iv, ctx->ZtZ);
    } else {
        ctools_matmul_atb(Z, Z, N, K_iv, K_iv, ctx->ZtZ);
    }

    /* Invert Z'Z */
    memcpy(ctx->ZtZ_inv, ctx->ZtZ, kiv_kiv_size);
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
        ctools_matmul_atdb(Z, ctx->X_all, weights, N, K_iv, K_total, ctx->ZtX);
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
        ctools_matmul_atb(Z, ctx->X_all, N, K_iv, K_total, ctx->ZtX);
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
    ctools_matmul_ab(ctx->ZtZ_inv, ctx->ZtX, K_iv, K_iv, K_total, ctx->temp_kiv_ktotal);

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
    if (!ZtZ_inv_Zty) return 920;
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
            ctools_matmul_atdb(ctx->X_all, ctx->X_all, ctx->weights, N, K_total, K_total, XtX);
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = ctx->X_all + i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * ctx->weights[k] * y[k];
                }
                Xty[i] = sum;
            }
        } else {
            ctools_matmul_atb(ctx->X_all, ctx->X_all, N, K_total, K_total, XtX);
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
    ST_double *resid,
    ST_double *XZWZX_inv_out
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
    ctools_matmul_ab(ZOmegaZ_inv, ctx->ZtX, K_iv, K_iv, K_total, WZtX);

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

    /* Output the GMM Hessian inverse (X'ZWZ'X)^-1 if requested */
    if (XZWZX_inv_out) {
        /* XZWZX_L contains the Cholesky factor, compute the inverse */
        ST_double *XZWZX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (XZWZX_inv) {
            /* Re-compute Cholesky and inverse from XZWZX */
            memcpy(XZWZX_inv, XZWZX, K_total * K_total * sizeof(ST_double));
            if (cholesky(XZWZX_inv, K_total) == 0) {
                invert_from_cholesky(XZWZX_inv, K_total, XZWZX_inv);
                memcpy(XZWZX_inv_out, XZWZX_inv, K_total * K_total * sizeof(ST_double));
            }
            free(XZWZX_inv);
        }
    }

    free(ZOmegaZ); free(ZOmegaZ_inv);
    free(WZtX); free(XZW); free(XZWZX); free(XZWZy);
    free(XZWZX_L); free(beta_temp);

    return STATA_OK;
}

/*
    Helper: Compute residuals given beta
*/
static void cue_compute_residuals(
    const ST_double *y,
    const ST_double *X_all,
    const ST_double *beta,
    ST_int N,
    ST_int K_total,
    ST_double *resid
)
{
    for (ST_int i = 0; i < N; i++) {
        ST_double pred = 0.0;
        for (ST_int k = 0; k < K_total; k++) {
            pred += X_all[k * N + i] * beta[k];
        }
        resid[i] = y[i] - pred;
    }
}

/*
    Helper: Compute CUE objective value Q(β) = g'W^{-1}g / N
    For homoskedastic: Q = e'Pz e / σ² where σ² = e'e/N
    For robust: Q = g'(Z'ΩZ)^{-1}g / N where Ω = diag(e²)
*/
static ST_double cue_objective(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *beta,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *work_resid  /* N-sized work buffer */
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    ST_int K_total = ctx->K_total;
    const ST_double *Z = ctx->Z;
    ST_int i, j, k;

    /* Compute residuals */
    cue_compute_residuals(y, ctx->X_all, beta, N, K_total, work_resid);

    /* Compute g = Z'e */
    ST_double *g = (ST_double *)calloc(K_iv, sizeof(ST_double));
    if (!g) return 1e30;  /* Return large value on allocation failure */

    for (i = 0; i < N; i++) {
        for (j = 0; j < K_iv; j++) {
            g[j] += Z[j * N + i] * work_resid[i];
        }
    }

    ST_double Q;

    if (vce_type == 0) {
        /* Homoskedastic: Q = e'Pz e / σ² = g'(Z'Z)^{-1}g / (e'e/N) */
        /* Compute g'(Z'Z)^{-1}g using pre-computed ZtZ_inv */
        ST_double *temp = (ST_double *)calloc(K_iv, sizeof(ST_double));
        if (!temp) {
            free(g);
            return 1e30;  /* Return large value on allocation failure */
        }
        for (i = 0; i < K_iv; i++) {
            for (j = 0; j < K_iv; j++) {
                temp[i] += ctx->ZtZ_inv[j * K_iv + i] * g[j];
            }
        }
        ST_double gWinvg = 0.0;
        for (i = 0; i < K_iv; i++) {
            gWinvg += g[i] * temp[i];
        }
        free(temp);

        /* Compute σ² = e'e / N */
        ST_double ete = 0.0;
        for (i = 0; i < N; i++) {
            ete += work_resid[i] * work_resid[i];
        }
        ST_double sigma2 = ete / N;

        Q = gWinvg / sigma2 / N;
    } else {
        /* Robust/Cluster: Q = g'W^{-1}g / N where W = Z'ΩZ */
        ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        ST_double *ZOmegaZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));

        if (!ZOmegaZ || !ZOmegaZ_inv) {
            free(g); free(ZOmegaZ); free(ZOmegaZ_inv);
            return 1e30;  /* Return large value on allocation failure */
        }

        if (cluster_ids && num_clusters > 0) {
            /* Cluster-robust */
            ST_double *cluster_ze = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
            if (!cluster_ze) {
                free(g); free(ZOmegaZ); free(ZOmegaZ_inv);
                return 1e30;  /* Return large value on allocation failure */
            }
            for (i = 0; i < N; i++) {
                ST_int c = cluster_ids[i] - 1;
                if (c < 0 || c >= num_clusters) continue;
                for (j = 0; j < K_iv; j++) {
                    cluster_ze[c * K_iv + j] += Z[j * N + i] * work_resid[i];
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
            /* Heteroskedastic robust */
            for (i = 0; i < N; i++) {
                ST_double e2 = work_resid[i] * work_resid[i];
                for (j = 0; j < K_iv; j++) {
                    ST_double z_j_e2 = Z[j * N + i] * e2;
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
            free(g); free(ZOmegaZ); free(ZOmegaZ_inv);
            return 1e30;  /* Singular, return large value */
        }
        invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv);

        /* Compute g'W^{-1}g */
        ST_double *temp = (ST_double *)calloc(K_iv, sizeof(ST_double));
        if (!temp) {
            free(g); free(ZOmegaZ); free(ZOmegaZ_inv);
            return 1e30;
        }
        for (i = 0; i < K_iv; i++) {
            for (j = 0; j < K_iv; j++) {
                temp[i] += ZOmegaZ_inv[j * K_iv + i] * g[j];
            }
        }
        ST_double gWinvg = 0.0;
        for (i = 0; i < K_iv; i++) {
            gWinvg += g[i] * temp[i];
        }

        Q = gWinvg / N;

        free(temp);
        free(ZOmegaZ);
        free(ZOmegaZ_inv);
    }

    free(g);
    return Q;
}

/*
    Compute CUE (Continuously Updated Estimator).

    CUE minimizes: Q(β) = g(β)' W(β)^-1 g(β)
    where g(β) = Z'(y - Xβ) and W(β) = Z'Ω(β)Z

    For homoskedastic (vce_type=0): Q = g'(Z'Z)^{-1}g / σ²
    For robust (vce_type=1): Q = g'(Z'diag(e²)Z)^{-1}g
    For cluster (vce_type>=2): Cluster-robust weighting

    Uses gradient descent with numerical differentiation to minimize Q(β).
*/
ST_retcode ivest_compute_cue(
    IVEstContext *ctx,
    const ST_double *y,
    const ST_double *initial_beta,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *beta,
    ST_double *resid,
    ST_int max_iter,
    ST_double tol,
    ST_double *XZWZX_inv_out
)
{
    ST_int N = ctx->N;
    ST_int K_iv = ctx->K_iv;
    ST_int K_total = ctx->K_total;
    const ST_double *Z = ctx->Z;
    ST_int i, j, k;
    ST_int use_robust_weighting = (vce_type > 0);

    /* Initialize beta from input */
    memcpy(beta, initial_beta, K_total * sizeof(ST_double));

    /* Allocate working arrays */
    ST_double *current_resid = (ST_double *)malloc(N * sizeof(ST_double));
    ST_double *beta_test = (ST_double *)malloc(K_total * sizeof(ST_double));
    ST_double *gradient = (ST_double *)calloc(K_total, sizeof(ST_double));

    if (!current_resid || !beta_test || !gradient) {
        free(current_resid); free(beta_test); free(gradient);
        return 920;
    }

    if (ctx->verbose) {
        SF_display("civreghdfe: Computing CUE (gradient descent)\n");
    }

    /* Compute initial objective */
    ST_double Q_current = cue_objective(ctx, y, beta, vce_type, cluster_ids, num_clusters, current_resid);

    if (ctx->verbose) {
        char buf[128];
        snprintf(buf, sizeof(buf), "civreghdfe: CUE initial Q = %.8e\n", Q_current);
        SF_display(buf);
    }

    ST_double eps = 1e-7;  /* Step for numerical gradient */
    ST_int iter;

    for (iter = 0; iter < max_iter; iter++) {
        /* Compute gradient using numerical differentiation */
        for (k = 0; k < K_total; k++) {
            memcpy(beta_test, beta, K_total * sizeof(ST_double));

            /* Forward step */
            beta_test[k] = beta[k] + eps;
            ST_double Q_plus = cue_objective(ctx, y, beta_test, vce_type, cluster_ids, num_clusters, current_resid);

            /* Backward step */
            beta_test[k] = beta[k] - eps;
            ST_double Q_minus = cue_objective(ctx, y, beta_test, vce_type, cluster_ids, num_clusters, current_resid);

            gradient[k] = (Q_plus - Q_minus) / (2.0 * eps);
        }

        /* Compute gradient norm */
        ST_double grad_norm = 0.0;
        for (k = 0; k < K_total; k++) {
            grad_norm += gradient[k] * gradient[k];
        }
        grad_norm = sqrt(grad_norm);

        if (grad_norm < tol) {
            if (ctx->verbose) {
                char buf[128];
                snprintf(buf, sizeof(buf), "civreghdfe: CUE converged in %d iterations (grad=%.2e)\n", iter, grad_norm);
                SF_display(buf);
            }
            break;
        }

        /* Line search along negative gradient direction */
        ST_double step = 1.0;
        ST_double armijo_c = 1e-4;
        ST_int line_search_iter = 0;
        ST_int max_ls_iter = 20;
        ST_double Q_new;

        while (line_search_iter < max_ls_iter) {
            /* Try step */
            for (k = 0; k < K_total; k++) {
                beta_test[k] = beta[k] - step * gradient[k];
            }

            Q_new = cue_objective(ctx, y, beta_test, vce_type, cluster_ids, num_clusters, current_resid);

            /* Armijo condition */
            if (Q_new <= Q_current - armijo_c * step * grad_norm * grad_norm) {
                break;
            }

            step *= 0.5;
            line_search_iter++;
        }

        if (line_search_iter >= max_ls_iter) {
            if (ctx->verbose) {
                SF_display("civreghdfe: CUE line search failed, stopping\n");
            }
            break;
        }

        /* Update beta */
        memcpy(beta, beta_test, K_total * sizeof(ST_double));
        Q_current = Q_new;
    }

    if (iter >= max_iter && ctx->verbose) {
        SF_display("civreghdfe: CUE reached max iterations\n");
    }

    free(gradient);
    free(beta_test);

    /* Continue with existing code for residuals and Hessian computation */
    /* Variables needed for the remaining code */
    ST_double *beta_old = (ST_double *)malloc(K_total * sizeof(ST_double));
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));

    if (!beta_old || !beta_temp) {
        free(current_resid); free(beta_old); free(beta_temp);
        return 920;
    }

    /* Compute final residuals */
    cue_compute_residuals(y, ctx->X_all, beta, N, K_total, current_resid);

    /* Now compute the final Hessian for VCE - use existing code structure */
    /* Compute optimal weighting matrix based on final residuals */
    ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    ST_double *ZOmegaZ_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));

    if (!ZOmegaZ || !ZOmegaZ_inv) {
        free(ZOmegaZ); free(ZOmegaZ_inv);
        free(current_resid); free(beta_old); free(beta_temp);
        return 920;
    }

    if (!use_robust_weighting) {
        /* Homoskedastic: W = Z'Z (use pre-computed ZtZ) */
        memcpy(ZOmegaZ, ctx->ZtZ, K_iv * K_iv * sizeof(ST_double));
    } else if (cluster_ids && num_clusters > 0) {
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
        /* Heteroskedastic robust */
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
        if (ctx->verbose) SF_display("civreghdfe: CUE final weighting matrix singular\n");
        free(ZOmegaZ); free(ZOmegaZ_inv);
        free(current_resid); free(beta_old); free(beta_temp);
        return 920;
    }
    invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv);

    /* Copy residuals to output if requested */
    if (resid) {
        memcpy(resid, current_resid, N * sizeof(ST_double));
    }

    /* Compute (X'ZWZ'X)^-1 for VCE using the already computed ZOmegaZ_inv */
    if (XZWZX_inv_out) {
        /* Compute X'ZWZ'X */
        ST_double *XZW_final = (ST_double *)calloc(K_total * K_iv, sizeof(ST_double));
        ST_double *XZWZX_final = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));

        if (XZW_final && XZWZX_final) {
            for (i = 0; i < K_total; i++) {
                for (j = 0; j < K_iv; j++) {
                    ST_double sum = 0.0;
                    for (k = 0; k < K_iv; k++) {
                        sum += ctx->ZtX[i * K_iv + k] * ZOmegaZ_inv[j * K_iv + k];
                    }
                    XZW_final[j * K_total + i] = sum;
                }
            }
            for (i = 0; i < K_total; i++) {
                for (j = 0; j < K_total; j++) {
                    ST_double sum = 0.0;
                    for (k = 0; k < K_iv; k++) {
                        sum += XZW_final[k * K_total + i] * ctx->ZtX[j * K_iv + k];
                    }
                    XZWZX_final[j * K_total + i] = sum;
                }
            }

            /* Invert X'ZWZ'X */
            if (cholesky(XZWZX_final, K_total) == 0) {
                invert_from_cholesky(XZWZX_final, K_total, XZWZX_inv_out);
            }
        }
        if (XZW_final) free(XZW_final);
        if (XZWZX_final) free(XZWZX_final);
    }

    free(ZOmegaZ);
    free(ZOmegaZ_inv);

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
/*
    Compute k-class IV estimation (includes 2SLS, LIML, Fuller, etc.)

    The k-class estimator is:
    beta = ((1-k)*X'X + k*X'P_Z X)^-1 ((1-k)*X'y + k*X'P_Z y)
    where P_Z = Z(Z'Z)^-1 Z' is the projection onto the instrument space

    When k=1, this reduces to 2SLS: beta = (X'P_Z X)^-1 X'P_Z y
    When k=lambda (min eigenvalue), this is LIML
    When k=lambda - alpha/(N-K_iv), this is Fuller LIML

    Parameters:
    - kclass: the k value (1.0 for 2SLS, lambda for LIML, etc.)
    - est_method: 0=2SLS, 1=LIML, 2=Fuller, 3=kclass, 4=GMM2S, 5=CUE
*/
ST_retcode ivest_compute_2sls(
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
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int df_a,
    ST_int nested_adj,
    ST_int verbose,
    ST_int est_method,
    ST_double kclass_user,
    ST_double fuller_alpha,
    ST_double *lambda_out,
    ST_int kernel_type,
    ST_int bw,
    ST_int kiefer,
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels,
    ST_int sdofminus,
    ST_int center,
    const ST_double *Z_original  /* Original (non-demeaned) Z for Kiefer VCE */
)
{
    ST_int K_total = K_exog + K_endog;  /* Total regressors */
    ST_int i, j, k;

    /* Suppress unused variable warnings for not-yet-implemented features */
    /* TODO: center is for centering HAC score vectors before outer product */
    (void)center;

    /* sdofminus is the FULL df_a (including all absorbed FE), used for test statistics DOF.
       This differs from the df_a parameter which is df_a_for_vce (excluding nested FE).
       For VCE computation, use df_a. For test statistics, use sdofminus. */
    ST_int df_a_full = (sdofminus > 0) ? sdofminus : df_a;

    if (verbose) {
        char buf[256];
        snprintf(buf, sizeof(buf), "civreghdfe: df_a=%d, sdofminus=%d, df_a_full=%d\n",
                 (int)df_a, (int)sdofminus, (int)df_a_full);
        SF_display(buf);
    }

    /* Determine k value based on estimation method */
    ST_double kclass = 1.0;  /* Default: 2SLS */
    ST_double lambda = 1.0;

    if (est_method == 1 || est_method == 2) {
        /* LIML or Fuller: compute lambda */
        lambda = compute_liml_lambda(y, X_endog, X_exog, Z, N, K_exog, K_endog, K_iv);

        /* For exactly identified models, lambda should be 1 */
        if (K_iv == K_total) {
            lambda = 1.0;
        }

        if (est_method == 1) {
            /* LIML */
            kclass = lambda;
        } else {
            /* Fuller: k = lambda - alpha/(N - K_iv) */
            if (fuller_alpha > (N - K_iv)) {
                SF_error("civreghdfe: Invalid Fuller parameter\n");
                return 198;
            }
            kclass = lambda - fuller_alpha / (ST_double)(N - K_iv);
        }

        if (lambda_out) *lambda_out = lambda;
    }
    else if (est_method == 3) {
        /* User-specified k-class */
        kclass = kclass_user;
    }
    /* est_method 0: k=1 (2SLS)
       est_method 4: GMM2S - will be computed after initial 2SLS
       est_method 5: CUE - not yet implemented */

    const char *method_name = "2SLS";
    if (est_method == 1) method_name = "LIML";
    else if (est_method == 2) method_name = "Fuller LIML";
    else if (est_method == 3) method_name = "k-class";
    else if (est_method == 4) method_name = "GMM2S";
    else if (est_method == 5) method_name = "CUE";

    if (verbose) {
        char buf[256];
        snprintf(buf, sizeof(buf), "civreghdfe: Computing %s estimation\n", method_name);
        SF_display(buf);
        snprintf(buf, sizeof(buf), "  N=%d, K_exog=%d, K_endog=%d, K_iv=%d\n",
                 (int)N, (int)K_exog, (int)K_endog, (int)K_iv);
        SF_display(buf);
        if (est_method == 1 || est_method == 2) {
            snprintf(buf, sizeof(buf), "  lambda=%.6f, k=%.6f\n", lambda, kclass);
            SF_display(buf);
        } else if (est_method == 3) {
            snprintf(buf, sizeof(buf), "  k=%.6f\n", kclass);
            SF_display(buf);
        }
    }

    /* Check identification */
    if (K_iv < K_total) {
        SF_error("civreghdfe: Model is underidentified\n");
        return 198;
    }

    /* Allocate work arrays (with overflow checks) */
    ST_double *ZtZ = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
    ST_double *ZtZ_inv = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_iv, sizeof(ST_double));
    ST_double *ZtX = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_total, sizeof(ST_double));
    ST_double *Zty = (ST_double *)ctools_safe_calloc2((size_t)K_iv, sizeof(ST_double));
    ST_double *XtPzX = (ST_double *)ctools_safe_calloc3((size_t)K_total, (size_t)K_total, sizeof(ST_double));
    ST_double *XtPzy = (ST_double *)ctools_safe_calloc2((size_t)K_total, sizeof(ST_double));
    ST_double *temp1 = (ST_double *)ctools_safe_calloc3((size_t)K_iv, (size_t)K_total, sizeof(ST_double));
    ST_double *X_all = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_total, sizeof(ST_double));
    ST_double *resid = (ST_double *)ctools_safe_malloc2((size_t)N, sizeof(ST_double));

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
        size_t exog_size;
        if (ctools_safe_alloc_size((size_t)N, (size_t)K_exog, sizeof(ST_double), &exog_size) != 0) {
            SF_error("civreghdfe: Size overflow in memcpy\n");
            free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
            free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
            return 920;
        }
        memcpy(X_all, X_exog, exog_size);
    }
    if (K_endog > 0 && X_endog) {
        size_t endog_size;
        if (ctools_safe_alloc_size((size_t)N, (size_t)K_endog, sizeof(ST_double), &endog_size) != 0) {
            SF_error("civreghdfe: Size overflow in memcpy\n");
            free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
            free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
            return 920;
        }
        memcpy(X_all + (size_t)N * (size_t)K_exog, X_endog, endog_size);
    }

    /* Step 1: Compute Z'Z (weighted if needed) */
    if (weights && weight_type != 0) {
        matmul_atdb(Z, Z, weights, N, K_iv, K_iv, ZtZ);
    } else {
        matmul_atb(Z, Z, N, K_iv, K_iv, ZtZ);
    }

    /* Step 2: Invert Z'Z */
    {
        size_t ztz_size;
        ctools_safe_alloc_size((size_t)K_iv, (size_t)K_iv, sizeof(ST_double), &ztz_size);
        memcpy(ZtZ_inv, ZtZ, ztz_size);
    }
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
    if (!ZtZ_inv_Zty) {
        SF_error("civreghdfe: Memory allocation failed for ZtZ_inv_Zty\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        return 920;
    }
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

    /* Step 6b: For k-class estimation (k != 1), compute X'X and X'y */
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_double *XkX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *Xky = (ST_double *)calloc(K_total, sizeof(ST_double));

    if (!XkX || !Xky) {
        SF_error("civreghdfe: Memory allocation failed for XkX/Xky\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XkX); free(Xky);
        return 920;
    }

    if (kclass != 1.0) {
        XtX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        Xty = (ST_double *)calloc(K_total, sizeof(ST_double));

        /* Compute X'X */
        if (weights && weight_type != 0) {
            matmul_atdb(X_all, X_all, weights, N, K_total, K_total, XtX);
        } else {
            matmul_atb(X_all, X_all, N, K_total, K_total, XtX);
        }

        /* Compute X'y */
        if (weights && weight_type != 0) {
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = X_all + i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * weights[k] * y[k];
                }
                Xty[i] = sum;
            }
        } else {
            for (i = 0; i < K_total; i++) {
                ST_double sum = 0.0;
                const ST_double *x_col = X_all + i * N;
                for (k = 0; k < N; k++) {
                    sum += x_col[k] * y[k];
                }
                Xty[i] = sum;
            }
        }

        /* Compute k-class weighted matrices:
           XkX = (1-k)*X'X + k*X'P_Z X
           Xky = (1-k)*X'y + k*X'P_Z y */
        ST_double one_minus_k = 1.0 - kclass;
        for (i = 0; i < K_total * K_total; i++) {
            XkX[i] = one_minus_k * XtX[i] + kclass * XtPzX[i];
        }
        for (i = 0; i < K_total; i++) {
            Xky[i] = one_minus_k * Xty[i] + kclass * XtPzy[i];
        }
    } else {
        /* k = 1 (2SLS): just use X'P_Z X and X'P_Z y */
        memcpy(XkX, XtPzX, K_total * K_total * sizeof(ST_double));
        memcpy(Xky, XtPzy, K_total * sizeof(ST_double));
    }

    /* Step 7: Solve XkX * beta = Xky */
    ST_double *XkX_copy = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
    if (!XkX_copy) {
        SF_error("civreghdfe: Memory allocation failed for XkX_copy\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty);
        if (XtX) free(XtX);
        if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 920;
    }
    memcpy(XkX_copy, XkX, K_total * K_total * sizeof(ST_double));

    /* Use Cholesky solve */
    if (cholesky(XkX_copy, K_total) != 0) {
        SF_error("civreghdfe: XkX matrix is singular\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XkX_copy);
        if (XtX) free(XtX);
        if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 198;
    }

    /* Forward substitution */
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));
    if (!beta_temp) {
        SF_error("civreghdfe: Memory allocation failed for beta_temp\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XkX_copy);
        if (XtX) free(XtX);
        if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 920;
    }
    for (i = 0; i < K_total; i++) {
        ST_double sum = Xky[i];
        for (j = 0; j < i; j++) {
            sum -= XkX_copy[i * K_total + j] * beta_temp[j];
        }
        beta_temp[i] = sum / XkX_copy[i * K_total + i];
    }

    /* Backward substitution */
    for (i = K_total - 1; i >= 0; i--) {
        ST_double sum = beta_temp[i];
        for (j = i + 1; j < K_total; j++) {
            sum -= XkX_copy[j * K_total + i] * beta[j];
        }
        beta[i] = sum / XkX_copy[i * K_total + i];
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
    /* OPTIMIZED: Parallelized with OpenMP */
    #pragma omp parallel for schedule(static) if(N > 10000)
    for (i = 0; i < N; i++) {
        ST_double pred = 0.0;
        #pragma omp simd reduction(+:pred)
        for (k = 0; k < K_total; k++) {
            pred += X_all[k * N + i] * beta[k];
        }
        resid[i] = y[i] - pred;
    }

    /* Allocate GMM/CUE Hessian inverse for VCE computation */
    ST_double *gmm_hessian_inv = NULL;

    /* Step 8b: For GMM2S, re-estimate with optimal weighting matrix */
    if (est_method == 4 && vce_type > 0) {
        /*
            Two-step efficient GMM using refactored helper from civreghdfe_estimate.c
            Step 1: 2SLS to get initial residuals (done above)
            Step 2: Compute optimal weighting matrix W = (Z'ΩZ)^-1 where Ω = diag(e²)
            Step 3: Re-estimate β = (X'ZWZ'X)^-1 X'ZWZ'y
        */

        /* Allocate storage for GMM Hessian inverse */
        gmm_hessian_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));

        /* Create estimation context for GMM2S computation */
        IVEstContext gmm_ctx;
        gmm_ctx.N = N;
        gmm_ctx.K_exog = K_exog;
        gmm_ctx.K_endog = K_endog;
        gmm_ctx.K_iv = K_iv;
        gmm_ctx.K_total = K_total;
        gmm_ctx.weights = weights;
        gmm_ctx.weight_type = weight_type;
        gmm_ctx.Z = Z;
        gmm_ctx.y = y;
        gmm_ctx.ZtZ = ZtZ;
        gmm_ctx.ZtZ_inv = ZtZ_inv;
        gmm_ctx.ZtX = ZtX;
        gmm_ctx.Zty = Zty;
        gmm_ctx.XtPzX = XtPzX;
        gmm_ctx.XtPzy = XtPzy;
        gmm_ctx.temp_kiv_ktotal = temp1;
        gmm_ctx.X_all = X_all;
        gmm_ctx.verbose = verbose;

        ST_retcode gmm_rc = ivest_compute_gmm2s(
            &gmm_ctx, y, resid, cluster_ids, num_clusters, beta, resid, gmm_hessian_inv
        );

        if (gmm_rc != STATA_OK) {
            if (verbose) SF_display("civreghdfe: GMM2S re-estimation failed, using 2SLS\n");
            if (gmm_hessian_inv) { free(gmm_hessian_inv); gmm_hessian_inv = NULL; }
        } else {
            if (verbose) {
                char buf[256];
                SF_display("civreghdfe: GMM2S estimates:\n");
                for (i = 0; i < K_total; i++) {
                    snprintf(buf, sizeof(buf), "  beta[%d] = %g\n", (int)i, beta[i]);
                    SF_display(buf);
                }
            }
        }

        /* Note: gmm_ctx uses pointers to existing arrays, no need to free context members */
    }

    /* Step 8c: CUE (Continuously Updated Estimator) using refactored helper */
    if (est_method == 5) {
        /*
            CUE using refactored helper from civreghdfe_estimate.c
            Iteratively re-weights until convergence.
        */

        /* Allocate storage for CUE Hessian inverse */
        gmm_hessian_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));

        /* Create estimation context for CUE computation */
        IVEstContext cue_ctx;
        cue_ctx.N = N;
        cue_ctx.K_exog = K_exog;
        cue_ctx.K_endog = K_endog;
        cue_ctx.K_iv = K_iv;
        cue_ctx.K_total = K_total;
        cue_ctx.weights = weights;
        cue_ctx.weight_type = weight_type;
        cue_ctx.Z = Z;
        cue_ctx.y = y;
        cue_ctx.ZtZ = ZtZ;
        cue_ctx.ZtZ_inv = ZtZ_inv;
        cue_ctx.ZtX = ZtX;
        cue_ctx.Zty = Zty;
        cue_ctx.XtPzX = XtPzX;
        cue_ctx.XtPzy = XtPzy;
        cue_ctx.temp_kiv_ktotal = temp1;
        cue_ctx.X_all = X_all;
        cue_ctx.verbose = verbose;

        const ST_int max_cue_iter = 100;
        const ST_double cue_tol = 1e-9;

        ST_retcode cue_rc = ivest_compute_cue(
            &cue_ctx, y, beta, vce_type, cluster_ids, num_clusters,
            beta, resid, max_cue_iter, cue_tol, gmm_hessian_inv
        );

        if (cue_rc != STATA_OK) {
            if (verbose) SF_display("civreghdfe: CUE computation failed, using 2SLS\n");
            if (gmm_hessian_inv) { free(gmm_hessian_inv); gmm_hessian_inv = NULL; }
        } else {
            if (verbose) {
                char buf[256];
                SF_display("civreghdfe: CUE final estimates:\n");
                for (i = 0; i < K_total; i++) {
                    snprintf(buf, sizeof(buf), "  beta[%d] = %g\n", (int)i, beta[i]);
                    SF_display(buf);
                }
            }
        }

        /* Note: cue_ctx uses pointers to existing arrays, no need to free context members */
    }

    /* Step 9: Compute VCE */
    /* For k-class estimators, the VCE is: sigma^2 * (XkX)^-1 */
    /* where sigma^2 = RSS / (N - K) using actual residuals */

    /* First, invert XkX (the k-class weighted matrix) */
    memcpy(XkX_copy, XkX, K_total * K_total * sizeof(ST_double));
    if (cholesky(XkX_copy, K_total) != 0) {
        SF_error("civreghdfe: Cannot compute VCE (XkX singular)\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XkX_copy); free(beta_temp);
        if (gmm_hessian_inv) free(gmm_hessian_inv);
        if (XtX) free(XtX);
        if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 198;
    }

    ST_double *XkX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!XkX_inv) {
        SF_error("civreghdfe: Memory allocation failed for XkX_inv\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XkX_copy); free(beta_temp);
        if (gmm_hessian_inv) free(gmm_hessian_inv);
        if (XtX) free(XtX);
        if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 198;
    }
    invert_from_cholesky(XkX_copy, K_total, XkX_inv);

    /* Compute RSS (needed for Sargan test; VCE helper computes its own internally) */
    ST_double rss = 0.0;
    if (weights && weight_type != 0) {
        for (i = 0; i < N; i++) rss += weights[i] * resid[i] * resid[i];
    } else {
        for (i = 0; i < N; i++) rss += resid[i] * resid[i];
    }

    if (est_method == 4 && gmm_hessian_inv != NULL) {
        /*
           Efficient GMM2S VCE:
           For two-step efficient GMM with optimal weights W = (Z'ΩZ)^-1,
           the asymptotic VCE is simply (X'ZWZ'X)^-1 regardless of whether
           robust or cluster is specified. The optimal weighting already
           accounts for heteroskedasticity/clustering in step 1.

           Apply small-sample DOF correction to match ivreghdfe behavior:
           - For robust: V = (X'ZWZ'X)^-1 * (N / df_r)
           - For cluster: V = (X'ZWZ'X)^-1 * (N/df_r) * (G/(G-1))
           where df_r = N - K_total - df_a
        */
        ST_int df_r_gmm = N - K_total - df_a;
        if (df_r_gmm <= 0) df_r_gmm = 1;
        ST_double dof_adj_gmm;
        if (vce_type == CIVREGHDFE_VCE_CLUSTER && num_clusters > 0) {
            /* Cluster DOF adjustment: N/df_r * G/(G-1) * (1 + K/N²)
               The small K/N² factor matches ivreghdfe's GMM2S cluster VCE */
            ST_double G = (ST_double)num_clusters;
            ST_double small_sample_corr = 1.0 + (ST_double)K_total / ((ST_double)N * (ST_double)N);
            dof_adj_gmm = ((ST_double)N / (ST_double)df_r_gmm) * (G / (G - 1.0)) * small_sample_corr;
        } else {
            /* Robust or unadjusted: N/df_r */
            dof_adj_gmm = (ST_double)N / (ST_double)df_r_gmm;
        }
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = gmm_hessian_inv[i] * dof_adj_gmm;
        }
    } else if (est_method == 5 && gmm_hessian_inv != NULL && vce_type > 0) {
        /*
           CUE VCE (robust/cluster case):
           For Continuously Updated Estimator with optimal weights,
           the asymptotic VCE is (X'ZWZ'X)^-1 where W is evaluated at final β.

           Apply same DOF correction as GMM2S.
           Note: For homoskedastic CUE (vce_type == 0), we fall through to
           the standard VCE computation which uses σ² * (X'PzX)^-1.
        */
        ST_int df_r_cue = N - K_total - df_a;
        if (df_r_cue <= 0) df_r_cue = 1;
        ST_double dof_adj_cue;
        if (vce_type == CIVREGHDFE_VCE_CLUSTER && num_clusters > 0) {
            /* Cluster DOF adjustment: N/df_r * G/(G-1) * (1 + K/N²) */
            ST_double G = (ST_double)num_clusters;
            ST_double small_sample_corr = 1.0 + (ST_double)K_total / ((ST_double)N * (ST_double)N);
            dof_adj_cue = ((ST_double)N / (ST_double)df_r_cue) * (G / (G - 1.0)) * small_sample_corr;
        } else {
            dof_adj_cue = (ST_double)N / (ST_double)df_r_cue;
        }
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = gmm_hessian_inv[i] * dof_adj_cue;
        }
    } else if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && cluster2_ids != NULL && num_clusters2 > 0) {
        /*
           Two-way clustered VCE using Cameron-Gelbach-Miller (2011) formula:
           V = V1 + V2 - V_intersection
        */
        ivvce_compute_twoway(
            Z, resid, temp1, XkX_inv,
            weights, weight_type,
            N, K_total, K_iv,
            cluster_ids, num_clusters,
            cluster2_ids, num_clusters2,
            df_a,
            V
        );
    } else if (kiefer && kernel_type > 0 && bw > 0 && cluster_ids != NULL && num_clusters > 0) {
        /*
           Kiefer (1980) VCE: homoskedastic within-panel autocorrelation.
           Unlike cluster-robust, Kiefer assumes homoskedasticity.
           Uses ivvce_compute_kiefer which implements the correct formula.
        */
        if (verbose) {
            char buf[256];
            snprintf(buf, sizeof(buf), "civreghdfe: Using Kiefer VCE (kernel=%d, bw=%d, panels=%d)\n",
                     (int)kernel_type, (int)bw, (int)num_clusters);
            SF_display(buf);
        }
        /* Use original (non-demeaned) Z for Kiefer VCE
           This ensures sum_t(Z'e) != 0, which is required for the Kiefer formula
           to produce non-zero results */
        const ST_double *Z_for_kiefer = Z_original ? Z_original : Z;
        /* Note: Kiefer VCE with time-clustering is currently not fully compatible
           with FE demeaning. The FE structure causes residuals to have zero mean
           within each panel, which makes sum_t(Z_it * e_it) small when summed
           across panels. This is a known limitation.
           TODO: Investigate ivreg2's approach for proper Kiefer + FE handling. */
        ivvce_compute_kiefer(
            Z_for_kiefer, resid, temp1, XkX_inv,
            weights, weight_type,
            N, K_total, K_iv,
            cluster_ids, num_clusters,
            df_a,
            V
        );
    } else {
        /*
           Standard VCE: Use refactored helper from civreghdfe_vce.c
           Handles unadjusted, robust (HC), HAC, and clustered VCE types.
        */
        /*
           Non-GMM2S VCE: Use refactored helper from civreghdfe_vce.c
           Handles unadjusted, robust (HC), HAC, and clustered VCE types.
        */
        ivvce_compute_full(
            Z, resid, temp1, XkX_inv,
            weights, weight_type,
            N, K_total, K_iv,
            vce_type, cluster_ids, num_clusters,
            df_a, nested_adj, kernel_type, bw,
            hac_panel_ids, num_hac_panels,
            V
        );
    }

    /* Step 10: Compute first-stage F statistics */
    /* For each endogenous variable, compute F-stat from first stage regression */
    /* The F tests whether excluded instruments are jointly significant,
       controlling for exogenous regressors. Uses partial R² approach:
       F = ((R²_full - R²_reduced) / L) / ((1 - R²_full) / df_resid) */
    if (first_stage_F) {
        for (ST_int e = 0; e < K_endog; e++) {
            const ST_double *X_e = X_endog + e * N;

            /* Compute X_e'X_e */
            ST_double xx = 0.0;
            for (i = 0; i < N; i++) {
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                xx += w * X_e[i] * X_e[i];
            }
            if (xx <= 0) xx = 1.0;

            /* Compute R²_full: R² from projecting X_e onto all instruments Z */
            ST_double xpx_full = 0.0;
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
                xpx_full += w * X_e[i] * pz_xe;
            }
            ST_double r2_full = xpx_full / xx;
            if (r2_full > 1.0) r2_full = 1.0;
            if (r2_full < 0.0) r2_full = 0.0;

            /* Compute R²_reduced: R² from projecting X_e onto only exogenous regressors */
            ST_double r2_reduced = 0.0;
            if (K_exog > 0) {
                /* Build Z_exog'Z_exog (K_exog x K_exog) from ZtZ */
                ST_double *ZeZe = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                ST_double *ZeZe_inv = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                ST_double *ZeXe = (ST_double *)calloc(K_exog, sizeof(ST_double));

                if (ZeZe && ZeZe_inv && ZeXe) {
                    /* Extract Z_exog'Z_exog from ZtZ */
                    for (ST_int l1 = 0; l1 < K_exog; l1++) {
                        for (ST_int l2 = 0; l2 < K_exog; l2++) {
                            ZeZe[l2 * K_exog + l1] = ZtZ[l2 * K_iv + l1];
                        }
                    }

                    /* Compute Z_exog'X_e */
                    for (ST_int l = 0; l < K_exog; l++) {
                        ST_double sum = 0.0;
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            sum += w * Z[l * N + i] * X_e[i];
                        }
                        ZeXe[l] = sum;
                    }

                    /* Invert Z_exog'Z_exog */
                    memcpy(ZeZe_inv, ZeZe, K_exog * K_exog * sizeof(ST_double));
                    if (cholesky(ZeZe_inv, K_exog) == 0 &&
                        invert_from_cholesky(ZeZe_inv, K_exog, ZeZe_inv) == 0) {

                        /* Compute X_e'P_exog X_e = ZeXe' * ZeZe_inv * ZeXe */
                        ST_double xpx_reduced = 0.0;
                        for (ST_int l1 = 0; l1 < K_exog; l1++) {
                            ST_double temp = 0.0;
                            for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                temp += ZeZe_inv[l2 * K_exog + l1] * ZeXe[l2];
                            }
                            xpx_reduced += ZeXe[l1] * temp;
                        }
                        r2_reduced = xpx_reduced / xx;
                        if (r2_reduced > 1.0) r2_reduced = 1.0;
                        if (r2_reduced < 0.0) r2_reduced = 0.0;
                    }
                }

                free(ZeZe);
                free(ZeZe_inv);
                free(ZeXe);
            }

            /* Partial R² = R²_full - R²_reduced */
            ST_double partial_r2 = r2_full - r2_reduced;
            if (partial_r2 < 0.0) partial_r2 = 0.0;

            /* F = (partial_R² / L) / ((1 - R²_full) / df_resid) */
            ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
            if (L <= 0) L = 1;
            ST_int denom_df = N - K_iv - df_a;
            if (denom_df <= 0) denom_df = 1;

            ST_double denom = (1.0 - r2_full) / (ST_double)denom_df;
            if (denom <= 0.0) denom = 1e-10;

            first_stage_F[e] = (partial_r2 / (ST_double)L) / denom;
            if (first_stage_F[e] < 0) first_stage_F[e] = 0;

            /* Save partial R² for ffirst display */
            if (!g_civreghdfe_noreturn) {
                char r2_name[64];
                snprintf(r2_name, sizeof(r2_name), "__civreghdfe_partial_r2_%d", (int)(e + 1));
                ctools_scal_save(r2_name, partial_r2);
            }

            if (verbose) {
                char buf[256];
                snprintf(buf, sizeof(buf), "  First-stage F[%d] = %g (partial R² = %g, R²_full = %g)\n",
                         (int)e, first_stage_F[e], partial_r2, r2_full);
                SF_display(buf);
            }
        }
    }

    /* Step 10b: Compute underidentification test and weak instrument stats */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double underid_stat = 0.0;
    ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
    ST_int underid_df = L;
    ST_double cd_f = 0.0;       /* Cragg-Donald Wald F (homoskedastic) */
    ST_double kp_f = 0.0;       /* Kleibergen-Paap rk Wald F (robust) */

    civreghdfe_compute_underid_test(
        X_endog, Z, ZtZ, ZtZ_inv, temp1, first_stage_F,
        weights, weight_type, N, K_exog, K_endog, K_iv, df_a_full,
        vce_type, cluster_ids, num_clusters,
        kernel_type, bw, kiefer,
        hac_panel_ids, num_hac_panels,
        &underid_stat, &underid_df, &cd_f, &kp_f
    );

    /* Store underidentification test result */
    if (!g_civreghdfe_noreturn) {
        ctools_scal_save("__civreghdfe_underid", underid_stat);
        ctools_scal_save("__civreghdfe_underid_df", (ST_double)underid_df);
    }

    /* Step 11: Compute Sargan/Hansen J overidentification test */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double sargan_stat = 0.0;
    ST_int overid_df = K_iv - K_total;

    civreghdfe_compute_sargan_j(
        resid, Z, ZtZ_inv, rss,
        weights, weight_type, N, K_exog, K_endog, K_iv,
        vce_type, cluster_ids, num_clusters,
        kernel_type, bw, kiefer,
        hac_panel_ids, num_hac_panels,
        &sargan_stat, &overid_df
    );

    /* Store diagnostic statistics as Stata scalars */
    if (!g_civreghdfe_noreturn) {
        ctools_scal_save("__civreghdfe_sargan", sargan_stat);
        ctools_scal_save("__civreghdfe_sargan_df", (ST_double)overid_df);
        ctools_scal_save("__civreghdfe_cd_f", cd_f);
        ctools_scal_save("__civreghdfe_kp_f", kp_f);
    }

    /* Step 12: Compute Durbin-Wu-Hausman endogeneity test */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double endog_chi2 = 0.0;
    ST_double endog_f = 0.0;
    ST_int endog_df = K_endog;

    civreghdfe_compute_dwh_test(
        y, X_exog, X_endog, Z, temp1,
        N, K_exog, K_endog, K_iv, df_a_full,
        &endog_chi2, &endog_f, &endog_df
    );

    if (!g_civreghdfe_noreturn) {
        ctools_scal_save("__civreghdfe_endog_chi2", endog_chi2);
        ctools_scal_save("__civreghdfe_endog_f", endog_f);
        ctools_scal_save("__civreghdfe_endog_df", (ST_double)endog_df);
    }

    /* Step 13: Compute optional diagnostic tests (orthog, endogtest, redundant) */
    ST_double dval_test;
    ST_int n_orthog = 0;
    if (SF_scal_use("__civreghdfe_n_orthog", &dval_test) == 0) {
        n_orthog = (ST_int)dval_test;
    }

    if (n_orthog > 0) {
        /* Read orthog indices */
        ST_int *orthog_indices = (ST_int *)malloc(n_orthog * sizeof(ST_int));
        if (orthog_indices) {
            char scal_name[64];
            for (ST_int oi = 0; oi < n_orthog; oi++) {
                snprintf(scal_name, sizeof(scal_name), "__civreghdfe_orthog_%d", (int)(oi + 1));
                if (SF_scal_use(scal_name, &dval_test) == 0) {
                    orthog_indices[oi] = (ST_int)dval_test;
                }
            }

            ST_double cstat = 0.0;
            ST_int cstat_df = n_orthog;

            civreghdfe_compute_cstat(
                y, X_exog, X_endog, Z, N, K_exog, K_endog, K_iv,
                orthog_indices, n_orthog,
                vce_type, weights, weight_type, cluster_ids, num_clusters,
                sargan_stat, rss, &cstat, &cstat_df
            );

            if (!g_civreghdfe_noreturn) {
                ctools_scal_save("__civreghdfe_cstat", cstat);
                ctools_scal_save("__civreghdfe_cstat_df", (ST_double)cstat_df);
            }

            free(orthog_indices);
        }
    }

    ST_int n_endogtest = 0;
    if (SF_scal_use("__civreghdfe_n_endogtest", &dval_test) == 0) {
        n_endogtest = (ST_int)dval_test;
    }

    if (n_endogtest > 0) {
        /* Read endogtest indices */
        ST_int *endogtest_indices = (ST_int *)malloc(n_endogtest * sizeof(ST_int));
        if (endogtest_indices) {
            char scal_name[64];
            for (ST_int ei = 0; ei < n_endogtest; ei++) {
                snprintf(scal_name, sizeof(scal_name), "__civreghdfe_endogtest_%d", (int)(ei + 1));
                if (SF_scal_use(scal_name, &dval_test) == 0) {
                    endogtest_indices[ei] = (ST_int)dval_test;
                }
            }

            ST_double endogtest_stat = 0.0;
            ST_int endogtest_df_out = n_endogtest;

            civreghdfe_compute_endogtest_subset(
                y, X_exog, X_endog, Z, resid, ZtZ_inv,
                N, K_exog, K_endog, K_iv,
                endogtest_indices, n_endogtest,
                &endogtest_stat, &endogtest_df_out
            );

            if (!g_civreghdfe_noreturn) {
                ctools_scal_save("__civreghdfe_endogtest_stat", endogtest_stat);
                ctools_scal_save("__civreghdfe_endogtest_df", (ST_double)endogtest_df_out);
            }

            free(endogtest_indices);
        }
    }

    ST_int n_redundant = 0;
    if (SF_scal_use("__civreghdfe_n_redundant", &dval_test) == 0) {
        n_redundant = (ST_int)dval_test;
    }

    if (n_redundant > 0) {
        /* Read redundant indices */
        ST_int *redund_indices = (ST_int *)malloc(n_redundant * sizeof(ST_int));
        if (redund_indices) {
            char scal_name[64];
            for (ST_int ri = 0; ri < n_redundant; ri++) {
                snprintf(scal_name, sizeof(scal_name), "__civreghdfe_redundant_%d", (int)(ri + 1));
                if (SF_scal_use(scal_name, &dval_test) == 0) {
                    redund_indices[ri] = (ST_int)dval_test;
                }
            }

            ST_double redund_stat = 0.0;
            ST_int redund_df = K_endog * n_redundant;

            civreghdfe_compute_redundant(
                X_endog, Z, N, K_exog, K_endog, K_iv,
                redund_indices, n_redundant,
                &redund_stat, &redund_df
            );

            if (!g_civreghdfe_noreturn) {
                ctools_scal_save("__civreghdfe_redund_stat", redund_stat);
                ctools_scal_save("__civreghdfe_redund_df", (ST_double)redund_df);
            }

            free(redund_indices);
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
    free(XkX_copy);
    free(beta_temp);
    free(XkX_inv);
    if (gmm_hessian_inv) free(gmm_hessian_inv);
    if (XtX) free(XtX);
    if (Xty) free(Xty);
    free(XkX);
    free(Xky);

    return STATA_OK;
}
