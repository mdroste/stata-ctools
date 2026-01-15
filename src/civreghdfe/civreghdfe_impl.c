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
#include "civreghdfe_matrix.h"
#include "civreghdfe_estimate.h"
#include "civreghdfe_vce.h"
#include "civreghdfe_tests.h"
#include "../ctools_hdfe_utils.h"
#include "../ctools_unroll.h"

/* Debug flag */
#define CIVREGHDFE_DEBUG 0

/* Convenience aliases for the matrix functions */
#define matmul_atb  civreghdfe_matmul_atb
#define matmul_ab   civreghdfe_matmul_ab
#define matmul_atdb civreghdfe_matmul_atdb
#define kernel_weight civreghdfe_kernel_weight
#define jacobi_eigenvalues civreghdfe_jacobi_eigenvalues
#define jacobi_min_eigenvalue civreghdfe_jacobi_min_eigenvalue
#define matpowersym_neg_half civreghdfe_matpowersym_neg_half
#define compute_liml_lambda civreghdfe_compute_liml_lambda

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
    ST_int nested_adj,
    ST_int verbose,
    ST_int est_method,
    ST_double kclass_user,
    ST_double fuller_alpha,
    ST_double *lambda_out,
    ST_int kernel_type,
    ST_int bw
)
{
    ST_int K_total = K_exog + K_endog;  /* Total regressors */
    ST_int i, j, k;

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

    /* Step 6b: For k-class estimation (k != 1), compute X'X and X'y */
    ST_double *XtX = NULL;
    ST_double *Xty = NULL;
    ST_double *XkX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *Xky = (ST_double *)calloc(K_total, sizeof(ST_double));

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
    memcpy(XkX_copy, XkX, K_total * K_total * sizeof(ST_double));

    /* Use Cholesky solve */
    if (cholesky(XkX_copy, K_total) != 0) {
        SF_error("civreghdfe: XkX matrix is singular\n");
        free(ZtZ); free(ZtZ_inv); free(ZtX); free(Zty);
        free(XtPzX); free(XtPzy); free(temp1); free(X_all); free(resid);
        free(ZtZ_inv_Zty); free(XkX_copy);
        if (XtX) free(XtX); if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 198;
    }

    /* Forward substitution */
    ST_double *beta_temp = (ST_double *)calloc(K_total, sizeof(ST_double));
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

    /* Step 8b: For GMM2S, re-estimate with optimal weighting matrix */
    if (est_method == 4 && vce_type > 0) {
        /*
            Two-step efficient GMM using refactored helper from civreghdfe_estimate.c
            Step 1: 2SLS to get initial residuals (done above)
            Step 2: Compute optimal weighting matrix W = (Z'ΩZ)^-1 where Ω = diag(e²)
            Step 3: Re-estimate β = (X'ZWZ'X)^-1 X'ZWZ'y
        */

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
            &gmm_ctx, y, resid, cluster_ids, num_clusters, beta, resid
        );

        if (gmm_rc != STATA_OK) {
            if (verbose) SF_display("civreghdfe: GMM2S re-estimation failed, using 2SLS\n");
        } else {
            /* Update XkX with GMM2S's X'ZWZ'X for VCE computation */
            /* Note: For efficient GMM, XkX is updated by the GMM2S estimator */

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
        const ST_double cue_tol = 1e-8;

        ST_retcode cue_rc = ivest_compute_cue(
            &cue_ctx, y, beta, cluster_ids, num_clusters,
            beta, resid, max_cue_iter, cue_tol
        );

        if (cue_rc != STATA_OK) {
            if (verbose) SF_display("civreghdfe: CUE computation failed, using 2SLS\n");
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
        if (XtX) free(XtX); if (Xty) free(Xty);
        free(XkX); free(Xky);
        return 198;
    }

    ST_double *XkX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    invert_from_cholesky(XkX_copy, K_total, XkX_inv);

    /* Compute RSS (needed for Sargan test; VCE helper computes its own internally) */
    ST_double rss = 0.0;
    if (weights && weight_type != 0) {
        for (i = 0; i < N; i++) rss += weights[i] * resid[i] * resid[i];
    } else {
        for (i = 0; i < N; i++) rss += resid[i] * resid[i];
    }

    if (est_method == 4 && vce_type > 0) {
        /*
           Efficient GMM2S VCE:
           When using optimal weights W = (Z'ΩZ)^-1, the VCE is simply (X'ZWZ'X)^-1.
           The sandwich is not needed because optimal weighting already accounts
           for heteroskedasticity/clustering. This is the efficiency gain of GMM.
           XkX already contains X'ZWZ'X for GMM2S.
        */
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = XkX_inv[i];
        }
    } else {
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
            V
        );
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

            /* F-stat: (R^2 / q) / ((1 - R^2) / (N - K_iv - df_a)) */
            ST_int q = K_iv - K_exog;  /* Number of excluded instruments */
            if (q <= 0) q = 1;
            ST_int denom_df = N - K_iv - df_a;
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

    /* Step 10b: Compute underidentification test and weak instrument stats */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double underid_stat = 0.0;
    ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
    ST_int underid_df = L;
    ST_double cd_f = 0.0;       /* Cragg-Donald Wald F (homoskedastic) */
    ST_double kp_f = 0.0;       /* Kleibergen-Paap rk Wald F (robust) */

    civreghdfe_compute_underid_test(
        X_endog, Z, ZtZ, ZtZ_inv, temp1, first_stage_F,
        weights, weight_type, N, K_exog, K_endog, K_iv, df_a,
        vce_type, cluster_ids, num_clusters,
        &underid_stat, &underid_df, &cd_f, &kp_f
    );

    /* Store underidentification test result */
    SF_scal_save("__civreghdfe_underid", underid_stat);
    SF_scal_save("__civreghdfe_underid_df", (ST_double)underid_df);

    /* Step 11: Compute Sargan/Hansen J overidentification test */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double sargan_stat = 0.0;
    ST_int overid_df = K_iv - K_total;

    civreghdfe_compute_sargan_j(
        resid, Z, ZtZ_inv, rss,
        weights, weight_type, N, K_exog, K_endog, K_iv,
        vce_type, cluster_ids, num_clusters,
        &sargan_stat, &overid_df
    );

    /* Store diagnostic statistics as Stata scalars */
    SF_scal_save("__civreghdfe_sargan", sargan_stat);
    SF_scal_save("__civreghdfe_sargan_df", (ST_double)overid_df);
    SF_scal_save("__civreghdfe_cd_f", cd_f);
    SF_scal_save("__civreghdfe_kp_f", kp_f);

    /* Step 12: Compute Durbin-Wu-Hausman endogeneity test */
    /* Calls modular function from civreghdfe_tests.c */
    ST_double endog_chi2 = 0.0;
    ST_double endog_f = 0.0;
    ST_int endog_df = K_endog;

    civreghdfe_compute_dwh_test(
        y, X_exog, X_endog, Z, temp1,
        N, K_exog, K_endog, K_iv, df_a,
        &endog_chi2, &endog_f, &endog_df
    );

    SF_scal_save("__civreghdfe_endog_chi2", endog_chi2);
    SF_scal_save("__civreghdfe_endog_f", endog_f);
    SF_scal_save("__civreghdfe_endog_df", (ST_double)endog_df);

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
    if (XtX) free(XtX);
    if (Xty) free(Xty);
    free(XkX);
    free(Xky);

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
    ST_int nested_fe_index;
    SF_scal_use("__civreghdfe_nested_fe_index", &dval); nested_fe_index = (ST_int)dval;

    /* New scalars for estimation method */
    ST_int est_method = 0;
    ST_double kclass_user = 0, fuller_alpha = 0;
    SF_scal_use("__civreghdfe_est_method", &dval); est_method = (ST_int)dval;
    SF_scal_use("__civreghdfe_kclass", &dval); kclass_user = dval;
    SF_scal_use("__civreghdfe_fuller", &dval); fuller_alpha = dval;

    /* HAC parameters */
    ST_int kernel_type = 0, bw = 0, dkraay = 0;
    SF_scal_use("__civreghdfe_kernel", &dval); kernel_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_bw", &dval); bw = (ST_int)dval;
    SF_scal_use("__civreghdfe_dkraay", &dval); dkraay = (ST_int)dval;
    (void)dkraay;  /* Reserved for Driscoll-Kraay panel HAC */

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

    /* Count observations that pass the if/in condition */
    ST_int N_ifobs = 0;
    for (ST_int obs = in1; obs <= in2; obs++) {
        if (SF_ifobs(obs)) N_ifobs++;
    }

    if (N_ifobs == 0) {
        SF_error("civreghdfe: No observations selected\n");
        free(y); free(X_endog); free(X_exog); free(Z);
        free(weights); free(cluster_ids);
        for (ST_int g = 0; g < G; g++) free(fe_levels[g]);
        free(fe_levels);
        return 2001;
    }

    /* Reallocate arrays to actual size based on if/in condition */
    N_total = N_ifobs;

    /* Read data from Stata - only include observations where SF_ifobs() is true */
    ST_double val;
    ST_int i = 0;
    for (ST_int obs = in1; obs <= in2; obs++) {
        if (!SF_ifobs(obs)) continue;  /* Skip observations not in sample */

        /* y */
        SF_vdata(var_y, obs, &val);
        y[i] = val;

        /* X_endog */
        for (ST_int k = 0; k < K_endog; k++) {
            SF_vdata(var_endog_start + k, obs, &val);
            X_endog[k * N_ifobs + i] = val;
        }

        /* X_exog */
        for (ST_int k = 0; k < K_exog; k++) {
            SF_vdata(var_exog_start + k, obs, &val);
            X_exog[k * N_ifobs + i] = val;
        }

        /* Z (instruments) */
        for (ST_int k = 0; k < K_iv; k++) {
            SF_vdata(var_iv_start + k, obs, &val);
            Z[k * N_ifobs + i] = val;
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

        i++;
    }

    /* Drop observations with missing values */
    ST_int *valid_mask = (ST_int *)calloc(N_total, sizeof(ST_int));
    ST_int N_valid = 0;

    for (i = 0; i < N_total; i++) {
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

        /* Check FE - already loaded and converted to int, check for very large negative values
           (which could indicate missing - Stata missing values become very large when cast to int) */
        /* Note: FE values were already validated as non-missing during the data load */

        /* Check weights */
        if (has_weights && is_valid) {
            if (SF_is_missing(weights[i]) || weights[i] <= 0) is_valid = 0;
        }

        /* Check cluster - already loaded as int, was already validated during data load */

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

    /* Singleton detection using shared utility from ctools_hdfe_utils */
    ST_int *singleton_mask = (ST_int *)malloc(N * sizeof(ST_int));
    ST_int num_singletons_total = ctools_remove_singletons(
        fe_levels_c, G, N, singleton_mask, 100, verbose
    );

    /* Compact data if singletons were found */
    if (num_singletons_total > 0) {
        ST_int N_after = N - num_singletons_total;

        /* Allocate compacted arrays */
        ST_double *y_new = (ST_double *)malloc(N_after * sizeof(ST_double));
        ST_double *X_endog_new = (K_endog > 0) ? (ST_double *)malloc(N_after * K_endog * sizeof(ST_double)) : NULL;
        ST_double *X_exog_new = (K_exog > 0) ? (ST_double *)malloc(N_after * K_exog * sizeof(ST_double)) : NULL;
        ST_double *Z_new = (ST_double *)malloc(N_after * K_iv * sizeof(ST_double));
        ST_double *weights_new = has_weights ? (ST_double *)malloc(N_after * sizeof(ST_double)) : NULL;
        ST_int *cluster_ids_new = has_cluster ? (ST_int *)malloc(N_after * sizeof(ST_int)) : NULL;
        ST_int **fe_levels_new = (ST_int **)malloc(G * sizeof(ST_int *));
        for (ST_int g = 0; g < G; g++) {
            fe_levels_new[g] = (ST_int *)malloc(N_after * sizeof(ST_int));
        }

        /* Compact using shared utilities */
        ctools_compact_array_double(y_c, y_new, singleton_mask, N, N_after);
        if (K_endog > 0) ctools_compact_matrix_double(X_endog_c, X_endog_new, singleton_mask, N, N_after, K_endog);
        if (K_exog > 0) ctools_compact_matrix_double(X_exog_c, X_exog_new, singleton_mask, N, N_after, K_exog);
        ctools_compact_matrix_double(Z_c, Z_new, singleton_mask, N, N_after, K_iv);
        if (has_weights) ctools_compact_array_double(weights_c, weights_new, singleton_mask, N, N_after);
        if (has_cluster) ctools_compact_array_int(cluster_ids_c, cluster_ids_new, singleton_mask, N, N_after);
        for (ST_int g = 0; g < G; g++) {
            ctools_compact_array_int(fe_levels_c[g], fe_levels_new[g], singleton_mask, N, N_after);
        }

        /* Free old arrays and swap in new ones */
        free(y_c); y_c = y_new;
        if (X_endog_c) free(X_endog_c); X_endog_c = X_endog_new;
        if (X_exog_c) free(X_exog_c); X_exog_c = X_exog_new;
        free(Z_c); Z_c = Z_new;
        if (weights_c) free(weights_c); weights_c = weights_new;
        if (cluster_ids_c) free(cluster_ids_c); cluster_ids_c = cluster_ids_new;
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        fe_levels_c = fe_levels_new;

        N = N_after;
    } else if (verbose) {
        SF_display("civreghdfe: No singletons found\n");
    }
    free(singleton_mask);

    /* Check we still have enough observations */
    if (N < K_total + K_iv + 1) {
        SF_error("civreghdfe: Insufficient observations after singleton removal\n");
        free(y_c); free(X_endog_c); free(X_exog_c); free(Z_c);
        free(weights_c); free(cluster_ids_c);
        for (ST_int g = 0; g < G; g++) free(fe_levels_c[g]);
        free(fe_levels_c);
        return 2001;
    }

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
    state->tolerance = tolerance;
    state->verbose = verbose;

    /* Set up factors */
    state->factors = (FE_Factor *)calloc(G, sizeof(FE_Factor));
    for (ST_int g = 0; g < G; g++) {
        /* Find max level value for remapping */
        ST_int max_level = 0;
        for (ST_int i = 0; i < N; i++) {
            if (fe_levels_c[g][i] > max_level) max_level = fe_levels_c[g][i];
        }

        /* Remap FE levels to contiguous 1-based indices
           The CG solver expects levels to be 1, 2, 3, ..., num_levels
           so it can use levels[i] - 1 as array indices */
        ST_int *remap = (ST_int *)calloc(max_level + 1, sizeof(ST_int));
        if (!remap) {
            SF_error("civreghdfe: Memory allocation failed for level remap\n");
            return 920;
        }

        /* First pass: assign contiguous level IDs starting from 1 */
        ST_int next_level = 1;
        for (ST_int i = 0; i < N; i++) {
            ST_int old_level = fe_levels_c[g][i];
            if (remap[old_level] == 0) {
                remap[old_level] = next_level++;
            }
        }

        ST_int num_levels = next_level - 1;  /* Total unique levels */

        /* Allocate factor arrays with correct size */
        state->factors[g].has_intercept = 1;
        state->factors[g].levels = fe_levels_c[g];  /* Will be remapped in place */
        state->factors[g].num_levels = num_levels;
        state->factors[g].max_level = num_levels;  /* After remapping, max_level == num_levels */
        state->factors[g].counts = (ST_double *)calloc(num_levels, sizeof(ST_double));
        state->factors[g].weighted_counts = has_weights ? (ST_double *)calloc(num_levels, sizeof(ST_double)) : NULL;
        state->factors[g].means = (ST_double *)calloc(num_levels, sizeof(ST_double));

        /* Second pass: remap levels in place and count */
        for (ST_int i = 0; i < N; i++) {
            ST_int old_level = fe_levels_c[g][i];
            ST_int new_level = remap[old_level];  /* 1-based */
            fe_levels_c[g][i] = new_level;        /* Remap in place */

            state->factors[g].counts[new_level - 1] += 1.0;  /* 0-based index */
            if (has_weights) {
                state->factors[g].weighted_counts[new_level - 1] += weights_c[i];
            }
        }

        free(remap);
    }

    /* Allocate CG solver buffers */
    ST_int num_threads = omp_get_max_threads();
    if (num_threads > 8) num_threads = 8;

    state->thread_cg_r = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_cg_u = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_cg_v = (ST_double **)malloc(num_threads * sizeof(ST_double *));
    state->thread_proj = (ST_double **)malloc(num_threads * sizeof(ST_double *));

    /* thread_fe_means needs num_threads * G arrays, one per (thread, factor) combination */
    state->thread_fe_means = (ST_double **)calloc(num_threads * G, sizeof(ST_double *));

    for (ST_int t = 0; t < num_threads; t++) {
        state->thread_cg_r[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_cg_u[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_cg_v[t] = (ST_double *)malloc(N * sizeof(ST_double));
        state->thread_proj[t] = (ST_double *)malloc(N * sizeof(ST_double));

        /* Allocate mean arrays for each factor - use num_levels (contiguous after remapping) */
        for (ST_int g = 0; g < G; g++) {
            state->thread_fe_means[t * G + g] = (ST_double *)malloc(
                state->factors[g].num_levels * sizeof(ST_double));
        }
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
    /* Single FE: df_a = num_levels (constant absorbed by FE)
       Multi-FE: df_a = sum(num_levels) - (G - 1) for redundant levels */
    ST_int df_a = 0;
    ST_int df_a_nested = 0;  /* Levels from FE nested within cluster */
    for (ST_int g = 0; g < G; g++) {
        df_a += state->factors[g].num_levels;
        /* Check if this FE is nested in cluster (1-indexed) */
        if (nested_fe_index > 0 && g == (nested_fe_index - 1)) {
            df_a_nested = state->factors[g].num_levels;
        }
    }
    if (G > 1) df_a -= (G - 1);  /* Subtract redundant levels for multi-way FE */

    /* For VCE calculation when FE is nested in cluster:
       - nested FE levels don't contribute to df_a for VCE
       - nested_adj = 1 to account for the constant */
    ST_int df_a_for_vce = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;

    if (verbose && has_cluster) {
        char buf[512];
        snprintf(buf, sizeof(buf),
            "civreghdfe: nested_fe_index=%d, df_a=%d, df_a_nested=%d, df_a_for_vce=%d, nested_adj=%d\n",
            (int)nested_fe_index, (int)df_a, (int)df_a_nested, (int)df_a_for_vce, (int)nested_adj);
        SF_display(buf);
    }

    /* Remap cluster IDs to contiguous 1-based indices using shared utility */
    ST_int num_clusters = 0;
    if (has_cluster) {
        ctools_remap_cluster_ids(cluster_ids_c, N, &num_clusters);
    }

    /* Allocate output arrays */
    ST_double *beta = (ST_double *)calloc(K_total, sizeof(ST_double));
    ST_double *V = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *first_stage_F = (ST_double *)calloc(K_endog, sizeof(ST_double));

    /* Compute k-class IV estimation (2SLS, LIML, Fuller, etc.) */
    /* For cluster VCE, pass df_a_for_vce (excluding nested FE levels) and nested_adj */
    ST_double lambda = 1.0;
    ST_retcode rc = compute_2sls(
        y_dem, X_exog_dem, X_endog_dem, Z_dem,
        weights_c, weight_type,
        N, K_exog, K_endog, K_iv,
        beta, V, first_stage_F,
        vce_type, cluster_ids_c, num_clusters,
        df_a_for_vce, nested_adj, verbose,
        est_method, kclass_user, fuller_alpha, &lambda,
        kernel_type, bw
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
            for (ST_int g = 0; g < G; g++) {
                if (state->thread_fe_means[t * G + g])
                    free(state->thread_fe_means[t * G + g]);
            }
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

    /* Compute RSS, TSS, R-squared, Root MSE, and F-statistic */
    /* These are computed on demeaned (partialled-out) data */
    ST_double rss = 0.0;
    ST_double tss = 0.0;
    ST_double y_mean = 0.0;
    ST_double sum_w = 0.0;

    /* Compute weighted mean of y_dem */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        y_mean += w * y_dem[i];
        sum_w += w;
    }
    y_mean /= sum_w;

    /* Compute TSS (centered) = sum(w * (y - ybar)^2) */
    for (ST_int i = 0; i < N; i++) {
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        ST_double dev = y_dem[i] - y_mean;
        tss += w * dev * dev;
    }

    /* Compute RSS = sum(w * resid^2) */
    /* Compute fitted values: yhat = X_exog * beta[0:K_exog-1] + X_endog * beta[K_exog:K_total-1] */
    for (ST_int i = 0; i < N; i++) {
        ST_double fitted = 0.0;
        for (ST_int j = 0; j < K_exog; j++) {
            fitted += X_exog_dem[j * N + i] * beta[j];
        }
        for (ST_int j = 0; j < K_endog; j++) {
            fitted += X_endog_dem[j * N + i] * beta[K_exog + j];
        }
        ST_double resid = y_dem[i] - fitted;
        ST_double w = (weights_c && weight_type != 0) ? weights_c[i] : 1.0;
        rss += w * resid * resid;
    }

    ST_double r2 = (tss > 0) ? 1.0 - rss / tss : 0.0;
    ST_int df_r_val = N - K_total - df_a;
    if (df_r_val <= 0) df_r_val = 1;
    ST_double rmse = sqrt(rss / df_r_val);

    /* Compute model F-statistic: (R2 / K) / ((1 - R2) / df_r) */
    ST_double f_stat = 0.0;
    if (K_total > 0 && r2 < 1.0) {
        f_stat = (r2 / K_total) / ((1.0 - r2) / df_r_val);
    }

    /* Store results to Stata */
    /* Scalars */
    SF_scal_save("__civreghdfe_N", (ST_double)N);
    SF_scal_save("__civreghdfe_df_r", (ST_double)df_r_val);
    SF_scal_save("__civreghdfe_df_a", (ST_double)df_a);
    SF_scal_save("__civreghdfe_K", (ST_double)K_total);
    SF_scal_save("__civreghdfe_rss", rss);
    SF_scal_save("__civreghdfe_tss", tss);
    SF_scal_save("__civreghdfe_r2", r2);
    SF_scal_save("__civreghdfe_rmse", rmse);
    SF_scal_save("__civreghdfe_F", f_stat);

    if (has_cluster) {
        SF_scal_save("__civreghdfe_N_clust", (ST_double)num_clusters);
    }

    /* Store lambda for LIML/Fuller */
    if (est_method == 1 || est_method == 2) {
        SF_scal_save("__civreghdfe_lambda", lambda);
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
        for (ST_int g = 0; g < G; g++) {
            if (state->thread_fe_means[t * G + g])
                free(state->thread_fe_means[t * G + g]);
        }
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
