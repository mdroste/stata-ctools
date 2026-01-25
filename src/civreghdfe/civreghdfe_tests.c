/*
    civreghdfe_tests.c
    IV Diagnostic Tests Implementation

    Implements:
    - First-stage F statistics
    - Underidentification tests (Anderson LM, Kleibergen-Paap rk LM)
    - Weak instrument tests (Cragg-Donald F, Kleibergen-Paap rk Wald F)
    - Overidentification tests (Sargan, Hansen J)
    - Endogeneity tests (Durbin-Wu-Hausman)
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "civreghdfe_tests.h"
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"

/* Shared OLS functions */
#include "../ctools_ols.h"
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/*
    Compute first-stage F statistics for each endogenous variable.
*/
void civreghdfe_compute_first_stage_F(
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *ZtZ_inv,
    const ST_double *ZtX,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_double *first_stage_F
)
{
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
        ST_double r2 = (xx > 0) ? xpx / xx : 0.0;
        if (r2 > 1.0) r2 = 1.0;
        if (r2 < 0.0) r2 = 0.0;

        /* F-stat: (R^2 / q) / ((1 - R^2) / (N - K_iv - df_a - 1))
           The -1 is for the constant partialled out by HDFE */
        ST_int q = K_iv - K_exog;  /* Number of excluded instruments */
        if (q <= 0) q = 1;
        ST_int denom_df = N - K_iv - df_a - 1;
        if (denom_df <= 0) denom_df = 1;

        first_stage_F[e] = (r2 / (ST_double)q) / ((1.0 - r2) / (ST_double)denom_df);
        if (first_stage_F[e] < 0) first_stage_F[e] = 0;
    }
}

/*
    Compute underidentification test (Anderson LM / Kleibergen-Paap rk LM).
*/
void civreghdfe_compute_underid_test(
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *ZtZ,
    const ST_double *ZtZ_inv,
    const ST_double *temp1,
    const ST_double *first_stage_F,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *underid_stat,
    ST_int *underid_df,
    ST_double *cd_f,
    ST_double *kp_f
)
{
    ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
    ST_int i, k;

    /* ZtZ_inv reserved for potential robust formula implementation */
    (void)ZtZ_inv;

    *underid_stat = 0.0;
    *underid_df = L;
    *cd_f = 0.0;
    *kp_f = 0.0;

    if (K_endog <= 0 || L <= 0) return;

    if (K_endog == 1 && first_stage_F) {
        /* Single endogenous variable: simple formulas */
        /* first_stage_F is already the partial F-stat from excluded instruments
           after controlling for exogenous regressors (computed in civreghdfe_estimate.c) */
        ST_double F = first_stage_F[0];
        ST_int df_resid = N - K_iv - df_a;
        if (df_resid <= 0) df_resid = 1;

        /* Recover partial_R² from F: partial_R² = (L * F) / (L * F + df_resid)
           This is the partial R² of excluded instruments after controlling for exogenous */
        ST_double partial_r2 = ((ST_double)L * F) / ((ST_double)L * F + (ST_double)df_resid);

        /* Anderson LM = N * partial_R² (based on canonical correlations)
           For single endogenous, this equals N * partial_R² */
        *underid_stat = (ST_double)N * partial_r2;

        /* Cragg-Donald F = partial first-stage F for single endogenous */
        *cd_f = F;

        /* Kleibergen-Paap for robust/cluster VCE */
        if (vce_type != 0) {
            /* Extract first-stage coefficients: pi = temp1[:, K_exog] */
            const ST_double *pi = temp1 + K_exog * K_iv;

            /* Compute first-stage residuals: v = X_endog - Z * pi */
            ST_double *v = (ST_double *)calloc(N, sizeof(ST_double));
            if (!v) return;

            for (i = 0; i < N; i++) {
                ST_double pred = 0.0;
                for (k = 0; k < K_iv; k++) {
                    pred += Z[k * N + i] * pi[k];
                }
                v[i] = X_endog[i] - pred;
            }

            /* Compute shat0 for Wald statistic */
            ST_double *shat0 = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
            ST_double *shat0_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
            if (!shat0 || !shat0_inv) {
                free(v); free(shat0); free(shat0_inv);
                return;
            }

            if (vce_type == 2 && cluster_ids && num_clusters > 0) {
                /* Cluster-robust with small-sample correction G/(G-1) */
                ST_double *cluster_Zv = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
                if (!cluster_Zv) {
                    free(v); free(shat0); free(shat0_inv);
                    return;
                }

                for (i = 0; i < N; i++) {
                    ST_int c = cluster_ids[i] - 1;
                    if (c < 0 || c >= num_clusters) continue;
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    for (k = 0; k < K_iv; k++) {
                        cluster_Zv[c * K_iv + k] += w * Z[k * N + i] * v[i];
                    }
                }

                for (ST_int c = 0; c < num_clusters; c++) {
                    for (ST_int ki = 0; ki < K_iv; ki++) {
                        for (ST_int kj = 0; kj < K_iv; kj++) {
                            shat0[kj * K_iv + ki] +=
                                cluster_Zv[c * K_iv + ki] * cluster_Zv[c * K_iv + kj];
                        }
                    }
                }
                free(cluster_Zv);

                /* Apply small-sample correction: G/(G-1) */
                if (num_clusters > 1) {
                    ST_double cluster_adj = (ST_double)num_clusters / (ST_double)(num_clusters - 1);
                    for (ST_int ki = 0; ki < K_iv * K_iv; ki++) {
                        shat0[ki] *= cluster_adj;
                    }
                }
            } else {
                /* HC robust */
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    ST_double v2w = w * v[i] * v[i];
                    for (ST_int ki = 0; ki < K_iv; ki++) {
                        for (ST_int kj = 0; kj < K_iv; kj++) {
                            shat0[kj * K_iv + ki] += v2w * Z[ki * N + i] * Z[kj * N + i];
                        }
                    }
                }
            }

            /* Invert shat0 */
            memcpy(shat0_inv, shat0, K_iv * K_iv * sizeof(ST_double));
            ST_int shat_ok = (cholesky(shat0_inv, K_iv) == 0);
            if (shat_ok) {
                shat_ok = (invert_from_cholesky(shat0_inv, K_iv, shat0_inv) == 0);
            }

            if (shat_ok) {
                /* DOF adjustment matches ivreghdfe: N - K_iv - df_a - 1
                   The -1 is for the constant partialled out by HDFE */
                ST_int df_adjust = N - K_iv - df_a - 1;
                if (df_adjust <= 0) df_adjust = 1;

                /* Compute ZtX_e = Z'X_endog (moment condition) */
                ST_double *ZtX_e = (ST_double *)calloc(K_iv, sizeof(ST_double));
                if (ZtX_e) {
                    for (ST_int ki = 0; ki < K_iv; ki++) {
                        ST_double sum = 0.0;
                        for (ST_int kj = 0; kj < K_iv; kj++) {
                            sum += ZtZ[kj * K_iv + ki] * pi[kj];
                        }
                        ZtX_e[ki] = sum;
                    }

                    /* Compute shat0_inv * ZtX_e */
                    ST_double *shat0_inv_ZtX = (ST_double *)calloc(K_iv, sizeof(ST_double));
                    if (shat0_inv_ZtX) {
                        for (ST_int ki = 0; ki < K_iv; ki++) {
                            ST_double sum = 0.0;
                            for (ST_int kj = 0; kj < K_iv; kj++) {
                                sum += shat0_inv[kj * K_iv + ki] * ZtX_e[kj];
                            }
                            shat0_inv_ZtX[ki] = sum;
                        }

                        /* Compute ZtZ * shat0_inv * ZtX */
                        ST_double *ZtZ_shat0_inv_ZtX = (ST_double *)calloc(K_iv, sizeof(ST_double));
                        if (ZtZ_shat0_inv_ZtX) {
                            for (ST_int ki = 0; ki < K_iv; ki++) {
                                ST_double sum = 0.0;
                                for (ST_int kj = 0; kj < K_iv; kj++) {
                                    sum += ZtZ[kj * K_iv + ki] * shat0_inv_ZtX[kj];
                                }
                                ZtZ_shat0_inv_ZtX[ki] = sum;
                            }

                            /* KP Wald = pi' * ZtZ * shat0_inv * ZtZ * pi */
                            ST_double kp_wald_raw = 0.0;
                            for (ST_int ki = 0; ki < K_iv; ki++) {
                                kp_wald_raw += pi[ki] * ZtZ_shat0_inv_ZtX[ki];
                            }

                            /* KP rk Wald F formula from ivreghdfe:
                               - Non-cluster: chi2 / N * (N - K_iv - sdofminus) / L
                               - Cluster: chi2 / (N-1) * (N - K_iv - sdofminus) * (G-1)/G / L
                               Note: shat0 has G/(G-1) for cluster, so shat0_inv has (G-1)/G.
                               df_adjust = N - K_iv - df_a - 1 accounts for partialled-out constant.
                               Non-cluster robust also needs N/(N-K_iv) adjustment to match ranktest. */
                            if (vce_type == 2 && num_clusters > 1) {
                                /* Cluster formula - (G-1)/G already in kp_wald_raw via shat0 */
                                ST_double kp_wald = kp_wald_raw / (ST_double)(N - 1);
                                kp_wald *= (ST_double)df_adjust;
                                *kp_f = kp_wald / (ST_double)L;
                            } else {
                                /* Non-cluster robust formula - includes N/(N-1) adjustment */
                                ST_double kp_wald = kp_wald_raw * (ST_double)df_adjust / (ST_double)N;
                                /* Apply small-sample adjustment N/(N-1) to match ranktest */
                                kp_wald *= (ST_double)N / (ST_double)(N - 1);
                                *kp_f = kp_wald / (ST_double)L;
                            }

                            free(ZtZ_shat0_inv_ZtX);
                        }
                        free(shat0_inv_ZtX);
                    }
                    free(ZtX_e);
                }

                /* Compute KP LM using X_endog directly */
                ST_double *shat0_lm = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
                if (shat0_lm) {
                    if (vce_type == 2 && cluster_ids && num_clusters > 0) {
                        ST_double *cluster_Zy = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
                        if (cluster_Zy) {
                            for (i = 0; i < N; i++) {
                                ST_int c = cluster_ids[i] - 1;
                                if (c < 0 || c >= num_clusters) continue;
                                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                for (k = 0; k < K_iv; k++) {
                                    cluster_Zy[c * K_iv + k] += w * Z[k * N + i] * X_endog[i];
                                }
                            }
                            for (ST_int c = 0; c < num_clusters; c++) {
                                for (ST_int ki = 0; ki < K_iv; ki++) {
                                    for (ST_int kj = 0; kj < K_iv; kj++) {
                                        shat0_lm[kj * K_iv + ki] +=
                                            cluster_Zy[c * K_iv + ki] * cluster_Zy[c * K_iv + kj];
                                    }
                                }
                            }
                            free(cluster_Zy);
                        }
                    } else {
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            ST_double y2w = w * X_endog[i] * X_endog[i];
                            for (ST_int ki = 0; ki < K_iv; ki++) {
                                for (ST_int kj = 0; kj < K_iv; kj++) {
                                    shat0_lm[kj * K_iv + ki] += y2w * Z[ki * N + i] * Z[kj * N + i];
                                }
                            }
                        }
                    }

                    ST_double *shat0_lm_inv = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
                    if (shat0_lm_inv) {
                        memcpy(shat0_lm_inv, shat0_lm, K_iv * K_iv * sizeof(ST_double));
                        if (cholesky(shat0_lm_inv, K_iv) == 0 &&
                            invert_from_cholesky(shat0_lm_inv, K_iv, shat0_lm_inv) == 0) {

                            /* Compute Z'X_endog */
                            ST_double *ZtY = (ST_double *)calloc(K_iv, sizeof(ST_double));
                            if (ZtY) {
                                for (k = 0; k < K_iv; k++) {
                                    ST_double sum = 0.0;
                                    for (i = 0; i < N; i++) {
                                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                        sum += w * Z[k * N + i] * X_endog[i];
                                    }
                                    ZtY[k] = sum;
                                }

                                ST_double *shat0_lm_inv_ZtY = (ST_double *)calloc(K_iv, sizeof(ST_double));
                                if (shat0_lm_inv_ZtY) {
                                    for (ST_int ki = 0; ki < K_iv; ki++) {
                                        ST_double sum = 0.0;
                                        for (ST_int kj = 0; kj < K_iv; kj++) {
                                            sum += shat0_lm_inv[kj * K_iv + ki] * ZtY[kj];
                                        }
                                        shat0_lm_inv_ZtY[ki] = sum;
                                    }

                                    ST_double quad_lm = 0.0;
                                    for (ST_int ki = 0; ki < K_iv; ki++) {
                                        quad_lm += ZtY[ki] * shat0_lm_inv_ZtY[ki];
                                    }

                                    *underid_stat = quad_lm;
                                    free(shat0_lm_inv_ZtY);
                                }
                                free(ZtY);
                            }
                        }
                        free(shat0_lm_inv);
                    }
                    free(shat0_lm);
                }
            }

            free(v);
            free(shat0);
            free(shat0_inv);
        }

    } else if (K_endog > 1) {
        /* Multiple endogenous variables: canonical correlations */
        /* Compute Y'Y where Y = X_endog */
        ST_double *YtY = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
        if (!YtY) return;

        for (ST_int e1 = 0; e1 < K_endog; e1++) {
            for (ST_int e2 = 0; e2 < K_endog; e2++) {
                ST_double sum = 0.0;
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    sum += w * X_endog[e1 * N + i] * X_endog[e2 * N + i];
                }
                YtY[e2 * K_endog + e1] = sum;
            }
        }

        /* Extract excluded instruments submatrix */
        ST_double *ZeZe = (ST_double *)calloc(L * L, sizeof(ST_double));
        ST_double *ZeY = (ST_double *)calloc(L * K_endog, sizeof(ST_double));

        if (ZeZe && ZeY) {
            /* Extract Z_excl'Z_excl from ZtZ */
            for (ST_int l1 = 0; l1 < L; l1++) {
                for (ST_int l2 = 0; l2 < L; l2++) {
                    ZeZe[l2 * L + l1] = ZtZ[(K_exog + l2) * K_iv + (K_exog + l1)];
                }
            }

            /* Compute Z_excl'Y */
            for (ST_int l = 0; l < L; l++) {
                for (ST_int e = 0; e < K_endog; e++) {
                    ST_double sum = 0.0;
                    for (i = 0; i < N; i++) {
                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                        sum += w * Z[(K_exog + l) * N + i] * X_endog[e * N + i];
                    }
                    ZeY[e * L + l] = sum;
                }
            }

            /* Invert matrices and compute canonical correlations */
            ST_double *YtY_inv = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
            ST_double *ZeZe_inv = (ST_double *)calloc(L * L, sizeof(ST_double));

            if (YtY_inv && ZeZe_inv) {
                memcpy(YtY_inv, YtY, K_endog * K_endog * sizeof(ST_double));
                memcpy(ZeZe_inv, ZeZe, L * L * sizeof(ST_double));

                ST_int yy_ok = (cholesky(YtY_inv, K_endog) == 0);
                if (yy_ok) yy_ok = (invert_from_cholesky(YtY_inv, K_endog, YtY_inv) == 0);

                ST_int zz_ok = (cholesky(ZeZe_inv, L) == 0);
                if (zz_ok) zz_ok = (invert_from_cholesky(ZeZe_inv, L, ZeZe_inv) == 0);

                if (yy_ok && zz_ok) {
                    /* Compute (Y'Y)^{-1/2} */
                    ST_double *YtY_inv_half = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
                    if (YtY_inv_half) {
                        civreghdfe_matpowersym_neg_half(YtY, K_endog, YtY_inv_half);

                        /* Compute temp1_mat = inv(ZeZe) * ZeY */
                        ST_double *temp1_mat = (ST_double *)calloc(L * K_endog, sizeof(ST_double));
                        if (temp1_mat) {
                            for (ST_int e = 0; e < K_endog; e++) {
                                for (ST_int l1 = 0; l1 < L; l1++) {
                                    ST_double sum = 0.0;
                                    for (ST_int l2 = 0; l2 < L; l2++) {
                                        sum += ZeZe_inv[l2 * L + l1] * ZeY[e * L + l2];
                                    }
                                    temp1_mat[e * L + l1] = sum;
                                }
                            }

                            /* temp2 = ZeY' * temp1 */
                            ST_double *temp2_mat = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
                            if (temp2_mat) {
                                for (ST_int e1 = 0; e1 < K_endog; e1++) {
                                    for (ST_int e2 = 0; e2 < K_endog; e2++) {
                                        ST_double sum = 0.0;
                                        for (ST_int l = 0; l < L; l++) {
                                            sum += ZeY[e1 * L + l] * temp1_mat[e2 * L + l];
                                        }
                                        temp2_mat[e2 * K_endog + e1] = sum;
                                    }
                                }

                                /* M = YtY_inv_half * temp2 * YtY_inv_half */
                                ST_double *temp3_mat = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
                                ST_double *M = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));

                                if (temp3_mat && M) {
                                    civreghdfe_matmul_ab(temp2_mat, YtY_inv_half, K_endog, K_endog, K_endog, temp3_mat);
                                    civreghdfe_matmul_ab(YtY_inv_half, temp3_mat, K_endog, K_endog, K_endog, M);

                                    /* Find minimum eigenvalue */
                                    ST_double *eigenvalues = (ST_double *)calloc(K_endog, sizeof(ST_double));
                                    if (eigenvalues) {
                                        civreghdfe_jacobi_eigenvalues(M, K_endog, eigenvalues);

                                        ST_double min_eval = eigenvalues[0];
                                        for (ST_int e = 1; e < K_endog; e++) {
                                            if (eigenvalues[e] < min_eval) min_eval = eigenvalues[e];
                                        }

                                        /* Anderson LM = N * min_eigenvalue */
                                        *underid_stat = (ST_double)N * min_eval;
                                        *underid_df = L - K_endog + 1;
                                        if (*underid_df < 1) *underid_df = 1;

                                        /* Cragg-Donald F */
                                        ST_int df_cd = N - df_a - K_iv;
                                        if (df_cd <= 0) df_cd = 1;
                                        if (min_eval > 0.0 && min_eval < 1.0) {
                                            ST_double cd = min_eval / (1.0 - min_eval);
                                            *cd_f = cd * ((ST_double)df_cd / (ST_double)L);
                                        }

                                        if (vce_type != 0) {
                                            *kp_f = *cd_f;  /* Placeholder */
                                        }

                                        free(eigenvalues);
                                    }
                                }
                                free(temp3_mat);
                                free(M);
                                free(temp2_mat);
                            }
                            free(temp1_mat);
                        }
                        free(YtY_inv_half);
                    }
                }
            }
            free(YtY_inv);
            free(ZeZe_inv);
        }

        free(YtY);
        free(ZeZe);
        free(ZeY);
    }
}

/*
    Compute Sargan/Hansen J overidentification test.
*/
void civreghdfe_compute_sargan_j(
    const ST_double *resid,
    const ST_double *Z,
    const ST_double *ZtZ_inv,
    ST_double rss,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double *sargan_stat,
    ST_int *overid_df
)
{
    ST_int K_total = K_exog + K_endog;
    ST_int i, k;

    *overid_df = K_iv - K_total;
    *sargan_stat = 0.0;

    if (*overid_df <= 0) return;

    /* Compute Z'r */
    ST_double *Ztr = (ST_double *)calloc(K_iv, sizeof(ST_double));
    if (!Ztr) return;

    for (k = 0; k < K_iv; k++) {
        ST_double sum = 0.0;
        const ST_double *z_col = Z + k * N;
        if (weights && weight_type != 0) {
            for (i = 0; i < N; i++) {
                sum += z_col[i] * weights[i] * resid[i];
            }
        } else {
            for (i = 0; i < N; i++) {
                sum += z_col[i] * resid[i];
            }
        }
        Ztr[k] = sum;
    }

    if (vce_type == 0) {
        /* Sargan statistic: N * r'Z(Z'Z)^-1 Z'r / rss */
        ST_double *ZtZ_inv_Ztr = (ST_double *)calloc(K_iv, sizeof(ST_double));
        if (ZtZ_inv_Ztr) {
            for (i = 0; i < K_iv; i++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_iv; k++) {
                    sum += ZtZ_inv[k * K_iv + i] * Ztr[k];
                }
                ZtZ_inv_Ztr[i] = sum;
            }

            ST_double quad_form = 0.0;
            for (i = 0; i < K_iv; i++) {
                quad_form += Ztr[i] * ZtZ_inv_Ztr[i];
            }

            *sargan_stat = (ST_double)N * quad_form / rss;
            free(ZtZ_inv_Ztr);
        }
    } else {
        /* Hansen J: r'Z(Z'ΩZ)^-1 Z'r */
        ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        if (!ZOmegaZ) {
            free(Ztr);
            return;
        }

        if (vce_type == 2 && cluster_ids && num_clusters > 0) {
            ivvce_compute_ZOmegaZ_cluster(Z, resid, weights, weight_type,
                                          cluster_ids, N, K_iv, num_clusters, ZOmegaZ);
        } else {
            ivvce_compute_ZOmegaZ_robust(Z, resid, weights, weight_type, N, K_iv, ZOmegaZ);
        }

        ST_double *ZOmegaZ_inv = (ST_double *)malloc(K_iv * K_iv * sizeof(ST_double));
        if (ZOmegaZ_inv) {
            memcpy(ZOmegaZ_inv, ZOmegaZ, K_iv * K_iv * sizeof(ST_double));

            if (cholesky(ZOmegaZ_inv, K_iv) == 0) {
                invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv);

                ST_double *temp_zr = (ST_double *)calloc(K_iv, sizeof(ST_double));
                if (temp_zr) {
                    for (i = 0; i < K_iv; i++) {
                        ST_double sum = 0.0;
                        for (k = 0; k < K_iv; k++) {
                            sum += ZOmegaZ_inv[k * K_iv + i] * Ztr[k];
                        }
                        temp_zr[i] = sum;
                    }

                    ST_double quad_form = 0.0;
                    for (i = 0; i < K_iv; i++) {
                        quad_form += Ztr[i] * temp_zr[i];
                    }

                    *sargan_stat = quad_form;
                    free(temp_zr);
                }
            }
            free(ZOmegaZ_inv);
        }
        free(ZOmegaZ);
    }

    free(Ztr);
}

/*
    Compute Durbin-Wu-Hausman endogeneity test.
*/
void civreghdfe_compute_dwh_test(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *temp1,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_double *endog_chi2,
    ST_double *endog_f,
    ST_int *endog_df
)
{
    ST_int K_total = K_exog + K_endog;
    ST_int i, j, k;

    *endog_chi2 = 0.0;
    *endog_f = 0.0;
    *endog_df = K_endog;

    if (K_endog <= 0) return;

    /* Compute first-stage residuals: v = X_endog - Z * pi */
    ST_double *v_resid = (ST_double *)calloc(N * K_endog, sizeof(ST_double));
    if (!v_resid) return;

    for (ST_int e = 0; e < K_endog; e++) {
        const ST_double *x_endog_col = X_endog + e * N;
        ST_double *v_col = v_resid + e * N;
        const ST_double *pi_col = temp1 + (K_exog + e) * K_iv;

        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_iv; k++) {
                pred += Z[k * N + i] * pi_col[k];
            }
            v_col[i] = x_endog_col[i] - pred;
        }
    }

    /* Build augmented design matrix: [X_exog, X_endog, v] */
    ST_int K_aug = K_total + K_endog;
    ST_double *X_aug = (ST_double *)calloc(N * K_aug, sizeof(ST_double));
    ST_double *XaXa = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *Xay = (ST_double *)calloc(K_aug, sizeof(ST_double));

    if (!X_aug || !XaXa || !Xay) {
        free(v_resid);
        free(X_aug); free(XaXa); free(Xay);
        return;
    }

    /* Copy X_exog, X_endog, v_resid */
    if (X_exog) {
        for (j = 0; j < K_exog; j++) {
            for (i = 0; i < N; i++) {
                X_aug[j * N + i] = X_exog[j * N + i];
            }
        }
    }
    if (X_endog) {
        for (j = 0; j < K_endog; j++) {
            for (i = 0; i < N; i++) {
                X_aug[(K_exog + j) * N + i] = X_endog[j * N + i];
            }
        }
    }
    for (j = 0; j < K_endog; j++) {
        for (i = 0; i < N; i++) {
            X_aug[(K_total + j) * N + i] = v_resid[j * N + i];
        }
    }

    /* Compute X_aug'X_aug */
    civreghdfe_matmul_atb(X_aug, X_aug, N, K_aug, K_aug, XaXa);

    /* Compute X_aug'y */
    for (j = 0; j < K_aug; j++) {
        ST_double sum = 0.0;
        const ST_double *xa_col = X_aug + j * N;
        for (i = 0; i < N; i++) {
            sum += xa_col[i] * y[i];
        }
        Xay[j] = sum;
    }

    /* Solve for augmented OLS coefficients */
    ST_double *XaXa_L = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *XaXa_inv = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *beta_aug = (ST_double *)calloc(K_aug, sizeof(ST_double));

    if (XaXa_L && XaXa_inv && beta_aug) {
        memcpy(XaXa_L, XaXa, K_aug * K_aug * sizeof(ST_double));
        if (cholesky(XaXa_L, K_aug) == 0) {
            invert_from_cholesky(XaXa_L, K_aug, XaXa_inv);

            /* beta_aug = XaXa_inv * Xay */
            for (i = 0; i < K_aug; i++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_aug; k++) {
                    sum += XaXa_inv[k * K_aug + i] * Xay[k];
                }
                beta_aug[i] = sum;
            }

            /* Compute residuals and sigma^2 */
            ST_double sse_aug = 0.0;
            for (i = 0; i < N; i++) {
                ST_double fitted = 0.0;
                for (j = 0; j < K_aug; j++) {
                    fitted += X_aug[j * N + i] * beta_aug[j];
                }
                ST_double r = y[i] - fitted;
                sse_aug += r * r;
            }

            ST_int df_aug = N - K_aug - df_a;
            if (df_aug <= 0) df_aug = 1;
            ST_double sigma2_aug = sse_aug / df_aug;

            /* gamma = beta_aug[K_total:K_aug-1] */
            ST_double *gamma = beta_aug + K_total;

            /* Extract v'v block from XaXa */
            ST_double *XaXa_vv = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
            if (XaXa_vv) {
                for (j = 0; j < K_endog; j++) {
                    for (i = 0; i < K_endog; i++) {
                        XaXa_vv[j * K_endog + i] = XaXa[(K_total + j) * K_aug + (K_total + i)];
                    }
                }

                /* chi2 = gamma' * XaXa_vv * gamma / sigma2_aug */
                ST_double quad = 0.0;
                for (j = 0; j < K_endog; j++) {
                    ST_double sum = 0.0;
                    for (i = 0; i < K_endog; i++) {
                        sum += XaXa_vv[j * K_endog + i] * gamma[i];
                    }
                    quad += gamma[j] * sum;
                }

                *endog_chi2 = quad / sigma2_aug;
                *endog_f = (*endog_chi2) / K_endog;

                free(XaXa_vv);
            }
        }
    }

    free(v_resid);
    free(X_aug);
    free(XaXa);
    free(Xay);
    free(XaXa_L);
    free(XaXa_inv);
    free(beta_aug);
}

/*
    Compute C-statistic for testing instrument orthogonality.

    C = J_full - J_restricted
    where J_restricted is computed using only the non-tested instruments.

    For efficiency, we use the algebraic form:
    C = e'Z_test (Z_test' M_{Z_rest} Z_test)^{-1} Z_test'e / sigma^2
    where M_{Z_rest} = I - Z_rest(Z_rest'Z_rest)^{-1}Z_rest'
*/
void civreghdfe_compute_cstat(
    const ST_double *resid,
    const ST_double *Z,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const ST_int *orthog_indices,
    ST_int n_orthog,
    ST_int vce_type,
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double sargan_full,
    ST_double *cstat,
    ST_int *cstat_df
)
{
    ST_int i, j, k;
    ST_int K_total = K_exog + K_endog;
    ST_int K_rest = K_iv - n_orthog;  /* Instruments not being tested */

    /* Suppress unused parameter warnings - these are reserved for future
       heteroskedasticity/cluster-robust C-stat implementation */
    (void)vce_type;
    (void)cluster_ids;
    (void)num_clusters;

    *cstat = 0.0;
    *cstat_df = n_orthog;

    if (n_orthog <= 0 || K_rest < K_total) return;

    /* Create mask for tested instruments (0 = keep, 1 = test) */
    ST_int *test_mask = (ST_int *)calloc(K_iv, sizeof(ST_int));
    if (!test_mask) return;

    for (i = 0; i < n_orthog; i++) {
        ST_int idx = orthog_indices[i] - 1;  /* Convert 1-based to 0-based */
        /* Excluded instruments start at position K_exog in Z */
        if (idx >= 0 && idx < K_iv - K_exog) {
            test_mask[K_exog + idx] = 1;
        }
    }

    /* Build Z_rest (instruments not being tested) - cast to size_t to prevent 32-bit overflow */
    ST_double *Z_rest = (ST_double *)malloc((size_t)N * K_rest * sizeof(ST_double));
    ST_double *Z_test = (ST_double *)malloc((size_t)N * n_orthog * sizeof(ST_double));
    if (!Z_rest || !Z_test) {
        free(test_mask);
        if (Z_rest) free(Z_rest);
        if (Z_test) free(Z_test);
        return;
    }

    ST_int rest_idx = 0, test_idx = 0;
    for (k = 0; k < K_iv; k++) {
        if (test_mask[k] == 0) {
            for (i = 0; i < N; i++) {
                Z_rest[rest_idx * N + i] = Z[k * N + i];
            }
            rest_idx++;
        } else {
            for (i = 0; i < N; i++) {
                Z_test[test_idx * N + i] = Z[k * N + i];
            }
            test_idx++;
        }
    }

    /* Compute Z_rest'Z_rest and its inverse */
    ST_double *ZrZr = (ST_double *)calloc(K_rest * K_rest, sizeof(ST_double));
    ST_double *ZrZr_inv = (ST_double *)calloc(K_rest * K_rest, sizeof(ST_double));
    if (!ZrZr || !ZrZr_inv) goto cleanup;

    civreghdfe_matmul_atb(Z_rest, Z_rest, N, K_rest, K_rest, ZrZr);

    memcpy(ZrZr_inv, ZrZr, K_rest * K_rest * sizeof(ST_double));
    if (cholesky(ZrZr_inv, K_rest) != 0) goto cleanup;
    if (invert_from_cholesky(ZrZr_inv, K_rest, ZrZr_inv) != 0) goto cleanup;

    /* Compute Z_rest'e */
    ST_double *Zr_e = (ST_double *)calloc(K_rest, sizeof(ST_double));
    if (!Zr_e) goto cleanup;

    for (k = 0; k < K_rest; k++) {
        ST_double sum = 0.0;
        for (i = 0; i < N; i++) {
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            sum += w * Z_rest[k * N + i] * resid[i];
        }
        Zr_e[k] = sum;
    }

    /* Compute restricted Sargan: e'Z_rest (Z_rest'Z_rest)^{-1} Z_rest'e */
    ST_double *temp_r = (ST_double *)calloc(K_rest, sizeof(ST_double));
    if (!temp_r) { free(Zr_e); goto cleanup; }

    for (k = 0; k < K_rest; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_rest; j++) {
            sum += ZrZr_inv[j * K_rest + k] * Zr_e[j];
        }
        temp_r[k] = sum;
    }

    ST_double quad_rest = 0.0;
    for (k = 0; k < K_rest; k++) {
        quad_rest += Zr_e[k] * temp_r[k];
    }

    /* Compute RSS for scaling */
    ST_double rss = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        rss += w * resid[i] * resid[i];
    }

    /* For homoskedastic case: restricted Sargan = N * quad / rss */
    ST_double sargan_rest = (ST_double)N * quad_rest / rss;

    /* C = J_full - J_restricted */
    *cstat = sargan_full - sargan_rest;
    if (*cstat < 0) *cstat = 0;  /* Can happen due to numerical issues */

    free(Zr_e);
    free(temp_r);

cleanup:
    free(test_mask);
    free(Z_rest);
    free(Z_test);
    if (ZrZr) free(ZrZr);
    if (ZrZr_inv) free(ZrZr_inv);
}

/*
    Compute endogeneity test for a subset of endogenous regressors.
    Uses DWH-style test: only tests specified variables, not all endogenous.
*/
void civreghdfe_compute_endogtest_subset(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *temp1,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    const ST_int *endogtest_indices,
    ST_int n_endogtest,
    ST_double *endogtest_stat,
    ST_int *endogtest_df
)
{
    ST_int K_total = K_exog + K_endog;
    ST_int i, j, k;

    *endogtest_stat = 0.0;
    *endogtest_df = n_endogtest;

    if (n_endogtest <= 0 || n_endogtest > K_endog) return;

    /* Compute first-stage residuals only for tested endogenous variables */
    ST_double *v_resid = (ST_double *)calloc(N * n_endogtest, sizeof(ST_double));
    if (!v_resid) return;

    for (ST_int t = 0; t < n_endogtest; t++) {
        ST_int e = endogtest_indices[t] - 1;  /* Convert 1-based to 0-based */
        if (e < 0 || e >= K_endog) {
            free(v_resid);
            return;
        }

        const ST_double *x_endog_col = X_endog + e * N;
        ST_double *v_col = v_resid + t * N;
        const ST_double *pi_col = temp1 + (K_exog + e) * K_iv;

        for (i = 0; i < N; i++) {
            ST_double pred = 0.0;
            for (k = 0; k < K_iv; k++) {
                pred += Z[k * N + i] * pi_col[k];
            }
            v_col[i] = x_endog_col[i] - pred;
        }
    }

    /* Build augmented design matrix: [X_exog, X_endog, v_test] */
    ST_int K_aug = K_total + n_endogtest;
    ST_double *X_aug = (ST_double *)calloc(N * K_aug, sizeof(ST_double));
    ST_double *XaXa = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *Xay = (ST_double *)calloc(K_aug, sizeof(ST_double));

    if (!X_aug || !XaXa || !Xay) {
        free(v_resid);
        free(X_aug); free(XaXa); free(Xay);
        return;
    }

    /* Copy X_exog, X_endog, v_resid */
    if (X_exog) {
        for (j = 0; j < K_exog; j++) {
            for (i = 0; i < N; i++) {
                X_aug[j * N + i] = X_exog[j * N + i];
            }
        }
    }
    for (j = 0; j < K_endog; j++) {
        for (i = 0; i < N; i++) {
            X_aug[(K_exog + j) * N + i] = X_endog[j * N + i];
        }
    }
    for (j = 0; j < n_endogtest; j++) {
        for (i = 0; i < N; i++) {
            X_aug[(K_total + j) * N + i] = v_resid[j * N + i];
        }
    }

    /* Compute X_aug'X_aug */
    civreghdfe_matmul_atb(X_aug, X_aug, N, K_aug, K_aug, XaXa);

    /* Compute X_aug'y */
    for (j = 0; j < K_aug; j++) {
        ST_double sum = 0.0;
        const ST_double *xa_col = X_aug + j * N;
        for (i = 0; i < N; i++) {
            sum += xa_col[i] * y[i];
        }
        Xay[j] = sum;
    }

    /* Solve for augmented OLS coefficients */
    ST_double *XaXa_L = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *XaXa_inv = (ST_double *)calloc(K_aug * K_aug, sizeof(ST_double));
    ST_double *beta_aug = (ST_double *)calloc(K_aug, sizeof(ST_double));

    if (XaXa_L && XaXa_inv && beta_aug) {
        memcpy(XaXa_L, XaXa, K_aug * K_aug * sizeof(ST_double));
        if (cholesky(XaXa_L, K_aug) == 0) {
            invert_from_cholesky(XaXa_L, K_aug, XaXa_inv);

            /* beta_aug = XaXa_inv * Xay */
            for (i = 0; i < K_aug; i++) {
                ST_double sum = 0.0;
                for (k = 0; k < K_aug; k++) {
                    sum += XaXa_inv[k * K_aug + i] * Xay[k];
                }
                beta_aug[i] = sum;
            }

            /* Compute residuals and sigma^2 */
            ST_double sse_aug = 0.0;
            for (i = 0; i < N; i++) {
                ST_double fitted = 0.0;
                for (j = 0; j < K_aug; j++) {
                    fitted += X_aug[j * N + i] * beta_aug[j];
                }
                ST_double r = y[i] - fitted;
                sse_aug += r * r;
            }

            ST_int df_aug = N - K_aug - df_a;
            if (df_aug <= 0) df_aug = 1;
            ST_double sigma2_aug = sse_aug / df_aug;

            /* gamma = beta_aug[K_total:K_aug-1] */
            ST_double *gamma = beta_aug + K_total;

            /* Extract v'v block from XaXa */
            ST_double *XaXa_vv = (ST_double *)calloc(n_endogtest * n_endogtest, sizeof(ST_double));
            if (XaXa_vv) {
                for (j = 0; j < n_endogtest; j++) {
                    for (i = 0; i < n_endogtest; i++) {
                        XaXa_vv[j * n_endogtest + i] = XaXa[(K_total + j) * K_aug + (K_total + i)];
                    }
                }

                /* chi2 = gamma' * XaXa_vv * gamma / sigma2_aug */
                ST_double quad = 0.0;
                for (j = 0; j < n_endogtest; j++) {
                    ST_double sum = 0.0;
                    for (i = 0; i < n_endogtest; i++) {
                        sum += XaXa_vv[j * n_endogtest + i] * gamma[i];
                    }
                    quad += gamma[j] * sum;
                }

                *endogtest_stat = quad / sigma2_aug;
                free(XaXa_vv);
            }
        }
    }

    free(v_resid);
    free(X_aug);
    free(XaXa);
    free(Xay);
    free(XaXa_L);
    free(XaXa_inv);
    free(beta_aug);
}

/*
    Compute instrument redundancy test.

    Tests if specified instruments add information to the first stage.
    LM statistic based on partial R^2 of tested instruments.
*/
void civreghdfe_compute_redundant(
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const ST_int *redund_indices,
    ST_int n_redund,
    ST_double *redund_stat,
    ST_int *redund_df
)
{
    ST_int i, j, k;
    ST_int K_rest = K_iv - n_redund;

    *redund_stat = 0.0;
    *redund_df = K_endog * n_redund;

    if (n_redund <= 0 || K_rest < K_exog + K_endog) return;

    /* Create mask for tested instruments */
    ST_int *test_mask = (ST_int *)calloc(K_iv, sizeof(ST_int));
    if (!test_mask) return;

    for (i = 0; i < n_redund; i++) {
        ST_int idx = redund_indices[i] - 1;  /* Convert 1-based to 0-based */
        if (idx >= 0 && idx < K_iv - K_exog) {
            test_mask[K_exog + idx] = 1;
        }
    }

    /* Build Z_rest and Z_test - cast to size_t to prevent 32-bit overflow */
    ST_double *Z_rest = (ST_double *)malloc((size_t)N * K_rest * sizeof(ST_double));
    ST_double *Z_test = (ST_double *)malloc((size_t)N * n_redund * sizeof(ST_double));
    if (!Z_rest || !Z_test) {
        free(test_mask);
        if (Z_rest) free(Z_rest);
        if (Z_test) free(Z_test);
        return;
    }

    ST_int rest_idx = 0, test_idx = 0;
    for (k = 0; k < K_iv; k++) {
        if (test_mask[k] == 0) {
            for (i = 0; i < N; i++) {
                Z_rest[rest_idx * N + i] = Z[k * N + i];
            }
            rest_idx++;
        } else {
            for (i = 0; i < N; i++) {
                Z_test[test_idx * N + i] = Z[k * N + i];
            }
            test_idx++;
        }
    }

    /* Compute partial correlation: X_endog on Z_test controlling for Z_rest */
    /* Residualize Z_test on Z_rest */
    ST_double *ZrZr = (ST_double *)calloc(K_rest * K_rest, sizeof(ST_double));
    ST_double *ZrZr_inv = (ST_double *)calloc(K_rest * K_rest, sizeof(ST_double));
    if (!ZrZr || !ZrZr_inv) {
        free(test_mask); free(Z_rest); free(Z_test);
        if (ZrZr) free(ZrZr);
        if (ZrZr_inv) free(ZrZr_inv);
        return;
    }

    civreghdfe_matmul_atb(Z_rest, Z_rest, N, K_rest, K_rest, ZrZr);
    memcpy(ZrZr_inv, ZrZr, K_rest * K_rest * sizeof(ST_double));

    if (cholesky(ZrZr_inv, K_rest) != 0 || invert_from_cholesky(ZrZr_inv, K_rest, ZrZr_inv) != 0) {
        free(test_mask); free(Z_rest); free(Z_test);
        free(ZrZr); free(ZrZr_inv);
        return;
    }

    /* Compute Z_rest'X_endog and Z_rest'Z_test */
    ST_double *ZrX = (ST_double *)calloc(K_rest * K_endog, sizeof(ST_double));
    ST_double *ZrZt = (ST_double *)calloc(K_rest * n_redund, sizeof(ST_double));

    if (ZrX && ZrZt) {
        civreghdfe_matmul_atb(Z_rest, X_endog, N, K_rest, K_endog, ZrX);
        civreghdfe_matmul_atb(Z_rest, Z_test, N, K_rest, n_redund, ZrZt);

        /* Residualize X_endog on Z_rest: M_r * X_endog */
        ST_double *X_resid = (ST_double *)calloc(N * K_endog, sizeof(ST_double));
        ST_double *Zt_resid = (ST_double *)calloc(N * n_redund, sizeof(ST_double));

        if (X_resid && Zt_resid) {
            /* X_resid = X_endog - Z_rest * (ZrZr_inv * ZrX) */
            ST_double *coef_x = (ST_double *)calloc(K_rest * K_endog, sizeof(ST_double));
            ST_double *coef_z = (ST_double *)calloc(K_rest * n_redund, sizeof(ST_double));

            if (coef_x && coef_z) {
                /* coef_x = ZrZr_inv * ZrX */
                for (j = 0; j < K_endog; j++) {
                    for (i = 0; i < K_rest; i++) {
                        ST_double sum = 0.0;
                        for (k = 0; k < K_rest; k++) {
                            sum += ZrZr_inv[k * K_rest + i] * ZrX[j * K_rest + k];
                        }
                        coef_x[j * K_rest + i] = sum;
                    }
                }

                /* coef_z = ZrZr_inv * ZrZt */
                for (j = 0; j < n_redund; j++) {
                    for (i = 0; i < K_rest; i++) {
                        ST_double sum = 0.0;
                        for (k = 0; k < K_rest; k++) {
                            sum += ZrZr_inv[k * K_rest + i] * ZrZt[j * K_rest + k];
                        }
                        coef_z[j * K_rest + i] = sum;
                    }
                }

                /* Compute residuals */
                for (j = 0; j < K_endog; j++) {
                    for (i = 0; i < N; i++) {
                        ST_double pred = 0.0;
                        for (k = 0; k < K_rest; k++) {
                            pred += Z_rest[k * N + i] * coef_x[j * K_rest + k];
                        }
                        X_resid[j * N + i] = X_endog[j * N + i] - pred;
                    }
                }

                for (j = 0; j < n_redund; j++) {
                    for (i = 0; i < N; i++) {
                        ST_double pred = 0.0;
                        for (k = 0; k < K_rest; k++) {
                            pred += Z_rest[k * N + i] * coef_z[j * K_rest + k];
                        }
                        Zt_resid[j * N + i] = Z_test[j * N + i] - pred;
                    }
                }

                /* Compute Zt_resid'X_resid and its squared norm relative to X'X */
                ST_double *ZtX = (ST_double *)calloc(n_redund * K_endog, sizeof(ST_double));
                ST_double *XtX = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
                ST_double *ZtZt = (ST_double *)calloc(n_redund * n_redund, sizeof(ST_double));

                if (ZtX && XtX && ZtZt) {
                    civreghdfe_matmul_atb(Zt_resid, X_resid, N, n_redund, K_endog, ZtX);
                    civreghdfe_matmul_atb(X_resid, X_resid, N, K_endog, K_endog, XtX);
                    civreghdfe_matmul_atb(Zt_resid, Zt_resid, N, n_redund, n_redund, ZtZt);

                    /* LM = N * tr((ZtX' * inv(ZtZt) * ZtX) * inv(XtX)) */
                    ST_double *ZtZt_inv = (ST_double *)calloc(n_redund * n_redund, sizeof(ST_double));
                    ST_double *XtX_inv = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));

                    if (ZtZt_inv && XtX_inv) {
                        memcpy(ZtZt_inv, ZtZt, n_redund * n_redund * sizeof(ST_double));
                        memcpy(XtX_inv, XtX, K_endog * K_endog * sizeof(ST_double));

                        if (cholesky(ZtZt_inv, n_redund) == 0 &&
                            invert_from_cholesky(ZtZt_inv, n_redund, ZtZt_inv) == 0 &&
                            cholesky(XtX_inv, K_endog) == 0 &&
                            invert_from_cholesky(XtX_inv, K_endog, XtX_inv) == 0) {

                            /* temp1 = ZtZt_inv * ZtX */
                            ST_double *temp1_lm = (ST_double *)calloc(n_redund * K_endog, sizeof(ST_double));
                            /* temp2 = ZtX' * temp1 = ZtX' * ZtZt_inv * ZtX */
                            ST_double *temp2_lm = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));

                            if (temp1_lm && temp2_lm) {
                                for (j = 0; j < K_endog; j++) {
                                    for (i = 0; i < n_redund; i++) {
                                        ST_double sum = 0.0;
                                        for (k = 0; k < n_redund; k++) {
                                            sum += ZtZt_inv[k * n_redund + i] * ZtX[j * n_redund + k];
                                        }
                                        temp1_lm[j * n_redund + i] = sum;
                                    }
                                }

                                for (ST_int j1 = 0; j1 < K_endog; j1++) {
                                    for (ST_int j2 = 0; j2 < K_endog; j2++) {
                                        ST_double sum = 0.0;
                                        for (i = 0; i < n_redund; i++) {
                                            sum += ZtX[j1 * n_redund + i] * temp1_lm[j2 * n_redund + i];
                                        }
                                        temp2_lm[j2 * K_endog + j1] = sum;
                                    }
                                }

                                /* Compute tr(temp2 * XtX_inv) */
                                ST_double trace = 0.0;
                                for (i = 0; i < K_endog; i++) {
                                    for (k = 0; k < K_endog; k++) {
                                        trace += temp2_lm[k * K_endog + i] * XtX_inv[i * K_endog + k];
                                    }
                                }

                                *redund_stat = (ST_double)N * trace;
                            }

                            free(temp1_lm);
                            free(temp2_lm);
                        }
                    }

                    free(ZtZt_inv);
                    free(XtX_inv);
                }

                free(ZtX);
                free(XtX);
                free(ZtZt);
            }

            free(coef_x);
            free(coef_z);
        }

        free(X_resid);
        free(Zt_resid);
    }

    free(ZrX);
    free(ZrZt);
    free(test_mask);
    free(Z_rest);
    free(Z_test);
    free(ZrZr);
    free(ZrZr_inv);
}
