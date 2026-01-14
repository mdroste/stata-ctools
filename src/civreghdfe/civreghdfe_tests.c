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

/* Forward declarations from creghdfe_solver */
extern int cholesky(double *A, int K);
extern int invert_from_cholesky(const double *L, int K, double *A_inv);

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

        /* F-stat: (R^2 / q) / ((1 - R^2) / (N - K_iv - df_a)) */
        ST_int q = K_iv - K_exog;  /* Number of excluded instruments */
        if (q <= 0) q = 1;
        ST_int denom_df = N - K_iv - df_a;
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
        ST_double F = first_stage_F[0];
        ST_int df_resid = N - K_iv - df_a;
        if (df_resid <= 0) df_resid = 1;

        /* partial_R² = (L * F) / (L * F + df_resid) */
        ST_double partial_r2 = ((ST_double)L * F) / ((ST_double)L * F + (ST_double)df_resid);

        /* Anderson LM = N * partial_R² */
        *underid_stat = (ST_double)N * partial_r2;

        /* Cragg-Donald F = first-stage F for single endogenous */
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
                /* Cluster-robust */
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
                ST_int Nminus_wald = N - df_a - 1;
                if (Nminus_wald <= 0) Nminus_wald = 1;

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

                            /* KP rk Wald F = (Nminus_wald / N) * kp_wald_raw / L */
                            ST_double kp_wald = kp_wald_raw * (ST_double)Nminus_wald / (ST_double)N;
                            *kp_f = kp_wald / (ST_double)L;

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
