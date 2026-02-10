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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "civreghdfe_tests.h"
#include "civreghdfe_matrix.h"
#include "civreghdfe_vce.h"
#include "../ctools_config.h"

/* Shared OLS functions */
#include "../ctools_ols.h"
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/*
    Helper: compute cluster meat for shat0 (one clustering dimension).
    shat[L x L] = sum_g (sum_i_in_g w_i * Z_i * v_i) * (sum_i_in_g w_i * Z_i * v_i)'
    Output is ADDED to shat (caller must zero-initialize).
*/
static void compute_cluster_shat0_dim(
    const ST_double *Z_excl,   /* L columns, N rows, column-major */
    const ST_double *v,        /* N x 1 residuals */
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int N,
    ST_int L,
    ST_double *shat
)
{
    ST_double *cluster_Zv = (ST_double *)calloc(num_clusters * L, sizeof(ST_double));
    if (!cluster_Zv) return;

    #pragma omp parallel if(N > 5000)
    {
        ST_double *local_cZv = (ST_double *)calloc(num_clusters * L, sizeof(ST_double));
        if (local_cZv) {
            #pragma omp for schedule(static)
            for (ST_int i = 0; i < N; i++) {
                ST_int c = cluster_ids[i] - 1;
                if (c < 0 || c >= num_clusters) continue;
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                for (ST_int k = 0; k < L; k++) {
                    local_cZv[c * L + k] += w * Z_excl[k * N + i] * v[i];
                }
            }
            #pragma omp critical
            {
                for (ST_int idx = 0; idx < num_clusters * L; idx++)
                    cluster_Zv[idx] += local_cZv[idx];
            }
            free(local_cZv);
        }
    }

    for (ST_int c = 0; c < num_clusters; c++) {
        for (ST_int ki = 0; ki < L; ki++) {
            for (ST_int kj = 0; kj < L; kj++) {
                shat[kj * L + ki] += cluster_Zv[c * L + ki] * cluster_Zv[c * L + kj];
            }
        }
    }
    free(cluster_Zv);
}

/*
    Helper: compute intersection cluster IDs for two-way clustering.
    Returns number of intersection clusters, or 0 on failure.
    Caller must free *intersection_ids if return > 0.
*/
static ST_int compute_intersection_ids(
    const ST_int *cluster1_ids,
    ST_int num_clusters1,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int N,
    ST_int **intersection_ids_out
)
{
    ST_int *intersection_ids = (ST_int *)calloc(N, sizeof(ST_int));
    if (!intersection_ids) return 0;

    size_t pair_array_size;
    if (ctools_safe_alloc_size((size_t)(num_clusters1 + 1), (size_t)(num_clusters2 + 1),
                               sizeof(ST_int), &pair_array_size) != 0) {
        free(intersection_ids);
        return 0;
    }
    ST_int *pair_to_int = (ST_int *)calloc(1, pair_array_size);
    if (!pair_to_int) {
        free(intersection_ids);
        return 0;
    }

    ST_int num_intersection = 0;
    for (ST_int i = 0; i < N; i++) {
        ST_int c1 = cluster1_ids[i];
        ST_int c2 = cluster2_ids[i];
        size_t pair_idx = (size_t)c1 * (size_t)(num_clusters2 + 1) + (size_t)c2;
        if (pair_to_int[pair_idx] == 0) {
            num_intersection++;
            pair_to_int[pair_idx] = num_intersection;
        }
        intersection_ids[i] = pair_to_int[pair_idx];
    }
    free(pair_to_int);

    *intersection_ids_out = intersection_ids;
    return num_intersection;
}

/*
    Helper: compute two-way CGM shat0.
    shat0[L x L] = shat1 + shat2 - shat_int (raw, no DOF correction)
*/
static void compute_twoway_shat0(
    const ST_double *Z_excl,
    const ST_double *v,
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster1_ids,
    ST_int num_clusters1,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int N,
    ST_int L,
    ST_double *shat0
)
{
    ST_double *shat1 = (ST_double *)calloc(L * L, sizeof(ST_double));
    ST_double *shat2 = (ST_double *)calloc(L * L, sizeof(ST_double));
    ST_double *shat_int = (ST_double *)calloc(L * L, sizeof(ST_double));
    if (!shat1 || !shat2 || !shat_int) {
        free(shat1); free(shat2); free(shat_int);
        return;
    }

    /* Compute intersection cluster IDs */
    ST_int *int_ids = NULL;
    ST_int num_int = compute_intersection_ids(cluster1_ids, num_clusters1,
                                               cluster2_ids, num_clusters2,
                                               N, &int_ids);
    if (num_int == 0 || !int_ids) {
        free(shat1); free(shat2); free(shat_int);
        return;
    }

    compute_cluster_shat0_dim(Z_excl, v, weights, weight_type, cluster1_ids, num_clusters1, N, L, shat1);
    compute_cluster_shat0_dim(Z_excl, v, weights, weight_type, cluster2_ids, num_clusters2, N, L, shat2);
    compute_cluster_shat0_dim(Z_excl, v, weights, weight_type, int_ids, num_int, N, L, shat_int);

    /* CGM: shat0 = shat1 + shat2 - shat_int */
    for (ST_int i = 0; i < L * L; i++) {
        shat0[i] = shat1[i] + shat2[i] - shat_int[i];
    }

    free(shat1);
    free(shat2);
    free(shat_int);
    free(int_ids);
}

/*
    Helper: compute homoskedastic HAC shat0 for Kiefer VCE (test statistics).

    Matches ivreg2/ranktest's m_omega homoskedastic path (livreg2.do lines 194-236):
    shat = sigma2 * Z'WZ + sum_{tau} kw(tau) * sigmahat_tau * (ZZhat_tau + ZZhat_tau')

    where sigma2 = resid'W*resid / (N_eff - dofminus)
          sigmahat_tau = quadcross(e[t], wv, e[t-tau]) / (N-dofminus)
          ZZhat_tau    = quadcross(Z[t], wv, Z[t-tau])

    Uses WITHIN-PANEL lag pairing, matching m_omega's use of Stata's L. operator.
    Data must be sorted by panel then time (tsset).

    shat is NOT divided by N at the end (= shat_accumulated from m_omega,
    before the final /N on line 235).

    Output is WRITTEN to shat (caller must provide zeroed L*L buffer).
*/
static void compute_kiefer_shat0_homoskedastic(
    const ST_double *Z,       /* L columns, N rows, column-major */
    const ST_double *resid,   /* N x 1 residuals */
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *panel_ids,  /* 1-based panel IDs for each obs */
    ST_int num_panels,
    ST_int N,
    ST_int N_eff,
    ST_int L,
    ST_int dofminus,          /* large-sample DOF correction */
    ST_int kernel_type,       /* kernel type for kernel weight */
    ST_int bw,                /* bandwidth (= T for Kiefer) */
    ST_double *shat
)
{
    ST_int i;

    /* Compute sigma2 = resid'W*resid / (N_eff - dofminus)
       Matches m_omega line 197: sigma2 = ee/(N-dofminus) */
    ST_double ee = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        ee += w * resid[i] * resid[i];
    }
    ST_int denom_sigma = N_eff - dofminus;
    if (denom_sigma <= 0) denom_sigma = 1;
    ST_double sigma2 = ee / (ST_double)denom_sigma;

    /* Base term: shat = sigma2 * Z'WZ (line 198) */
    for (ST_int ki = 0; ki < L; ki++) {
        for (ST_int kj = 0; kj <= ki; kj++) {
            ST_double sum = 0.0;
            for (i = 0; i < N; i++) {
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                sum += w * Z[ki * N + i] * Z[kj * N + i];
            }
            shat[ki * L + kj] = sigma2 * sum;
            if (ki != kj) shat[kj * L + ki] = sigma2 * sum;
        }
    }

    /* Build panel start/end indices for within-panel lag pairing.
       Data is sorted by panel then time. panel_ids are 1-based. */
    ST_int *panel_start = (ST_int *)calloc(num_panels, sizeof(ST_int));
    ST_int *panel_end   = (ST_int *)calloc(num_panels, sizeof(ST_int));
    ST_double *ZZhat_tau = (ST_double *)calloc(L * L, sizeof(ST_double));
    if (!panel_start || !panel_end || !ZZhat_tau) {
        free(panel_start); free(panel_end); free(ZZhat_tau);
        return;
    }

    /* Initialize panel boundaries */
    for (i = 0; i < num_panels; i++) {
        panel_start[i] = N;  /* sentinel */
        panel_end[i] = -1;
    }
    for (i = 0; i < N; i++) {
        ST_int p = panel_ids[i] - 1;
        if (p < 0 || p >= num_panels) continue;
        if (i < panel_start[p]) panel_start[p] = i;
        if (i > panel_end[p])   panel_end[p] = i;
    }

    /* Cross-lag terms: within-panel lag pairing.
       Matches m_omega lines 207-232 using Stata's L. operator:
         For each tau, pair obs i with obs i+tau WITHIN the same panel.
         sigmahat = sum_{within-panel pairs} wv * e_i * e_{i+tau} / (N-dofminus)
         ZZhat    = sum_{within-panel pairs} wv * Z_i * Z_{i+tau}'
         ghat = sigmahat # ZZhat  (scalar Kronecker for K1=1)
         shat += kw * (ghat + ghat') */
    for (ST_int tau = 1; tau < bw; tau++) {
        ST_double kw = civreghdfe_kernel_weight(kernel_type, tau, bw);
        if (fabs(kw) < 1e-10) continue;

        ST_double sigmahat_tau = 0.0;
        memset(ZZhat_tau, 0, L * L * sizeof(ST_double));

        /* Iterate over panels, pairing obs i with obs i+tau within each panel */
        for (ST_int p = 0; p < num_panels; p++) {
            ST_int pstart = panel_start[p];
            ST_int pend   = panel_end[p];
            if (pstart >= N || pend < 0) continue;
            ST_int T_p = pend - pstart + 1;
            if (tau >= T_p) continue;  /* No valid pairs in this panel */

            for (i = pstart; i <= pend - tau; i++) {
                ST_int i_lag = i + tau;
                ST_double w_i   = (weights && weight_type != 0) ? weights[i] : 1.0;
                ST_double w_lag = (weights && weight_type != 0) ? weights[i_lag] : 1.0;
                ST_double wv = w_i * w_lag;

                sigmahat_tau += wv * resid[i] * resid[i_lag];

                for (ST_int ki = 0; ki < L; ki++) {
                    for (ST_int kj = 0; kj < L; kj++) {
                        ZZhat_tau[ki * L + kj] +=
                            wv * Z[ki * N + i] * Z[kj * N + i_lag];
                    }
                }
            }
        }

        /* Large-sample DOF correction (m_omega line 226) */
        sigmahat_tau /= (ST_double)denom_sigma;

        /* shat += kw * sigmahat_tau * (ZZhat_tau + ZZhat_tau') */
        for (ST_int ki = 0; ki < L; ki++) {
            for (ST_int kj = 0; kj < L; kj++) {
                shat[kj * L + ki] += kw * sigmahat_tau *
                    (ZZhat_tau[ki * L + kj] + ZZhat_tau[kj * L + ki]);
            }
        }
    }

    free(panel_start);
    free(panel_end);
    free(ZZhat_tau);
}

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
    ST_int N_eff,
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
        ST_int denom_df = N_eff - K_iv - df_a - 1;
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
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int kernel_type,
    ST_int bw,
    ST_int kiefer,
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels,
    ST_double *underid_stat,
    ST_int *underid_df,
    ST_double *cd_f,
    ST_double *kp_f
)
{
    ST_int L = K_iv - K_exog;  /* Number of excluded instruments */
    ST_int i, k;

    /* ivreg2/ranktest always uses Bartlett kernel for test statistics,
       regardless of the kernel specified for VCE */
    const ST_int test_kernel = (kernel_type > 0) ? 1 : 0;

    /* ZtZ, ZtZ_inv reserved for potential robust formula implementation */
    (void)ZtZ;
    (void)ZtZ_inv;
    (void)temp1;

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
        ST_int df_resid = N_eff - K_iv - df_a;
        if (df_resid <= 0) df_resid = 1;

        /* Recover partial_R² from F: partial_R² = (L * F) / (L * F + df_resid)
           This is the partial R² of excluded instruments after controlling for exogenous */
        ST_double partial_r2 = ((ST_double)L * F) / ((ST_double)L * F + (ST_double)df_resid);

        /* Anderson LM = N_eff * partial_R² (based on canonical correlations)
           For single endogenous, this equals N_eff * partial_R² */
        *underid_stat = (ST_double)N_eff * partial_r2;

        /* Cragg-Donald F = partial first-stage F for single endogenous */
        *cd_f = F;

        /* Kleibergen-Paap for robust/cluster VCE */
        if (vce_type != 0) {
            /* For the KP test, we need to partial out included exogenous variables
               from both X_endog and Z_excl before computing the test statistic.
               This is required because the test is on the EXCLUDED instruments' coefficients
               after controlling for exogenous variables. */

            const ST_double *Z_exog = Z;                /* First K_exog columns */
            const ST_double *Z_excl_raw = Z + K_exog * N;  /* Last L columns */

            /* Allocate partialled-out versions */
            ST_double *X_endog_partial = (ST_double *)malloc(N * sizeof(ST_double));
            ST_double *Z_excl_partial = (ST_double *)malloc(N * L * sizeof(ST_double));
            ST_double *v = (ST_double *)calloc(N, sizeof(ST_double));
            if (!X_endog_partial || !Z_excl_partial || !v) {
                free(X_endog_partial); free(Z_excl_partial); free(v);
                return;
            }

            /* Copy raw values */
            memcpy(X_endog_partial, X_endog, N * sizeof(ST_double));
            memcpy(Z_excl_partial, Z_excl_raw, N * L * sizeof(ST_double));

            /* Partial out exogenous variables if K_exog > 0 */
            if (K_exog > 0) {
                /* Compute (Z_exog' * Z_exog)^-1 */
                ST_double *ZeZe = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                ST_double *ZeZe_inv = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                if (ZeZe && ZeZe_inv) {
                    for (ST_int l1 = 0; l1 < K_exog; l1++) {
                        for (ST_int l2 = 0; l2 < K_exog; l2++) {
                            ST_double sum = 0.0;
                            for (i = 0; i < N; i++) {
                                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                sum += w * Z_exog[l1 * N + i] * Z_exog[l2 * N + i];
                            }
                            ZeZe[l2 * K_exog + l1] = sum;
                        }
                    }

                    memcpy(ZeZe_inv, ZeZe, K_exog * K_exog * sizeof(ST_double));
                    if (cholesky(ZeZe_inv, K_exog) == 0) {
                        invert_from_cholesky(ZeZe_inv, K_exog, ZeZe_inv);

                        /* Partial out Z_exog from X_endog */
                        ST_double *ZeX = (ST_double *)calloc(K_exog, sizeof(ST_double));
                        if (ZeX) {
                            for (ST_int l = 0; l < K_exog; l++) {
                                ST_double sum = 0.0;
                                for (i = 0; i < N; i++) {
                                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                    sum += w * Z_exog[l * N + i] * X_endog[i];
                                }
                                ZeX[l] = sum;
                            }
                            /* gamma = (Z_exog'Z_exog)^-1 * Z_exog'X_endog */
                            ST_double *gamma = (ST_double *)calloc(K_exog, sizeof(ST_double));
                            if (gamma) {
                                for (ST_int l1 = 0; l1 < K_exog; l1++) {
                                    ST_double sum = 0.0;
                                    for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                        sum += ZeZe_inv[l2 * K_exog + l1] * ZeX[l2];
                                    }
                                    gamma[l1] = sum;
                                }
                                /* X_endog_partial = X_endog - Z_exog * gamma */
                                for (i = 0; i < N; i++) {
                                    ST_double pred = 0.0;
                                    for (ST_int l = 0; l < K_exog; l++) {
                                        pred += Z_exog[l * N + i] * gamma[l];
                                    }
                                    X_endog_partial[i] -= pred;
                                }
                                free(gamma);
                            }
                            free(ZeX);
                        }

                        /* Partial out Z_exog from each column of Z_excl */
                        for (ST_int col = 0; col < L; col++) {
                            ST_double *ZeZcol = (ST_double *)calloc(K_exog, sizeof(ST_double));
                            if (ZeZcol) {
                                for (ST_int l = 0; l < K_exog; l++) {
                                    ST_double sum = 0.0;
                                    for (i = 0; i < N; i++) {
                                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                        sum += w * Z_exog[l * N + i] * Z_excl_raw[col * N + i];
                                    }
                                    ZeZcol[l] = sum;
                                }
                                ST_double *gamma_z = (ST_double *)calloc(K_exog, sizeof(ST_double));
                                if (gamma_z) {
                                    for (ST_int l1 = 0; l1 < K_exog; l1++) {
                                        ST_double sum = 0.0;
                                        for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                            sum += ZeZe_inv[l2 * K_exog + l1] * ZeZcol[l2];
                                        }
                                        gamma_z[l1] = sum;
                                    }
                                    for (i = 0; i < N; i++) {
                                        ST_double pred = 0.0;
                                        for (ST_int l = 0; l < K_exog; l++) {
                                            pred += Z_exog[l * N + i] * gamma_z[l];
                                        }
                                        Z_excl_partial[col * N + i] -= pred;
                                    }
                                    free(gamma_z);
                                }
                                free(ZeZcol);
                            }
                        }
                    }
                }
                free(ZeZe);
                free(ZeZe_inv);
            }

            /* Now compute first-stage on partialled data: X_endog_partial = Z_excl_partial * pi_excl + v */
            /* First compute (Z_excl_partial' * Z_excl_partial)^-1 * Z_excl_partial' * X_endog_partial */
            ST_double *ZtZ_excl = (ST_double *)calloc(L * L, sizeof(ST_double));
            ST_double *ZtX = (ST_double *)calloc(L, sizeof(ST_double));
            ST_double *pi_excl = (ST_double *)calloc(L, sizeof(ST_double));

            if (ZtZ_excl && ZtX && pi_excl) {
                for (ST_int l1 = 0; l1 < L; l1++) {
                    for (ST_int l2 = 0; l2 < L; l2++) {
                        ST_double sum = 0.0;
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            sum += w * Z_excl_partial[l1 * N + i] * Z_excl_partial[l2 * N + i];
                        }
                        ZtZ_excl[l2 * L + l1] = sum;
                    }
                    ST_double sum = 0.0;
                    for (i = 0; i < N; i++) {
                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                        sum += w * Z_excl_partial[l1 * N + i] * X_endog_partial[i];
                    }
                    ZtX[l1] = sum;
                }

                /* Solve for pi_excl = ZtZ_excl^-1 * ZtX */
                ST_double *ZtZ_excl_inv = (ST_double *)calloc(L * L, sizeof(ST_double));
                if (ZtZ_excl_inv) {
                    memcpy(ZtZ_excl_inv, ZtZ_excl, L * L * sizeof(ST_double));
                    if (cholesky(ZtZ_excl_inv, L) == 0 && invert_from_cholesky(ZtZ_excl_inv, L, ZtZ_excl_inv) == 0) {
                        for (ST_int l1 = 0; l1 < L; l1++) {
                            ST_double sum = 0.0;
                            for (ST_int l2 = 0; l2 < L; l2++) {
                                sum += ZtZ_excl_inv[l2 * L + l1] * ZtX[l2];
                            }
                            pi_excl[l1] = sum;
                        }
                    }
                    free(ZtZ_excl_inv);
                }

                /* Compute residuals: v = X_endog_partial - Z_excl_partial * pi_excl */
                for (i = 0; i < N; i++) {
                    ST_double pred = 0.0;
                    for (k = 0; k < L; k++) {
                        pred += Z_excl_partial[k * N + i] * pi_excl[k];
                    }
                    v[i] = X_endog_partial[i] - pred;
                }
            }

            /* Use partialled Z_excl for shat0 computation */
            const ST_double *Z_excl = Z_excl_partial;

            ST_double *shat0 = (ST_double *)calloc(L * L, sizeof(ST_double));
            ST_double *shat0_inv = (ST_double *)calloc(L * L, sizeof(ST_double));
            if (!shat0 || !shat0_inv) {
                free(X_endog_partial); free(Z_excl_partial);
                free(v); free(shat0); free(shat0_inv);
                if (ZtZ_excl) free(ZtZ_excl);
                if (ZtX) free(ZtX);
                if (pi_excl) free(pi_excl);
                return;
            }

            if (kiefer) {
                /* Kiefer VCE: ivreg2 calls ranktest WITHOUT bw/kernel for
                   test statistics (bwopt/kernopt are empty because bw is set
                   after ivparse returns). So KP tests use IID shat0, not HAC.
                   IID shat0 = sigma2 * Z'Z where sigma2 = v'v/N */
                ST_double ee_kf = 0.0;
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    ee_kf += w * v[i] * v[i];
                }
                ST_double sigma2_kf = ee_kf / (ST_double)N_eff;
                for (ST_int ki = 0; ki < L; ki++) {
                    for (ST_int kj = 0; kj <= ki; kj++) {
                        ST_double sum = 0.0;
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            sum += w * Z_excl[ki * N + i] * Z_excl[kj * N + i];
                        }
                        shat0[ki * L + kj] = sigma2_kf * sum;
                        if (ki != kj) shat0[kj * L + ki] = sigma2_kf * sum;
                    }
                }
            } else if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && cluster_ids && cluster2_ids && num_clusters > 0 && num_clusters2 > 0) {
                /* Two-way cluster-robust: CGM shat0 = shat1 + shat2 - shat_int */
                compute_twoway_shat0(Z_excl, v, weights, weight_type,
                                     cluster_ids, num_clusters,
                                     cluster2_ids, num_clusters2,
                                     N, L, shat0);
            } else if ((vce_type == CIVREGHDFE_VCE_CLUSTER || vce_type == CIVREGHDFE_VCE_DKRAAY) && cluster_ids && num_clusters > 0) {
                /* Cluster-robust or HAC (Driscoll-Kraay) */
                ST_double *cluster_Zv = (ST_double *)calloc(num_clusters * L, sizeof(ST_double));
                if (!cluster_Zv) {
                    free(X_endog_partial); free(Z_excl_partial);
                    free(v); free(shat0); free(shat0_inv);
                    if (ZtZ_excl) free(ZtZ_excl);
                    if (ZtX) free(ZtX);
                    if (pi_excl) free(pi_excl);
                    return;
                }

                /* Sum Z_excl * v within each cluster (time period for dkraay) */
                #pragma omp parallel if(N > 5000)
                {
                    ST_double *local_cZv = (ST_double *)calloc(num_clusters * L, sizeof(ST_double));
                    if (local_cZv) {
                        #pragma omp for schedule(static)
                        for (i = 0; i < N; i++) {
                            ST_int c = cluster_ids[i] - 1;
                            if (c < 0 || c >= num_clusters) continue;
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            for (ST_int kk = 0; kk < L; kk++) {
                                local_cZv[c * L + kk] += w * Z_excl[kk * N + i] * v[i];
                            }
                        }
                        #pragma omp critical
                        {
                            for (ST_int idx = 0; idx < num_clusters * L; idx++)
                                cluster_Zv[idx] += local_cZv[idx];
                        }
                        free(local_cZv);
                    }
                }

                /* Check if HAC kernel should be applied (dkraay with kernel) */
                if (vce_type == CIVREGHDFE_VCE_DKRAAY && kernel_type > 0 && bw > 0) {
                    /* HAC: shat0 = sum_t sum_s kernel_weight(|t-s|) * cluster_Zv_t' * cluster_Zv_s
                       For dkraay, clusters are time periods (1, 2, 3, ..., T) */
                    for (ST_int t = 0; t < num_clusters; t++) {
                        for (ST_int s = 0; s < num_clusters; s++) {
                            ST_int lag = (t > s) ? (t - s) : (s - t);
                            ST_double kw = civreghdfe_kernel_weight(test_kernel, lag, bw);
                            if (fabs(kw) < 1e-10) continue;  /* Skip if no weight */

                            for (ST_int ki = 0; ki < L; ki++) {
                                for (ST_int kj = 0; kj < L; kj++) {
                                    shat0[kj * L + ki] += kw *
                                        cluster_Zv[t * L + ki] * cluster_Zv[s * L + kj];
                                }
                            }
                        }
                    }
                } else {
                    /* Standard cluster-robust: sum of outer products */
                    for (ST_int c = 0; c < num_clusters; c++) {
                        for (ST_int ki = 0; ki < L; ki++) {
                            for (ST_int kj = 0; kj < L; kj++) {
                                shat0[kj * L + ki] +=
                                    cluster_Zv[c * L + ki] * cluster_Zv[c * L + kj];
                            }
                        }
                    }
                }
                free(cluster_Zv);
            } else {
                /* HC robust */
                #pragma omp parallel if(N > 5000)
                {
                    ST_double *local_shat = (ST_double *)calloc(L * L, sizeof(ST_double));
                    if (local_shat) {
                        #pragma omp for schedule(static)
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            ST_double v2w = (weight_type == 1 || weight_type == 3) ?
                                w * w * v[i] * v[i] : w * v[i] * v[i];
                            for (ST_int ki = 0; ki < L; ki++) {
                                for (ST_int kj = 0; kj < L; kj++) {
                                    local_shat[kj * L + ki] += v2w * Z_excl[ki * N + i] * Z_excl[kj * N + i];
                                }
                            }
                        }
                        #pragma omp critical
                        {
                            for (ST_int idx = 0; idx < L * L; idx++)
                                shat0[idx] += local_shat[idx];
                        }
                        free(local_shat);
                    }
                }
            }

            /* Invert shat0 (L x L matrix).
             * Try Cholesky first; fall back to pseudo-inverse for rank-deficient
             * matrices (e.g., when G <= L clusters with cluster constraint). */
            memcpy(shat0_inv, shat0, L * L * sizeof(ST_double));
            ST_int shat_ok = (cholesky(shat0_inv, L) == 0);
            if (shat_ok) {
                /* Check condition: if min diagonal of L is too small relative to
                 * max, the matrix is near-singular and Cholesky inverse unreliable */
                ST_double chol_max = 0.0, chol_min = 1e308;
                for (ST_int ci = 0; ci < L; ci++) {
                    ST_double d = fabs(shat0_inv[ci * L + ci]);
                    if (d > chol_max) chol_max = d;
                    if (d < chol_min) chol_min = d;
                }
                if (chol_max > 0 && chol_min / chol_max < 1e-8) {
                    shat_ok = 0;  /* Near-singular: fall back to pseudo-inverse */
                } else {
                    shat_ok = (invert_from_cholesky(shat0_inv, L, shat0_inv) == 0);
                }
            }
            if (!shat_ok) {
                /* Sweep-based generalized inverse matching Stata's invsym():
                 * pivots on largest diagonal, zeros out rank-deficient columns.
                 * This matches ranktest's use of invsym() for the KP statistic. */
                shat_ok = (ctools_invsym(shat0, L, shat0_inv, NULL) == 0);
            }

            /* Handle standard HAC (kernel+bw without cluster/dkraay) */
            if (shat_ok && vce_type == CIVREGHDFE_VCE_ROBUST && kernel_type > 0 && bw > 0) {
                /* Recompute shat0 with HAC kernel weights */
                memset(shat0, 0, L * L * sizeof(ST_double));

                if (hac_panel_ids && num_hac_panels > 0) {
                    /* Panel-aware HAC: apply kernel weights within each panel */
                    ST_int *panel_counts = (ST_int *)calloc(num_hac_panels, sizeof(ST_int));
                    ST_int *panel_starts = (ST_int *)calloc(num_hac_panels + 1, sizeof(ST_int));
                    ST_int *obs_by_panel = (ST_int *)malloc(N * sizeof(ST_int));

                    if (panel_counts && panel_starts && obs_by_panel) {
                        /* Count observations per panel */
                        for (i = 0; i < N; i++) {
                            ST_int p = hac_panel_ids[i] - 1;
                            if (p >= 0 && p < num_hac_panels) {
                                panel_counts[p]++;
                            }
                        }

                        /* Compute start indices */
                        panel_starts[0] = 0;
                        for (ST_int p = 0; p < num_hac_panels; p++) {
                            panel_starts[p + 1] = panel_starts[p] + panel_counts[p];
                        }

                        /* Reset counts for filling */
                        memset(panel_counts, 0, num_hac_panels * sizeof(ST_int));

                        /* Fill obs_by_panel */
                        for (i = 0; i < N; i++) {
                            ST_int p = hac_panel_ids[i] - 1;
                            if (p >= 0 && p < num_hac_panels) {
                                ST_int idx = panel_starts[p] + panel_counts[p];
                                obs_by_panel[idx] = i;
                                panel_counts[p]++;
                            }
                        }

                        /* Compute HAC within each panel */
                        for (ST_int p = 0; p < num_hac_panels; p++) {
                            ST_int start = panel_starts[p];
                            ST_int end = panel_starts[p + 1];
                            ST_int T_p = end - start;

                            for (ST_int t1 = 0; t1 < T_p; t1++) {
                                ST_int i1 = obs_by_panel[start + t1];
                                for (ST_int t2 = 0; t2 < T_p; t2++) {
                                    ST_int i2 = obs_by_panel[start + t2];
                                    ST_int time_lag = (t1 > t2) ? (t1 - t2) : (t2 - t1);

                                    ST_double kw = civreghdfe_kernel_weight(test_kernel, time_lag, bw);
                                    if (fabs(kw) < 1e-10) continue;

                                    ST_double w1 = (weights && weight_type != 0) ? weights[i1] : 1.0;
                                    ST_double w2 = (weights && weight_type != 0) ? weights[i2] : 1.0;
                                    ST_double v_prod = w1 * v[i1] * w2 * v[i2];

                                    for (ST_int ki = 0; ki < L; ki++) {
                                        for (ST_int kj = 0; kj < L; kj++) {
                                            shat0[kj * L + ki] += kw * v_prod *
                                                Z_excl[ki * N + i1] * Z_excl[kj * N + i2];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (panel_counts) free(panel_counts);
                    if (panel_starts) free(panel_starts);
                    if (obs_by_panel) free(obs_by_panel);
                } else {
                    /* Non-panel HAC: use raw observation indices */
                    for (ST_int t = 0; t < N; t++) {
                        for (ST_int s = 0; s < N; s++) {
                            ST_int lag = (t > s) ? (t - s) : (s - t);
                            ST_double kw = civreghdfe_kernel_weight(test_kernel, lag, bw);
                            if (fabs(kw) < 1e-10) continue;

                            ST_double wt = (weights && weight_type != 0) ? weights[t] : 1.0;
                            ST_double ws = (weights && weight_type != 0) ? weights[s] : 1.0;
                            ST_double v_prod = wt * v[t] * ws * v[s];

                            for (ST_int ki = 0; ki < L; ki++) {
                                for (ST_int kj = 0; kj < L; kj++) {
                                    shat0[kj * L + ki] += kw * v_prod *
                                        Z_excl[ki * N + t] * Z_excl[kj * N + s];
                                }
                            }
                        }
                    }
                }

                /* Re-invert shat0 (with pseudo-inverse fallback) */
                memcpy(shat0_inv, shat0, L * L * sizeof(ST_double));
                shat_ok = (cholesky(shat0_inv, L) == 0);
                if (shat_ok) {
                    ST_double chol_max = 0.0, chol_min = 1e308;
                    for (ST_int ci = 0; ci < L; ci++) {
                        ST_double d = fabs(shat0_inv[ci * L + ci]);
                        if (d > chol_max) chol_max = d;
                        if (d < chol_min) chol_min = d;
                    }
                    if (chol_max > 0 && chol_min / chol_max < 1e-8) {
                        shat_ok = 0;
                    } else {
                        shat_ok = (invert_from_cholesky(shat0_inv, L, shat0_inv) == 0);
                    }
                }
                if (!shat_ok) {
                    shat_ok = (ctools_invsym(shat0, L, shat0_inv, NULL) == 0);
                }
            }

            if (shat_ok) {
                /* DOF adjustment: N_eff - K_iv - df_a - 1
                   The -1 accounts for degrees of freedom adjustment */
                ST_int df_adjust = N_eff - K_iv - df_a - 1;
                if (df_adjust <= 0) df_adjust = 1;

                /* Compute Z_excl'Z_excl (L x L matrix) */
                ST_double *ZtZ_excl = (ST_double *)calloc(L * L, sizeof(ST_double));
                if (ZtZ_excl) {
                    for (ST_int ki = 0; ki < L; ki++) {
                        for (ST_int kj = 0; kj < L; kj++) {
                            ST_double sum = 0.0;
                            for (i = 0; i < N; i++) {
                                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                sum += w * Z_excl[ki * N + i] * Z_excl[kj * N + i];
                            }
                            ZtZ_excl[kj * L + ki] = sum;
                        }
                    }

                    /* Compute ZtZ_excl * pi_excl */
                    ST_double *ZtX_e = (ST_double *)calloc(L, sizeof(ST_double));
                    if (ZtX_e) {
                        for (ST_int ki = 0; ki < L; ki++) {
                            ST_double sum = 0.0;
                            for (ST_int kj = 0; kj < L; kj++) {
                                sum += ZtZ_excl[kj * L + ki] * pi_excl[kj];
                            }
                            ZtX_e[ki] = sum;
                        }

                        /* Compute shat0_inv * ZtX_e */
                        ST_double *shat0_inv_ZtX = (ST_double *)calloc(L, sizeof(ST_double));
                        if (shat0_inv_ZtX) {
                            for (ST_int ki = 0; ki < L; ki++) {
                                ST_double sum = 0.0;
                                for (ST_int kj = 0; kj < L; kj++) {
                                    sum += shat0_inv[kj * L + ki] * ZtX_e[kj];
                                }
                                shat0_inv_ZtX[ki] = sum;
                            }

                            /* Compute ZtZ_excl * shat0_inv * ZtX */
                            ST_double *ZtZ_shat0_inv_ZtX = (ST_double *)calloc(L, sizeof(ST_double));
                            if (ZtZ_shat0_inv_ZtX) {
                                for (ST_int ki = 0; ki < L; ki++) {
                                    ST_double sum = 0.0;
                                    for (ST_int kj = 0; kj < L; kj++) {
                                        sum += ZtZ_excl[kj * L + ki] * shat0_inv_ZtX[kj];
                                    }
                                    ZtZ_shat0_inv_ZtX[ki] = sum;
                                }

                                /* KP Wald = pi_excl' * ZtZ_excl * shat0_inv * ZtZ_excl * pi_excl */
                                ST_double kp_wald_raw = 0.0;
                                for (ST_int ki = 0; ki < L; ki++) {
                                    kp_wald_raw += pi_excl[ki] * ZtZ_shat0_inv_ZtX[ki];
                                }

                            /* KP rk Wald F formula from ivreghdfe/ranktest:

                               In ranktest, the chi2 statistic is computed as:
                               gbar = (1/N) * Z' * X_endog
                               chi2 = gbar' * shat0_inv * gbar * N = (Z'X)' * shat0_inv * Z'X / N

                               In civreghdfe:
                               kp_wald_raw = pi' * ZtZ * shat0_inv * ZtZ * pi
                                           = (Z'X)' * shat0_inv * Z'X  (no 1/N division)
                                           = N * chi2

                               ivreghdfe formula: widstat = chi2/N * dof/L = (kp_wald_raw/N)/N * dof/L

                               So: kp_f = kp_wald_raw / N^2 * dof / L

                               For cluster, shat0 has G/(G-1) correction applied. */

                            /* KP F-stat formula differs between robust and cluster VCE:
                               - Robust: F = chi2/N * (N - K_iv - df_a) / L
                               - Cluster: F = chi2/(N-1) * (N - K_iv - sdofminus) * (G-1)/G / L
                               See ivreghdfe.ado lines 1829 and 1832-1835.

                               For HAC (dkraay), the HAC S matrix computation differs from
                               ranktest's m_omega, requiring an extra (G-1)/G factor. */
                            ST_int dof;
                            ST_double denom;
                            ST_double cluster_adj;
                            if ((vce_type == CIVREGHDFE_VCE_CLUSTER || vce_type == CIVREGHDFE_VCE_CLUSTER2 ||
                                vce_type == CIVREGHDFE_VCE_DKRAAY) && !kiefer) {
                                /* Cluster/dkraay (non-Kiefer): (G-1)/G adjustment
                                   ivreghdfe formula (lines 1832-1835):
                                   rkf = chi2/(N_eff-1) * (N_eff - K - sdofminus - dofminus) * (G-1)/G / L
                                   where sdofminus = partial_ct (absorbed FE count + constant)
                                   For two-way: G = min(G1, G2), matching ivreg2's N_clust convention */
                                dof = N_eff - K_iv - df_a;
                                denom = (ST_double)(N_eff - 1);
                                ST_int G = num_clusters;
                                if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && num_clusters2 > 0) {
                                    G = (num_clusters < num_clusters2) ? num_clusters : num_clusters2;
                                }
                                cluster_adj = (ST_double)(G - 1) / (ST_double)G;
                            } else {
                                /* Robust: ivreghdfe formula is chi2/N_eff * (N_eff - K_iv - sdofminus) / L
                                   where sdofminus = partial_ct + constant_flag
                                   For absorb(), sdofminus = G (# FE variables) + 1 (constant) = df_a
                                   (since df_a = absorbed FE levels + absorbed constant) */
                                dof = N_eff - K_iv - df_a;
                                denom = (ST_double)N_eff;
                                cluster_adj = 1.0;
                            }
                            if (dof <= 0) dof = 1;

                            /* ivreghdfe formula:
                               rkf = r(chi2)/(N-1) * (N - iv1_ct - sdofminus) * (G-1)/G / L
                               r(chi2) from ranktest already has scaling adjustments.
                               Our kp_wald_raw needs similar scaling. */
                            *kp_f = kp_wald_raw / denom * (ST_double)dof * cluster_adj / (ST_double)L;

                                free(ZtZ_shat0_inv_ZtX);
                            }
                            free(shat0_inv_ZtX);
                        }
                        free(ZtX_e);
                    }
                    free(ZtZ_excl);
                }

                /* Compute KP LM using excluded instruments and residualized X_endog
                   The correct formula uses:
                   1. Only excluded instruments (Z_excl, L columns)
                   2. X_endog residualized by exogenous instruments (not raw X_endog)

                   This matches ranktest's partial correlation approach. */

                /* First, compute X_endog residualized by exogenous instruments:
                   y_resid = X_endog - Z_exog * (Z_exog'Z_exog)^-1 * Z_exog'X_endog */
                ST_double *y_resid = (ST_double *)calloc(N, sizeof(ST_double));
                if (y_resid) {
                    if (K_exog > 0) {
                        /* Compute Z_exog'Z_exog and Z_exog'X_endog */
                        ST_double *ZeZe = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
                        ST_double *ZeY = (ST_double *)calloc(K_exog, sizeof(ST_double));

                        if (ZeZe && ZeY) {
                            /* Z_exog is first K_exog columns of Z */
                            for (ST_int l1 = 0; l1 < K_exog; l1++) {
                                for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                    ST_double sum = 0.0;
                                    for (i = 0; i < N; i++) {
                                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                        sum += w * Z[l1 * N + i] * Z[l2 * N + i];
                                    }
                                    ZeZe[l2 * K_exog + l1] = sum;
                                }
                                ST_double sum = 0.0;
                                for (i = 0; i < N; i++) {
                                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                    sum += w * Z[l1 * N + i] * X_endog[i];
                                }
                                ZeY[l1] = sum;
                            }

                            /* Solve for gamma = (Z_exog'Z_exog)^-1 * Z_exog'X_endog */
                            ST_double *gamma = (ST_double *)calloc(K_exog, sizeof(ST_double));
                            if (gamma && cholesky(ZeZe, K_exog) == 0) {
                                /* Forward substitution */
                                for (ST_int l = 0; l < K_exog; l++) {
                                    ST_double sum = ZeY[l];
                                    for (ST_int m = 0; m < l; m++) {
                                        sum -= ZeZe[l * K_exog + m] * gamma[m];
                                    }
                                    gamma[l] = sum / ZeZe[l * K_exog + l];
                                }
                                /* Backward substitution */
                                ST_double *gamma_temp = (ST_double *)calloc(K_exog, sizeof(ST_double));
                                if (gamma_temp) {
                                    memcpy(gamma_temp, gamma, K_exog * sizeof(ST_double));
                                    for (ST_int l = K_exog - 1; l >= 0; l--) {
                                        ST_double sum = gamma_temp[l];
                                        for (ST_int m = l + 1; m < K_exog; m++) {
                                            sum -= ZeZe[m * K_exog + l] * gamma[m];
                                        }
                                        gamma[l] = sum / ZeZe[l * K_exog + l];
                                    }
                                    free(gamma_temp);
                                }

                                /* Compute y_resid = X_endog - Z_exog * gamma */
                                for (i = 0; i < N; i++) {
                                    ST_double pred = 0.0;
                                    for (ST_int l = 0; l < K_exog; l++) {
                                        pred += Z[l * N + i] * gamma[l];
                                    }
                                    y_resid[i] = X_endog[i] - pred;
                                }
                            } else {
                                /* If Cholesky fails, use raw X_endog */
                                memcpy(y_resid, X_endog, N * sizeof(ST_double));
                            }
                            if (gamma) free(gamma);
                        }
                        if (ZeZe) free(ZeZe);
                        if (ZeY) free(ZeY);
                    } else {
                        /* No exogenous regressors, use raw X_endog */
                        memcpy(y_resid, X_endog, N * sizeof(ST_double));
                    }

                    /* Now compute KP LM using Z_excl and y_resid */
                    ST_double *shat0_lm = (ST_double *)calloc(L * L, sizeof(ST_double));
                    if (shat0_lm) {
                        if (kiefer) {
                            /* Kiefer VCE: IID shat0 for LM test (see Wald block comment).
                               sigma2 = y_resid'y_resid/N for LM test */
                            ST_double ee_lm = 0.0;
                            for (i = 0; i < N; i++) {
                                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                ee_lm += w * y_resid[i] * y_resid[i];
                            }
                            ST_double sigma2_lm = ee_lm / (ST_double)N_eff;
                            for (ST_int ki = 0; ki < L; ki++) {
                                for (ST_int kj = 0; kj <= ki; kj++) {
                                    ST_double sum = 0.0;
                                    for (i = 0; i < N; i++) {
                                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                        sum += w * Z_excl[ki * N + i] * Z_excl[kj * N + i];
                                    }
                                    shat0_lm[ki * L + kj] = sigma2_lm * sum;
                                    if (ki != kj) shat0_lm[kj * L + ki] = sigma2_lm * sum;
                                }
                            }
                        } else if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && cluster_ids && cluster2_ids && num_clusters > 0 && num_clusters2 > 0) {
                            /* Two-way cluster-robust: CGM shat0_lm = shat1 + shat2 - shat_int */
                            compute_twoway_shat0(Z_excl, y_resid, weights, weight_type,
                                                 cluster_ids, num_clusters,
                                                 cluster2_ids, num_clusters2,
                                                 N, L, shat0_lm);
                        } else if ((vce_type == CIVREGHDFE_VCE_CLUSTER || vce_type == CIVREGHDFE_VCE_DKRAAY) && cluster_ids && num_clusters > 0) {
                            ST_double *cluster_Zy = (ST_double *)calloc(num_clusters * L, sizeof(ST_double));
                            if (cluster_Zy) {
                                /* Sum Z_excl * y_resid within each cluster */
                                for (i = 0; i < N; i++) {
                                    ST_int c = cluster_ids[i] - 1;
                                    if (c < 0 || c >= num_clusters) continue;
                                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                    for (k = 0; k < L; k++) {
                                        cluster_Zy[c * L + k] += w * Z_excl[k * N + i] * y_resid[i];
                                    }
                                }

                                /* Check if HAC kernel should be applied (dkraay only, kiefer handled above) */
                                if (vce_type == CIVREGHDFE_VCE_DKRAAY && kernel_type > 0 && bw > 0) {
                                    /* HAC with kernel weights */
                                    for (ST_int t = 0; t < num_clusters; t++) {
                                        for (ST_int s = 0; s < num_clusters; s++) {
                                            ST_int lag = (t > s) ? (t - s) : (s - t);
                                            ST_double kw = civreghdfe_kernel_weight(test_kernel, lag, bw);
                                            if (fabs(kw) < 1e-10) continue;

                                            for (ST_int ki = 0; ki < L; ki++) {
                                                for (ST_int kj = 0; kj < L; kj++) {
                                                    shat0_lm[kj * L + ki] += kw *
                                                        cluster_Zy[t * L + ki] * cluster_Zy[s * L + kj];
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    /* Standard cluster-robust */
                                    for (ST_int c = 0; c < num_clusters; c++) {
                                        for (ST_int ki = 0; ki < L; ki++) {
                                            for (ST_int kj = 0; kj < L; kj++) {
                                                shat0_lm[kj * L + ki] +=
                                                    cluster_Zy[c * L + ki] * cluster_Zy[c * L + kj];
                                            }
                                        }
                                    }
                                }
                                free(cluster_Zy);
                            }
                        } else if (hac_panel_ids && num_hac_panels > 0 && kernel_type > 0 && bw > 0) {
                            /* Panel-aware HAC: hac_panel_ids contains panel IDs
                               Apply kernel weights only within panels */
                            ST_int *panel_counts = (ST_int *)calloc(num_hac_panels, sizeof(ST_int));
                            ST_int *panel_starts = (ST_int *)calloc(num_hac_panels + 1, sizeof(ST_int));
                            ST_int *obs_by_panel = (ST_int *)calloc(N, sizeof(ST_int));

                            if (panel_counts && panel_starts && obs_by_panel) {
                                /* Count observations per panel */
                                for (i = 0; i < N; i++) {
                                    ST_int p = hac_panel_ids[i] - 1;
                                    if (p >= 0 && p < num_hac_panels) {
                                        panel_counts[p]++;
                                    }
                                }

                                /* Compute start indices */
                                panel_starts[0] = 0;
                                for (ST_int p = 0; p < num_hac_panels; p++) {
                                    panel_starts[p + 1] = panel_starts[p] + panel_counts[p];
                                }

                                /* Reset counts */
                                memset(panel_counts, 0, num_hac_panels * sizeof(ST_int));

                                /* Fill obs_by_panel */
                                for (i = 0; i < N; i++) {
                                    ST_int p = hac_panel_ids[i] - 1;
                                    if (p >= 0 && p < num_hac_panels) {
                                        ST_int idx = panel_starts[p] + panel_counts[p];
                                        obs_by_panel[idx] = i;
                                        panel_counts[p]++;
                                    }
                                }

                                /* Compute HAC within each panel */
                                for (ST_int p = 0; p < num_hac_panels; p++) {
                                    ST_int start = panel_starts[p];
                                    ST_int end = panel_starts[p + 1];
                                    ST_int T_p = end - start;

                                    for (ST_int t1 = 0; t1 < T_p; t1++) {
                                        ST_int i1 = obs_by_panel[start + t1];
                                        for (ST_int t2 = 0; t2 < T_p; t2++) {
                                            ST_int i2 = obs_by_panel[start + t2];
                                            ST_int time_lag = (t1 > t2) ? (t1 - t2) : (t2 - t1);

                                            ST_double kw = civreghdfe_kernel_weight(test_kernel, time_lag, bw);
                                            if (fabs(kw) < 1e-10) continue;

                                            ST_double w1 = (weights && weight_type != 0) ? weights[i1] : 1.0;
                                            ST_double w2 = (weights && weight_type != 0) ? weights[i2] : 1.0;
                                            ST_double y_prod = w1 * y_resid[i1] * w2 * y_resid[i2];

                                            for (ST_int ki = 0; ki < L; ki++) {
                                                for (ST_int kj = 0; kj < L; kj++) {
                                                    shat0_lm[kj * L + ki] += kw * y_prod *
                                                        Z_excl[ki * N + i1] * Z_excl[kj * N + i2];
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            if (panel_counts) free(panel_counts);
                            if (panel_starts) free(panel_starts);
                            if (obs_by_panel) free(obs_by_panel);
                        } else if (vce_type == CIVREGHDFE_VCE_ROBUST && kernel_type > 0 && bw > 0) {
                            /* Standard HAC (kernel + bw without panel info)
                               Apply kernel weights at observation level */
                            for (ST_int t = 0; t < N; t++) {
                                for (ST_int s = 0; s < N; s++) {
                                    ST_int lag = (t > s) ? (t - s) : (s - t);
                                    ST_double kw = civreghdfe_kernel_weight(test_kernel, lag, bw);
                                    if (fabs(kw) < 1e-10) continue;

                                    ST_double wt = (weights && weight_type != 0) ? weights[t] : 1.0;
                                    ST_double ws = (weights && weight_type != 0) ? weights[s] : 1.0;
                                    ST_double y_prod = wt * y_resid[t] * ws * y_resid[s];

                                    for (ST_int ki = 0; ki < L; ki++) {
                                        for (ST_int kj = 0; kj < L; kj++) {
                                            shat0_lm[kj * L + ki] += kw * y_prod *
                                                Z_excl[ki * N + t] * Z_excl[kj * N + s];
                                        }
                                    }
                                }
                            }
                        } else {
                            /* HC robust (no kernel) */
                            for (i = 0; i < N; i++) {
                                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                /* For aw/pw: w^2 * y^2; for fw/iw: w * y^2 */
                                ST_double y2w = (weight_type == 1 || weight_type == 3) ?
                                    w * w * y_resid[i] * y_resid[i] : w * y_resid[i] * y_resid[i];
                                for (ST_int ki = 0; ki < L; ki++) {
                                    for (ST_int kj = 0; kj < L; kj++) {
                                        shat0_lm[kj * L + ki] += y2w * Z_excl[ki * N + i] * Z_excl[kj * N + i];
                                    }
                                }
                            }
                        }

                        ST_double *shat0_lm_inv = (ST_double *)calloc(L * L, sizeof(ST_double));
                        if (shat0_lm_inv) {
                            memcpy(shat0_lm_inv, shat0_lm, L * L * sizeof(ST_double));
                            ST_int lm_inv_ok = (cholesky(shat0_lm_inv, L) == 0);
                            if (lm_inv_ok) {
                                ST_double chol_max = 0.0, chol_min = 1e308;
                                for (ST_int ci = 0; ci < L; ci++) {
                                    ST_double d = fabs(shat0_lm_inv[ci * L + ci]);
                                    if (d > chol_max) chol_max = d;
                                    if (d < chol_min) chol_min = d;
                                }
                                if (chol_max > 0 && chol_min / chol_max < 1e-8) {
                                    lm_inv_ok = 0;
                                } else {
                                    lm_inv_ok = (invert_from_cholesky(shat0_lm_inv, L, shat0_lm_inv) == 0);
                                }
                            }
                            if (!lm_inv_ok) {
                                lm_inv_ok = (ctools_invsym(shat0_lm, L, shat0_lm_inv, NULL) == 0);
                            }
                            if (lm_inv_ok) {

                                /* Compute Z_excl'y_resid */
                                ST_double *ZtY = (ST_double *)calloc(L, sizeof(ST_double));
                                if (ZtY) {
                                    for (k = 0; k < L; k++) {
                                        ST_double sum = 0.0;
                                        for (i = 0; i < N; i++) {
                                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                            sum += w * Z_excl[k * N + i] * y_resid[i];
                                        }
                                        ZtY[k] = sum;
                                    }

                                    ST_double *shat0_lm_inv_ZtY = (ST_double *)calloc(L, sizeof(ST_double));
                                    if (shat0_lm_inv_ZtY) {
                                        for (ST_int ki = 0; ki < L; ki++) {
                                            ST_double sum = 0.0;
                                            for (ST_int kj = 0; kj < L; kj++) {
                                                sum += shat0_lm_inv[kj * L + ki] * ZtY[kj];
                                            }
                                            shat0_lm_inv_ZtY[ki] = sum;
                                        }

                                        ST_double quad_lm = 0.0;
                                        for (ST_int ki = 0; ki < L; ki++) {
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
                    free(y_resid);
                }
            }

            free(v);
            free(shat0);
            free(shat0_inv);
            free(X_endog_partial);
            free(Z_excl_partial);
            free(ZtZ_excl);
            free(ZtX);
            free(pi_excl);
        }

    } else if (K_endog > 1) {
        /* Multiple endogenous variables: canonical correlations */
        /* Must partial out exogenous instruments (FWL) before computing
           canonical correlations, matching the K_endog==1 path */

        const ST_double *Z_exog = Z;                /* First K_exog columns */
        const ST_double *Z_excl_raw = Z + K_exog * N;  /* Last L columns */

        /* Allocate partialled-out copies of X_endog and Z_excl */
        ST_double *Y_partial = (ST_double *)malloc(N * K_endog * sizeof(ST_double));
        ST_double *Ze_partial = (ST_double *)malloc(N * L * sizeof(ST_double));
        if (!Y_partial || !Ze_partial) {
            free(Y_partial); free(Ze_partial);
            return;
        }

        /* Copy raw values */
        memcpy(Y_partial, X_endog, N * K_endog * sizeof(ST_double));
        memcpy(Ze_partial, Z_excl_raw, N * L * sizeof(ST_double));

        /* Partial out exogenous instruments if K_exog > 0 */
        if (K_exog > 0) {
            ST_double *ZxZx = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
            ST_double *ZxZx_inv = (ST_double *)calloc(K_exog * K_exog, sizeof(ST_double));
            if (ZxZx && ZxZx_inv) {
                /* Compute Z_exog'Z_exog */
                for (ST_int l1 = 0; l1 < K_exog; l1++) {
                    for (ST_int l2 = 0; l2 < K_exog; l2++) {
                        ST_double sum = 0.0;
                        for (i = 0; i < N; i++) {
                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                            sum += w * Z_exog[l1 * N + i] * Z_exog[l2 * N + i];
                        }
                        ZxZx[l2 * K_exog + l1] = sum;
                    }
                }
                memcpy(ZxZx_inv, ZxZx, K_exog * K_exog * sizeof(ST_double));
                if (cholesky(ZxZx_inv, K_exog) == 0 &&
                    invert_from_cholesky(ZxZx_inv, K_exog, ZxZx_inv) == 0) {

                    /* Partial out Z_exog from each column of X_endog */
                    for (ST_int e = 0; e < K_endog; e++) {
                        ST_double *ZxY = (ST_double *)calloc(K_exog, sizeof(ST_double));
                        if (ZxY) {
                            for (ST_int l = 0; l < K_exog; l++) {
                                ST_double sum = 0.0;
                                for (i = 0; i < N; i++) {
                                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                    sum += w * Z_exog[l * N + i] * X_endog[e * N + i];
                                }
                                ZxY[l] = sum;
                            }
                            for (ST_int l1 = 0; l1 < K_exog; l1++) {
                                ST_double gamma = 0.0;
                                for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                    gamma += ZxZx_inv[l2 * K_exog + l1] * ZxY[l2];
                                }
                                for (i = 0; i < N; i++) {
                                    Y_partial[e * N + i] -= Z_exog[l1 * N + i] * gamma;
                                }
                            }
                            free(ZxY);
                        }
                    }

                    /* Partial out Z_exog from each column of Z_excl */
                    for (ST_int col = 0; col < L; col++) {
                        ST_double *ZxZcol = (ST_double *)calloc(K_exog, sizeof(ST_double));
                        if (ZxZcol) {
                            for (ST_int l = 0; l < K_exog; l++) {
                                ST_double sum = 0.0;
                                for (i = 0; i < N; i++) {
                                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                    sum += w * Z_exog[l * N + i] * Z_excl_raw[col * N + i];
                                }
                                ZxZcol[l] = sum;
                            }
                            for (ST_int l1 = 0; l1 < K_exog; l1++) {
                                ST_double gamma_z = 0.0;
                                for (ST_int l2 = 0; l2 < K_exog; l2++) {
                                    gamma_z += ZxZx_inv[l2 * K_exog + l1] * ZxZcol[l2];
                                }
                                for (i = 0; i < N; i++) {
                                    Ze_partial[col * N + i] -= Z_exog[l1 * N + i] * gamma_z;
                                }
                            }
                            free(ZxZcol);
                        }
                    }
                }
            }
            free(ZxZx);
            free(ZxZx_inv);
        }

        /* Now compute cross-products on partialled data */
        /* Y'Y where Y = Y_partial */
        ST_double *YtY = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
        /* Ze'Ze where Ze = Ze_partial */
        ST_double *ZeZe = (ST_double *)calloc(L * L, sizeof(ST_double));
        /* Ze'Y */
        ST_double *ZeY = (ST_double *)calloc(L * K_endog, sizeof(ST_double));

        if (!YtY || !ZeZe || !ZeY) {
            free(Y_partial); free(Ze_partial);
            free(YtY); free(ZeZe); free(ZeY);
            return;
        }

        for (ST_int e1 = 0; e1 < K_endog; e1++) {
            for (ST_int e2 = 0; e2 < K_endog; e2++) {
                ST_double sum = 0.0;
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    sum += w * Y_partial[e1 * N + i] * Y_partial[e2 * N + i];
                }
                YtY[e2 * K_endog + e1] = sum;
            }
        }

        for (ST_int l1 = 0; l1 < L; l1++) {
            for (ST_int l2 = 0; l2 < L; l2++) {
                ST_double sum = 0.0;
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    sum += w * Ze_partial[l1 * N + i] * Ze_partial[l2 * N + i];
                }
                ZeZe[l2 * L + l1] = sum;
            }
        }

        for (ST_int l = 0; l < L; l++) {
            for (ST_int e = 0; e < K_endog; e++) {
                ST_double sum = 0.0;
                for (i = 0; i < N; i++) {
                    ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                    sum += w * Ze_partial[l * N + i] * Y_partial[e * N + i];
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

                        /* temp2 = ZeY' * temp1 = ZeY' * inv(ZeZe) * ZeY */
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
                                ctools_matmul_ab(temp2_mat, YtY_inv_half, K_endog, K_endog, K_endog, temp3_mat);
                                ctools_matmul_ab(YtY_inv_half, temp3_mat, K_endog, K_endog, K_endog, M);

                                /* Find minimum eigenvalue */
                                ST_double *eigenvalues = (ST_double *)calloc(K_endog, sizeof(ST_double));
                                if (eigenvalues) {
                                    civreghdfe_jacobi_eigenvalues(M, K_endog, eigenvalues);

                                    ST_double min_eval = eigenvalues[0];
                                    for (ST_int e = 1; e < K_endog; e++) {
                                        if (eigenvalues[e] < min_eval) min_eval = eigenvalues[e];
                                    }

                                    /* Anderson LM = N_eff * min_eigenvalue */
                                    *underid_stat = (ST_double)N_eff * min_eval;
                                    *underid_df = L - K_endog + 1;
                                    if (*underid_df < 1) *underid_df = 1;

                                    /* Cragg-Donald F */
                                    ST_int df_cd = N_eff - df_a - K_iv;
                                    if (df_cd <= 0) df_cd = 1;
                                    if (min_eval > 0.0 && min_eval < 1.0) {
                                        ST_double cd = min_eval / (1.0 - min_eval);
                                        *cd_f = cd * ((ST_double)df_cd / (ST_double)L);
                                    }

                                    /* KP Wald F for robust/cluster VCE */
                                    if (vce_type != 0) {
                                        /* Compute first-stage coefficient matrix:
                                           Pi_excl = (Ze'Ze)^-1 * Ze'Y  (L x K_endog)
                                           Already have this in temp1_mat */

                                        /* Compute first-stage residuals: V = Y_partial - Ze_partial * Pi_excl
                                           V is N x K_endog */
                                        ST_double *V_resid = (ST_double *)calloc(N * K_endog, sizeof(ST_double));
                                        if (V_resid) {
                                            for (ST_int e = 0; e < K_endog; e++) {
                                                for (i = 0; i < N; i++) {
                                                    ST_double pred = 0.0;
                                                    for (k = 0; k < L; k++) {
                                                        pred += Ze_partial[k * N + i] * temp1_mat[e * L + k];
                                                    }
                                                    V_resid[e * N + i] = Y_partial[e * N + i] - pred;
                                                }
                                            }

                                            /* Build shat0: L*K_endog x L*K_endog "vec" covariance
                                               For the KP test we need shat0 of vec(gbar) where
                                               gbar = (1/N) * Ze' * V_resid

                                               We compute the Kronecker-structured form:
                                               For each pair of endogenous vars (e1,e2) and instrument
                                               pair (l1,l2), build the (L*K_endog x L*K_endog) matrix.

                                               Index mapping: row = e1*L + l1, col = e2*L + l2 */
                                            ST_int shat_dim = L * K_endog;
                                            ST_double *shat0_kp = (ST_double *)calloc(shat_dim * shat_dim, sizeof(ST_double));

                                            if (shat0_kp) {
                                                if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && cluster_ids && cluster2_ids &&
                                                    num_clusters > 0 && num_clusters2 > 0) {
                                                    /* Two-way cluster-robust shat0_kp: CGM = shat1 + shat2 - shat_int */
                                                    /* Use compute_twoway_shat0 for each (l1*K_endog+e1, l2*K_endog+e2) block.
                                                       Treat V_resid and Ze_partial as a combined shat_dim-column residual vector. */
                                                    /* Compute intersection IDs */
                                                    ST_int *int_ids_kp = NULL;
                                                    ST_int num_int_kp = compute_intersection_ids(cluster_ids, num_clusters,
                                                                                                 cluster2_ids, num_clusters2,
                                                                                                 N, &int_ids_kp);
                                                    if (num_int_kp > 0 && int_ids_kp) {
                                                        /* Build cluster moments for each dimension and combine */
                                                        /* Cluster moments: sum_i_in_g w_i * (ze_l * v_e) for each (e,l) */
                                                        ST_double *cm1 = (ST_double *)calloc(num_clusters * shat_dim, sizeof(ST_double));
                                                        ST_double *cm2 = (ST_double *)calloc(num_clusters2 * shat_dim, sizeof(ST_double));
                                                        ST_double *cm_int = (ST_double *)calloc(num_int_kp * shat_dim, sizeof(ST_double));
                                                        if (cm1 && cm2 && cm_int) {
                                                            for (i = 0; i < N; i++) {
                                                                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                                                ST_int c1 = cluster_ids[i] - 1;
                                                                ST_int c2 = cluster2_ids[i] - 1;
                                                                ST_int ci = int_ids_kp[i] - 1;
                                                                for (ST_int e = 0; e < K_endog; e++) {
                                                                    for (k = 0; k < L; k++) {
                                                                        ST_double moment = w * Ze_partial[k * N + i] * V_resid[e * N + i];
                                                                        ST_int idx = e * L + k;
                                                                        if (c1 >= 0 && c1 < num_clusters) cm1[c1 * shat_dim + idx] += moment;
                                                                        if (c2 >= 0 && c2 < num_clusters2) cm2[c2 * shat_dim + idx] += moment;
                                                                        if (ci >= 0 && ci < num_int_kp) cm_int[ci * shat_dim + idx] += moment;
                                                                    }
                                                                }
                                                            }
                                                            /* shat = shat1 + shat2 - shat_int */
                                                            for (ST_int c = 0; c < num_clusters; c++) {
                                                                for (ST_int r = 0; r < shat_dim; r++) {
                                                                    for (ST_int s = 0; s < shat_dim; s++) {
                                                                        shat0_kp[s * shat_dim + r] += cm1[c * shat_dim + r] * cm1[c * shat_dim + s];
                                                                    }
                                                                }
                                                            }
                                                            for (ST_int c = 0; c < num_clusters2; c++) {
                                                                for (ST_int r = 0; r < shat_dim; r++) {
                                                                    for (ST_int s = 0; s < shat_dim; s++) {
                                                                        shat0_kp[s * shat_dim + r] += cm2[c * shat_dim + r] * cm2[c * shat_dim + s];
                                                                    }
                                                                }
                                                            }
                                                            for (ST_int c = 0; c < num_int_kp; c++) {
                                                                for (ST_int r = 0; r < shat_dim; r++) {
                                                                    for (ST_int s = 0; s < shat_dim; s++) {
                                                                        shat0_kp[s * shat_dim + r] -= cm_int[c * shat_dim + r] * cm_int[c * shat_dim + s];
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        free(cm1); free(cm2); free(cm_int);
                                                        free(int_ids_kp);
                                                    }
                                                } else if ((vce_type == CIVREGHDFE_VCE_CLUSTER || vce_type == CIVREGHDFE_VCE_DKRAAY) &&
                                                    cluster_ids && num_clusters > 0) {
                                                    /* Cluster-robust shat0 */
                                                    ST_double *cluster_moments = (ST_double *)calloc(num_clusters * shat_dim, sizeof(ST_double));
                                                    if (cluster_moments) {
                                                        for (i = 0; i < N; i++) {
                                                            ST_int c = cluster_ids[i] - 1;
                                                            if (c < 0 || c >= num_clusters) continue;
                                                            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                                            for (ST_int e = 0; e < K_endog; e++) {
                                                                for (k = 0; k < L; k++) {
                                                                    cluster_moments[c * shat_dim + e * L + k] +=
                                                                        w * Ze_partial[k * N + i] * V_resid[e * N + i];
                                                                }
                                                            }
                                                        }
                                                        for (ST_int c = 0; c < num_clusters; c++) {
                                                            for (ST_int r = 0; r < shat_dim; r++) {
                                                                for (ST_int s = 0; s < shat_dim; s++) {
                                                                    shat0_kp[s * shat_dim + r] +=
                                                                        cluster_moments[c * shat_dim + r] *
                                                                        cluster_moments[c * shat_dim + s];
                                                                }
                                                            }
                                                        }
                                                        free(cluster_moments);
                                                    }
                                                } else {
                                                    /* HC robust shat0 */
                                                    for (i = 0; i < N; i++) {
                                                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                                                        ST_double w2 = (weight_type == 1 || weight_type == 3) ? w * w : w;
                                                        for (ST_int e1 = 0; e1 < K_endog; e1++) {
                                                            for (ST_int l1 = 0; l1 < L; l1++) {
                                                                ST_double zv1 = Ze_partial[l1 * N + i] * V_resid[e1 * N + i];
                                                                ST_int r = e1 * L + l1;
                                                                for (ST_int e2 = 0; e2 < K_endog; e2++) {
                                                                    for (ST_int l2 = 0; l2 < L; l2++) {
                                                                        ST_double zv2 = Ze_partial[l2 * N + i] * V_resid[e2 * N + i];
                                                                        ST_int s = e2 * L + l2;
                                                                        shat0_kp[s * shat_dim + r] += w2 * zv1 * zv2;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }

                                                /* Verify shat0_kp is positive definite */
                                                {
                                                    ST_double *shat0_kp_chk = (ST_double *)calloc(shat_dim * shat_dim, sizeof(ST_double));
                                                    ST_int kp_ok = 0;
                                                    if (shat0_kp_chk) {
                                                        memcpy(shat0_kp_chk, shat0_kp, shat_dim * shat_dim * sizeof(ST_double));
                                                        kp_ok = (cholesky(shat0_kp_chk, shat_dim) == 0);
                                                        free(shat0_kp_chk);
                                                    }

                                                    if (kp_ok) {
                                                        /* Compute (Ze'Ze)^{-1/2} */
                                                        ST_double *ZeZe_inv_half = (ST_double *)calloc(L * L, sizeof(ST_double));

                                                        if (ZeZe_inv_half) {
                                                            civreghdfe_matpowersym_neg_half(ZeZe, L, ZeZe_inv_half);
                                                            {

                                                                /* Build theta = ZeZe_inv_half * ZeY * YtY_inv_half
                                                                   theta is L x K_endog */
                                                                ST_double *theta_temp = (ST_double *)calloc(L * K_endog, sizeof(ST_double));
                                                                ST_double *theta = (ST_double *)calloc(L * K_endog, sizeof(ST_double));

                                                                if (theta_temp && theta) {
                                                                    /* theta_temp = ZeZe_inv_half * ZeY (L x K_endog) */
                                                                    for (ST_int e = 0; e < K_endog; e++) {
                                                                        for (ST_int l1 = 0; l1 < L; l1++) {
                                                                            ST_double sum = 0.0;
                                                                            for (ST_int l2 = 0; l2 < L; l2++) {
                                                                                sum += ZeZe_inv_half[l2 * L + l1] * ZeY[e * L + l2];
                                                                            }
                                                                            theta_temp[e * L + l1] = sum;
                                                                        }
                                                                    }
                                                                    /* theta = theta_temp * YtY_inv_half (L x K_endog) */
                                                                    for (ST_int l = 0; l < L; l++) {
                                                                        for (ST_int e1 = 0; e1 < K_endog; e1++) {
                                                                            ST_double sum = 0.0;
                                                                            for (ST_int e2 = 0; e2 < K_endog; e2++) {
                                                                                sum += theta_temp[e2 * L + l] * YtY_inv_half[e2 * K_endog + e1];
                                                                            }
                                                                            theta[e1 * L + l] = sum;
                                                                        }
                                                                    }

                                                                    /* Build kronecker-transformed variance:
                                                                       vhat = (YtY_inv_half' kron ZeZe_inv_half) * shat0_kp * (YtY_inv_half kron ZeZe_inv_half')
                                                                       This is shat_dim x shat_dim */
                                                                    ST_double *kron_left = (ST_double *)calloc(shat_dim * shat_dim, sizeof(ST_double));
                                                                    if (kron_left) {
                                                                        /* Compute (YtY_inv_half kron ZeZe_inv_half) */
                                                                        for (ST_int e1 = 0; e1 < K_endog; e1++) {
                                                                            for (ST_int e2 = 0; e2 < K_endog; e2++) {
                                                                                for (ST_int l1 = 0; l1 < L; l1++) {
                                                                                    for (ST_int l2 = 0; l2 < L; l2++) {
                                                                                        ST_int r = e1 * L + l1;
                                                                                        ST_int c = e2 * L + l2;
                                                                                        kron_left[c * shat_dim + r] =
                                                                                            YtY_inv_half[e2 * K_endog + e1] *
                                                                                            ZeZe_inv_half[l2 * L + l1];
                                                                                    }
                                                                                }
                                                                            }
                                                                        }

                                                                        /* vhat = kron_left' * shat0_kp * kron_left */
                                                                        ST_double *vhat_temp = (ST_double *)calloc(shat_dim * shat_dim, sizeof(ST_double));
                                                                        ST_double *vhat = (ST_double *)calloc(shat_dim * shat_dim, sizeof(ST_double));

                                                                        if (vhat_temp && vhat) {
                                                                            /* vhat_temp = shat0_kp * kron_left */
                                                                            for (ST_int r = 0; r < shat_dim; r++) {
                                                                                for (ST_int c = 0; c < shat_dim; c++) {
                                                                                    ST_double sum = 0.0;
                                                                                    for (ST_int m = 0; m < shat_dim; m++) {
                                                                                        sum += shat0_kp[m * shat_dim + r] * kron_left[c * shat_dim + m];
                                                                                    }
                                                                                    vhat_temp[c * shat_dim + r] = sum;
                                                                                }
                                                                            }
                                                                            /* vhat = kron_left' * vhat_temp */
                                                                            for (ST_int r = 0; r < shat_dim; r++) {
                                                                                for (ST_int c = 0; c < shat_dim; c++) {
                                                                                    ST_double sum = 0.0;
                                                                                    for (ST_int m = 0; m < shat_dim; m++) {
                                                                                        sum += kron_left[r * shat_dim + m] * vhat_temp[c * shat_dim + m];
                                                                                    }
                                                                                    vhat[c * shat_dim + r] = sum;
                                                                                }
                                                                            }

                                                                            /* Invert vhat */
                                                                            ST_double *vhat_inv = (ST_double *)calloc(shat_dim * shat_dim, sizeof(ST_double));
                                                                            if (vhat_inv) {
                                                                                memcpy(vhat_inv, vhat, shat_dim * shat_dim * sizeof(ST_double));
                                                                                ST_int vhat_ok = (cholesky(vhat_inv, shat_dim) == 0);
                                                                                if (vhat_ok) vhat_ok = (invert_from_cholesky(vhat_inv, shat_dim, vhat_inv) == 0);

                                                                                if (vhat_ok) {
                                                                                    /* Compute theta' * vhat_inv * theta (K_endog x K_endog)
                                                                                       theta is stored as L x K_endog, with theta[e*L + l] */
                                                                                    /* First: vec(theta) as shat_dim vector */
                                                                                    /* vhat_inv * vec(theta) for each "column" of theta
                                                                                       Actually we need: theta' * vhat_inv * theta where
                                                                                       theta is treated as vec */

                                                                                    /* Build W = theta' * vhat_inv * theta (K_endog x K_endog)
                                                                                       theta[:,e] occupies positions [e*L .. e*L+L-1] in vec space.
                                                                                       W[e1,e2] = sum_{l1,l2} theta[e1*L+l1] * vhat_inv[e1*L+l1, e2*L+l2] * theta[e2*L+l2] */
                                                                                    ST_double *W = (ST_double *)calloc(K_endog * K_endog, sizeof(ST_double));
                                                                                    if (W) {
                                                                                        for (ST_int e1 = 0; e1 < K_endog; e1++) {
                                                                                            for (ST_int e2 = 0; e2 < K_endog; e2++) {
                                                                                                ST_double sum = 0.0;
                                                                                                for (ST_int l1 = 0; l1 < L; l1++) {
                                                                                                    ST_int r = e1 * L + l1;
                                                                                                    for (ST_int l2 = 0; l2 < L; l2++) {
                                                                                                        ST_int s = e2 * L + l2;
                                                                                                        sum += theta[e1 * L + l1] *
                                                                                                               vhat_inv[s * shat_dim + r] *
                                                                                                               theta[e2 * L + l2];
                                                                                                    }
                                                                                                }
                                                                                                W[e2 * K_endog + e1] = sum;
                                                                                            }
                                                                                        }

                                                                                        /* Find minimum eigenvalue of W */
                                                                                        ST_double *W_evals = (ST_double *)calloc(K_endog, sizeof(ST_double));
                                                                                        if (W_evals) {
                                                                                            civreghdfe_jacobi_eigenvalues(W, K_endog, W_evals);
                                                                                            ST_double kp_min_eval = W_evals[0];
                                                                                            for (ST_int e = 1; e < K_endog; e++) {
                                                                                                if (W_evals[e] < kp_min_eval) kp_min_eval = W_evals[e];
                                                                                            }

                                                                                            /* Convert to F statistic */
                                                                                            ST_int kp_dof;
                                                                                            ST_double kp_denom;
                                                                                            ST_double kp_cluster_adj;
                                                                                            if ((vce_type == CIVREGHDFE_VCE_CLUSTER || vce_type == CIVREGHDFE_VCE_CLUSTER2 ||
                                                                                                vce_type == CIVREGHDFE_VCE_DKRAAY) && !kiefer) {
                                                                                                kp_dof = N_eff - K_iv - df_a;
                                                                                                kp_denom = (ST_double)(N_eff - 1);
                                                                                                ST_int G_kp = num_clusters;
                                                                                                if (vce_type == CIVREGHDFE_VCE_CLUSTER2 && num_clusters2 > 0) {
                                                                                                    G_kp = (num_clusters < num_clusters2) ? num_clusters : num_clusters2;
                                                                                                }
                                                                                                kp_cluster_adj = (ST_double)(G_kp - 1) / (ST_double)G_kp;
                                                                                            } else {
                                                                                                kp_dof = N_eff - K_iv - df_a;
                                                                                                kp_denom = (ST_double)N_eff;
                                                                                                kp_cluster_adj = 1.0;
                                                                                            }
                                                                                            if (kp_dof <= 0) kp_dof = 1;

                                                                                            *kp_f = kp_min_eval / kp_denom * (ST_double)kp_dof * kp_cluster_adj / (ST_double)L;

                                                                                            free(W_evals);
                                                                                        }
                                                                                        free(W);
                                                                                    }
                                                                                }
                                                                                free(vhat_inv);
                                                                            }
                                                                        }
                                                                        free(vhat_temp);
                                                                        free(vhat);
                                                                        free(kron_left);
                                                                    }
                                                                }
                                                                free(theta_temp);
                                                                free(theta);
                                                            }
                                                            free(ZeZe_inv_half);
                                                        }
                                                    }
                                                }
                                                    free(V_resid);
                                                }
                                                free(shat0_kp);
                                            }
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

        free(Y_partial);
        free(Ze_partial);
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
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int kernel_type,
    ST_int bw,
    ST_int kiefer,
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels,
    const ST_double *y,
    const ST_double *X_all,
    const ST_double *ZtX,
    const ST_double *Zty,
    ST_double *sargan_stat,
    ST_int *overid_df
)
{
    ST_int K_total = K_exog + K_endog;
    ST_int i, j, k;

    /* ivreg2/ranktest always uses Bartlett kernel for test statistics */
    const ST_int test_kernel = (kernel_type > 0) ? 1 : 0;

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

            *sargan_stat = (ST_double)N_eff * quad_form / rss;
            free(ZtZ_inv_Ztr);
        }
    } else {
        /* Hansen J: r'Z(Z'ΩZ)^-1 Z'r */
        ST_double *ZOmegaZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        if (!ZOmegaZ) {
            free(Ztr);
            return;
        }

        if (kiefer && cluster_ids && num_clusters > 0 && kernel_type > 0 && bw > 0) {
            /* Kiefer VCE: homoskedastic HAC with Kronecker structure */
            compute_kiefer_shat0_homoskedastic(
                Z, resid, weights, weight_type,
                cluster_ids, num_clusters,
                N, N_eff, K_iv, 0, kernel_type, bw, ZOmegaZ);
        } else if ((vce_type == CIVREGHDFE_VCE_CLUSTER || vce_type == CIVREGHDFE_VCE_DKRAAY) && cluster_ids && num_clusters > 0) {
            /* Cluster-robust or Driscoll-Kraay Z'ΩZ */
            if (vce_type == CIVREGHDFE_VCE_DKRAAY && kernel_type > 0 && bw > 0) {
                /* HAC with kernel weights for dkraay */
                ST_double *cluster_Zr = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
                if (cluster_Zr) {
                    /* Sum Z * resid within each time period */
                    for (i = 0; i < N; i++) {
                        ST_int c = cluster_ids[i] - 1;
                        if (c < 0 || c >= num_clusters) continue;
                        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                        for (k = 0; k < K_iv; k++) {
                            cluster_Zr[c * K_iv + k] += w * Z[k * N + i] * resid[i];
                        }
                    }

                    /* HAC: sum_t sum_s kernel_weight(|t-s|) * Zr_t * Zr_s' */
                    for (ST_int t = 0; t < num_clusters; t++) {
                        for (ST_int s = 0; s < num_clusters; s++) {
                            ST_int lag = (t > s) ? (t - s) : (s - t);
                            ST_double kw = civreghdfe_kernel_weight(test_kernel, lag, bw);
                            if (fabs(kw) < 1e-10) continue;

                            for (ST_int ki = 0; ki < K_iv; ki++) {
                                for (ST_int kj = 0; kj < K_iv; kj++) {
                                    ZOmegaZ[kj * K_iv + ki] += kw *
                                        cluster_Zr[t * K_iv + ki] * cluster_Zr[s * K_iv + kj];
                                }
                            }
                        }
                    }
                    free(cluster_Zr);
                }
            } else {
                /* Standard cluster-robust */
                ivvce_compute_ZOmegaZ_cluster(Z, resid, weights, weight_type,
                                              cluster_ids, N, K_iv, num_clusters, ZOmegaZ);
            }
        } else if (vce_type == CIVREGHDFE_VCE_ROBUST && kernel_type > 0 && bw > 0) {
            /* HAC (bw + robust without cluster) - apply kernel weights */
            if (hac_panel_ids && num_hac_panels > 0) {
                /* Panel-aware HAC: apply kernel weights within each panel */
                ST_int *panel_counts = (ST_int *)calloc(num_hac_panels, sizeof(ST_int));
                ST_int *panel_starts = (ST_int *)calloc(num_hac_panels + 1, sizeof(ST_int));
                ST_int *obs_by_panel = (ST_int *)malloc(N * sizeof(ST_int));

                if (panel_counts && panel_starts && obs_by_panel) {
                    /* Count observations per panel */
                    for (i = 0; i < N; i++) {
                        ST_int p = hac_panel_ids[i] - 1;
                        if (p >= 0 && p < num_hac_panels) {
                            panel_counts[p]++;
                        }
                    }

                    /* Compute start indices */
                    panel_starts[0] = 0;
                    for (ST_int p = 0; p < num_hac_panels; p++) {
                        panel_starts[p + 1] = panel_starts[p] + panel_counts[p];
                    }

                    /* Reset counts for filling */
                    memset(panel_counts, 0, num_hac_panels * sizeof(ST_int));

                    /* Fill obs_by_panel */
                    for (i = 0; i < N; i++) {
                        ST_int p = hac_panel_ids[i] - 1;
                        if (p >= 0 && p < num_hac_panels) {
                            ST_int idx = panel_starts[p] + panel_counts[p];
                            obs_by_panel[idx] = i;
                            panel_counts[p]++;
                        }
                    }

                    /* Compute HAC within each panel */
                    for (ST_int p = 0; p < num_hac_panels; p++) {
                        ST_int start = panel_starts[p];
                        ST_int end = panel_starts[p + 1];
                        ST_int T_p = end - start;

                        for (ST_int t1 = 0; t1 < T_p; t1++) {
                            ST_int i1 = obs_by_panel[start + t1];
                            for (ST_int t2 = 0; t2 < T_p; t2++) {
                                ST_int i2 = obs_by_panel[start + t2];
                                ST_int time_lag = (t1 > t2) ? (t1 - t2) : (t2 - t1);

                                ST_double kw = civreghdfe_kernel_weight(test_kernel, time_lag, bw);
                                if (fabs(kw) < 1e-10) continue;

                                ST_double w1 = (weights && weight_type != 0) ? weights[i1] : 1.0;
                                ST_double w2 = (weights && weight_type != 0) ? weights[i2] : 1.0;
                                ST_double r_prod = w1 * resid[i1] * w2 * resid[i2];

                                for (j = 0; j < K_iv; j++) {
                                    for (k = 0; k < K_iv; k++) {
                                        ZOmegaZ[k * K_iv + j] += kw * r_prod * Z[j * N + i1] * Z[k * N + i2];
                                    }
                                }
                            }
                        }
                    }
                }
                if (panel_counts) free(panel_counts);
                if (panel_starts) free(panel_starts);
                if (obs_by_panel) free(obs_by_panel);
            } else {
                /* Non-panel HAC: use raw observation indices */
                for (ST_int t = 0; t < N; t++) {
                    for (ST_int s = 0; s < N; s++) {
                        ST_int lag = (t > s) ? (t - s) : (s - t);
                        ST_double kw = civreghdfe_kernel_weight(test_kernel, lag, bw);
                        if (fabs(kw) < 1e-10) continue;

                        ST_double wt = (weights && weight_type != 0) ? weights[t] : 1.0;
                        ST_double ws = (weights && weight_type != 0) ? weights[s] : 1.0;
                        ST_double r_prod = wt * resid[t] * ws * resid[s];

                        for (j = 0; j < K_iv; j++) {
                            for (k = 0; k < K_iv; k++) {
                                ZOmegaZ[k * K_iv + j] += kw * r_prod * Z[j * N + t] * Z[k * N + s];
                            }
                        }
                    }
                }
            }
        } else {
            /* HC robust (no kernel) */
            ivvce_compute_ZOmegaZ_robust(Z, resid, weights, weight_type, N, K_iv, ZOmegaZ);
        }

        ST_double *ZOmegaZ_inv = (ST_double *)malloc(K_iv * K_iv * sizeof(ST_double));
        if (ZOmegaZ_inv) {
            memcpy(ZOmegaZ_inv, ZOmegaZ, K_iv * K_iv * sizeof(ST_double));

            if (cholesky(ZOmegaZ_inv, K_iv) == 0) {
                invert_from_cholesky(ZOmegaZ_inv, K_iv, ZOmegaZ_inv);

                /* ivreg2's Hansen J requires efficient GMM residuals, not 2SLS residuals.
                   ivreg2's s_iegmm re-estimates beta using W = (Z'ΩZ)^{-1} as optimal weights:
                     beta_eff = (X'Z * (Z'ΩZ)^{-1} * Z'X)^{-1} * X'Z * (Z'ΩZ)^{-1} * Z'y
                   Then J = (Z'e_eff)' * (Z'ΩZ)^{-1} * (Z'e_eff) using the efficient residuals.
                   For homoskedastic VCE, efficient GMM = 2SLS so this only matters for
                   robust/cluster/HAC cases. */
                if (y && X_all && ZtX && Zty) {
                    /* Step 1: temp_sinv_ztx = ZOmegaZ_inv * ZtX  (K_iv × K_total) */
                    ST_double *temp_sinv_ztx = (ST_double *)malloc(K_iv * K_total * sizeof(ST_double));
                    ST_double *XZSZX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
                    ST_double *beta_eff = (ST_double *)calloc(K_total, sizeof(ST_double));
                    ST_double *rhs = (ST_double *)calloc(K_total, sizeof(ST_double));

                    if (temp_sinv_ztx && XZSZX && beta_eff && rhs) {
                        for (j = 0; j < K_total; j++) {
                            for (i = 0; i < K_iv; i++) {
                                ST_double sum = 0.0;
                                for (k = 0; k < K_iv; k++) {
                                    sum += ZOmegaZ_inv[k * K_iv + i] * ZtX[j * K_iv + k];
                                }
                                temp_sinv_ztx[j * K_iv + i] = sum;
                            }
                        }

                        /* Step 2: XZSZX = ZtX' * temp_sinv_ztx  (K_total × K_total) */
                        for (j = 0; j < K_total; j++) {
                            for (i = 0; i < K_total; i++) {
                                ST_double sum = 0.0;
                                for (k = 0; k < K_iv; k++) {
                                    sum += ZtX[i * K_iv + k] * temp_sinv_ztx[j * K_iv + k];
                                }
                                XZSZX[j * K_total + i] = sum;
                            }
                        }

                        /* Step 3: sinv_zty = ZOmegaZ_inv * Zty  (K_iv × 1) */
                        ST_double *sinv_zty = (ST_double *)calloc(K_iv, sizeof(ST_double));
                        if (sinv_zty) {
                            for (i = 0; i < K_iv; i++) {
                                ST_double sum = 0.0;
                                for (k = 0; k < K_iv; k++) {
                                    sum += ZOmegaZ_inv[k * K_iv + i] * Zty[k];
                                }
                                sinv_zty[i] = sum;
                            }

                            /* Step 4: rhs = ZtX' * sinv_zty  (K_total × 1) */
                            for (i = 0; i < K_total; i++) {
                                ST_double sum = 0.0;
                                for (k = 0; k < K_iv; k++) {
                                    sum += ZtX[i * K_iv + k] * sinv_zty[k];
                                }
                                rhs[i] = sum;
                            }
                            free(sinv_zty);
                        }

                        /* Step 5: Solve XZSZX * beta_eff = rhs via Cholesky */
                        ST_double *XZSZX_L = (ST_double *)malloc(K_total * K_total * sizeof(ST_double));
                        if (XZSZX_L) {
                            memcpy(XZSZX_L, XZSZX, K_total * K_total * sizeof(ST_double));
                            if (cholesky(XZSZX_L, K_total) == 0) {
                                /* Forward substitution: L * z = rhs */
                                for (i = 0; i < K_total; i++) {
                                    ST_double sum = rhs[i];
                                    for (j = 0; j < i; j++) {
                                        sum -= XZSZX_L[i * K_total + j] * beta_eff[j];
                                    }
                                    beta_eff[i] = sum / XZSZX_L[i * K_total + i];
                                }
                                /* Back substitution: L' * beta = z */
                                for (i = K_total - 1; i >= 0; i--) {
                                    ST_double sum = beta_eff[i];
                                    for (j = i + 1; j < K_total; j++) {
                                        sum -= XZSZX_L[j * K_total + i] * beta_eff[j];
                                    }
                                    beta_eff[i] = sum / XZSZX_L[i * K_total + i];
                                }

                                /* Step 6: Recompute Ztr using efficient GMM residuals */
                                for (k = 0; k < K_iv; k++) {
                                    ST_double sum = 0.0;
                                    const ST_double *z_col = Z + k * N;
                                    if (weights && weight_type != 0) {
                                        for (i = 0; i < N; i++) {
                                            ST_double e_eff = y[i];
                                            for (j = 0; j < K_total; j++) {
                                                e_eff -= X_all[j * N + i] * beta_eff[j];
                                            }
                                            sum += z_col[i] * weights[i] * e_eff;
                                        }
                                    } else {
                                        for (i = 0; i < N; i++) {
                                            ST_double e_eff = y[i];
                                            for (j = 0; j < K_total; j++) {
                                                e_eff -= X_all[j * N + i] * beta_eff[j];
                                            }
                                            sum += z_col[i] * e_eff;
                                        }
                                    }
                                    Ztr[k] = sum;
                                }
                            }
                            free(XZSZX_L);
                        }
                    }

                    if (temp_sinv_ztx) free(temp_sinv_ztx);
                    if (XZSZX) free(XZSZX);
                    if (beta_eff) free(beta_eff);
                    if (rhs) free(rhs);
                }

                /* Compute J = Ztr' * ZOmegaZ_inv * Ztr */
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
    ctools_matmul_atb(X_aug, X_aug, N, K_aug, K_aug, XaXa);

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
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    ST_int N_eff,
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
    ST_double rss_full,
    ST_double *cstat,
    ST_int *cstat_df
)
{
    ST_int i, j, k;
    ST_int K_total = K_exog + K_endog;
    ST_int K_rest = K_iv - n_orthog;  /* Instruments not being tested */
    int use_weights = (weights != NULL && weight_type != 0);

    /* Variables freed in cleanup — must be initialized before any goto */
    ST_double *ZrZr = NULL, *ZrZr_inv = NULL;
    ST_double *ZrX = NULL, *Zry = NULL;
    ST_double *PZrX = NULL, *PZry = NULL, *temp_rv = NULL;
    ST_double *PZrX_L = NULL, *PZrX_inv = NULL, *beta_rest = NULL;

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
        ST_int idx = orthog_indices[i] - 1;
        if (idx >= 0 && idx < K_iv - K_exog) {
            test_mask[K_exog + idx] = 1;
        }
    }

    /* Build Z_rest (instruments not being tested) */
    ST_double *Z_rest = (ST_double *)malloc((size_t)N * K_rest * sizeof(ST_double));
    if (!Z_rest) { free(test_mask); return; }

    ST_int rest_idx = 0;
    for (k = 0; k < K_iv; k++) {
        if (test_mask[k] == 0) {
            for (i = 0; i < N; i++) {
                Z_rest[rest_idx * N + i] = Z[k * N + i];
            }
            rest_idx++;
        }
    }
    free(test_mask);

    /* Build X_all = [X_exog, X_endog] (column-major) */
    ST_double *X_all = (ST_double *)malloc((size_t)N * K_total * sizeof(ST_double));
    if (!X_all) { free(Z_rest); return; }

    if (X_exog && K_exog > 0) {
        memcpy(X_all, X_exog, (size_t)N * K_exog * sizeof(ST_double));
    }
    if (X_endog && K_endog > 0) {
        memcpy(X_all + (size_t)N * K_exog, X_endog, (size_t)N * K_endog * sizeof(ST_double));
    }

    /* Re-estimate 2SLS with restricted instruments Z_rest.
       C-stat = J_full - J_restricted where J_restricted uses re-estimated
       coefficients and residuals from the restricted instrument set. */

    /* Step 1: Z_rest'Z_rest and its inverse */
    ZrZr = (ST_double *)calloc(K_rest * K_rest, sizeof(ST_double));
    ZrZr_inv = (ST_double *)calloc(K_rest * K_rest, sizeof(ST_double));
    if (!ZrZr || !ZrZr_inv) goto cleanup;

    if (use_weights)
        ctools_matmul_atdb(Z_rest, Z_rest, weights, N, K_rest, K_rest, ZrZr);
    else
        ctools_matmul_atb(Z_rest, Z_rest, N, K_rest, K_rest, ZrZr);
    memcpy(ZrZr_inv, ZrZr, K_rest * K_rest * sizeof(ST_double));
    if (cholesky(ZrZr_inv, K_rest) != 0) goto cleanup;
    if (invert_from_cholesky(ZrZr_inv, K_rest, ZrZr_inv) != 0) goto cleanup;

    /* Step 2: Z_rest'X and Z_rest'y */
    ZrX = (ST_double *)calloc(K_rest * K_total, sizeof(ST_double));
    Zry = (ST_double *)calloc(K_rest, sizeof(ST_double));
    if (!ZrX || !Zry) goto cleanup;

    if (use_weights)
        ctools_matmul_atdb(Z_rest, X_all, weights, N, K_rest, K_total, ZrX);
    else
        ctools_matmul_atb(Z_rest, X_all, N, K_rest, K_total, ZrX);
    for (k = 0; k < K_rest; k++) {
        ST_double sum = 0.0;
        if (use_weights) {
            for (i = 0; i < N; i++) sum += Z_rest[k * N + i] * weights[i] * y[i];
        } else {
            for (i = 0; i < N; i++) sum += Z_rest[k * N + i] * y[i];
        }
        Zry[k] = sum;
    }

    /* Step 3: 2SLS with Z_rest
       X_hat = Z_rest (Z_rest'Z_rest)^{-1} Z_rest'X
       beta_rest = (X_hat'X_hat)^{-1} X_hat'y
       which simplifies to: beta_rest = (X'P_Zr X)^{-1} X'P_Zr y
       where P_Zr = Z_rest (Z_rest'Z_rest)^{-1} Z_rest'
       i.e., beta_rest = (ZrX' ZrZr_inv ZrX)^{-1} ZrX' ZrZr_inv Zry */
    PZrX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));  /* X'P_Zr X */
    PZry = (ST_double *)calloc(K_total, sizeof(ST_double));           /* X'P_Zr y */
    temp_rv = (ST_double *)calloc(K_rest * K_total, sizeof(ST_double)); /* ZrZr_inv * ZrX */
    if (!PZrX || !PZry || !temp_rv) goto cleanup;

    /* temp_rv = ZrZr_inv * ZrX  (K_rest x K_total) */
    ctools_matmul_ab(ZrZr_inv, ZrX, K_rest, K_rest, K_total, temp_rv);

    /* PZrX = ZrX' * temp_rv  (K_total x K_total) */
    ctools_matmul_atb(ZrX, temp_rv, K_rest, K_total, K_total, PZrX);

    /* PZry = ZrX' * ZrZr_inv * Zry  (K_total x 1) */
    ST_double *temp_ry = (ST_double *)calloc(K_rest, sizeof(ST_double));
    if (!temp_ry) goto cleanup;

    for (k = 0; k < K_rest; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_rest; j++) sum += ZrZr_inv[j * K_rest + k] * Zry[j];
        temp_ry[k] = sum;
    }
    /* PZry = ZrX' * temp_ry where ZrX is K_rest x K_total (col-major) */
    for (k = 0; k < K_total; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_rest; j++) sum += ZrX[k * K_rest + j] * temp_ry[j];
        PZry[k] = sum;
    }
    free(temp_ry);

    /* Solve beta_rest = PZrX^{-1} * PZry */
    PZrX_L = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    PZrX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    beta_rest = (ST_double *)calloc(K_total, sizeof(ST_double));
    if (!PZrX_L || !PZrX_inv || !beta_rest) goto cleanup;

    memcpy(PZrX_L, PZrX, K_total * K_total * sizeof(ST_double));
    if (cholesky(PZrX_L, K_total) != 0) goto cleanup;
    if (invert_from_cholesky(PZrX_L, K_total, PZrX_inv) != 0) goto cleanup;

    for (k = 0; k < K_total; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_total; j++) sum += PZrX_inv[j * K_total + k] * PZry[j];
        beta_rest[k] = sum;
    }

    /* Step 4: Compute restricted residuals: e_rest = y - X_all * beta_rest */
    ST_double *resid_rest = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    if (!resid_rest) goto cleanup;

    for (i = 0; i < N; i++) {
        ST_double fitted = 0.0;
        for (k = 0; k < K_total; k++) {
            fitted += X_all[k * N + i] * beta_rest[k];
        }
        resid_rest[i] = y[i] - fitted;
    }

    /* Step 5: Compute restricted Sargan = N_eff * e_rest' P_Zr e_rest / rss_rest */
    ST_double *Zr_e = (ST_double *)calloc(K_rest, sizeof(ST_double));
    if (!Zr_e) { free(resid_rest); goto cleanup; }

    for (k = 0; k < K_rest; k++) {
        ST_double sum = 0.0;
        if (use_weights) {
            for (i = 0; i < N; i++) sum += Z_rest[k * N + i] * weights[i] * resid_rest[i];
        } else {
            for (i = 0; i < N; i++) sum += Z_rest[k * N + i] * resid_rest[i];
        }
        Zr_e[k] = sum;
    }

    ST_double *temp_r = (ST_double *)calloc(K_rest, sizeof(ST_double));
    if (!temp_r) { free(resid_rest); free(Zr_e); goto cleanup; }

    for (k = 0; k < K_rest; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_rest; j++) sum += ZrZr_inv[j * K_rest + k] * Zr_e[j];
        temp_r[k] = sum;
    }

    ST_double quad_rest = 0.0;
    for (k = 0; k < K_rest; k++) quad_rest += Zr_e[k] * temp_r[k];

    ST_double sargan_rest = (rss_full > 0) ? (ST_double)N_eff * quad_rest / rss_full : 0.0;

    /* C = J_full - J_restricted */
    *cstat = sargan_full - sargan_rest;
    if (*cstat < 0) *cstat = 0;

    free(resid_rest);
    free(Zr_e);
    free(temp_r);

cleanup:
    free(Z_rest);
    free(X_all);
    free(ZrZr);
    free(ZrZr_inv);
    free(ZrX);
    free(Zry);
    free(PZrX);
    free(PZry);
    free(temp_rv);
    free(PZrX_L);
    free(PZrX_inv);
    free(beta_rest);
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
    const ST_double *resid,
    const ST_double *ZtZ_inv,
    ST_int N,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const ST_int *endogtest_indices,
    ST_int n_endogtest,
    ST_double *endogtest_stat,
    ST_int *endogtest_df
)
{
    /*
        Difference-in-Sargan C-statistic approach (matches ivreg2):
        1. Build Z_aug = [Z, tested_endogenous_vars] (add tested vars to instruments)
        2. Run 2SLS with Z_aug: beta_exog, residuals, rss_exog
        3. Compute J_exog = N_eff * e_exog'P_{Z_aug}e_exog / rss_exog
        4. Compute J_orig = N_eff * e_orig'P_Z e_orig / rss_exog  (uses rss_exog!)
        5. C = J_exog - J_orig
    */
    ST_int K_total = K_exog + K_endog;
    ST_int K_iv_aug = K_iv + n_endogtest;
    ST_int i, j, k;
    int use_weights = (weights != NULL && weight_type != 0);

    *endogtest_stat = 0.0;
    *endogtest_df = n_endogtest;

    if (n_endogtest <= 0 || n_endogtest > K_endog) return;

    /* Allocate all working memory */
    ST_double *Z_aug = NULL, *X_all = NULL;
    ST_double *ZaZa = NULL, *ZaZa_inv = NULL;
    ST_double *ZaX = NULL, *Zay = NULL;
    ST_double *temp_rv = NULL, *PZaX = NULL, *PZay = NULL, *temp_ry = NULL;
    ST_double *PZaX_L = NULL, *PZaX_inv = NULL, *beta_exog = NULL;
    ST_double *resid_exog = NULL, *Za_e = NULL, *temp_a = NULL;
    ST_double *Z_e_orig = NULL, *temp_o = NULL;

    /* Step 1: Build Z_aug = [Z, tested_endogenous_vars] */
    Z_aug = (ST_double *)malloc((size_t)N * K_iv_aug * sizeof(ST_double));
    if (!Z_aug) return;

    /* Copy original Z */
    memcpy(Z_aug, Z, (size_t)N * K_iv * sizeof(ST_double));

    /* Append tested endogenous variables */
    for (ST_int t = 0; t < n_endogtest; t++) {
        ST_int idx = endogtest_indices[t] - 1;
        if (idx < 0 || idx >= K_endog) goto endogtest_cleanup;
        memcpy(Z_aug + (size_t)(K_iv + t) * N, X_endog + (size_t)idx * N,
               (size_t)N * sizeof(ST_double));
    }

    /* Build X_all = [X_exog, X_endog] */
    X_all = (ST_double *)malloc((size_t)N * K_total * sizeof(ST_double));
    if (!X_all) goto endogtest_cleanup;

    if (X_exog && K_exog > 0)
        memcpy(X_all, X_exog, (size_t)N * K_exog * sizeof(ST_double));
    if (X_endog && K_endog > 0)
        memcpy(X_all + (size_t)N * K_exog, X_endog, (size_t)N * K_endog * sizeof(ST_double));

    /* Step 2: 2SLS with Z_aug */
    /* Z_aug'Z_aug and its inverse */
    ZaZa = (ST_double *)calloc(K_iv_aug * K_iv_aug, sizeof(ST_double));
    ZaZa_inv = (ST_double *)calloc(K_iv_aug * K_iv_aug, sizeof(ST_double));
    if (!ZaZa || !ZaZa_inv) goto endogtest_cleanup;

    if (use_weights)
        ctools_matmul_atdb(Z_aug, Z_aug, weights, N, K_iv_aug, K_iv_aug, ZaZa);
    else
        ctools_matmul_atb(Z_aug, Z_aug, N, K_iv_aug, K_iv_aug, ZaZa);
    memcpy(ZaZa_inv, ZaZa, (size_t)K_iv_aug * K_iv_aug * sizeof(ST_double));
    if (cholesky(ZaZa_inv, K_iv_aug) != 0) goto endogtest_cleanup;
    if (invert_from_cholesky(ZaZa_inv, K_iv_aug, ZaZa_inv) != 0) goto endogtest_cleanup;

    /* Z_aug'X and Z_aug'y */
    ZaX = (ST_double *)calloc(K_iv_aug * K_total, sizeof(ST_double));
    Zay = (ST_double *)calloc(K_iv_aug, sizeof(ST_double));
    if (!ZaX || !Zay) goto endogtest_cleanup;

    if (use_weights)
        ctools_matmul_atdb(Z_aug, X_all, weights, N, K_iv_aug, K_total, ZaX);
    else
        ctools_matmul_atb(Z_aug, X_all, N, K_iv_aug, K_total, ZaX);
    for (k = 0; k < K_iv_aug; k++) {
        ST_double sum = 0.0;
        if (use_weights) {
            for (i = 0; i < N; i++) sum += Z_aug[(size_t)k * N + i] * weights[i] * y[i];
        } else {
            for (i = 0; i < N; i++) sum += Z_aug[(size_t)k * N + i] * y[i];
        }
        Zay[k] = sum;
    }

    /* beta_exog = (X'P_{Z_aug}X)^{-1} X'P_{Z_aug}y */
    /* temp_rv = ZaZa_inv * ZaX  (K_iv_aug x K_total) */
    temp_rv = (ST_double *)calloc(K_iv_aug * K_total, sizeof(ST_double));
    PZaX = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!temp_rv || !PZaX) goto endogtest_cleanup;

    ctools_matmul_ab(ZaZa_inv, ZaX, K_iv_aug, K_iv_aug, K_total, temp_rv);
    ctools_matmul_atb(ZaX, temp_rv, K_iv_aug, K_total, K_total, PZaX);

    /* PZay = ZaX' * ZaZa_inv * Zay */
    temp_ry = (ST_double *)calloc(K_iv_aug, sizeof(ST_double));
    PZay = (ST_double *)calloc(K_total, sizeof(ST_double));
    if (!temp_ry || !PZay) goto endogtest_cleanup;

    for (k = 0; k < K_iv_aug; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_iv_aug; j++) sum += ZaZa_inv[(size_t)j * K_iv_aug + k] * Zay[j];
        temp_ry[k] = sum;
    }
    for (k = 0; k < K_total; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_iv_aug; j++) sum += ZaX[(size_t)k * K_iv_aug + j] * temp_ry[j];
        PZay[k] = sum;
    }

    /* Solve beta_exog = PZaX^{-1} * PZay */
    PZaX_L = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    PZaX_inv = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    beta_exog = (ST_double *)calloc(K_total, sizeof(ST_double));
    if (!PZaX_L || !PZaX_inv || !beta_exog) goto endogtest_cleanup;

    memcpy(PZaX_L, PZaX, (size_t)K_total * K_total * sizeof(ST_double));
    if (cholesky(PZaX_L, K_total) != 0) goto endogtest_cleanup;
    if (invert_from_cholesky(PZaX_L, K_total, PZaX_inv) != 0) goto endogtest_cleanup;

    for (k = 0; k < K_total; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_total; j++) sum += PZaX_inv[(size_t)j * K_total + k] * PZay[j];
        beta_exog[k] = sum;
    }

    /* Step 3: Compute exogenous model residuals and RSS */
    resid_exog = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    if (!resid_exog) goto endogtest_cleanup;

    ST_double rss_exog = 0.0;
    for (i = 0; i < N; i++) {
        ST_double fitted = 0.0;
        for (k = 0; k < K_total; k++)
            fitted += X_all[(size_t)k * N + i] * beta_exog[k];
        resid_exog[i] = y[i] - fitted;
        if (use_weights)
            rss_exog += weights[i] * resid_exog[i] * resid_exog[i];
        else
            rss_exog += resid_exog[i] * resid_exog[i];
    }

    if (rss_exog <= 0.0) goto endogtest_cleanup;

    /* Step 4: Compute J_exog = N * e_exog' Z_aug (Z_aug'Z_aug)^{-1} Z_aug' e_exog / rss_exog */
    Za_e = (ST_double *)calloc(K_iv_aug, sizeof(ST_double));
    temp_a = (ST_double *)calloc(K_iv_aug, sizeof(ST_double));
    if (!Za_e || !temp_a) goto endogtest_cleanup;

    for (k = 0; k < K_iv_aug; k++) {
        ST_double sum = 0.0;
        if (use_weights) {
            for (i = 0; i < N; i++) sum += Z_aug[(size_t)k * N + i] * weights[i] * resid_exog[i];
        } else {
            for (i = 0; i < N; i++) sum += Z_aug[(size_t)k * N + i] * resid_exog[i];
        }
        Za_e[k] = sum;
    }
    for (k = 0; k < K_iv_aug; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_iv_aug; j++) sum += ZaZa_inv[(size_t)j * K_iv_aug + k] * Za_e[j];
        temp_a[k] = sum;
    }
    ST_double quad_exog = 0.0;
    for (k = 0; k < K_iv_aug; k++) quad_exog += Za_e[k] * temp_a[k];

    ST_double J_exog = (ST_double)N_eff * quad_exog / rss_exog;

    /* Step 5: Compute J_orig = N * e_orig' Z (Z'Z)^{-1} Z' e_orig / rss_exog
       (uses rss_exog, NOT rss_orig — this matches ivreg2's smatrix approach) */
    Z_e_orig = (ST_double *)calloc(K_iv, sizeof(ST_double));
    temp_o = (ST_double *)calloc(K_iv, sizeof(ST_double));
    if (!Z_e_orig || !temp_o) goto endogtest_cleanup;

    for (k = 0; k < K_iv; k++) {
        ST_double sum = 0.0;
        if (use_weights) {
            for (i = 0; i < N; i++) sum += Z[(size_t)k * N + i] * weights[i] * resid[i];
        } else {
            for (i = 0; i < N; i++) sum += Z[(size_t)k * N + i] * resid[i];
        }
        Z_e_orig[k] = sum;
    }
    for (k = 0; k < K_iv; k++) {
        ST_double sum = 0.0;
        for (j = 0; j < K_iv; j++) sum += ZtZ_inv[(size_t)j * K_iv + k] * Z_e_orig[j];
        temp_o[k] = sum;
    }
    ST_double quad_orig = 0.0;
    for (k = 0; k < K_iv; k++) quad_orig += Z_e_orig[k] * temp_o[k];

    ST_double J_orig = (ST_double)N_eff * quad_orig / rss_exog;

    /* C = J_exog - J_orig */
    *endogtest_stat = J_exog - J_orig;
    if (*endogtest_stat < 0) *endogtest_stat = 0;

endogtest_cleanup:
    free(Z_aug);
    free(X_all);
    free(ZaZa);
    free(ZaZa_inv);
    free(ZaX);
    free(Zay);
    free(temp_rv);
    free(PZaX);
    free(PZay);
    free(temp_ry);
    free(PZaX_L);
    free(PZaX_inv);
    free(beta_exog);
    free(resid_exog);
    free(Za_e);
    free(temp_a);
    free(Z_e_orig);
    free(temp_o);
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
    const ST_double *weights,
    ST_int weight_type,
    ST_int N_eff,
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
    int use_weights = (weights != NULL && weight_type != 0);

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

    if (use_weights)
        ctools_matmul_atdb(Z_rest, Z_rest, weights, N, K_rest, K_rest, ZrZr);
    else
        ctools_matmul_atb(Z_rest, Z_rest, N, K_rest, K_rest, ZrZr);
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
        if (use_weights) {
            ctools_matmul_atdb(Z_rest, X_endog, weights, N, K_rest, K_endog, ZrX);
            ctools_matmul_atdb(Z_rest, Z_test, weights, N, K_rest, n_redund, ZrZt);
        } else {
            ctools_matmul_atb(Z_rest, X_endog, N, K_rest, K_endog, ZrX);
            ctools_matmul_atb(Z_rest, Z_test, N, K_rest, n_redund, ZrZt);
        }

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
                    if (use_weights) {
                        ctools_matmul_atdb(Zt_resid, X_resid, weights, N, n_redund, K_endog, ZtX);
                        ctools_matmul_atdb(X_resid, X_resid, weights, N, K_endog, K_endog, XtX);
                        ctools_matmul_atdb(Zt_resid, Zt_resid, weights, N, n_redund, n_redund, ZtZt);
                    } else {
                        ctools_matmul_atb(Zt_resid, X_resid, N, n_redund, K_endog, ZtX);
                        ctools_matmul_atb(X_resid, X_resid, N, K_endog, K_endog, XtX);
                        ctools_matmul_atb(Zt_resid, Zt_resid, N, n_redund, n_redund, ZtZt);
                    }

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

                                *redund_stat = (ST_double)N_eff * trace;
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
