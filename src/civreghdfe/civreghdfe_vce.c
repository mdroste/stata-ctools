/*
    civreghdfe_vce.c
    Variance-Covariance Estimation for IV Regression

    Implements unadjusted, robust (HC), clustered, and HAC VCE
    for k-class IV estimators.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "civreghdfe_vce.h"
#include "../ctools_config.h"
#include "../ctools_matrix.h"
#include "civreghdfe_matrix.h"

/* Shared OLS functions */
#include "../ctools_ols.h"
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/*
    Compute unadjusted VCE for k-class estimator.
*/
static void ivvce_compute_unadjusted(
    const ST_double *XkX_inv,
    ST_double rss,
    ST_int df_r,
    ST_int K_total,
    ST_double *V
)
{
    ST_int i;
    ST_double sigma2 = rss / df_r;

    for (i = 0; i < K_total * K_total; i++) {
        V[i] = sigma2 * XkX_inv[i];
    }
}

/*
    Compute efficient GMM2S VCE.
*/
static void ivvce_compute_gmm2s(
    const ST_double *XZWZX_inv,
    ST_int K_total,
    ST_double *V
)
{
    ST_int i;
    /* For efficient GMM, V = (X'ZWZ'X)^-1 */
    for (i = 0; i < K_total * K_total; i++) {
        V[i] = XZWZX_inv[i];
    }
}

/*
    Compute Z'ΩZ matrix for heteroskedastic case.
*/
void ivvce_compute_ZOmegaZ_robust(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_iv,
    ST_double *ZOmegaZ
)
{
    ST_int i;

    memset(ZOmegaZ, 0, K_iv * K_iv * sizeof(ST_double));

    #pragma omp parallel if(N > 5000)
    {
        ST_double *local_ZOZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        if (local_ZOZ) {
            #pragma omp for schedule(static)
            for (i = 0; i < N; i++) {
                ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
                ST_double e2 = w * resid[i] * resid[i];
                for (ST_int jj = 0; jj < K_iv; jj++) {
                    for (ST_int kk = 0; kk <= jj; kk++) {
                        ST_double contrib = Z[jj * N + i] * e2 * Z[kk * N + i];
                        local_ZOZ[jj * K_iv + kk] += contrib;
                        if (kk != jj) local_ZOZ[kk * K_iv + jj] += contrib;
                    }
                }
            }
            #pragma omp critical
            {
                for (ST_int idx = 0; idx < K_iv * K_iv; idx++)
                    ZOmegaZ[idx] += local_ZOZ[idx];
            }
            free(local_ZOZ);
        }
    }
}

/*
    Compute Z'ΩZ matrix for cluster case.
*/
void ivvce_compute_ZOmegaZ_cluster(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster_ids,
    ST_int N,
    ST_int K_iv,
    ST_int num_clusters,
    ST_double *ZOmegaZ
)
{
    ST_int j, k, c;

    memset(ZOmegaZ, 0, K_iv * K_iv * sizeof(ST_double));

    /* Allocate per-cluster Z'e sums */
    ST_double *cluster_ze = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
    if (!cluster_ze) return;

    /* Sum Z_i * e_i within each cluster */
    #pragma omp parallel if(N > 5000)
    {
        ST_double *local_cze = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
        if (local_cze) {
            #pragma omp for schedule(static)
            for (ST_int ii = 0; ii < N; ii++) {
                ST_int cc = cluster_ids[ii] - 1;
                if (cc < 0 || cc >= num_clusters) continue;
                ST_double w = (weights && weight_type != 0) ? weights[ii] : 1.0;
                ST_double we = w * resid[ii];
                for (ST_int jj = 0; jj < K_iv; jj++) {
                    local_cze[cc * K_iv + jj] += Z[jj * N + ii] * we;
                }
            }
            #pragma omp critical
            {
                for (ST_int idx = 0; idx < num_clusters * K_iv; idx++)
                    cluster_ze[idx] += local_cze[idx];
            }
            free(local_cze);
        }
    }

    /* Compute outer products: Z'ΩZ = sum_c (ze_c)(ze_c)' */
    for (c = 0; c < num_clusters; c++) {
        for (j = 0; j < K_iv; j++) {
            for (k = 0; k <= j; k++) {
                ST_double contrib = cluster_ze[c * K_iv + j] * cluster_ze[c * K_iv + k];
                ZOmegaZ[j * K_iv + k] += contrib;
                if (k != j) ZOmegaZ[k * K_iv + j] += contrib;
            }
        }
    }

    free(cluster_ze);
}

/*
    Helper function to compute one-way cluster meat matrix.
    Returns meat = sum_c (sum_i (PzX_i * e_i))' (sum_i (PzX_i * e_i))
*/
static void compute_cluster_meat(
    const ST_double *PzX,
    const ST_double *resid,
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster_ids,
    ST_int N,
    ST_int K_total,
    ST_int num_clusters,
    ST_double *meat
)
{
    ST_int j, k, c;

    memset(meat, 0, K_total * K_total * sizeof(ST_double));

    /* Allocate per-cluster sums */
    ST_double *cluster_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));
    if (!cluster_sums) return;

    /* Sum (PzX_i * e_i) within each cluster */
    #pragma omp parallel if(N > 5000)
    {
        ST_double *local_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));
        if (local_sums) {
            #pragma omp for schedule(static)
            for (ST_int ii = 0; ii < N; ii++) {
                ST_int cc = cluster_ids[ii] - 1;
                if (cc < 0 || cc >= num_clusters) continue;
                ST_double w = (weights && weight_type != 0) ? weights[ii] : 1.0;
                ST_double we = w * resid[ii];
                for (ST_int jj = 0; jj < K_total; jj++) {
                    local_sums[cc * K_total + jj] += PzX[jj * N + ii] * we;
                }
            }
            #pragma omp critical
            {
                for (ST_int idx = 0; idx < num_clusters * K_total; idx++)
                    cluster_sums[idx] += local_sums[idx];
            }
            free(local_sums);
        }
    }

    /* Compute meat: sum over clusters of outer products */
    for (c = 0; c < num_clusters; c++) {
        for (j = 0; j < K_total; j++) {
            for (k = 0; k <= j; k++) {
                ST_double contrib = cluster_sums[c * K_total + j] * cluster_sums[c * K_total + k];
                meat[j * K_total + k] += contrib;
                if (k != j) meat[k * K_total + j] += contrib;
            }
        }
    }

    free(cluster_sums);
}

/*
    Compute two-way clustered VCE using Cameron-Gelbach-Miller (2011) formula.
    V = V1 + V2 - V_intersection
*/
void ivvce_compute_twoway(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
    ST_int K_total,
    ST_int K_iv,
    const ST_int *cluster1_ids,
    ST_int num_clusters1,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int df_a,
    ST_double *V
)
{
    ST_int i;

    /* Compute P_Z X = Z * (Z'Z)^-1 Z'X */
    ST_double *PzX = (ST_double *)calloc(N * K_total, sizeof(ST_double));
    if (!PzX) return;

    #pragma omp parallel for schedule(static) if(N > 1000)
    for (i = 0; i < N; i++) {
        for (ST_int jj = 0; jj < K_total; jj++) {
            ST_double sum = 0.0;
            for (ST_int kk = 0; kk < K_iv; kk++) {
                sum += Z[kk * N + i] * temp_kiv_ktotal[jj * K_iv + kk];
            }
            PzX[jj * N + i] = sum;
        }
    }

    /* Allocate meat matrices for each clustering */
    ST_double *meat1 = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *meat2 = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *meat_int = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *V1 = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *V2 = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *V_int = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));

    if (!meat1 || !meat2 || !meat_int || !V1 || !V2 || !V_int || !temp_v) {
        free(PzX);
        if (meat1) free(meat1);
        if (meat2) free(meat2);
        if (meat_int) free(meat_int);
        if (V1) free(V1);
        if (V2) free(V2);
        if (V_int) free(V_int);
        if (temp_v) free(temp_v);
        return;
    }

    /* Create intersection cluster IDs */
    /* Each unique (cluster1, cluster2) pair becomes a single cluster */
    ST_int *intersection_ids = (ST_int *)calloc(N, sizeof(ST_int));
    if (!intersection_ids) {
        free(PzX); free(meat1); free(meat2); free(meat_int);
        free(V1); free(V2); free(V_int); free(temp_v);
        return;
    }

    /* Simple approach: map (c1, c2) to unique integer */
    /* Use hash: intersection_id = c1 * max_c2 + c2 */
    /* Then compact to 1..num_intersection */
    /* Check for overflow: (num_clusters1 + 1) * (num_clusters2 + 1) */
    size_t pair_array_size;
    if (ctools_safe_alloc_size((size_t)(num_clusters1 + 1), (size_t)(num_clusters2 + 1),
                               sizeof(ST_int), &pair_array_size) != 0) {
        /* Overflow - fall back to direct intersection computation */
        free(PzX); free(meat1); free(meat2); free(meat_int);
        free(V1); free(V2); free(V_int); free(temp_v);
        free(intersection_ids);
        return;
    }
    ST_int *pair_to_int = (ST_int *)calloc(1, pair_array_size);
    if (!pair_to_int) {
        free(PzX); free(meat1); free(meat2); free(meat_int);
        free(V1); free(V2); free(V_int); free(temp_v);
        free(intersection_ids);
        return;
    }

    ST_int num_intersection = 0;
    for (i = 0; i < N; i++) {
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

    /* Compute meat matrices for each clustering dimension */
    compute_cluster_meat(PzX, resid, weights, weight_type, cluster1_ids, N, K_total, num_clusters1, meat1);
    compute_cluster_meat(PzX, resid, weights, weight_type, cluster2_ids, N, K_total, num_clusters2, meat2);
    compute_cluster_meat(PzX, resid, weights, weight_type, intersection_ids, N, K_total, num_intersection, meat_int);

    /* CGM: combine raw meats first, then apply single DOF adjustment.
       This matches ivreg2's approach:
       1. shat = shat1 + shat2 - shat3 (combine in Z-space, no per-dimension DOF)
       2. V = sandwich(shat_combined) * (N-1)/(N-K) * N_clust/(N_clust-1)
       where N_clust = min(G1, G2). */

    /* meat_combined = meat1 + meat2 - meat_int */
    for (i = 0; i < K_total * K_total; i++) {
        meat1[i] = meat1[i] + meat2[i] - meat_int[i];
    }

    /* V = XkX_inv * meat_combined * XkX_inv */
    ctools_matmul_ab(XkX_inv, meat1, K_total, K_total, K_total, temp_v);
    ctools_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);

    /* Single DOF adjustment: (N_eff-1)/(N_eff-K-df_a) * N_clust/(N_clust-1)
       where N_clust = min(G1, G2), matching ivreg2's small-sample correction */
    ST_int df_r = N_eff - K_total - df_a;
    if (df_r <= 0) df_r = 1;
    ST_double N_clust = (ST_double)((num_clusters1 < num_clusters2) ? num_clusters1 : num_clusters2);
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (ST_double)df_r) * (N_clust / (N_clust - 1.0));
    for (i = 0; i < K_total * K_total; i++) {
        V[i] *= dof_adj;
    }

    /* Clean up */
    free(PzX);
    free(meat1);
    free(meat2);
    free(meat_int);
    free(V1);
    free(V2);
    free(V_int);
    free(temp_v);
    free(intersection_ids);
}

/*
    Compute Kiefer VCE (within-panel autocorrelation robust).

    Kiefer (1980) VCE: homoskedastic HAC with truncated kernel, bandwidth = T.
    Unlike cluster-robust VCE which uses individual Z_i*e_i score vectors,
    Kiefer uses a Kronecker product structure separating the residual
    autocovariance from the instrument cross-products at each lag.

    Following ivreg2's m_omega (homoskedastic + kernel path, lines 194-236):
    1. Base: shat = sigma2 * Z'WZ  where sigma2 = e'We / N_eff
    2. For each lag tau = 1..T-1 (truncated kernel, kw=1):
       sigmahat_tau = sum of w_t * e_t * w_{t-tau} * e_{t-tau} / N_eff
       ZZhat_tau    = sum of w_t * Z_t' * w_{t-tau} * Z_{t-tau}
       ghat = sigmahat_tau * ZZhat_tau   (Kronecker for scalar)
       shat += kw * (ghat + ghat')
    3. shat /= N_eff
    4. Transform to X-space: meat = A' * shat * A
    5. V = XkX_inv * meat * XkX_inv * (N_eff / df_r)
*/
void ivvce_compute_kiefer(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
    ST_int K_total,
    ST_int K_iv,
    const ST_int *panel_ids,
    ST_int num_panels,
    ST_int df_a,
    ST_double *V
)
{
    ST_int i, p;

    ST_int df_r = N_eff - K_total - df_a;
    if (df_r <= 0) df_r = 1;

    /* Build panel membership lists */
    ST_int *panel_counts = (ST_int *)calloc(num_panels, sizeof(ST_int));
    ST_int *panel_starts = (ST_int *)calloc(num_panels + 1, sizeof(ST_int));
    ST_int *obs_by_panel = (ST_int *)calloc(N, sizeof(ST_int));

    if (!panel_counts || !panel_starts || !obs_by_panel) {
        if (panel_counts) free(panel_counts);
        if (panel_starts) free(panel_starts);
        if (obs_by_panel) free(obs_by_panel);
        return;
    }

    for (i = 0; i < N; i++) {
        p = panel_ids[i] - 1;
        if (p >= 0 && p < num_panels) {
            panel_counts[p]++;
        }
    }

    panel_starts[0] = 0;
    for (p = 0; p < num_panels; p++) {
        panel_starts[p + 1] = panel_starts[p] + panel_counts[p];
    }

    memset(panel_counts, 0, num_panels * sizeof(ST_int));
    for (i = 0; i < N; i++) {
        p = panel_ids[i] - 1;
        if (p >= 0 && p < num_panels) {
            obs_by_panel[panel_starts[p] + panel_counts[p]++] = i;
        }
    }

    /* Find T = max panel length */
    ST_int T_max = 0;
    for (p = 0; p < num_panels; p++) {
        ST_int T_p = panel_starts[p + 1] - panel_starts[p];
        if (T_p > T_max) T_max = T_p;
    }

    if (T_max == 0) {
        free(panel_counts); free(panel_starts); free(obs_by_panel);
        return;
    }

    /* === Homoskedastic HAC (Kiefer) ===
       Base term: shat_ZZ = sigma2 * Z'WZ
       sigma2 = e'We / N_eff  (large-sample, dofminus=0 matching ivreghdfe->ivreg2) */

    /* Compute sigma2 = e'We / N_eff */
    ST_double ee = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ee) if(N > 10000)
    for (i = 0; i < N; i++) {
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        ee += w * resid[i] * resid[i];
    }
    ST_double sigma2 = ee / (ST_double)N_eff;

    /* Compute shat_ZZ = sigma2 * Z'WZ */
    ST_double *shat_ZZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    if (!shat_ZZ) {
        free(panel_counts); free(panel_starts); free(obs_by_panel);
        return;
    }

    if (weights && weight_type != 0)
        ctools_matmul_atdb(Z, Z, weights, N, K_iv, K_iv, shat_ZZ);
    else
        ctools_matmul_atb(Z, Z, N, K_iv, K_iv, shat_ZZ);

    for (i = 0; i < K_iv * K_iv; i++)
        shat_ZZ[i] *= sigma2;

    /* Cross-lag terms: for tau=1..T_max-1, truncated kernel (kw=1 for all lags)
       For each lag tau, accumulate within-panel (t, t-tau) pairs:
         sigmahat_tau = Σ_{panels} Σ_{t} wv * e[t] * e[t-tau] / N_eff
         ZZhat_tau    = Σ_{panels} Σ_{t} wv * Z[t]' Z[t-tau]
       Then: shat_ZZ += ghat + ghat'  where ghat = sigmahat_tau * ZZhat_tau */

    ST_double *ZZhat_tau = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
    if (!ZZhat_tau) {
        free(shat_ZZ);
        free(panel_counts); free(panel_starts); free(obs_by_panel);
        return;
    }

    for (ST_int tau = 1; tau < T_max; tau++) {
        ST_double sigmahat_tau = 0.0;
        memset(ZZhat_tau, 0, K_iv * K_iv * sizeof(ST_double));

        /* Panels are independent — parallelize with thread-local reduction */
        #pragma omp parallel if(num_panels > 4)
        {
            ST_double local_sigma = 0.0;
            ST_double *local_ZZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
            if (local_ZZ) {
                #pragma omp for schedule(dynamic)
                for (p = 0; p < num_panels; p++) {
                    ST_int start = panel_starts[p];
                    ST_int T_p = panel_starts[p + 1] - start;

                    for (ST_int t = tau; t < T_p; t++) {
                        ST_int obs_t     = obs_by_panel[start + t];
                        ST_int obs_t_tau = obs_by_panel[start + t - tau];

                        ST_double w_t     = (weights && weight_type != 0) ? weights[obs_t] : 1.0;
                        ST_double w_t_tau = (weights && weight_type != 0) ? weights[obs_t_tau] : 1.0;
                        ST_double wv = w_t * w_t_tau;

                        local_sigma += wv * resid[obs_t] * resid[obs_t_tau];

                        for (ST_int l = 0; l < K_iv; l++) {
                            for (ST_int m = 0; m < K_iv; m++) {
                                local_ZZ[l * K_iv + m] +=
                                    wv * Z[l * N + obs_t] * Z[m * N + obs_t_tau];
                            }
                        }
                    }
                }
                #pragma omp critical
                {
                    sigmahat_tau += local_sigma;
                    for (ST_int idx = 0; idx < K_iv * K_iv; idx++)
                        ZZhat_tau[idx] += local_ZZ[idx];
                }
                free(local_ZZ);
            }
        }

        sigmahat_tau /= (ST_double)N_eff;  /* large-sample dof correction */

        /* ghat = sigmahat_tau * ZZhat_tau (Kronecker for scalar)
           shat_ZZ += ghat + ghat' */
        for (ST_int l = 0; l < K_iv; l++) {
            for (ST_int m = 0; m < K_iv; m++) {
                shat_ZZ[l * K_iv + m] +=
                    sigmahat_tau * (ZZhat_tau[l * K_iv + m] + ZZhat_tau[m * K_iv + l]);
            }
        }
    }

    free(ZZhat_tau);
    free(panel_counts); free(panel_starts); free(obs_by_panel);

    /* Note: shat_ZZ is NOT divided by N here.  In ivreg2, m_omega divides by N
       (line 235) but s_iegmm multiplies by N in the sandwich (line 5557).
       Our convention matches the robust/HAC path in ivvce_compute_full:
       shat_ZZ stores the raw (N * omega) and dof_adj = N/df_r handles scaling. */

    /* Transform shat_ZZ to X-space:
       meat = A' * shat_ZZ * A where A = temp_kiv_ktotal = (Z'Z)^{-1} * Z'X */
    ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *temp_kk = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));

    if (!meat || !temp_kk) {
        free(shat_ZZ);
        if (meat) free(meat);
        if (temp_kk) free(temp_kk);
        return;
    }

    ctools_matmul_ab(shat_ZZ, temp_kiv_ktotal, K_iv, K_iv, K_total, temp_kk);
    ctools_matmul_atb(temp_kiv_ktotal, temp_kk, K_iv, K_total, K_total, meat);

    free(shat_ZZ);
    free(temp_kk);

    /* Small-sample DOF adjustment: V *= (N_eff - dofminus) / (N_eff - K - dofminus - sdofminus)
       With dofminus=0 (ivreghdfe convention): V *= N_eff / (N_eff - K_total - df_a) */
    ST_double dof_adj = (ST_double)N_eff / (ST_double)df_r;

    /* V = XkX_inv * meat * XkX_inv * dof_adj */
    ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (temp_v) {
        ctools_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
        ctools_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
        for (i = 0; i < K_total * K_total; i++) {
            V[i] *= dof_adj;
        }
        free(temp_v);
    }

    free(meat);
}

/*
    Full VCE computation with P_Z X calculation.

    This is the main entry point for VCE computation when Z is available.

    Parameters:
    - Z: Instruments (N x K_iv)
    - resid: Residuals (N x 1)
    - temp_kiv_ktotal: (Z'Z)^-1 Z'X (K_iv x K_total)
    - XkX_inv: Inverse of k-class matrix (K_total x K_total)
    - weights, weight_type: Weighting
    - N, K_total, K_iv: Dimensions
    - vce_type: CIVREGHDFE_VCE_* constant
    - cluster_ids, num_clusters: Clustering info
    - df_a, nested_adj: DOF adjustments
    - kernel_type, bw: HAC parameters
    - V: Output VCE matrix (K_total x K_total)
*/
void ivvce_compute_full(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,  /* Effective sample size: sum(fw) for fweights, N otherwise */
    ST_int K_total,
    ST_int K_iv,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_int kernel_type,
    ST_int bw,
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels,
    ST_double *V
)
{
    ST_int i, j, k;

    /* Compute residual sum of squares */
    ST_double rss = 0.0;
    if (weights && weight_type != 0) {
        #pragma omp parallel for schedule(static) reduction(+:rss) if(N > 10000)
        for (i = 0; i < N; i++) {
            rss += weights[i] * resid[i] * resid[i];
        }
    } else {
        #pragma omp parallel for schedule(static) reduction(+:rss) if(N > 10000)
        for (i = 0; i < N; i++) {
            rss += resid[i] * resid[i];
        }
    }

    ST_int df_r = N_eff - K_total - df_a;
    if (df_r <= 0) df_r = 1;

    if (vce_type == CIVREGHDFE_VCE_UNADJUSTED) {
        /* Unadjusted: V = sigma^2 * XkX_inv */
        ST_double sigma2 = rss / df_r;
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = sigma2 * XkX_inv[i];
        }
        return;
    }

    /* Compute P_Z X = Z * (Z'Z)^-1 Z'X = Z * temp_kiv_ktotal */
    ST_double *PzX = (ST_double *)calloc(N * K_total, sizeof(ST_double));
    if (!PzX) return;

    #pragma omp parallel for schedule(static) if(N > 1000)
    for (i = 0; i < N; i++) {
        for (ST_int jj = 0; jj < K_total; jj++) {
            ST_double sum = 0.0;
            for (ST_int kk = 0; kk < K_iv; kk++) {
                sum += Z[kk * N + i] * temp_kiv_ktotal[jj * K_iv + kk];
            }
            PzX[jj * N + i] = sum;
        }
    }

    /* Compute meat matrix based on VCE type */
    ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!meat) {
        free(PzX);
        return;
    }

    if (vce_type == CIVREGHDFE_VCE_DKRAAY && cluster_ids && num_clusters > 0 && kernel_type > 0 && bw > 0) {
        /*
           Driscoll-Kraay VCE: Cross-sectional averaging + HAC on time series
           1. Compute sum of moment contributions at each time cluster:
              u_t = sum_{i in time t} (PzX_i * e_i)
           2. Apply HAC kernel to the time series of sums
        */
        ST_double *cluster_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));
        if (!cluster_sums) {
            free(PzX); free(meat);
            return;
        }

        /* Sum (PzX_i * e_i) at each time period (cluster) */
        #pragma omp parallel if(N > 5000)
        {
            ST_double *local_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));
            if (local_sums) {
                #pragma omp for schedule(static)
                for (ST_int ii = 0; ii < N; ii++) {
                    ST_int t = cluster_ids[ii] - 1;
                    if (t < 0 || t >= num_clusters) continue;
                    ST_double w = (weights && weight_type != 0) ? weights[ii] : 1.0;
                    ST_double we = w * resid[ii];
                    for (ST_int jj = 0; jj < K_total; jj++) {
                        local_sums[t * K_total + jj] += PzX[jj * N + ii] * we;
                    }
                }
                #pragma omp critical
                {
                    for (ST_int idx = 0; idx < num_clusters * K_total; idx++)
                        cluster_sums[idx] += local_sums[idx];
                }
                free(local_sums);
            }
        }

        /* Apply HAC kernel to time series of cluster sums
           For spectral windows (QS=3, Danielle=6, Tent=7), loop through all lags up to T-1.
           For lag windows (Bartlett=1, Parzen=2, Truncated=4, Tukey-Hanning=5), loop through bw.
           Following ivreg2's m_omega approach. */
        ST_int is_spectral = (kernel_type == 3);  /* QS is spectral */
        ST_int max_lag = is_spectral ? (num_clusters - 1) : bw;
        if (max_lag >= num_clusters) max_lag = num_clusters - 1;

        for (ST_int lag = 0; lag <= max_lag; lag++) {
            ST_double kw = civreghdfe_kernel_weight(kernel_type, lag, bw);
            /* Skip if weight is essentially zero. Note: spectral kernels can have
               negative weights, so we use fabs() to check magnitude, not value. */
            if (fabs(kw) < 1e-10) continue;

            if (lag == 0) {
                /* Diagonal: sum_t u_t * u_t' */
                for (ST_int t = 0; t < num_clusters; t++) {
                    for (j = 0; j < K_total; j++) {
                        for (k = 0; k <= j; k++) {
                            ST_double contrib = kw * cluster_sums[t * K_total + j] * cluster_sums[t * K_total + k];
                            meat[j * K_total + k] += contrib;
                            if (k != j) meat[k * K_total + j] += contrib;
                        }
                    }
                }
            } else {
                /* Off-diagonal: sum_t u_t * u_{t+lag}' + u_{t+lag} * u_t' */
                for (ST_int t = 0; t < num_clusters - lag; t++) {
                    for (j = 0; j < K_total; j++) {
                        for (k = 0; k < K_total; k++) {
                            ST_double contrib = kw * (cluster_sums[t * K_total + j] * cluster_sums[(t + lag) * K_total + k] +
                                                      cluster_sums[(t + lag) * K_total + j] * cluster_sums[t * K_total + k]);
                            meat[k * K_total + j] += contrib;
                        }
                    }
                }
            }
        }

        free(cluster_sums);

        /* Driscoll-Kraay DOF adjustment: (N-1)/(N-K-df_a) * T/(T-1)
           This matches ivreg2's cluster-robust small-sample adjustment.
           The (N-1)/(N-K-df_a) factor accounts for partialled-out fixed effects. */
        ST_double T = (ST_double)num_clusters;
        ST_int effective_df_a = df_a;
        if (df_a == 0 && nested_adj == 1) effective_df_a = 1;
        ST_double denom = (ST_double)(N_eff - K_total - effective_df_a);
        if (denom <= 0) denom = 1.0;
        ST_double dof_adj = ((ST_double)(N_eff - 1) / denom) * (T / (T - 1.0));

        /* V = XkX_inv * meat * XkX_inv * dof_adj */
        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (temp_v) {
            ctools_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
            ctools_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
            for (i = 0; i < K_total * K_total; i++) {
                V[i] *= dof_adj;
            }
            free(temp_v);
        }

    } else if (vce_type == CIVREGHDFE_VCE_CLUSTER && cluster_ids && num_clusters > 0) {
        /* Clustered VCE — delegate to shared VCE engine */
        /* Convert 1-based cluster IDs to 0-based for shared function */
        ST_int *cluster_ids_0 = (ST_int *)malloc(N * sizeof(ST_int));
        if (!cluster_ids_0) {
            free(PzX); free(meat);
            return;
        }
        for (i = 0; i < N; i++) {
            cluster_ids_0[i] = cluster_ids[i] - 1;
        }

        /* DOF adjustment */
        ST_int effective_df_a = df_a;
        if (df_a == 0 && nested_adj == 1) effective_df_a = 1;
        ST_double denom = (ST_double)(N_eff - K_total - effective_df_a);
        if (denom <= 0) denom = 1.0;
        ST_double G = (ST_double)num_clusters;
        ST_double dof_adj = ((ST_double)(N_eff - 1) / denom) * (G / (G - 1.0));

        ctools_vce_data d;
        d.X_eff = PzX;
        d.D = XkX_inv;
        d.resid = resid;
        d.weights = weights;
        d.weight_type = weight_type;
        d.N = N;
        d.K = K_total;
        d.normalize_weights = 1;  /* Normalize aw/pw weights for correct meat */

        ctools_vce_cluster(&d, cluster_ids_0, num_clusters, dof_adj, V);
        free(cluster_ids_0);

    } else if (kernel_type > 0 && bw > 0) {
        /* HAC VCE - Compute in Z-space like ivreg2, then transform to X-space
           Following livreg2.do m_omega() exactly:
           1. Score vectors in Z-space: uZ[l,i] = Z[l,i] * e[i]
           2. Compute shat_ZZ using HAC formula with kernel weights
           3. Transform: meat_X = A' * shat_ZZ * A where A = (Z'Z)^{-1} * Z'X
           4. V = (X'PzX)^{-1} * meat_X * (X'PzX)^{-1} * (N / df_r)
         */

        /* Allocate shat_ZZ (K_iv x K_iv) and score vectors */
        ST_double *shat_ZZ = (ST_double *)calloc(K_iv * K_iv, sizeof(ST_double));
        ST_double *uZ = (ST_double *)calloc(N * K_iv, sizeof(ST_double));
        if (!shat_ZZ || !uZ) {
            if (shat_ZZ) free(shat_ZZ);
            if (uZ) free(uZ);
            free(PzX); free(meat);
            return;
        }

        /* Compute Z-space score vectors: uZ[l,i] = Z[l,i] * e[i] */
        #pragma omp parallel for schedule(static) if(N > 5000)
        for (i = 0; i < N; i++) {
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            ST_double we = w * resid[i];
            for (ST_int l = 0; l < K_iv; l++) {
                uZ[l * N + i] = Z[l * N + i] * we;
            }
        }

        /* Compute HAC meat in Z-space: shat_ZZ = Σ kw(|t-s|) * uZ_t * uZ_s' */
        if (hac_panel_ids && num_hac_panels > 0) {
            /* Panel-aware HAC in Z-space */
            ST_int *panel_counts = (ST_int *)calloc(num_hac_panels, sizeof(ST_int));
            ST_int *panel_starts = (ST_int *)calloc(num_hac_panels + 1, sizeof(ST_int));
            ST_int *obs_by_panel = (ST_int *)calloc(N, sizeof(ST_int));

            if (!panel_counts || !panel_starts || !obs_by_panel) {
                free(uZ); free(shat_ZZ);
                if (panel_counts) free(panel_counts);
                if (panel_starts) free(panel_starts);
                if (obs_by_panel) free(obs_by_panel);
                free(PzX); free(meat);
                return;
            }

            for (i = 0; i < N; i++) {
                ST_int p = hac_panel_ids[i] - 1;
                if (p >= 0 && p < num_hac_panels) panel_counts[p]++;
            }
            panel_starts[0] = 0;
            for (ST_int p = 0; p < num_hac_panels; p++) {
                panel_starts[p + 1] = panel_starts[p] + panel_counts[p];
            }
            memset(panel_counts, 0, num_hac_panels * sizeof(ST_int));
            for (i = 0; i < N; i++) {
                ST_int p = hac_panel_ids[i] - 1;
                if (p >= 0 && p < num_hac_panels) {
                    obs_by_panel[panel_starts[p] + panel_counts[p]++] = i;
                }
            }

            for (ST_int p = 0; p < num_hac_panels; p++) {
                ST_int start = panel_starts[p];
                ST_int end = panel_starts[p + 1];
                ST_int T_p = end - start;

                for (ST_int t1 = 0; t1 < T_p; t1++) {
                    ST_int i1 = obs_by_panel[start + t1];
                    for (ST_int t2 = 0; t2 < T_p; t2++) {
                        ST_int i2 = obs_by_panel[start + t2];
                        ST_int time_lag = (t1 > t2) ? (t1 - t2) : (t2 - t1);
                        ST_double kw = civreghdfe_kernel_weight(kernel_type, time_lag, bw);
                        if (fabs(kw) < 1e-10) continue;  /* Use fabs for spectral kernels with negative weights */

                        /* Use full (l, m) loop without premature symmetrization.
                           Each (i1, i2) contribution uZ[l,i1]*uZ[m,i2] is NOT symmetric in l,m.
                           The overall shat_ZZ will be symmetric because we sum over all (t1,t2)
                           pairs, and Gamma(tau) + Gamma(-tau) is symmetric. */
                        for (ST_int l = 0; l < K_iv; l++) {
                            for (ST_int m = 0; m < K_iv; m++) {
                                ST_double contrib = kw * uZ[l * N + i1] * uZ[m * N + i2];
                                shat_ZZ[l * K_iv + m] += contrib;
                            }
                        }
                    }
                }
            }
            free(panel_counts); free(panel_starts); free(obs_by_panel);

        } else {
            /* Standard (non-panel) HAC in Z-space
               For spectral windows (QS), loop through all lags up to N-1.
               For lag windows, loop through bw. */
            ST_int is_spectral = (kernel_type == 3);  /* QS is spectral */
            ST_int max_lag = is_spectral ? (N - 1) : bw;
            if (max_lag >= N) max_lag = N - 1;

            for (ST_int lag = 0; lag <= max_lag; lag++) {
                ST_double kw = civreghdfe_kernel_weight(kernel_type, lag, bw);
                if (fabs(kw) < 1e-10) continue;

                if (lag == 0) {
                    /* Gamma(0) = Σ_i uZ_i * uZ_i' */
                    for (i = 0; i < N; i++) {
                        for (ST_int l = 0; l < K_iv; l++) {
                            for (ST_int m = 0; m <= l; m++) {
                                ST_double contrib = kw * uZ[l * N + i] * uZ[m * N + i];
                                shat_ZZ[l * K_iv + m] += contrib;
                                if (m != l) shat_ZZ[m * K_iv + l] += contrib;
                            }
                        }
                    }
                } else {
                    /* Gamma(lag) + Gamma(lag)': add contributions for both directions */
                    for (i = 0; i < N - lag; i++) {
                        for (ST_int l = 0; l < K_iv; l++) {
                            for (ST_int m = 0; m < K_iv; m++) {
                                /* ghat[l,m] = uZ[l,i] * uZ[m,i+lag]
                                   (ghat + ghat')[l,m] = uZ[l,i]*uZ[m,i+lag] + uZ[l,i+lag]*uZ[m,i] */
                                ST_double contrib = kw * (uZ[l * N + i] * uZ[m * N + (i + lag)] +
                                                          uZ[l * N + (i + lag)] * uZ[m * N + i]);
                                shat_ZZ[m * K_iv + l] += contrib;
                            }
                        }
                    }
                }
            }
        }

        free(uZ);

        /* Transform shat_ZZ to X-space:
           meat_X = A' * shat_ZZ * A where A = temp_kiv_ktotal = (Z'Z)^{-1} * Z'X
           This is equivalent to: (Z'X)' * (Z'Z)^{-1} * shat_ZZ * (Z'Z)^{-1} * (Z'X)
         */
        ST_double *temp_kk = (ST_double *)calloc(K_iv * K_total, sizeof(ST_double));
        if (temp_kk) {
            /* temp_kk = shat_ZZ * A  (K_iv x K_total) */
            ctools_matmul_ab(shat_ZZ, temp_kiv_ktotal, K_iv, K_iv, K_total, temp_kk);
            /* meat = A' * temp_kk  (K_total x K_total) */
            ctools_matmul_atb(temp_kiv_ktotal, temp_kk, K_iv, K_total, K_total, meat);
            free(temp_kk);
        }
        free(shat_ZZ);

        /* HAC VCE sandwich: V = XkX_inv * meat * XkX_inv * dof_adj
           ivreghdfe applies small-sample adjustment: V = V * N / (N - K - G)
           where K = number of regressors and G = absorbed FE count.
           This is equivalent to V = V * N / df_r.
         */
        ST_double dof_adj = (ST_double)N_eff / (ST_double)df_r;

        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (temp_v) {
            ctools_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
            ctools_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
            /* Small-sample adjustment */
            for (i = 0; i < K_total * K_total; i++) {
                V[i] *= dof_adj;
            }
            free(temp_v);
        }

    } else {
        /* Standard HC robust — delegate to shared VCE engine */
        ST_double dof_adj = (ST_double)N_eff / (ST_double)df_r;

        ctools_vce_data d;
        d.X_eff = PzX;
        d.D = XkX_inv;
        d.resid = resid;
        d.weights = weights;
        d.weight_type = weight_type;
        d.N = N;
        d.K = K_total;
        d.normalize_weights = 1;  /* Normalize aw/pw weights for correct meat (w^2*e^2) */

        ctools_vce_robust(&d, dof_adj, V);
    }

    free(PzX);
    free(meat);
}
