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

/*
    Sandwich VCE from a Z-space moment covariance S (K_iv x K_iv, raw, from
    civ_omega_build): V = XkX_inv (A'SA) XkX_inv * dof_adj with
    A = temp_kiv_ktotal = (Z'Z)^-1 Z'X, ivreg2's s_iegmm / s_liml sandwich.
    Used for the kernel VCEs (HAC, AC, Driscoll-Kraay, Kiefer).
*/
ST_retcode ivvce_compute_sandwich(
    const ST_double *temp_kiv_ktotal,
    const ST_double *S,
    const ST_double *XkX_inv,
    ST_int K_total,
    ST_int K_iv,
    ST_double dof_adj,
    ST_double *V
)
{
    ST_double *SA = (ST_double *)calloc((size_t)K_iv * K_total, sizeof(ST_double));
    ST_double *meat = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
    ST_double *temp_v = (ST_double *)calloc((size_t)K_total * K_total, sizeof(ST_double));
    if (!SA || !meat || !temp_v) {
        free(SA); free(meat); free(temp_v);
        return 920;
    }

    /* meat = A' S A */
    ctools_matmul_ab(S, temp_kiv_ktotal, K_iv, K_iv, K_total, SA);
    ctools_matmul_atb(temp_kiv_ktotal, SA, K_iv, K_total, K_total, meat);

    /* V = XkX_inv * meat * XkX_inv * dof_adj */
    ctools_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
    ctools_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
    for (ST_int i = 0; i < K_total * K_total; i++) V[i] *= dof_adj;

    free(SA); free(meat); free(temp_v);
    return 0;
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
    ST_double *cluster_sums = (ST_double *)calloc((size_t)num_clusters * K_total, sizeof(ST_double));
    if (!cluster_sums) return;

    /* Sum (PzX_i * e_i) within each cluster */
    #pragma omp parallel if(N > 5000)
    {
        ST_double *local_sums = (ST_double *)calloc((size_t)num_clusters * K_total, sizeof(ST_double));
        if (local_sums) {
            #pragma omp for schedule(static)
            for (ST_int ii = 0; ii < N; ii++) {
                ST_int cc = cluster_ids[ii] - 1;
                if (cc < 0 || cc >= num_clusters) continue;
                ST_double w = (weights && weight_type != 0) ? weights[ii] : 1.0;
                ST_double we = w * resid[ii];
                for (ST_int jj = 0; jj < K_total; jj++) {
                    local_sums[(size_t)cc * K_total + jj] += PzX[(size_t)jj * N + ii] * we;
                }
            }
            #pragma omp critical
            {
                for (size_t idx = 0; idx < (size_t)num_clusters * K_total; idx++)
                    cluster_sums[idx] += local_sums[idx];
            }
            free(local_sums);
        }
    }

    /* Compute meat: sum over clusters of outer products */
    for (c = 0; c < num_clusters; c++) {
        for (j = 0; j < K_total; j++) {
            for (k = 0; k <= j; k++) {
                ST_double contrib = cluster_sums[(size_t)c * K_total + j] * cluster_sums[(size_t)c * K_total + k];
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
    ST_double *PzX = (ST_double *)calloc((size_t)N * K_total, sizeof(ST_double));
    if (!PzX) return;

    #pragma omp parallel for schedule(static) if(N > 1000)
    for (i = 0; i < N; i++) {
        for (ST_int jj = 0; jj < K_total; jj++) {
            ST_double sum = 0.0;
            for (ST_int kk = 0; kk < K_iv; kk++) {
                sum += Z[(size_t)kk * N + i] * temp_kiv_ktotal[jj * K_iv + kk];
            }
            PzX[(size_t)jj * N + i] = sum;
        }
    }

    /* Allocate meat matrices for each clustering */
    ST_double *meat1 = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *meat2 = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *meat_int = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));

    if (!meat1 || !meat2 || !meat_int || !temp_v) {
        free(PzX);
        if (meat1) free(meat1);
        if (meat2) free(meat2);
        if (meat_int) free(meat_int);
        if (temp_v) free(temp_v);
        return;
    }

    /* Create intersection cluster IDs */
    /* Each unique (cluster1, cluster2) pair becomes a single cluster */
    ST_int *intersection_ids = (ST_int *)calloc(N, sizeof(ST_int));
    if (!intersection_ids) {
        free(PzX); free(meat1); free(meat2); free(meat_int);
        free(temp_v);
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
        free(temp_v);
        free(intersection_ids);
        return;
    }
    ST_int *pair_to_int = (ST_int *)calloc(1, pair_array_size);
    if (!pair_to_int) {
        free(PzX); free(meat1); free(meat2); free(meat_int);
        free(temp_v);
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
    free(temp_v);
    free(intersection_ids);
}

/*
    Full VCE computation with P_Z X calculation.

    This is the main entry point for the unadjusted, robust (HC) and
    one-way cluster VCEs; the kernel VCEs use ivvce_compute_sandwich.

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
    - V: Output VCE matrix (K_total x K_total)
*/
ST_retcode ivvce_compute_full(
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
    ST_double *V
)
{
    ST_int i;
    ST_retcode rc = 0;

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
        return 0;
    }

    /* Compute P_Z X = Z * (Z'Z)^-1 Z'X = Z * temp_kiv_ktotal */
    ST_double *PzX = (ST_double *)calloc((size_t)N * K_total, sizeof(ST_double));
    if (!PzX) return 920;

    #pragma omp parallel for schedule(static) if(N > 1000)
    for (i = 0; i < N; i++) {
        for (ST_int jj = 0; jj < K_total; jj++) {
            ST_double sum = 0.0;
            for (ST_int kk = 0; kk < K_iv; kk++) {
                sum += Z[(size_t)kk * N + i] * temp_kiv_ktotal[jj * K_iv + kk];
            }
            PzX[(size_t)jj * N + i] = sum;
        }
    }

    if (vce_type == CIVREGHDFE_VCE_CLUSTER && cluster_ids && num_clusters > 0) {
        /* Clustered VCE — delegate to shared VCE engine */
        /* Convert 1-based cluster IDs to 0-based for shared function */
        ST_int *cluster_ids_0 = (ST_int *)malloc(N * sizeof(ST_int));
        if (!cluster_ids_0) {
            free(PzX);
            return 920;
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

        rc = ctools_vce_cluster(&d, cluster_ids_0, num_clusters, dof_adj, V);
        free(cluster_ids_0);

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

        rc = ctools_vce_robust(&d, dof_adj, V);
    }

    free(PzX);
    return rc;
}
