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
#include "../ctools_hdfe_utils.h"
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
    P_Z X = Z (Z'Z)^-1 Z'X = Z * temp_kiv_ktotal, N x K_total column-major.
    Returns NULL on allocation failure; the caller frees the result.
*/
static ST_double *ivvce_project_x(const ST_double *Z, const ST_double *temp_kiv_ktotal,
                                  ST_int N, ST_int K_total, ST_int K_iv)
{
    ST_double *PzX = (ST_double *)ctools_safe_malloc3((size_t)N, (size_t)K_total, sizeof(ST_double));
    if (!PzX) return NULL;

    /* Each element sums Z[i,k] * temp[k,j] over k in order */
    ctools_matmul_xb(Z, N, K_iv, temp_kiv_ktotal, K_total, PzX);
    return PzX;
}

/*
    Compute two-way clustered VCE using Cameron-Gelbach-Miller (2011) formula.
    V = V1 + V2 - V_intersection
    Returns 0, 920 on allocation failure or the shared engine's error code;
    V is written only on success.
*/
ST_retcode ivvce_compute_twoway(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights (may exceed INT_MAX), N otherwise */
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
    ST_retcode rc = 0;
    const size_t KK = (size_t)K_total * K_total;
    ST_double *PzX = ivvce_project_x(Z, temp_kiv_ktotal, N, K_total, K_iv);
    ST_double *V_sum = (ST_double *)calloc(KK, sizeof(ST_double));
    ST_double *V_term = (ST_double *)calloc(KK, sizeof(ST_double));
    /* 0-based ids for the shared engine: one dimension, then the pairs */
    ST_int *ids = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
    ST_int *pair_ids = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
    ST_int num_pairs = 0;

    if (!PzX || !V_sum || !V_term || !ids || !pair_ids ||
        ctools_intersect_cluster_ids(cluster1_ids, cluster2_ids, N, pair_ids, &num_pairs) != 0) {
        rc = 920;
        goto twoway_done;
    }

    /* As ivreg2: the three terms share one small-sample factor,
       (N_eff-1)/(N_eff-K-df_a) * G/(G-1) with G = min(G1, G2), which equals
       combining the raw meats first (shat1 + shat2 - shat12) */
    ST_double df_r = N_eff - K_total - df_a;
    if (df_r <= 0) df_r = 1;
    ST_double N_clust = (ST_double)((num_clusters1 < num_clusters2) ? num_clusters1 : num_clusters2);
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (ST_double)df_r) * (N_clust / (N_clust - 1.0));

    ctools_vce_data d;
    d.X_eff = PzX;
    d.D = XkX_inv;
    d.resid = resid;
    d.weights = (weights && weight_type != 0) ? weights : NULL;
    d.weight_type = weight_type;
    d.N = N;
    d.K = K_total;
    d.normalize_weights = 0;  /* scores are w * e with the estimator's weights */

    for (ST_int term = 0; term < 3 && !rc; term++) {
        const ST_int *term_ids = pair_ids;
        ST_int term_clusters = num_pairs;
        if (term < 2) {
            const ST_int *src = term == 0 ? cluster1_ids : cluster2_ids;
            for (i = 0; i < N; i++) ids[i] = src[i] - 1;  /* 1-based on input */
            term_ids = ids;
            term_clusters = term == 0 ? num_clusters1 : num_clusters2;
        }
        rc = ctools_vce_cluster(&d, term_ids, term_clusters, dof_adj, V_term);
        ST_double sign = term < 2 ? 1.0 : -1.0;
        for (size_t e = 0; e < KK && !rc; e++) V_sum[e] += sign * V_term[e];
    }
    if (!rc) memcpy(V, V_sum, KK * sizeof(ST_double));

twoway_done:
    free(PzX);
    free(V_sum);
    free(V_term);
    free(ids);
    free(pair_ids);
    return rc;
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
    ST_double N_eff,  /* Effective sample size: sum(fw) for fweights (may exceed INT_MAX), N otherwise */
    ST_int K_total,
    ST_int K_iv,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_int dofminus,
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

    /* ivreg2: sigma^2 = e'e/(N-dofminus), V scaled by
       (N-dofminus)/(N-K-dofminus-sdofminus) without clusters */
    ST_double df_r = N_eff - K_total - df_a - dofminus;
    if (df_r <= 0) df_r = 1;

    if (vce_type == CIVREGHDFE_VCE_UNADJUSTED) {
        /* Unadjusted: V = sigma^2 * XkX_inv */
        ST_double sigma2 = rss / df_r;
        for (i = 0; i < K_total * K_total; i++) {
            V[i] = sigma2 * XkX_inv[i];
        }
        return 0;
    }

    ST_double *PzX = ivvce_project_x(Z, temp_kiv_ktotal, N, K_total, K_iv);
    if (!PzX) return 920;

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
        ST_double dof_adj = (ST_double)(N_eff - dofminus) / (ST_double)df_r;

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
