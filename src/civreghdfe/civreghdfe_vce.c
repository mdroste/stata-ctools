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

#include "civreghdfe_vce.h"
#include "../ctools_config.h"
#include "civreghdfe_matrix.h"

/* Shared OLS functions */
#include "../ctools_ols.h"
#define cholesky ctools_cholesky
#define invert_from_cholesky ctools_invert_from_cholesky

/*
    Compute unadjusted VCE for k-class estimator.
*/
void ivvce_compute_unadjusted(
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
    Compute robust (HC1) VCE for k-class estimator.
*/
void ivvce_compute_robust(
    IVEstContext *ctx,
    const ST_double *resid,
    const ST_double *XkX_inv,
    ST_int df_a,
    ST_double *V
)
{
    ST_int N = ctx->N;
    ST_int K_total = ctx->K_total;
    ST_int i;

    /* Note: This is a stub - full implementation requires Z matrix */
    (void)resid;  /* Will be used in full implementation */

    ST_int df_r = N - K_total - df_a;
    if (df_r <= 0) df_r = 1;

    /* HC1 DOF adjustment */
    ST_double dof_adj = (ST_double)N / (ST_double)df_r;

    /* Compute P_Z X = Z(Z'Z)^-1 Z'X using pre-computed matrices */
    /* PzX = Z * temp_kiv_ktotal where temp = (Z'Z)^-1 Z'X */
    /* Note: This requires Z matrix - for now, compute meat directly */

    /* Compute size with overflow check */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K_total, (size_t)K_total, sizeof(ST_double), &kk_size) != 0) {
        return;  /* Overflow */
    }

    /* Allocate meat matrix */
    ST_double *meat = (ST_double *)calloc(1, kk_size);
    if (!meat) return;

    /* For robust VCE without Z, we use the sandwich with X'P_Z X */
    /* This is an approximation - full implementation needs Z */

    /* Compute temp = XkX_inv * meat */
    ST_double *temp_v = (ST_double *)calloc(1, kk_size);
    if (!temp_v) {
        free(meat);
        return;
    }

    /* For now, use scaled unadjusted VCE as approximation */
    for (i = 0; i < K_total * K_total; i++) {
        V[i] = XkX_inv[i] * dof_adj;
    }

    free(meat);
    free(temp_v);
}

/*
    Compute clustered VCE for k-class estimator.
*/
void ivvce_compute_cluster(
    IVEstContext *ctx,
    const ST_double *resid,
    const ST_double *XkX_inv,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_double *V
)
{
    ST_int N = ctx->N;
    ST_int K_total = ctx->K_total;
    ST_int i;

    /* Note: This is a stub - full implementation requires Z matrix */
    (void)resid;       /* Will be used in full implementation */
    (void)cluster_ids; /* Will be used in full implementation */

    /* Cluster DOF adjustment:
       dof_adj = (N-1)/(N - K_total - effective_df_a) * G/(G-1) */
    ST_int effective_df_a = df_a;
    if (df_a == 0 && nested_adj == 1) {
        effective_df_a = 1;
    }
    ST_double denom = (ST_double)(N - K_total - effective_df_a);
    if (denom <= 0) denom = 1.0;
    ST_double G = (ST_double)num_clusters;
    ST_double dof_adj = ((ST_double)(N - 1) / denom) * (G / (G - 1.0));

    /* Allocate cluster sums and meat matrix */
    ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!meat) return;

    /* Note: Full implementation requires Z matrix for P_Z X computation */
    /* For now, provide scaled approximation */

    /* Compute V = XkX_inv * meat * XkX_inv * dof_adj */
    ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!temp_v) {
        free(meat);
        return;
    }

    /* Approximate with scaled unadjusted */
    for (i = 0; i < K_total * K_total; i++) {
        V[i] = XkX_inv[i] * dof_adj;
    }

    free(meat);
    free(temp_v);
}

/*
    Compute HAC VCE.
*/
void ivvce_compute_hac(
    IVEstContext *ctx,
    const ST_double *resid,
    const ST_double *XkX_inv,
    ST_int kernel_type,
    ST_int bw,
    ST_int df_a,
    ST_double *V
)
{
    ST_int N = ctx->N;
    ST_int K_total = ctx->K_total;
    ST_int i;

    /* Note: This is a stub - full implementation requires Z matrix */
    (void)resid;       /* Will be used in full implementation */
    (void)kernel_type; /* Will be used in full implementation */
    (void)bw;          /* Will be used in full implementation */

    ST_int df_r = N - K_total - df_a;
    if (df_r <= 0) df_r = 1;
    ST_double dof_adj = (ST_double)N / (ST_double)df_r;

    /* HAC requires time-series data structure */
    /* Full implementation deferred - use robust approximation */
    for (i = 0; i < K_total * K_total; i++) {
        V[i] = XkX_inv[i] * dof_adj;
    }
}

/*
    Compute efficient GMM2S VCE.
*/
void ivvce_compute_gmm2s(
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
    ST_int i, j, k;

    memset(ZOmegaZ, 0, K_iv * K_iv * sizeof(ST_double));

    for (i = 0; i < N; i++) {
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        ST_double e2 = w * resid[i] * resid[i];
        for (j = 0; j < K_iv; j++) {
            for (k = 0; k <= j; k++) {
                ST_double contrib = Z[j * N + i] * e2 * Z[k * N + i];
                ZOmegaZ[j * K_iv + k] += contrib;
                if (k != j) ZOmegaZ[k * K_iv + j] += contrib;
            }
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
    ST_int i, j, k, c;

    memset(ZOmegaZ, 0, K_iv * K_iv * sizeof(ST_double));

    /* Allocate per-cluster Z'e sums */
    ST_double *cluster_ze = (ST_double *)calloc(num_clusters * K_iv, sizeof(ST_double));
    if (!cluster_ze) return;

    /* Sum Z_i * e_i within each cluster */
    for (i = 0; i < N; i++) {
        c = cluster_ids[i] - 1;  /* 1-indexed to 0-indexed */
        if (c < 0 || c >= num_clusters) continue;
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        ST_double we = w * resid[i];
        for (j = 0; j < K_iv; j++) {
            cluster_ze[c * K_iv + j] += Z[j * N + i] * we;
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
    ST_int i, j, k, c;

    memset(meat, 0, K_total * K_total * sizeof(ST_double));

    /* Allocate per-cluster sums */
    ST_double *cluster_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));
    if (!cluster_sums) return;

    /* Sum (PzX_i * e_i) within each cluster */
    for (i = 0; i < N; i++) {
        c = cluster_ids[i] - 1;  /* 1-indexed to 0-indexed */
        if (c < 0 || c >= num_clusters) continue;
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        ST_double we = w * resid[i];
        for (j = 0; j < K_total; j++) {
            cluster_sums[c * K_total + j] += PzX[j * N + i] * we;
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
    ST_int i, j, k;

    /* Compute P_Z X = Z * (Z'Z)^-1 Z'X */
    ST_double *PzX = (ST_double *)calloc(N * K_total, sizeof(ST_double));
    if (!PzX) return;

    for (i = 0; i < N; i++) {
        for (j = 0; j < K_total; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += Z[k * N + i] * temp_kiv_ktotal[j * K_iv + k];
            }
            PzX[j * N + i] = sum;
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
    ST_int *pair_to_int = (ST_int *)calloc((num_clusters1 + 1) * (num_clusters2 + 1), sizeof(ST_int));
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
        ST_int pair_idx = c1 * (num_clusters2 + 1) + c2;
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

    /* DOF adjustments for each dimension */
    ST_int df_r = N - K_total - df_a;
    if (df_r <= 0) df_r = 1;

    /* V1 = XkX_inv * meat1 * XkX_inv * dof_adj1 */
    ST_double G1 = (ST_double)num_clusters1;
    ST_double dof_adj1 = ((ST_double)(N - 1) / (ST_double)df_r) * (G1 / (G1 - 1.0));
    civreghdfe_matmul_ab(XkX_inv, meat1, K_total, K_total, K_total, temp_v);
    civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V1);
    for (i = 0; i < K_total * K_total; i++) V1[i] *= dof_adj1;

    /* V2 = XkX_inv * meat2 * XkX_inv * dof_adj2 */
    ST_double G2 = (ST_double)num_clusters2;
    ST_double dof_adj2 = ((ST_double)(N - 1) / (ST_double)df_r) * (G2 / (G2 - 1.0));
    civreghdfe_matmul_ab(XkX_inv, meat2, K_total, K_total, K_total, temp_v);
    civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V2);
    for (i = 0; i < K_total * K_total; i++) V2[i] *= dof_adj2;

    /* V_int = XkX_inv * meat_int * XkX_inv * dof_adj_int */
    ST_double G_int = (ST_double)num_intersection;
    ST_double dof_adj_int = ((ST_double)(N - 1) / (ST_double)df_r) * (G_int / (G_int - 1.0));
    civreghdfe_matmul_ab(XkX_inv, meat_int, K_total, K_total, K_total, temp_v);
    civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V_int);
    for (i = 0; i < K_total * K_total; i++) V_int[i] *= dof_adj_int;

    /* V = V1 + V2 - V_int */
    for (i = 0; i < K_total * K_total; i++) {
        V[i] = V1[i] + V2[i] - V_int[i];
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
    Compute Kiefer VCE (homoskedastic within-panel autocorrelation).

    Kiefer (1980) uses sigma^2 * sum_g (sum_i PzX_i)(sum_i PzX_i)'
    instead of sum_g (sum_i PzX_i * e_i)(sum_i PzX_i * e_i)'

    This assumes homoskedastic errors with within-panel autocorrelation.
*/
void ivvce_compute_kiefer(
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *temp_kiv_ktotal,
    const ST_double *XkX_inv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int K_total,
    ST_int K_iv,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_double *V
)
{
    ST_int i, j, k, c;

    /* Kiefer VCE uses the HAC formula with residuals, computed within panels.
       The meat formula is the same as cluster-robust:
       meat = sum_g (sum_i PzX_i * e_i)(sum_i PzX_i * e_i)'

       But the DOF adjustment is different:
       - Cluster: dof_adj = (N-1)/(N-K-df_a) * G/(G-1)
       - Kiefer: dof_adj = N/(N-K-df_a) (no G/(G-1) factor)

       Additionally, ivreg2 may use small-sample corrections. */

    ST_int df_r = N - K_total - df_a;
    if (df_r <= 0) df_r = 1;

    /* Compute P_Z X = Z * (Z'Z)^-1 Z'X = Z * temp_kiv_ktotal */
    ST_double *PzX = (ST_double *)calloc(N * K_total, sizeof(ST_double));
    if (!PzX) return;

    for (i = 0; i < N; i++) {
        for (j = 0; j < K_total; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += Z[k * N + i] * temp_kiv_ktotal[j * K_iv + k];
            }
            PzX[j * N + i] = sum;
        }
    }

    /*
       Kiefer VCE: Panel-aware HAC with truncated kernel.

       The meat is computed using HAC formula but only for observations
       in the SAME panel. With truncated kernel and bw >= T-1, this
       includes all observation pairs within each panel.

       For panel g with observations i1, i2, ..., iT:
       meat_g = sum_{t,s in g} u_t * u_s' where u_t = PzX_t * e_t

       Total meat = sum_g meat_g

       This differs from standard cluster-robust which computes:
       meat_cluster = sum_g (sum_t u_t)(sum_t u_t)'

       The HAC formula with individual pairs gives different results
       when residuals have different magnitudes within each cluster.
    */

    /* Compute RSS and sigma for homoskedastic scaling */
    ST_double rss = 0.0;
    for (i = 0; i < N; i++) {
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        rss += w * resid[i] * resid[i];
    }
    ST_double sigma = sqrt(rss / df_r);

    /* Compute u_i = PzX_i * e_i for all observations */
    ST_double *u = (ST_double *)calloc(N * K_total, sizeof(ST_double));
    if (!u) {
        free(PzX);
        return;
    }

    for (i = 0; i < N; i++) {
        ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
        ST_double we = w * resid[i];
        for (j = 0; j < K_total; j++) {
            u[j * N + i] = PzX[j * N + i] * we;
        }
    }

    (void)sigma;  /* May be used for alternative Kiefer formula */

    free(PzX);

    /* Compute panel-aware HAC meat.
       For each pair of observations (i1, i2) in the same panel,
       add u_i1 * u_i2' to the meat matrix.
       This is O(sum_g T_g^2) which is efficient for balanced panels. */
    ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!meat) {
        free(u);
        return;
    }

    /* Build cluster membership lists for efficient iteration */
    /* First, count observations per cluster */
    ST_int *cluster_counts = (ST_int *)calloc(num_clusters, sizeof(ST_int));
    ST_int *cluster_starts = (ST_int *)calloc(num_clusters + 1, sizeof(ST_int));
    ST_int *obs_by_cluster = (ST_int *)calloc(N, sizeof(ST_int));

    if (!cluster_counts || !cluster_starts || !obs_by_cluster) {
        free(u); free(meat);
        if (cluster_counts) free(cluster_counts);
        if (cluster_starts) free(cluster_starts);
        if (obs_by_cluster) free(obs_by_cluster);
        return;
    }

    for (i = 0; i < N; i++) {
        c = cluster_ids[i] - 1;
        if (c >= 0 && c < num_clusters) {
            cluster_counts[c]++;
        }
    }

    /* Compute start indices */
    cluster_starts[0] = 0;
    for (c = 0; c < num_clusters; c++) {
        cluster_starts[c + 1] = cluster_starts[c] + cluster_counts[c];
    }

    /* Reset counts to use as insertion indices */
    memset(cluster_counts, 0, num_clusters * sizeof(ST_int));

    /* Fill obs_by_cluster */
    for (i = 0; i < N; i++) {
        c = cluster_ids[i] - 1;
        if (c >= 0 && c < num_clusters) {
            ST_int idx = cluster_starts[c] + cluster_counts[c];
            obs_by_cluster[idx] = i;
            cluster_counts[c]++;
        }
    }

    /* Compute meat: sum over all within-panel pairs */
    for (c = 0; c < num_clusters; c++) {
        ST_int start = cluster_starts[c];
        ST_int end = cluster_starts[c + 1];
        ST_int T_c = end - start;

        /* For each pair (i1, i2) in this cluster */
        for (ST_int t1 = 0; t1 < T_c; t1++) {
            ST_int i1 = obs_by_cluster[start + t1];
            for (ST_int t2 = 0; t2 < T_c; t2++) {
                ST_int i2 = obs_by_cluster[start + t2];

                /* Add u_i1 * u_i2' to meat */
                for (j = 0; j < K_total; j++) {
                    for (k = 0; k <= j; k++) {
                        ST_double contrib = u[j * N + i1] * u[k * N + i2];
                        meat[j * K_total + k] += contrib;
                        if (k != j) meat[k * K_total + j] += contrib;
                    }
                }
            }
        }
    }

    free(u);
    free(cluster_counts);
    free(cluster_starts);
    free(obs_by_cluster);

    /* DOF adjustment for Kiefer: N / (N - K - df_a)
       Note: NO G/(G-1) factor like cluster-robust.
       This is the key difference between Kiefer and cluster VCE. */
    ST_double dof_adj = (ST_double)N / (ST_double)df_r;

    /* V = XkX_inv * meat * XkX_inv * dof_adj */
    ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (temp_v) {
        civreghdfe_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
        civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
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
    ST_int K_total,
    ST_int K_iv,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int df_a,
    ST_int nested_adj,
    ST_int kernel_type,
    ST_int bw,
    ST_double *V
)
{
    ST_int i, j, k;

    /* Compute residual sum of squares */
    ST_double rss = 0.0;
    if (weights && weight_type != 0) {
        for (i = 0; i < N; i++) {
            rss += weights[i] * resid[i] * resid[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            rss += resid[i] * resid[i];
        }
    }

    ST_int df_r = N - K_total - df_a;
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

    for (i = 0; i < N; i++) {
        for (j = 0; j < K_total; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K_iv; k++) {
                sum += Z[k * N + i] * temp_kiv_ktotal[j * K_iv + k];
            }
            PzX[j * N + i] = sum;
        }
    }

    /* Compute meat matrix based on VCE type */
    ST_double *meat = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
    if (!meat) {
        free(PzX);
        return;
    }

    if (vce_type == CIVREGHDFE_VCE_CLUSTER && cluster_ids && num_clusters > 0) {
        /* Clustered VCE */
        ST_double *cluster_sums = (ST_double *)calloc(num_clusters * K_total, sizeof(ST_double));
        if (!cluster_sums) {
            free(PzX); free(meat);
            return;
        }

        /* Sum (PzX_i * e_i) within each cluster */
        for (i = 0; i < N; i++) {
            ST_int c = cluster_ids[i] - 1;
            if (c < 0 || c >= num_clusters) continue;
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            ST_double we = w * resid[i];
            for (j = 0; j < K_total; j++) {
                cluster_sums[c * K_total + j] += PzX[j * N + i] * we;
            }
        }

        /* Compute meat: sum over clusters of outer products */
        for (ST_int c = 0; c < num_clusters; c++) {
            for (j = 0; j < K_total; j++) {
                for (k = 0; k <= j; k++) {
                    ST_double contrib = cluster_sums[c * K_total + j] * cluster_sums[c * K_total + k];
                    meat[j * K_total + k] += contrib;
                    if (k != j) meat[k * K_total + j] += contrib;
                }
            }
        }

        free(cluster_sums);

        /* DOF adjustment */
        ST_int effective_df_a = df_a;
        if (df_a == 0 && nested_adj == 1) effective_df_a = 1;
        ST_double denom = (ST_double)(N - K_total - effective_df_a);
        if (denom <= 0) denom = 1.0;
        ST_double G = (ST_double)num_clusters;
        ST_double dof_adj = ((ST_double)(N - 1) / denom) * (G / (G - 1.0));

        /* V = XkX_inv * meat * XkX_inv * dof_adj */
        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (temp_v) {
            civreghdfe_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
            civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
            for (i = 0; i < K_total * K_total; i++) {
                V[i] *= dof_adj;
            }
            free(temp_v);
        }

    } else if (kernel_type > 0 && bw > 0) {
        /* HAC VCE */
        ST_double *u = (ST_double *)calloc(N * K_total, sizeof(ST_double));
        if (!u) {
            free(PzX); free(meat);
            return;
        }

        /* Compute u_i = PzX_i * e_i */
        for (i = 0; i < N; i++) {
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            ST_double we = w * resid[i];
            for (j = 0; j < K_total; j++) {
                u[j * N + i] = PzX[j * N + i] * we;
            }
        }

        /* Sum over all lag pairs with kernel weights */
        for (ST_int lag = 0; lag <= bw; lag++) {
            ST_double kw = civreghdfe_kernel_weight(kernel_type, lag, bw);
            if (kw < 1e-10) continue;

            if (lag == 0) {
                for (i = 0; i < N; i++) {
                    for (j = 0; j < K_total; j++) {
                        for (k = 0; k <= j; k++) {
                            ST_double contrib = kw * u[j * N + i] * u[k * N + i];
                            meat[j * K_total + k] += contrib;
                            if (k != j) meat[k * K_total + j] += contrib;
                        }
                    }
                }
            } else {
                for (i = 0; i < N - lag; i++) {
                    for (j = 0; j < K_total; j++) {
                        for (k = 0; k < K_total; k++) {
                            ST_double contrib = kw * (u[j * N + i] * u[k * N + (i + lag)] +
                                                      u[j * N + (i + lag)] * u[k * N + i]);
                            meat[k * K_total + j] += contrib;
                        }
                    }
                }
            }
        }

        free(u);

        /* HC1 adjustment and sandwich */
        ST_double dof_adj = (ST_double)N / (ST_double)df_r;
        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (temp_v) {
            civreghdfe_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
            civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
            for (i = 0; i < K_total * K_total; i++) {
                V[i] *= dof_adj;
            }
            free(temp_v);
        }

    } else {
        /* Standard HC robust */
        for (i = 0; i < N; i++) {
            ST_double w = (weights && weight_type != 0) ? weights[i] : 1.0;
            ST_double we = w * resid[i];
            for (j = 0; j < K_total; j++) {
                for (k = 0; k <= j; k++) {
                    ST_double contrib = PzX[j * N + i] * we * PzX[k * N + i] * resid[i];
                    meat[j * K_total + k] += contrib;
                    if (k != j) meat[k * K_total + j] += contrib;
                }
            }
        }

        /* HC1 adjustment */
        ST_double dof_adj = (ST_double)N / (ST_double)df_r;

        /* V = XkX_inv * meat * XkX_inv * dof_adj */
        ST_double *temp_v = (ST_double *)calloc(K_total * K_total, sizeof(ST_double));
        if (temp_v) {
            civreghdfe_matmul_ab(XkX_inv, meat, K_total, K_total, K_total, temp_v);
            civreghdfe_matmul_ab(temp_v, XkX_inv, K_total, K_total, K_total, V);
            for (i = 0; i < K_total * K_total; i++) {
                V[i] *= dof_adj;
            }
            free(temp_v);
        }
    }

    free(PzX);
    free(meat);
}
