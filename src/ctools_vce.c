/*
 * ctools_vce.c
 *
 * Shared sandwich VCE computation for creghdfe and civreghdfe.
 * Computes V = D * meat * D * dof_adj for both OLS and IV estimators.
 * Part of the ctools Stata plugin suite.
 */

#include <stdlib.h>
#include <string.h>
#include "ctools_vce.h"
#include "ctools_matrix.h"

/*
 * Sort observations by cluster_id using counting sort.
 * cluster_ids are 0-indexed and contiguous.
 * O(N + num_clusters) time.
 */
static int vce_sort_by_cluster(const ST_int *cluster_ids, ST_int N, ST_int num_clusters,
                               ST_int *sort_perm, ST_int *boundaries)
{
    ST_int i;

    if (N <= 0 || !cluster_ids || !sort_perm || !boundaries || num_clusters <= 0) {
        return -1;
    }

    /* Count observations per cluster */
    ST_int *counts = (ST_int *)calloc(num_clusters, sizeof(ST_int));
    if (!counts) return -1;

    for (i = 0; i < N; i++) {
        counts[cluster_ids[i]]++;
    }

    /* Compute prefix sums to get starting positions */
    boundaries[0] = 0;
    for (ST_int c = 0; c < num_clusters; c++) {
        boundaries[c + 1] = boundaries[c] + counts[c];
    }

    /* Reset counts to use as insertion pointers */
    memset(counts, 0, num_clusters * sizeof(ST_int));

    /* Fill sort_perm */
    for (i = 0; i < N; i++) {
        ST_int c = cluster_ids[i];
        ST_int pos = boundaries[c] + counts[c];
        sort_perm[pos] = i;
        counts[c]++;
    }

    free(counts);
    return 0;
}

/*
 * Compute robust (HC1) sandwich VCE.
 */
void ctools_vce_robust(const ctools_vce_data *d, ST_double dof_adj, ST_double *V)
{
    ST_int idx, j, k;
    ST_int N = d->N;
    ST_int K = d->K;

    /* Normalize weights for aweight/pweight */
    ST_double *w_norm = NULL;
    if (d->weights != NULL && (d->weight_type == 1 || d->weight_type == 3)) {
        ST_double sum_w = 0.0;
        w_norm = (ST_double *)malloc(N * sizeof(ST_double));
        if (!w_norm) return;
        for (idx = 0; idx < N; idx++) sum_w += d->weights[idx];
        ST_double scale = (ST_double)N / sum_w;
        for (idx = 0; idx < N; idx++) w_norm[idx] = d->weights[idx] * scale;
    }

    /* Allocate meat and temp matrices */
    ST_double *meat = (ST_double *)calloc(K * K, sizeof(ST_double));
    ST_double *temp = (ST_double *)malloc(K * K * sizeof(ST_double));

    if (!meat || !temp) {
        if (meat) free(meat);
        if (temp) free(temp);
        if (w_norm) free(w_norm);
        return;
    }

    /* Compute meat = X_eff' * diag(w_i) * X_eff
     * where w_i depends on weight type:
     * - No weights:     w_i = e_i^2
     * - fweight:        w_i = e_i^2 * fw_i
     * - aweight/pweight: w_i = (e_i * w_norm_i)^2 */
    for (idx = 0; idx < N; idx++) {
        ST_double w_i;
        ST_double e = d->resid[idx];

        if (d->weights == NULL) {
            w_i = e * e;
        } else if (d->weight_type == 2) {
            w_i = e * e * d->weights[idx];
        } else {
            ST_double ew = e * w_norm[idx];
            w_i = ew * ew;
        }

        for (j = 0; j < K; j++) {
            ST_double xj = d->X_eff[j * N + idx];
            for (k = 0; k < K; k++) {
                ST_double xk = d->X_eff[k * N + idx];
                meat[j * K + k] += w_i * xj * xk;
            }
        }
    }

    if (w_norm) free(w_norm);

    /* V = D * meat * D * dof_adj */
    ctools_matmul_ab(d->D, meat, K, K, K, temp);
    ctools_matmul_ab(temp, d->D, K, K, K, V);
    for (ST_int i = 0; i < K * K; i++) {
        V[i] *= dof_adj;
    }

    free(meat);
    free(temp);
}

/*
 * Compute clustered sandwich VCE.
 */
void ctools_vce_cluster(const ctools_vce_data *d,
                        const ST_int *cluster_ids, ST_int num_clusters,
                        ST_double dof_adj, ST_double *V)
{
    ST_int i, k, idx;
    ST_int N = d->N;
    ST_int K = d->K;

    /* Normalize weights for aweight/pweight */
    ST_double *w_norm = NULL;
    if (d->weights != NULL && (d->weight_type == 1 || d->weight_type == 3)) {
        ST_double sum_w = 0.0;
        w_norm = (ST_double *)malloc(N * sizeof(ST_double));
        if (!w_norm) return;
        for (idx = 0; idx < N; idx++) sum_w += d->weights[idx];
        ST_double scale = (ST_double)N / sum_w;
        for (idx = 0; idx < N; idx++) w_norm[idx] = d->weights[idx] * scale;
    }

    /* Allocate buffers */
    ST_double *meat = (ST_double *)calloc(K * K, sizeof(ST_double));
    ST_double *temp = (ST_double *)malloc(K * K * sizeof(ST_double));
    ST_double *ecX = (ST_double *)malloc(K * sizeof(ST_double));
    ST_int *sort_perm = (ST_int *)malloc(N * sizeof(ST_int));
    ST_int *boundaries = (ST_int *)malloc((num_clusters + 1) * sizeof(ST_int));

    if (!meat || !temp || !ecX || !sort_perm || !boundaries) {
        if (meat) free(meat);
        if (temp) free(temp);
        if (ecX) free(ecX);
        if (sort_perm) free(sort_perm);
        if (boundaries) free(boundaries);
        if (w_norm) free(w_norm);
        return;
    }

    /* Sort observations by cluster using counting sort */
    if (vce_sort_by_cluster(cluster_ids, N, num_clusters, sort_perm, boundaries) != 0) {
        free(meat); free(temp); free(ecX);
        free(sort_perm); free(boundaries);
        if (w_norm) free(w_norm);
        return;
    }

    /* Stream through clusters, accumulating ecX then outer product into meat */
    for (ST_int c = 0; c < num_clusters; c++) {
        ST_int start = boundaries[c];
        ST_int end = boundaries[c + 1];

        /* Reset ecX buffer */
        for (k = 0; k < K; k++) {
            ecX[k] = 0.0;
        }

        /* Accumulate e'X for this cluster */
        for (ST_int pos = start; pos < end; pos++) {
            idx = sort_perm[pos];
            ST_double e_w;

            if (d->weights == NULL) {
                e_w = d->resid[idx];
            } else if (d->weight_type == 2) {
                e_w = d->resid[idx] * d->weights[idx];
            } else {
                e_w = d->resid[idx] * w_norm[idx];
            }

            for (k = 0; k < K; k++) {
                ecX[k] += e_w * d->X_eff[k * N + idx];
            }
        }

        /* Accumulate outer product ecX * ecX' into meat */
        for (ST_int jj = 0; jj < K; jj++) {
            for (ST_int kk = 0; kk < K; kk++) {
                meat[jj * K + kk] += ecX[jj] * ecX[kk];
            }
        }
    }

    free(ecX);
    free(sort_perm);
    free(boundaries);
    if (w_norm) free(w_norm);

    /* V = D * meat * D * dof_adj */
    ctools_matmul_ab(d->D, meat, K, K, K, temp);
    ctools_matmul_ab(temp, d->D, K, K, K, V);
    for (i = 0; i < K * K; i++) {
        V[i] *= dof_adj;
    }

    free(meat);
    free(temp);
}
