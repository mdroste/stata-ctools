/*
 * ctools_matrix.c
 *
 * Shared matrix operations and sandwich VCE for creghdfe and civreghdfe.
 * Provides optimized matrix multiplication with OpenMP and SIMD,
 * plus robust and clustered variance-covariance estimation.
 * Part of the ctools Stata plugin suite.
 */

#include <stdlib.h>
#include <string.h>

#include "ctools_matrix.h"
#include "ctools_config.h"
#include "ctools_unroll.h"
#include "ctools_simd.h"

/* ============================================================================
 * Matrix Multiplication
 * ============================================================================ */

/*
 * Matrix multiply C = A' * B
 * A is N x K1, B is N x K2, Result C is K1 x K2
 *
 * OPTIMIZED: Uses K-way unrolled dot products and OpenMP parallelization
 */
void ctools_matmul_atb(const ST_double * restrict A, const ST_double * restrict B,
                       ST_int N, ST_int K1, ST_int K2,
                       ST_double * restrict C)
{
    ST_int j;

    memset(C, 0, K1 * K2 * sizeof(ST_double));

    /* Parallelize over output columns - i must be private to avoid race condition */
    #pragma omp parallel for schedule(static) if(K1 * K2 > 4)
    for (j = 0; j < K2; j++) {
        ST_int i;  /* Declare inside parallel region to make it private */
        const ST_double *b_col = B + j * N;
        for (i = 0; i < K1; i++) {
            const ST_double *a_col = A + i * N;
            C[j * K1 + i] = ctools_dot_unrolled(a_col, b_col, N);
        }
    }
}

/*
 * Matrix multiply C = A * B
 * A is K1 x K2, B is K2 x K3, Result C is K1 x K3
 *
 * OPTIMIZED: Uses cache-friendly loop order and OpenMP parallelization
 */
void ctools_matmul_ab(const ST_double * restrict A, const ST_double * restrict B,
                      ST_int K1, ST_int K2, ST_int K3,
                      ST_double * restrict C)
{
    ST_int j;

    memset(C, 0, K1 * K3 * sizeof(ST_double));

    /* C[i,j] = sum_k A[i,k] * B[k,j] */
    /* A is K1 x K2, stored column-major: A[i,k] = A[k*K1 + i] */
    /* B is K2 x K3, stored column-major: B[k,j] = B[j*K2 + k] */
    /* C is K1 x K3, stored column-major: C[i,j] = C[j*K1 + i] */

    /* Reorder loops for better cache locality: k-j-i instead of j-i-k */
    /* This allows sequential access to A columns and B columns */
    #pragma omp parallel for schedule(static) if(K1 * K3 > 16)
    for (j = 0; j < K3; j++) {
        ST_int i, k;  /* Private to each thread */
        const ST_double *b_col = B + j * K2;
        ST_double *c_col = C + j * K1;
        for (k = 0; k < K2; k++) {
            ST_double b_kj = b_col[k];
            const ST_double *a_col = A + k * K1;
            #pragma omp simd
            for (i = 0; i < K1; i++) {
                c_col[i] += a_col[i] * b_kj;
            }
        }
    }
}

/*
 * Weighted matrix multiply C = A' * diag(w) * B
 *
 * OPTIMIZED: Uses SIMD-accelerated weighted dot products (AVX2/NEON with FMA)
 * and OpenMP parallelization
 */
void ctools_matmul_atdb(const ST_double * restrict A, const ST_double * restrict B,
                        const ST_double * restrict w, ST_int N, ST_int K1, ST_int K2,
                        ST_double * restrict C)
{
    ST_int j;

    memset(C, 0, K1 * K2 * sizeof(ST_double));

    if (w) {
        /* Weighted case: use SIMD-accelerated weighted dot product */
        #pragma omp parallel for schedule(static) if(K1 * K2 > 4)
        for (j = 0; j < K2; j++) {
            ST_int i;  /* Private to each thread */
            const ST_double *b_col = B + j * N;
            for (i = 0; i < K1; i++) {
                const ST_double *a_col = A + i * N;
                C[j * K1 + i] = ctools_simd_dot_weighted(w, a_col, b_col, (size_t)N);
            }
        }
    } else {
        /* Unweighted case: use SIMD-accelerated dot product */
        #pragma omp parallel for schedule(static) if(K1 * K2 > 4)
        for (j = 0; j < K2; j++) {
            ST_int i;  /* Private to each thread */
            const ST_double *b_col = B + j * N;
            for (i = 0; i < K1; i++) {
                const ST_double *a_col = A + i * N;
                C[j * K1 + i] = ctools_simd_dot(a_col, b_col, (size_t)N);
            }
        }
    }
}

/* ============================================================================
 * Sandwich VCE
 * ============================================================================ */

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
        ST_int c = cluster_ids[i];
        if (c < 0 || c >= num_clusters) { free(counts); return -1; }
        counts[c]++;
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

    /* Zero V on entry so any early return leaves a known state (zero SE) */
    memset(V, 0, K * K * sizeof(ST_double));

    /* Normalize weights for aweight/pweight (OLS convention) */
    ST_double *w_norm = NULL;
    if (d->normalize_weights && d->weights != NULL && (d->weight_type == 1 || d->weight_type == 3)) {
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
     * where w_i depends on weight type and normalization:
     * - No weights:      w_i = e_i^2
     * - fweight:          w_i = e_i^2 * fw_i
     * - aw/pw normalized: w_i = (e_i * w_norm_i)^2
     * - aw/pw raw (IV):   w_i = (w_i * e_i)^2    [but actually w_i * e_i^2]
     * Symmetric: accumulate upper triangle, then mirror. */
    for (idx = 0; idx < N; idx++) {
        ST_double w_i;
        ST_double e = d->resid[idx];

        if (d->weights == NULL) {
            w_i = e * e;
        } else if (d->weight_type == 2) {
            w_i = e * e * d->weights[idx];
        } else if (w_norm) {
            ST_double ew = e * w_norm[idx];
            w_i = ew * ew;
        } else {
            /* Raw weights (IV path): w * e^2 */
            w_i = d->weights[idx] * e * e;
        }

        for (j = 0; j < K; j++) {
            ST_double wxj = w_i * d->X_eff[j * N + idx];
            meat[j * K + j] += wxj * d->X_eff[j * N + idx];
            for (k = j + 1; k < K; k++) {
                meat[j * K + k] += wxj * d->X_eff[k * N + idx];
            }
        }
    }
    /* Mirror upper triangle to lower */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            meat[k * K + j] = meat[j * K + k];
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

    /* Zero V on entry so any early return leaves a known state (zero SE) */
    memset(V, 0, K * K * sizeof(ST_double));

    /* Normalize weights for aweight/pweight (OLS convention) */
    ST_double *w_norm = NULL;
    if (d->normalize_weights && d->weights != NULL && (d->weight_type == 1 || d->weight_type == 3)) {
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
            } else if (w_norm) {
                e_w = d->resid[idx] * w_norm[idx];
            } else {
                e_w = d->resid[idx] * d->weights[idx];
            }

            for (k = 0; k < K; k++) {
                ecX[k] += e_w * d->X_eff[k * N + idx];
            }
        }

        /* Accumulate outer product ecX * ecX' into meat (upper triangle) */
        for (ST_int jj = 0; jj < K; jj++) {
            ST_double ecXj = ecX[jj];
            meat[jj * K + jj] += ecXj * ecXj;
            for (ST_int kk = jj + 1; kk < K; kk++) {
                meat[jj * K + kk] += ecXj * ecX[kk];
            }
        }
    }

    free(ecX);
    free(sort_perm);
    free(boundaries);
    if (w_norm) free(w_norm);

    /* Mirror upper triangle to lower */
    for (ST_int jj = 0; jj < K; jj++) {
        for (ST_int kk = jj + 1; kk < K; kk++) {
            meat[kk * K + jj] = meat[jj * K + kk];
        }
    }

    /* V = D * meat * D * dof_adj */
    ctools_matmul_ab(d->D, meat, K, K, K, temp);
    ctools_matmul_ab(temp, d->D, K, K, K, V);
    for (i = 0; i < K * K; i++) {
        V[i] *= dof_adj;
    }

    free(meat);
    free(temp);
}
