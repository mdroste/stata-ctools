/*
 * ctools_matrix.c
 *
 * Shared matrix operations and sandwich VCE for creghdfe and civreghdfe.
 * Provides optimized matrix multiplication with OpenMP and SIMD,
 * plus robust and clustered variance-covariance estimation.
 * Part of the ctools Stata plugin suite.
 */

#include <stdlib.h>
#include <math.h>
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

    memset(C, 0, (size_t)K1 * K2 * sizeof(ST_double));

    /* Parallelize over output columns - i must be private to avoid race condition */
    #pragma omp parallel for schedule(static) if((size_t)K1 * K2 > 4)
    for (j = 0; j < K2; j++) {
        ST_int i;  /* Declare inside parallel region to make it private */
        const ST_double *b_col = B + (size_t)j * N;
        for (i = 0; i < K1; i++) {
            const ST_double *a_col = A + (size_t)i * N;
            C[(size_t)j * K1 + i] = ctools_dot_unrolled(a_col, b_col, N);
        }
    }
}

/*
 * Gram matrix C = A' * A (upper triangle computed, then mirrored)
 */
void ctools_matmul_ata(const ST_double *A, ST_int N, ST_int K, ST_double *C)
{
    ST_int j;

    #pragma omp parallel for schedule(dynamic) if(K > 2)
    for (j = 0; j < K; j++) {
        ST_int i;
        const ST_double *a_j = A + (size_t)j * N;
        for (i = 0; i <= j; i++) {
            C[(size_t)j * K + i] = ctools_dot_unrolled(A + (size_t)i * N, a_j, N);
        }
    }
    for (j = 0; j < K; j++) {
        for (ST_int i = 0; i < j; i++) {
            C[(size_t)i * K + j] = C[(size_t)j * K + i];
        }
    }
}

/* Rows per block of the row-blocked kernels: a block of a few columns
 * stays in the L2 cache while the K passes over it run */
#define CTOOLS_ROW_BLOCK 2048

/*
 * C = X * B (and y += X * b): row blocks in parallel, columns in order
 */
static void matmul_xb_blocks(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                             ST_double *C, int zero_first)
{
    /* N + CTOOLS_ROW_BLOCK - 1 would overflow int for N within a block of
     * INT_MAX, leaving every row uncomputed */
    const ST_int num_blocks = N / CTOOLS_ROW_BLOCK + (N % CTOOLS_ROW_BLOCK != 0);
    ST_int blk;

    #pragma omp parallel for schedule(static) if((size_t)N * (size_t)K * (size_t)M > 400000)
    for (blk = 0; blk < num_blocks; blk++) {
        const ST_int i0 = blk * CTOOLS_ROW_BLOCK;
        const ST_int i1 = (N - i0 < CTOOLS_ROW_BLOCK) ? N : i0 + CTOOLS_ROW_BLOCK;
        for (ST_int m = 0; m < M; m++) {
            ST_double *c = C + (size_t)m * N;
            const ST_double *b = B + (size_t)m * K;
            if (zero_first) {
                for (ST_int i = i0; i < i1; i++) c[i] = 0.0;
            }
            for (ST_int k = 0; k < K; k++) {
                const ST_double *x = X + (size_t)k * N;
                const ST_double bk = b[k];
                for (ST_int i = i0; i < i1; i++) c[i] += x[i] * bk;
            }
        }
    }
}

void ctools_matvec_add(const ST_double *X, ST_int N, ST_int K, const ST_double *b, ST_double *y)
{
    matmul_xb_blocks(X, N, K, b, 1, y, 0);
}

void ctools_matmul_xb(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                      ST_double *C)
{
    matmul_xb_blocks(X, N, K, B, M, C, 1);
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

    memset(C, 0, (size_t)K1 * K3 * sizeof(ST_double));

    /* C[i,j] = sum_k A[i,k] * B[k,j] */
    /* A is K1 x K2, stored column-major: A[i,k] = A[k*K1 + i] */
    /* B is K2 x K3, stored column-major: B[k,j] = B[j*K2 + k] */
    /* C is K1 x K3, stored column-major: C[i,j] = C[j*K1 + i] */

    /* Reorder loops for better cache locality: k-j-i instead of j-i-k */
    /* This allows sequential access to A columns and B columns */
    #pragma omp parallel for schedule(static) if((size_t)K1 * K3 > 16)
    for (j = 0; j < K3; j++) {
        ST_int i, k;  /* Private to each thread */
        const ST_double *b_col = B + (size_t)j * K2;
        ST_double *c_col = C + (size_t)j * K1;
        for (k = 0; k < K2; k++) {
            ST_double b_kj = b_col[k];
            const ST_double *a_col = A + (size_t)k * K1;
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

    memset(C, 0, (size_t)K1 * K2 * sizeof(ST_double));

    if (w) {
        /* Weighted case: use SIMD-accelerated weighted dot product */
        #pragma omp parallel for schedule(static) if((size_t)K1 * K2 > 4)
        for (j = 0; j < K2; j++) {
            ST_int i;  /* Private to each thread */
            const ST_double *b_col = B + (size_t)j * N;
            for (i = 0; i < K1; i++) {
                const ST_double *a_col = A + (size_t)i * N;
                C[(size_t)j * K1 + i] = ctools_simd_dot_weighted(w, a_col, b_col, (size_t)N);
            }
        }
    } else {
        /* Unweighted case: use SIMD-accelerated dot product */
        #pragma omp parallel for schedule(static) if((size_t)K1 * K2 > 4)
        for (j = 0; j < K2; j++) {
            ST_int i;  /* Private to each thread */
            const ST_double *b_col = B + (size_t)j * N;
            for (i = 0; i < K1; i++) {
                const ST_double *a_col = A + (size_t)i * N;
                C[(size_t)j * K1 + i] = ctools_simd_dot(a_col, b_col, (size_t)N);
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
ST_retcode ctools_vce_robust(const ctools_vce_data *d, ST_double dof_adj, ST_double *V)
{
    ST_int idx, i, j, k;
    ST_int N = d->N;
    ST_int K = d->K;

    if (N < 1 || K < 1 || !isfinite(dof_adj) || dof_adj <= 0)
        return 498;

    /* Normalize weights for aweight/pweight (OLS convention) */
    ST_double *w_norm = NULL;
    if (d->normalize_weights && d->weights != NULL && (d->weight_type == 1 || d->weight_type == 3)) {
        ST_double sum_w = 0.0;
        w_norm = (ST_double *)malloc(N * sizeof(ST_double));
        if (!w_norm) return 920;
        for (idx = 0; idx < N; idx++) sum_w += d->weights[idx];
        if (!isfinite(sum_w) || sum_w <= 0) { free(w_norm); return 498; }
        ST_double scale = (ST_double)N / sum_w;
        for (idx = 0; idx < N; idx++) w_norm[idx] = d->weights[idx] * scale;
    }

    /* Allocate high-precision meat and temp matrices */
    long double *meat = (long double *)calloc((size_t)K * K, sizeof(long double));
    long double *temp = (long double *)malloc((size_t)K * K * sizeof(long double));

    if (!meat || !temp) {
        if (meat) free(meat);
        if (temp) free(temp);
        if (w_norm) free(w_norm);
        return 920;
    }

    /* Compute meat = X_eff' * diag(w_i) * X_eff
     * where w_i depends on weight type and normalization:
     * - No weights:      w_i = e_i^2
     * - fweight:          w_i = e_i^2 * fw_i
     * - aw/pw normalized: w_i = (e_i * w_norm_i)^2
     * - aw/pw raw (IV):   w_i = (w_i * e_i)^2    [but actually w_i * e_i^2]
     * Symmetric: accumulate upper triangle, then mirror. */
    ST_double *obs_w = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    if (!obs_w) {
        free(meat); free(temp);
        if (w_norm) free(w_norm);
        return 920;
    }
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
        obs_w[idx] = w_i;
    }

    /* meat[j,k] = sum_i (w_i * x_ij) * x_ik, each element accumulated over
     * the observations in order. Column pairs are independent, so the rows
     * of the upper triangle run in parallel, a few columns k per pass. */
    #pragma omp parallel for schedule(dynamic) if(K > 1 && (size_t)N * (size_t)K > 100000)
    for (j = 0; j < K; j++) {
        const ST_double *xj = d->X_eff + (size_t)j * N;
        ST_int kk = j;
        for (; kk + 3 < K; kk += 4) {
            const ST_double *x0 = d->X_eff + (size_t)kk * N;
            const ST_double *x1 = x0 + N, *x2 = x1 + N, *x3 = x2 + N;
            long double a0 = 0.0L, a1 = 0.0L, a2 = 0.0L, a3 = 0.0L;
            for (ST_int ii = 0; ii < N; ii++) {
                long double wxj = (long double)obs_w[ii] * (long double)xj[ii];
                a0 += wxj * (long double)x0[ii];
                a1 += wxj * (long double)x1[ii];
                a2 += wxj * (long double)x2[ii];
                a3 += wxj * (long double)x3[ii];
            }
            meat[(size_t)j * K + kk] = a0;
            meat[(size_t)j * K + kk + 1] = a1;
            meat[(size_t)j * K + kk + 2] = a2;
            meat[(size_t)j * K + kk + 3] = a3;
        }
        for (; kk < K; kk++) {
            const ST_double *xk = d->X_eff + (size_t)kk * N;
            long double acc = 0.0L;
            for (ST_int ii = 0; ii < N; ii++) {
                long double wxj = (long double)obs_w[ii] * (long double)xj[ii];
                acc += wxj * (long double)xk[ii];
            }
            meat[(size_t)j * K + kk] = acc;
        }
    }
    free(obs_w);
    /* Mirror upper triangle to lower */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            meat[(size_t)k * K + j] = meat[(size_t)j * K + k];
        }
    }

    if (w_norm) free(w_norm);

    /* High-precision V = D * meat * D * dof_adj (column-major). */
    for (j = 0; j < K; j++) {
        for (i = 0; i < K; i++) {
            long double acc = 0.0L;
            for (k = 0; k < K; k++) {
                acc += (long double)d->D[(size_t)k * K + i] * meat[(size_t)j * K + k];
            }
            temp[(size_t)j * K + i] = acc;
        }
    }
    for (j = 0; j < K; j++) {
        for (i = 0; i < K; i++) {
            long double acc = 0.0L;
            for (k = 0; k < K; k++) {
                acc += temp[(size_t)k * K + i] * (long double)d->D[(size_t)j * K + k];
            }
            ST_double value = (ST_double)(acc * (long double)dof_adj);
            if (!isfinite(value)) { free(meat); free(temp); return 498; }
            meat[(size_t)j * K + i] = value;
        }
    }

    for (size_t entry = 0; entry < (size_t)K * K; entry++) V[entry] = (ST_double)meat[entry];
    free(meat);
    free(temp);
    return 0;
}

/*
 * Compute clustered sandwich VCE.
 */
ST_retcode ctools_vce_cluster(const ctools_vce_data *d,
                        const ST_int *cluster_ids, ST_int num_clusters,
                        ST_double dof_adj, ST_double *V)
{
    ST_int k, idx;
    ST_int N = d->N;
    ST_int K = d->K;

    if (N < 1 || K < 1 || !isfinite(dof_adj) || dof_adj <= 0)
        return 498;

    if (!cluster_ids || num_clusters <= 1) return 498;
    for (idx = 0; idx < N; idx++) {
        if (cluster_ids[idx] < 0 || cluster_ids[idx] >= num_clusters) return 498;
    }

    /* Normalize weights for aweight/pweight (OLS convention) */
    ST_double *w_norm = NULL;
    if (d->normalize_weights && d->weights != NULL && (d->weight_type == 1 || d->weight_type == 3)) {
        ST_double sum_w = 0.0;
        w_norm = (ST_double *)malloc(N * sizeof(ST_double));
        if (!w_norm) return 920;
        for (idx = 0; idx < N; idx++) sum_w += d->weights[idx];
        if (!isfinite(sum_w) || sum_w <= 0) { free(w_norm); return 498; }
        ST_double scale = (ST_double)N / sum_w;
        for (idx = 0; idx < N; idx++) w_norm[idx] = d->weights[idx] * scale;
    }

    /* Allocate high-precision buffers */
    long double *meat = (long double *)calloc((size_t)K * K, sizeof(long double));
    long double *temp = (long double *)malloc((size_t)K * K * sizeof(long double));
    long double *ecX = (long double *)malloc(K * sizeof(long double));
    ST_int *sort_perm = (ST_int *)malloc(N * sizeof(ST_int));
    ST_int *boundaries = (ST_int *)malloc(((size_t)num_clusters + 1) * sizeof(ST_int));

    if (!meat || !temp || !ecX || !sort_perm || !boundaries) {
        if (meat) free(meat);
        if (temp) free(temp);
        if (ecX) free(ecX);
        if (sort_perm) free(sort_perm);
        if (boundaries) free(boundaries);
        if (w_norm) free(w_norm);
        return 920;
    }

    /* Sort observations by cluster using counting sort */
    if (vce_sort_by_cluster(cluster_ids, N, num_clusters, sort_perm, boundaries) != 0) {
        free(meat); free(temp); free(ecX);
        free(sort_perm); free(boundaries);
        if (w_norm) free(w_norm);
        return 920;
    }

    /* Weighted residual of every observation */
    ST_double *e_w = (ST_double *)malloc((size_t)N * sizeof(ST_double));
    /* Clusters are processed in blocks: the e'X sums of a block are formed
     * one regressor per thread (each sum over its cluster's observations in
     * sorted order), then the outer products are added in cluster order. */
    ST_int block_clusters = 65536 / (K > 0 ? K : 1);
    if (block_clusters < 1) block_clusters = 1;
    if (block_clusters > num_clusters) block_clusters = num_clusters;
    long double *block_sums = (long double *)malloc((size_t)block_clusters * (size_t)K * sizeof(long double));
    if (!e_w || !block_sums) {
        free(e_w); free(block_sums);
        free(meat); free(temp); free(ecX);
        free(sort_perm); free(boundaries);
        if (w_norm) free(w_norm);
        return 920;
    }
    for (idx = 0; idx < N; idx++) {
        if (d->weights == NULL) {
            e_w[idx] = d->resid[idx];
        } else if (d->weight_type == 2) {
            e_w[idx] = d->resid[idx] * d->weights[idx];
        } else if (w_norm) {
            e_w[idx] = d->resid[idx] * w_norm[idx];
        } else {
            e_w[idx] = d->resid[idx] * d->weights[idx];
        }
    }

    /* c0 advances to c1, which never passes num_clusters (c0 += block_clusters
     * could overflow int on the last block) */
    for (ST_int c0 = 0, c1; c0 < num_clusters; c0 = c1) {
        c1 = (num_clusters - c0 < block_clusters) ? num_clusters : c0 + block_clusters;
        const size_t block_obs = (size_t)(boundaries[c1] - boundaries[c0]);

        #pragma omp parallel for schedule(static) if(K > 1 && block_obs * (size_t)K > 100000)
        for (k = 0; k < K; k++) {
            const ST_double *xk = d->X_eff + (size_t)k * N;
            for (ST_int c = c0; c < c1; c++) {
                long double acc = 0.0L;
                for (ST_int pos = boundaries[c]; pos < boundaries[c + 1]; pos++) {
                    ST_int obs = sort_perm[pos];
                    acc += (long double)e_w[obs] * (long double)xk[obs];
                }
                block_sums[(size_t)(c - c0) * K + k] = acc;
            }
        }

        /* Accumulate outer product ecX * ecX' into meat (upper triangle) */
        for (ST_int c = c0; c < c1; c++) {
            const long double *sums = block_sums + (size_t)(c - c0) * K;
            for (ST_int jj = 0; jj < K; jj++) {
                long double ecXj = sums[jj];
                meat[(size_t)jj * K + jj] += ecXj * ecXj;
                for (ST_int kk = jj + 1; kk < K; kk++) {
                    meat[(size_t)jj * K + kk] += ecXj * sums[kk];
                }
            }
        }
    }
    free(e_w);
    free(block_sums);

    free(ecX);
    free(sort_perm);
    free(boundaries);
    if (w_norm) free(w_norm);

    /* Mirror upper triangle to lower */
    for (ST_int jj = 0; jj < K; jj++) {
        for (ST_int kk = jj + 1; kk < K; kk++) {
            meat[(size_t)kk * K + jj] = meat[(size_t)jj * K + kk];
        }
    }

    /* High-precision V = D * meat * D * dof_adj (column-major). */
    for (ST_int jj = 0; jj < K; jj++) {
        for (ST_int ii = 0; ii < K; ii++) {
            long double acc = 0.0L;
            for (ST_int kk = 0; kk < K; kk++) {
                acc += (long double)d->D[(size_t)kk * K + ii] * meat[(size_t)jj * K + kk];
            }
            temp[(size_t)jj * K + ii] = acc;
        }
    }
    for (ST_int jj = 0; jj < K; jj++) {
        for (ST_int ii = 0; ii < K; ii++) {
            long double acc = 0.0L;
            for (ST_int kk = 0; kk < K; kk++) {
                acc += temp[(size_t)kk * K + ii] * (long double)d->D[(size_t)jj * K + kk];
            }
            ST_double value = (ST_double)(acc * (long double)dof_adj);
            if (!isfinite(value)) { free(meat); free(temp); return 498; }
            meat[(size_t)jj * K + ii] = value;
        }
    }

    for (size_t entry = 0; entry < (size_t)K * K; entry++) V[entry] = (ST_double)meat[entry];
    free(meat);
    free(temp);
    return 0;
}
