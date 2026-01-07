/*
 * cqreg_vce.c
 *
 * Variance-covariance estimation for quantile regression.
 * Part of the ctools suite.
 */

#include "cqreg_vce.h"
#include "cqreg_linalg.h"
#include "cqreg_sparsity.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================================
 * Helper Functions
 * ============================================================================ */

ST_int cqreg_compute_xtx_inv(ST_double *XtX_inv,
                             const ST_double *X,
                             ST_int N, ST_int K)
{
    ST_int j, k, i;
    ST_double *XtX = NULL;
    ST_double *L = NULL;
    ST_int rc = 0;

    /* Allocate temporary storage */
    XtX = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    L = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);

    if (XtX == NULL || L == NULL) {
        rc = -1;
        goto cleanup;
    }

    /* Compute X'X */
    memset(XtX, 0, K * K * sizeof(ST_double));

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) if(K > 2)
    #endif
    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[j * N];

        /* Diagonal */
        XtX[j * K + j] = cqreg_dot_self(Xj, N);

        /* Off-diagonal */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[k * N];
            ST_double dot = cqreg_dot(Xj, Xk, N);
            XtX[j * K + k] = dot;
            XtX[k * K + j] = dot;
        }
    }

    /* Cholesky decomposition */
    memcpy(L, XtX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(L, K) != 0) {
        /* Try with regularization */
        memcpy(L, XtX, K * K * sizeof(ST_double));
        cqreg_add_regularization(L, K, 1e-10);
        if (cqreg_cholesky(L, K) != 0) {
            rc = -1;
            goto cleanup;
        }
    }

    /* Compute inverse using Cholesky factor */
    cqreg_invert_cholesky(XtX_inv, L, K);

cleanup:
    cqreg_aligned_free(XtX);
    cqreg_aligned_free(L);

    return rc;
}

void cqreg_sandwich_product(ST_double *V,
                            const ST_double *A,
                            const ST_double *B,
                            ST_int K)
{
    ST_int i, j, k, l;
    ST_double *AB = NULL;

    /* Allocate temporary for A * B */
    AB = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (AB == NULL) {
        memset(V, 0, K * K * sizeof(ST_double));
        return;
    }

    /* Compute AB = A * B */
    memset(AB, 0, K * K * sizeof(ST_double));
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += A[i * K + k] * B[k * K + j];
            }
            AB[i * K + j] = sum;
        }
    }

    /* Compute V = AB * A' */
    memset(V, 0, K * K * sizeof(ST_double));
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += AB[i * K + k] * A[j * K + k];  /* A' means A[j,k] = A[j*K + k] */
            }
            V[i * K + j] = sum;
        }
    }

    cqreg_aligned_free(AB);
}

ST_int cqreg_map_clusters(const ST_int *cluster_ids,
                          ST_int N,
                          ST_int *cluster_map,
                          ST_int *num_clusters)
{
    /* Simple approach: find unique values and assign indices */
    /* For large N, could use hash table for O(N) instead of O(N*G) */

    ST_int *unique = NULL;
    ST_int n_unique = 0;
    ST_int i, j;
    ST_int capacity = 1000;

    unique = (ST_int *)malloc(capacity * sizeof(ST_int));
    if (unique == NULL) return -1;

    for (i = 0; i < N; i++) {
        ST_int cid = cluster_ids[i];
        ST_int found = -1;

        /* Search for existing cluster */
        for (j = 0; j < n_unique; j++) {
            if (unique[j] == cid) {
                found = j;
                break;
            }
        }

        if (found < 0) {
            /* New cluster */
            if (n_unique >= capacity) {
                capacity *= 2;
                ST_int *tmp = (ST_int *)realloc(unique, capacity * sizeof(ST_int));
                if (tmp == NULL) {
                    free(unique);
                    return -1;
                }
                unique = tmp;
            }
            unique[n_unique] = cid;
            found = n_unique;
            n_unique++;
        }

        cluster_map[i] = found;
    }

    *num_clusters = n_unique;
    free(unique);

    return 0;
}

/* ============================================================================
 * IID VCE
 * ============================================================================ */

ST_int cqreg_vce_iid(ST_double *V,
                     const ST_double *X,
                     const ST_double *residuals,
                     ST_int N, ST_int K,
                     ST_double q,
                     ST_double sparsity)
{
    ST_int i, j;
    ST_double *XtX_inv = NULL;


    /* Allocate (X'X)^{-1} */
    XtX_inv = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (XtX_inv == NULL) {
        return -1;
    }

    /* Compute (X'X)^{-1} */
    if (cqreg_compute_xtx_inv(XtX_inv, X, N, K) != 0) {
        cqreg_aligned_free(XtX_inv);
        return -1;
    }

    /* V = sparsity^2 * q*(1-q) * (X'X)^{-1}
     * Note: The asymptotic variance formula is:
     * Var(β̂) = τ(1-τ) / f(0)² * (X'X)^{-1}
     * where sparsity = 1/f(0)
     * So V = q*(1-q) * sparsity² * (X'X)^{-1}
     * WITHOUT the (1/n) factor that some sources incorrectly include.
     */
    ST_double scale = sparsity * sparsity * q * (1.0 - q);

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            V[i * K + j] = scale * XtX_inv[i * K + j];
        }
    }

    cqreg_aligned_free(XtX_inv);

    return 0;
}

/* ============================================================================
 * Robust (Sandwich) VCE
 * ============================================================================ */

ST_int cqreg_vce_robust(ST_double *V,
                        const ST_double *X,
                        const ST_double *residuals,
                        ST_int N, ST_int K,
                        ST_double q,
                        ST_double bandwidth)
{
    ST_int i, j, k;
    ST_double *XtX_inv = NULL;
    ST_double *M = NULL;
    ST_int rc = 0;

    /* Allocate matrices */
    XtX_inv = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    M = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);

    if (XtX_inv == NULL || M == NULL) {
        rc = -1;
        goto cleanup;
    }

    /* Compute (X'X)^{-1} */
    if (cqreg_compute_xtx_inv(XtX_inv, X, N, K) != 0) {
        rc = -1;
        goto cleanup;
    }

    /* Compute M = (1/nh) * sum_i K(r_i/h) * X_i * X_i'
     * where K is the kernel function */
    memset(M, 0, K * K * sizeof(ST_double));

    ST_double h_inv = 1.0 / bandwidth;
    ST_double scale = h_inv / (ST_double)N;

    #ifdef _OPENMP
    /* Use thread-local M matrices and reduce */
    #pragma omp parallel
    {
        ST_double *M_local = (ST_double *)calloc(K * K, sizeof(ST_double));

        #pragma omp for schedule(static)
        for (i = 0; i < N; i++) {
            ST_double u = residuals[i] * h_inv;
            ST_double kernel_val = cqreg_kernel_epanechnikov(u);

            if (kernel_val > 0) {
                for (j = 0; j < K; j++) {
                    ST_double Xij = X[j * N + i];
                    for (k = j; k < K; k++) {
                        ST_double Xik = X[k * N + i];
                        ST_double contrib = kernel_val * Xij * Xik;
                        M_local[j * K + k] += contrib;
                        if (k != j) {
                            M_local[k * K + j] += contrib;
                        }
                    }
                }
            }
        }

        #pragma omp critical
        {
            for (j = 0; j < K * K; j++) {
                M[j] += M_local[j];
            }
        }

        free(M_local);
    }
    #else
    for (i = 0; i < N; i++) {
        ST_double u = residuals[i] * h_inv;
        ST_double kernel_val = cqreg_kernel_epanechnikov(u);

        if (kernel_val > 0) {
            for (j = 0; j < K; j++) {
                ST_double Xij = X[j * N + i];
                for (k = j; k < K; k++) {
                    ST_double Xik = X[k * N + i];
                    ST_double contrib = kernel_val * Xij * Xik;
                    M[j * K + k] += contrib;
                    if (k != j) {
                        M[k * K + j] += contrib;
                    }
                }
            }
        }
    }
    #endif

    /* Scale M */
    for (j = 0; j < K * K; j++) {
        M[j] *= scale;
    }

    /* Compute sandwich: V = XtX_inv * M * XtX_inv */
    /* But for QR, the formula is actually:
     * V = q*(1-q) * (X'X)^{-1} * (X' * diag(f_i) * X) * (X'X)^{-1}
     * where f_i is kernel density at residual i
     *
     * Simpler: V = (X'X)^{-1} * M_robust * (X'X)^{-1}
     * where M_robust = q*(1-q) * X' * X (for sandwich form)
     */

    /* For robust VCE, we use:
     * V = (X'X)^{-1} * Omega * (X'X)^{-1}
     * where Omega = X' * W * X
     * and W[i,i] = q*(1-q) (under heteroskedasticity)
     */

    /* Scale M by q*(1-q) for the proper sandwich form */
    ST_double q_scale = q * (1.0 - q);
    for (j = 0; j < K * K; j++) {
        M[j] *= q_scale;
    }

    /* Sandwich product */
    cqreg_sandwich_product(V, XtX_inv, M, K);

cleanup:
    cqreg_aligned_free(XtX_inv);
    cqreg_aligned_free(M);

    return rc;
}

/* ============================================================================
 * Cluster-Robust VCE
 * ============================================================================ */

ST_int cqreg_vce_cluster(ST_double *V,
                         const ST_double *X,
                         const ST_double *residuals,
                         const ST_int *cluster_ids,
                         ST_int num_clusters,
                         ST_int N, ST_int K,
                         ST_double q,
                         ST_double bandwidth)
{
    ST_int i, j, k, g;
    ST_double *XtX_inv = NULL;
    ST_double *M = NULL;
    ST_double *score_g = NULL;
    ST_int *cluster_map = NULL;
    ST_int rc = 0;

    /* Allocate matrices */
    XtX_inv = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    M = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    score_g = (ST_double *)cqreg_aligned_alloc(K * sizeof(ST_double), CQREG_CACHE_LINE);
    cluster_map = (ST_int *)malloc(N * sizeof(ST_int));

    if (XtX_inv == NULL || M == NULL || score_g == NULL || cluster_map == NULL) {
        rc = -1;
        goto cleanup;
    }

    /* Map cluster IDs to 0-indexed */
    ST_int actual_num_clusters = num_clusters;
    if (cqreg_map_clusters(cluster_ids, N, cluster_map, &actual_num_clusters) != 0) {
        rc = -1;
        goto cleanup;
    }

    /* Compute (X'X)^{-1} */
    if (cqreg_compute_xtx_inv(XtX_inv, X, N, K) != 0) {
        rc = -1;
        goto cleanup;
    }

    /* Compute cluster-robust meat:
     * M = sum_g (sum_i in g: psi_i * X_i) * (sum_i in g: psi_i * X_i)'
     *
     * where psi_i = q - I(resid_i < 0) is the influence function score
     */

    memset(M, 0, K * K * sizeof(ST_double));

    /* Allocate storage for cluster scores */
    ST_double **cluster_scores = (ST_double **)calloc(actual_num_clusters, sizeof(ST_double *));
    if (cluster_scores == NULL) {
        rc = -1;
        goto cleanup;
    }

    for (g = 0; g < actual_num_clusters; g++) {
        cluster_scores[g] = (ST_double *)calloc(K, sizeof(ST_double));
        if (cluster_scores[g] == NULL) {
            for (i = 0; i < g; i++) free(cluster_scores[i]);
            free(cluster_scores);
            rc = -1;
            goto cleanup;
        }
    }

    /* Accumulate scores by cluster */
    for (i = 0; i < N; i++) {
        ST_int g_idx = cluster_map[i];
        ST_double psi = (residuals[i] < 0) ? (q - 1.0) : q;

        for (j = 0; j < K; j++) {
            cluster_scores[g_idx][j] += psi * X[j * N + i];
        }
    }

    /* Compute outer products and sum */
    for (g = 0; g < actual_num_clusters; g++) {
        for (j = 0; j < K; j++) {
            for (k = 0; k < K; k++) {
                M[j * K + k] += cluster_scores[g][j] * cluster_scores[g][k];
            }
        }
    }

    /* Free cluster scores */
    for (g = 0; g < actual_num_clusters; g++) {
        free(cluster_scores[g]);
    }
    free(cluster_scores);

    /* Small-sample adjustment:
     * Multiply by G/(G-1) * (N-1)/(N-K)
     * where G = number of clusters
     */
    ST_double G = (ST_double)actual_num_clusters;
    ST_double adj = (G / (G - 1.0)) * ((ST_double)(N - 1) / (ST_double)(N - K));

    for (j = 0; j < K * K; j++) {
        M[j] *= adj;
    }

    /* Sandwich product */
    cqreg_sandwich_product(V, XtX_inv, M, K);

cleanup:
    cqreg_aligned_free(XtX_inv);
    cqreg_aligned_free(M);
    cqreg_aligned_free(score_g);
    free(cluster_map);

    return rc;
}

/* ============================================================================
 * Main Dispatcher
 * ============================================================================ */

ST_int cqreg_compute_vce(cqreg_state *state, const ST_double *X)
{
    if (state == NULL || X == NULL) {
        return -1;
    }

    ST_int rc;

    switch (state->vce_type) {
        case CQREG_VCE_IID:
            rc = cqreg_vce_iid(state->V, X, state->residuals,
                              state->N, state->K,
                              state->quantile, state->sparsity);
            break;

        case CQREG_VCE_ROBUST:
            rc = cqreg_vce_robust(state->V, X, state->residuals,
                                 state->N, state->K,
                                 state->quantile, state->bandwidth);
            break;

        case CQREG_VCE_CLUSTER:
            if (state->cluster_ids == NULL) {
                return -1;
            }
            rc = cqreg_vce_cluster(state->V, X, state->residuals,
                                  state->cluster_ids, state->num_clusters,
                                  state->N, state->K,
                                  state->quantile, state->bandwidth);
            break;

        default:
            rc = -1;
    }

    return rc;
}
