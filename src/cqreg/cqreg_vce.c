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
#include <stdarg.h>

/* Debug logging */
#define VCE_DEBUG 0

#if VCE_DEBUG
static FILE *vce_debug_file = NULL;

static void vce_debug_open(void) {
    if (vce_debug_file == NULL) {
        vce_debug_file = fopen("/tmp/cqreg_vce_debug.log", "a");
    }
}

static void vce_debug_log(const char *fmt, ...) {
    if (vce_debug_file) {
        va_list args;
        va_start(args, fmt);
        vfprintf(vce_debug_file, fmt, args);
        va_end(args);
        fflush(vce_debug_file);
    }
}

static void vce_debug_close(void) {
    if (vce_debug_file) {
        fclose(vce_debug_file);
        vce_debug_file = NULL;
    }
}
#else
#define vce_debug_open()
#define vce_debug_log(...)
#define vce_debug_close()
#endif

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

    vce_debug_log("cqreg_compute_xtx_inv: N=%d K=%d\n", N, K);

    /* Allocate temporary storage */
    XtX = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    L = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);

    vce_debug_log("  XtX=%p, L=%p\n", (void*)XtX, (void*)L);

    if (XtX == NULL || L == NULL) {
        vce_debug_log("  ERROR: allocation failed\n");
        rc = -1;
        goto cleanup;
    }

    vce_debug_log("  Computing X'X...\n");
    /* Compute X'X */
    memset(XtX, 0, K * K * sizeof(ST_double));

    /* NOTE: OpenMP is INTENTIONALLY DISABLED here.
     * When enabled, it causes heap corruption during cqreg_aligned_free()
     * later in this function. This appears to be an interaction between
     * OpenMP and the Stata plugin's memory allocation (posix_memalign).
     * Since K is typically small (3-20), parallelization provides minimal
     * benefit anyway. Keep this serial for stability.
     */
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
    vce_debug_log("  X'X computed, diagonal[0]=%.4f\n", XtX[0]);

    /* Cholesky decomposition */
    vce_debug_log("  Cholesky decomposition...\n");
    memcpy(L, XtX, K * K * sizeof(ST_double));
    if (cqreg_cholesky(L, K) != 0) {
        vce_debug_log("  Cholesky failed, trying with regularization...\n");
        /* Try with regularization */
        memcpy(L, XtX, K * K * sizeof(ST_double));
        cqreg_add_regularization(L, K, 1e-10);
        if (cqreg_cholesky(L, K) != 0) {
            vce_debug_log("  ERROR: Cholesky still failed\n");
            rc = -1;
            goto cleanup;
        }
    }
    vce_debug_log("  Cholesky done, L[0]=%.4f\n", L[0]);

    /* Compute inverse using Cholesky factor */
    vce_debug_log("  Inverting via Cholesky...\n");
    cqreg_invert_cholesky(XtX_inv, L, K);
    vce_debug_log("  Inverse computed, XtX_inv[0]=%.6e\n", XtX_inv[0]);

cleanup:
    vce_debug_log("  Cleanup: freeing XtX=%p...\n", (void*)XtX);
    cqreg_aligned_free(XtX);
    vce_debug_log("  XtX freed. Now freeing L=%p...\n", (void*)L);
    cqreg_aligned_free(L);
    vce_debug_log("  L freed. cqreg_compute_xtx_inv: returning %d\n", rc);

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

    vce_debug_open();
    vce_debug_log("cqreg_vce_iid: ENTRY N=%d K=%d q=%.4f sparsity=%.4f\n", N, K, q, sparsity);
    vce_debug_log("  V=%p, X=%p\n", (void*)V, (void*)X);

    /* Allocate (X'X)^{-1} */
    XtX_inv = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    vce_debug_log("  XtX_inv=%p\n", (void*)XtX_inv);
    if (XtX_inv == NULL) {
        vce_debug_log("  ERROR: XtX_inv alloc failed\n");
        vce_debug_close();
        return -1;
    }

    vce_debug_log("  Calling cqreg_compute_xtx_inv...\n");
    /* Compute (X'X)^{-1} */
    if (cqreg_compute_xtx_inv(XtX_inv, X, N, K) != 0) {
        vce_debug_log("  ERROR: cqreg_compute_xtx_inv failed\n");
        cqreg_aligned_free(XtX_inv);
        vce_debug_close();
        return -1;
    }
    vce_debug_log("  cqreg_compute_xtx_inv returned successfully\n");

    /* V = sparsity^2 * q*(1-q) * (X'X)^{-1}
     * Note: The asymptotic variance formula is:
     * Var(β̂) = τ(1-τ) / f(0)² * (X'X)^{-1}
     * where sparsity = 1/f(0)
     * So V = q*(1-q) * sparsity² * (X'X)^{-1}
     * WITHOUT the (1/n) factor that some sources incorrectly include.
     */
    ST_double scale = sparsity * sparsity * q * (1.0 - q);
    vce_debug_log("  scale=%.6e, computing V...\n", scale);

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            V[i * K + j] = scale * XtX_inv[i * K + j];
        }
    }
    vce_debug_log("  V computed, V[0]=%.6e\n", V[0]);

    vce_debug_log("  Freeing XtX_inv...\n");
    cqreg_aligned_free(XtX_inv);

    vce_debug_log("cqreg_vce_iid: EXIT\n");
    vce_debug_close();
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
                        ST_double bandwidth,
                        ST_double sparsity)
{
    /*
     * Robust sandwich VCE for quantile regression.
     *
     * Formula: V = (X'DX)^{-1} * Omega * (X'DX)^{-1}
     *
     * where:
     *   D = diag(f_i) with f_i = density estimate at residual i
     *   Omega = sum_i psi_i^2 * X_i * X_i'
     *   psi_i = tau - I(r_i < 0)  (influence function)
     *
     * We use the simplified approach with a constant density estimate:
     *   D = (1/sparsity) * I_n
     *
     * This gives: V = sparsity^2 * (X'X)^{-1} * Omega * (X'X)^{-1}
     *
     * Note: This may differ from Stata's qreg vce(robust) which uses
     * a different "fitted" density estimation method.
     */

    vce_debug_open();
    vce_debug_log("cqreg_vce_robust: ENTRY N=%d K=%d q=%.4f bw=%.6f sparsity=%.4f\n",
                  N, K, q, bandwidth, sparsity);

    ST_int i, j, k;
    ST_double *XtX_inv = NULL;  /* (X'X)^{-1} */
    ST_double *Omega = NULL;    /* Score variance matrix */
    ST_int rc = 0;

    /* Allocate matrices */
    XtX_inv = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    Omega = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);

    if (XtX_inv == NULL || Omega == NULL) {
        vce_debug_log("  ERROR: allocation failed\n");
        rc = -1;
        goto cleanup;
    }

    /* Compute (X'X)^{-1} */
    if (cqreg_compute_xtx_inv(XtX_inv, X, N, K) != 0) {
        vce_debug_log("  ERROR: cqreg_compute_xtx_inv failed\n");
        rc = -1;
        goto cleanup;
    }

    /*
     * Compute Omega = sum_i psi_i^2 * X_i * X_i'
     * where psi_i = tau - I(r_i < 0)
     *
     * For r_i >= 0: psi_i = tau, psi_i^2 = tau^2
     * For r_i < 0:  psi_i = tau - 1, psi_i^2 = (1-tau)^2
     */
    memset(Omega, 0, K * K * sizeof(ST_double));

    ST_double psi_pos_sq = q * q;           /* tau^2 */
    ST_double psi_neg_sq = (1.0 - q) * (1.0 - q);  /* (1-tau)^2 */

    for (i = 0; i < N; i++) {
        ST_double psi_sq = (residuals[i] >= 0) ? psi_pos_sq : psi_neg_sq;

        for (j = 0; j < K; j++) {
            ST_double Xij = X[j * N + i];
            for (k = j; k < K; k++) {
                ST_double Xik = X[k * N + i];
                ST_double contrib = psi_sq * Xij * Xik;
                Omega[j * K + k] += contrib;
                if (k != j) {
                    Omega[k * K + j] += contrib;
                }
            }
        }
    }

    vce_debug_log("  Omega[0,0] = %.6e\n", Omega[0]);

    /*
     * Compute V = sparsity^2 * (X'X)^{-1} * Omega * (X'X)^{-1}
     */
    ST_double scale = sparsity * sparsity;
    vce_debug_log("  sparsity^2 = %.6e\n", scale);

    /* First compute temp = Omega * (X'X)^{-1} */
    ST_double *temp = (ST_double *)cqreg_aligned_alloc(K * K * sizeof(ST_double), CQREG_CACHE_LINE);
    if (temp == NULL) {
        vce_debug_log("  ERROR: temp allocation failed\n");
        rc = -1;
        goto cleanup;
    }

    memset(temp, 0, K * K * sizeof(ST_double));
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += Omega[i * K + k] * XtX_inv[k * K + j];
            }
            temp[i * K + j] = sum;
        }
    }

    /* Now compute V = scale * (X'X)^{-1} * temp */
    memset(V, 0, K * K * sizeof(ST_double));
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += XtX_inv[i * K + k] * temp[k * K + j];
            }
            V[i * K + j] = scale * sum;
        }
    }

    cqreg_aligned_free(temp);

    vce_debug_log("  V[0,0] = %.6e, SE[0] = %.4f\n", V[0], sqrt(V[0]));

cleanup:
    vce_debug_log("cqreg_vce_robust: EXIT rc=%d\n", rc);
    vce_debug_close();
    cqreg_aligned_free(XtX_inv);
    cqreg_aligned_free(Omega);

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
                         ST_double sparsity)
{
    /*
     * Cluster-robust VCE for quantile regression.
     *
     * Formula: V = sparsity^2 * (X'X)^{-1} * M * (X'X)^{-1}
     *
     * where M = sum_g (sum_i in g: psi_i * X_i)(sum_i in g: psi_i * X_i)'
     * and psi_i = tau - I(r_i < 0) is the influence function score.
     *
     * The sparsity^2 factor converts from the score-space variance
     * to the coefficient-space variance.
     */

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

    /* Sandwich product: temp = (X'X)^{-1} * M * (X'X)^{-1} */
    cqreg_sandwich_product(V, XtX_inv, M, K);

    /* Scale by sparsity^2 to convert to coefficient variance */
    ST_double sparsity_sq = sparsity * sparsity;
    for (j = 0; j < K * K; j++) {
        V[j] *= sparsity_sq;
    }

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
                                 state->quantile, state->bandwidth,
                                 state->sparsity);
            break;

        case CQREG_VCE_CLUSTER:
            if (state->cluster_ids == NULL) {
                return -1;
            }
            rc = cqreg_vce_cluster(state->V, X, state->residuals,
                                  state->cluster_ids, state->num_clusters,
                                  state->N, state->K,
                                  state->quantile, state->sparsity);
            break;

        default:
            rc = -1;
    }

    return rc;
}
