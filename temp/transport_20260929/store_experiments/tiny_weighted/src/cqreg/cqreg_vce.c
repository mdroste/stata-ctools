/*
 * cqreg_vce.c
 *
 * Variance-covariance estimation for quantile regression.
 * Part of the ctools suite.
 */

#include "cqreg_vce.h"
#include "cqreg_linalg.h"
#include "cqreg_blas.h"
#include "cqreg_sparsity.h"
#include "../ctools_config.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>

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

/* Compute X'X (K x K, full symmetric)
 * OPTIMIZED: Uses cqreg_dot and cqreg_dot_self which have 8x unrolling.
 * NOTE: OpenMP is INTENTIONALLY DISABLED here.
 * When enabled, it causes heap corruption during ctools_aligned_free()
 * later in cqreg_compute_xtx_inv. This appears to be an interaction between
 * OpenMP and the Stata plugin's memory allocation (posix_memalign).
 * Since K is typically small (3-20), parallelization provides minimal
 * benefit anyway. Keep this serial for stability.
 */
static void cqreg_compute_xtx(ST_double *XtX, const ST_double *X,
                              ST_int N, ST_int K)
{
    ST_int j, k;

    memset(XtX, 0, (size_t)K * K * sizeof(ST_double));

    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[(size_t)j * N];

        /* Diagonal - cqreg_dot_self already uses 8x unrolling */
        XtX[j * K + j] = cqreg_dot_self(Xj, N);

        /* Off-diagonal - cqreg_dot already uses 8x unrolling */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[(size_t)k * N];
            ST_double dot = cqreg_dot(Xj, Xk, N);
            XtX[j * K + k] = dot;
            XtX[k * K + j] = dot;
        }
    }
}

static ST_int cqreg_compute_xtx_inv(ST_double *XtX_inv,
                                    const ST_double *X,
                                    ST_int N, ST_int K)
{
    ST_double *XtX = NULL;
    ST_double *L = NULL;
    ST_int rc = 0;

    vce_debug_log("cqreg_compute_xtx_inv: N=%d K=%d\n", N, K);

    /* Allocate temporary storage (with overflow check) */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        vce_debug_log("  ERROR: size overflow K*K*sizeof\n");
        return -1;
    }
    XtX = (ST_double *)ctools_cacheline_alloc(kk_size);
    L = (ST_double *)ctools_cacheline_alloc(kk_size);

    vce_debug_log("  XtX=%p, L=%p\n", (void*)XtX, (void*)L);

    if (XtX == NULL || L == NULL) {
        vce_debug_log("  ERROR: allocation failed\n");
        rc = -1;
        goto cleanup;
    }

    vce_debug_log("  Computing X'X...\n");
    cqreg_compute_xtx(XtX, X, N, K);
    vce_debug_log("  X'X computed, diagonal[0]=%.4f\n", XtX[0]);

    /* Cholesky decomposition */
    vce_debug_log("  Cholesky decomposition...\n");
    memcpy(L, XtX, kk_size);
    if (cqreg_cholesky(L, K) != 0) {
        vce_debug_log("  Cholesky failed, trying with regularization...\n");
        /* Try with regularization */
        memcpy(L, XtX, kk_size);
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
    ctools_aligned_free(XtX);
    vce_debug_log("  XtX freed. Now freeing L=%p...\n", (void*)L);
    ctools_aligned_free(L);
    vce_debug_log("  L freed. cqreg_compute_xtx_inv: returning %d\n", rc);

    return rc;
}

static void cqreg_sandwich_product(ST_double *V,
                                   const ST_double *A,
                                   const ST_double *B,
                                   ST_int K)
{
    ST_int i, j, k;
    ST_double *AB = NULL;

    /* Compute size with overflow check */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        memset(V, 0, (size_t)K * (size_t)K * sizeof(ST_double));
        return;
    }

    /* Allocate temporary for A * B */
    AB = (ST_double *)ctools_cacheline_alloc(kk_size);
    if (AB == NULL) {
        memset(V, 0, kk_size);
        return;
    }

    /* Compute AB = A * B */
    memset(AB, 0, kk_size);
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
    memset(V, 0, kk_size);
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += AB[i * K + k] * A[j * K + k];  /* A' means A[j,k] = A[j*K + k] */
            }
            V[i * K + j] = sum;
        }
    }

    ctools_aligned_free(AB);
}

static ST_int cqreg_map_clusters(const ST_int *cluster_ids,
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
                /* Check for overflow in capacity doubling */
                if (capacity > INT32_MAX / 2) {
                    free(unique);
                    return -1;
                }
                ST_int new_capacity = capacity * 2;
                /* Check for overflow in realloc size */
                if ((size_t)new_capacity > SIZE_MAX / sizeof(ST_int)) {
                    free(unique);
                    return -1;
                }
                ST_int *tmp = (ST_int *)realloc(unique, (size_t)new_capacity * sizeof(ST_int));
                if (tmp == NULL) {
                    free(unique);
                    return -1;
                }
                unique = tmp;
                capacity = new_capacity;
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

static ST_int cqreg_vce_iid(ST_double *V,
                            const ST_double *X,
                            const ST_double *residuals,
                            ST_int N, ST_int K,
                            ST_double q,
                            ST_double sparsity)
{
    (void)residuals;  /* Unused - IID VCE uses only sparsity estimate */
    ST_int i, j;
    ST_double *XtX_inv = NULL;

    vce_debug_open();
    vce_debug_log("cqreg_vce_iid: ENTRY N=%d K=%d q=%.4f sparsity=%.4f\n", N, K, q, sparsity);
    vce_debug_log("  V=%p, X=%p\n", (void*)V, (void*)X);

    /* Compute size with overflow check */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        vce_debug_log("  ERROR: size overflow\n");
        vce_debug_close();
        return -1;
    }

    /* Allocate (X'X)^{-1} */
    XtX_inv = (ST_double *)ctools_cacheline_alloc(kk_size);
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
        ctools_aligned_free(XtX_inv);
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
    ctools_aligned_free(XtX_inv);

    vce_debug_log("cqreg_vce_iid: EXIT\n");
    vce_debug_close();
    return 0;
}

/* ============================================================================
 * Robust and Cluster VCE for the Kernel Method (qreg's VCE_kernel)
 * ============================================================================ */

/*
 * Generalized inverse of a symmetric matrix as Stata's invsym(): sweep the
 * pivots in order; a pivot that is not positive or is below 1e-13 times its
 * original diagonal element is skipped and its row and column are zero in
 * the result. With all pivots swept this is the ordinary inverse.
 * A and inv are K x K; they may not overlap.
 */
static ST_int cqreg_invsym(const ST_double *A, ST_int K, ST_double *inv)
{
    ST_int i, j, p;
    unsigned char *swept = (unsigned char *)calloc((size_t)K, 1);
    if (swept == NULL) return -1;

    memcpy(inv, A, (size_t)K * K * sizeof(ST_double));
    for (p = 0; p < K; p++) {
        ST_double d = inv[p * K + p];
        if (!(d > 0.0) || !(d > 1e-13 * A[p * K + p])) continue;
        /* Rank-one update of the other rows and columns */
        for (i = 0; i < K; i++) {
            if (i == p) continue;
            ST_double f = inv[i * K + p] / d;
            if (f == 0.0) continue;
            for (j = 0; j < K; j++) {
                if (j != p) inv[i * K + j] -= f * inv[p * K + j];
            }
        }
        /* Row p scaled by 1/d, column p by -1/d, pivot 1/d */
        for (j = 0; j < K; j++) {
            if (j == p) continue;
            inv[p * K + j] /= d;
            inv[j * K + p] = -inv[j * K + p] / d;
        }
        inv[p * K + p] = 1.0 / d;
        swept[p] = 1;
    }
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            if (!swept[i] || !swept[j]) inv[i * K + j] = 0.0;
        }
    }
    /* Symmetrize the rounding differences of the two triangles */
    for (i = 0; i < K; i++) {
        for (j = i + 1; j < K; j++) {
            ST_double v = 0.5 * (inv[i * K + j] + inv[j * K + i]);
            inv[i * K + j] = v;
            inv[j * K + i] = v;
        }
    }
    free(swept);
    return 0;
}

/*
 * H = invsym(X' diag(kval) X / kbwidth), as qreg's mat accum [iw=kval] and
 * invsym(): the sweep generalized inverse zeroes singular directions (e.g. a
 * group of observations without kernel weight) instead of failing.
 */
static ST_int cqreg_kernel_hessian_inv(ST_double *H, const ST_double *X,
                                       const ST_double *kval,
                                       ST_int N, ST_int K, ST_double kbwidth)
{
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        return -1;
    }
    ST_double *XKX = (ST_double *)ctools_cacheline_alloc(kk_size);
    if (XKX == NULL) return -1;

    blas_xtdx(XKX, X, N, K, kval);
    for (ST_int j = 0; j < K * K; j++) {
        XKX[j] /= kbwidth;
    }
    ST_int rc = cqreg_invsym(XKX, K, H);

    ctools_aligned_free(XKX);
    return rc;
}

/*
 * qreg vce(robust, kernel): V = tau(1-tau) * H * X'X * H with the kernel
 * Hessian H above.
 */
static ST_int cqreg_vce_robust_kernel(ST_double *V,
                                      const ST_double *X,
                                      const ST_double *kval,
                                      ST_int N, ST_int K,
                                      ST_double q,
                                      ST_double kbwidth)
{
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        return -1;
    }
    ST_double *H = (ST_double *)ctools_cacheline_alloc(kk_size);
    ST_double *XtX = (ST_double *)ctools_cacheline_alloc(kk_size);
    ST_int rc = -1;

    if (H == NULL || XtX == NULL) goto cleanup;
    if (cqreg_kernel_hessian_inv(H, X, kval, N, K, kbwidth) != 0) goto cleanup;
    cqreg_compute_xtx(XtX, X, N, K);

    cqreg_sandwich_product(V, H, XtX, K);
    for (ST_int j = 0; j < K * K; j++) {
        V[j] *= q * (1.0 - q);
    }
    rc = 0;

cleanup:
    ctools_aligned_free(H);
    ctools_aligned_free(XtX);
    return rc;
}

/*
 * qreg vce(excluster, kernel) (GetVCE): _robust with scores
 * -tau (r > 0) or 1-tau (r < 0), omitting |r| < 1e-10, and the kernel
 * Hessian H as the bread: V = G/(G-1) * H * (sum_g u_g u_g') * H, where
 * u_g sums score_i * x_i within cluster g and G counts the clusters with a
 * used observation (stored in *num_clusters).
 */
static ST_int cqreg_vce_cluster_kernel(ST_double *V,
                                       const ST_double *X,
                                       const ST_double *residuals,
                                       const ST_double *kval,
                                       const ST_int *cluster_ids,
                                       ST_int *num_clusters,
                                       ST_int N, ST_int K,
                                       ST_double q,
                                       ST_double kbwidth)
{
    ST_int i, j, k, g;
    ST_int rc = -1;
    ST_int *cluster_map = NULL;
    ST_double *H = NULL, *M = NULL, *u = NULL;
    unsigned char *used = NULL;

    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        return -1;
    }
    H = (ST_double *)ctools_cacheline_alloc(kk_size);
    M = (ST_double *)ctools_cacheline_alloc(kk_size);
    cluster_map = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));
    if (H == NULL || M == NULL || cluster_map == NULL) goto cleanup;

    /* Dense 0-based cluster indices: the ado passes egen group() ids 1..G;
     * other ids fall back to the general mapping */
    ST_int G_all = 0;
    for (i = 0; i < N; i++) {
        if (cluster_ids[i] < 1 || cluster_ids[i] > N) break;
        cluster_map[i] = cluster_ids[i] - 1;
        if (cluster_ids[i] > G_all) G_all = cluster_ids[i];
    }
    if (i < N && cqreg_map_clusters(cluster_ids, N, cluster_map, &G_all) != 0) goto cleanup;

    u = (ST_double *)calloc((size_t)G_all * K, sizeof(ST_double));
    used = (unsigned char *)calloc((size_t)G_all, 1);
    if (u == NULL || used == NULL) goto cleanup;

    for (i = 0; i < N; i++) {
        if (fabs(residuals[i]) < 1e-10) continue;
        ST_double score = (residuals[i] > 0) ? -q : 1.0 - q;
        g = cluster_map[i];
        used[g] = 1;
        for (j = 0; j < K; j++) {
            u[(size_t)g * K + j] += score * X[(size_t)j * N + i];
        }
    }

    ST_int G = 0;
    memset(M, 0, kk_size);
    for (g = 0; g < G_all; g++) {
        if (!used[g]) continue;
        G++;
        const ST_double *ug = &u[(size_t)g * K];
        for (j = 0; j < K; j++) {
            for (k = 0; k < K; k++) {
                M[j * K + k] += ug[j] * ug[k];
            }
        }
    }
    *num_clusters = G;
    if (G < 2) goto cleanup;

    if (cqreg_kernel_hessian_inv(H, X, kval, N, K, kbwidth) != 0) goto cleanup;
    cqreg_sandwich_product(V, H, M, K);
    ST_double adj = (ST_double)G / (ST_double)(G - 1);
    for (j = 0; j < K * K; j++) {
        V[j] *= adj;
    }
    rc = 0;

cleanup:
    ctools_aligned_free(H);
    ctools_aligned_free(M);
    free(cluster_map);
    free(u);
    free(used);
    return rc;
}

/* ============================================================================
 * Robust VCE with Per-Observation Densities (Fitted Method)
 * ============================================================================ */

static ST_int cqreg_vce_robust_fitted(ST_double *V,
                                      const ST_double *X,
                                      const ST_double *obs_density,
                                      ST_int N, ST_int K,
                                      ST_double q)
{
    /*
     * Powell sandwich VCE with per-observation densities (Stata's "fitted" method).
     *
     * Formula: V = (X'DX)^{-1} * tau(1-tau) * (X'X) * (X'DX)^{-1}
     *
     * where D = diag(f_i(0)) and f_i(0) is the density estimate at observation i.
     *
     * This uses the asymptotic variance formula from Powell (1991):
     * J = E[f_i(0) X_i X_i']  (approximate Hessian)
     * I = tau(1-tau) E[X_i X_i']  (score variance)
     * V = J^{-1} I J^{-1}
     */

    vce_debug_open();
    vce_debug_log("cqreg_vce_robust_fitted: ENTRY N=%d K=%d q=%.4f\n", N, K, q);

    ST_int i, j, k;
    ST_double *XDX = NULL;      /* X' D X where D = diag(obs_density) */
    ST_double *XDX_inv = NULL;  /* (X'DX)^{-1} */
    ST_double *XtX = NULL;      /* X' X */
    ST_double *L = NULL;        /* Cholesky factor */
    ST_int rc = 0;

    /* Compute size with overflow check */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        vce_debug_log("  ERROR: size overflow\n");
        return -1;
    }

    /* Allocate matrices */
    XDX = (ST_double *)ctools_cacheline_alloc(kk_size);
    XDX_inv = (ST_double *)ctools_cacheline_alloc(kk_size);
    XtX = (ST_double *)ctools_cacheline_alloc(kk_size);
    L = (ST_double *)ctools_cacheline_alloc(kk_size);

    if (XDX == NULL || XDX_inv == NULL || XtX == NULL || L == NULL) {
        vce_debug_log("  ERROR: allocation failed\n");
        rc = -1;
        goto cleanup;
    }

    /* Compute X'DX = sum_i f_i * X_i * X_i'
     * OPTIMIZED: Column-major traversal with 8x loop unrolling */
    memset(XDX, 0, kk_size);
    ST_int N8 = N - (N & 7);

    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[(size_t)j * N];

        /* Diagonal element with 8x unrolling */
        ST_double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
        ST_double d4 = 0.0, d5 = 0.0, d6 = 0.0, d7 = 0.0;
        for (i = 0; i < N8; i += 8) {
            d0 += obs_density[i]     * Xj[i]     * Xj[i];
            d1 += obs_density[i + 1] * Xj[i + 1] * Xj[i + 1];
            d2 += obs_density[i + 2] * Xj[i + 2] * Xj[i + 2];
            d3 += obs_density[i + 3] * Xj[i + 3] * Xj[i + 3];
            d4 += obs_density[i + 4] * Xj[i + 4] * Xj[i + 4];
            d5 += obs_density[i + 5] * Xj[i + 5] * Xj[i + 5];
            d6 += obs_density[i + 6] * Xj[i + 6] * Xj[i + 6];
            d7 += obs_density[i + 7] * Xj[i + 7] * Xj[i + 7];
        }
        for (; i < N; i++) {
            d0 += obs_density[i] * Xj[i] * Xj[i];
        }
        XDX[j * K + j] = ((d0 + d4) + (d1 + d5)) + ((d2 + d6) + (d3 + d7));

        /* Off-diagonal elements with 8x unrolling */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[(size_t)k * N];
            ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
            for (i = 0; i < N8; i += 8) {
                s0 += obs_density[i]     * Xj[i]     * Xk[i];
                s1 += obs_density[i + 1] * Xj[i + 1] * Xk[i + 1];
                s2 += obs_density[i + 2] * Xj[i + 2] * Xk[i + 2];
                s3 += obs_density[i + 3] * Xj[i + 3] * Xk[i + 3];
                s4 += obs_density[i + 4] * Xj[i + 4] * Xk[i + 4];
                s5 += obs_density[i + 5] * Xj[i + 5] * Xk[i + 5];
                s6 += obs_density[i + 6] * Xj[i + 6] * Xk[i + 6];
                s7 += obs_density[i + 7] * Xj[i + 7] * Xk[i + 7];
            }
            for (; i < N; i++) {
                s0 += obs_density[i] * Xj[i] * Xk[i];
            }
            ST_double total = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
            XDX[j * K + k] = total;
            XDX[k * K + j] = total;
        }
    }

    vce_debug_log("  XDX[0,0] = %.6e\n", XDX[0]);

    /* Compute X'X
     * OPTIMIZED: 8x loop unrolling */
    memset(XtX, 0, kk_size);

    for (j = 0; j < K; j++) {
        const ST_double *Xj = &X[(size_t)j * N];

        /* Diagonal with 8x unrolling */
        ST_double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
        ST_double d4 = 0.0, d5 = 0.0, d6 = 0.0, d7 = 0.0;
        for (i = 0; i < N8; i += 8) {
            d0 += Xj[i]     * Xj[i];
            d1 += Xj[i + 1] * Xj[i + 1];
            d2 += Xj[i + 2] * Xj[i + 2];
            d3 += Xj[i + 3] * Xj[i + 3];
            d4 += Xj[i + 4] * Xj[i + 4];
            d5 += Xj[i + 5] * Xj[i + 5];
            d6 += Xj[i + 6] * Xj[i + 6];
            d7 += Xj[i + 7] * Xj[i + 7];
        }
        for (; i < N; i++) {
            d0 += Xj[i] * Xj[i];
        }
        XtX[j * K + j] = ((d0 + d4) + (d1 + d5)) + ((d2 + d6) + (d3 + d7));

        /* Off-diagonal with 8x unrolling */
        for (k = j + 1; k < K; k++) {
            const ST_double *Xk = &X[(size_t)k * N];
            ST_double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            ST_double s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
            for (i = 0; i < N8; i += 8) {
                s0 += Xj[i]     * Xk[i];
                s1 += Xj[i + 1] * Xk[i + 1];
                s2 += Xj[i + 2] * Xk[i + 2];
                s3 += Xj[i + 3] * Xk[i + 3];
                s4 += Xj[i + 4] * Xk[i + 4];
                s5 += Xj[i + 5] * Xk[i + 5];
                s6 += Xj[i + 6] * Xk[i + 6];
                s7 += Xj[i + 7] * Xk[i + 7];
            }
            for (; i < N; i++) {
                s0 += Xj[i] * Xk[i];
            }
            ST_double total = ((s0 + s4) + (s1 + s5)) + ((s2 + s6) + (s3 + s7));
            XtX[j * K + k] = total;
            XtX[k * K + j] = total;
        }
    }

    vce_debug_log("  XtX[0,0] = %.6e\n", XtX[0]);

    /* Compute (X'DX)^{-1} using Cholesky */
    memcpy(L, XDX, kk_size);

    if (cqreg_cholesky(L, K) != 0) {
        /* Try again with regularization for numerical stability */
        vce_debug_log("  Cholesky failed, retrying with regularization...\n");
        memcpy(L, XDX, kk_size);
        for (j = 0; j < K; j++) {
            L[j * K + j] += 1e-10;
        }
        if (cqreg_cholesky(L, K) != 0) {
            vce_debug_log("  ERROR: Cholesky of XDX failed even with regularization\n");
            rc = -1;
            goto cleanup;
        }
    }

    cqreg_invert_cholesky(XDX_inv, L, K);

    vce_debug_log("  XDX_inv[0,0] = %.6e\n", XDX_inv[0]);

    /* Compute V = tau(1-tau) * (X'DX)^{-1} * X'X * (X'DX)^{-1} */
    ST_double scale = q * (1.0 - q);
    vce_debug_log("  scale (q*(1-q)) = %.6e\n", scale);

    /* temp = X'X * (X'DX)^{-1} */
    ST_double *temp = (ST_double *)ctools_cacheline_alloc(kk_size);
    if (temp == NULL) {
        rc = -1;
        goto cleanup;
    }

    memset(temp, 0, kk_size);
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += XtX[i * K + k] * XDX_inv[k * K + j];
            }
            temp[i * K + j] = sum;
        }
    }

    /* V = scale * (X'DX)^{-1} * temp */
    memset(V, 0, kk_size);
    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            ST_double sum = 0.0;
            for (k = 0; k < K; k++) {
                sum += XDX_inv[i * K + k] * temp[k * K + j];
            }
            V[i * K + j] = scale * sum;
        }
    }

    ctools_aligned_free(temp);

    vce_debug_log("  V[0,0] = %.6e, SE[0] = %.4f\n", V[0], sqrt(V[0]));

cleanup:
    vce_debug_log("cqreg_vce_robust_fitted: EXIT rc=%d\n", rc);
    vce_debug_close();
    ctools_aligned_free(XDX);
    ctools_aligned_free(XDX_inv);
    ctools_aligned_free(XtX);
    ctools_aligned_free(L);

    return rc;
}

/* ============================================================================
 * Cluster-Robust VCE
 * ============================================================================ */

static ST_int cqreg_vce_cluster(ST_double *V,
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

    /* Compute size with overflow check */
    size_t kk_size;
    if (ctools_safe_alloc_size((size_t)K, (size_t)K, sizeof(ST_double), &kk_size) != 0) {
        return -1;
    }

    /* Allocate matrices */
    XtX_inv = (ST_double *)ctools_cacheline_alloc(kk_size);
    M = (ST_double *)ctools_cacheline_alloc(kk_size);
    score_g = (ST_double *)ctools_safe_malloc2((size_t)K, sizeof(ST_double));
    cluster_map = (ST_int *)ctools_safe_malloc2((size_t)N, sizeof(ST_int));

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

    memset(M, 0, kk_size);

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
            cluster_scores[g_idx][j] += psi * X[(size_t)j * N + i];
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
    ctools_aligned_free(XtX_inv);
    ctools_aligned_free(M);
    free(score_g);
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
            /*
             * IID VCE: Always uses scalar sparsity, just different estimation methods.
             * - Residual: difference quotient on residuals (qreg vce(iid, residual))
             * - Fitted: Siddiqui difference quotient at mean(X) (qreg default)
             * - Kernel: kband / mean(K(r/kband)) (qreg vce(iid, kernel()))
             * All use the simple IID formula: V = q(1-q) * s^2 * (X'X)^{-1}
             */
            rc = cqreg_vce_iid(state->V, X, state->residuals,
                              state->N, state->K,
                              state->quantile, state->sparsity);
            break;

        case CQREG_VCE_ROBUST:
            /*
             * Robust VCE: Powell sandwich with per-observation density weights
             * (fitted) or kernel weights (kernel). qreg defines no residual
             * robust VCE; the ado and plugin reject it (r(184)).
             */
            if (state->density_method == CQREG_DENSITY_FITTED) {
                /* Fitted method: use per-observation densities */
                rc = cqreg_vce_robust_fitted(state->V, X, state->obs_density,
                                            state->N, state->K,
                                            state->quantile);
            } else if (state->density_method == CQREG_DENSITY_KERNEL) {
                /* Kernel method: obs_density holds the kernel weights */
                rc = cqreg_vce_robust_kernel(state->V, X, state->obs_density,
                                            state->N, state->K,
                                            state->quantile, state->kbwidth);
            } else {
                rc = -1;
            }
            break;

        case CQREG_VCE_CLUSTER:
            if (state->cluster_ids == NULL) {
                return -1;
            }
            if (state->density_method == CQREG_DENSITY_KERNEL) {
                /* qreg's cluster VCE with the kernel Hessian */
                rc = cqreg_vce_cluster_kernel(state->V, X, state->residuals,
                                             state->obs_density,
                                             state->cluster_ids,
                                             &state->num_clusters,
                                             state->N, state->K,
                                             state->quantile, state->kbwidth);
                break;
            }
            /* Cluster VCE uses scalar sparsity */
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
