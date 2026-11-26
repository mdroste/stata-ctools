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
#include "ctools_ols.h"

/* ============================================================================
 * Matrix Multiplication
 * ============================================================================ */

/*
 * Cross products C = A' diag(w) B (w == NULL: unit weights).
 *
 * Every element is one fixed sum, whatever the thread count and whichever of
 * the two loop arrangements below computes it. The rows are cut into blocks
 * of CTOOLS_CROSS_BLOCK; a block's products (w[i] * a[i]) * b[i] are added by
 * cross_dot(); the blocks of a chunk are added in order, and the chunk sums
 * are added in order. Where the chunks begin depends only on N, so an element
 * depends only on its two columns, the weights and N.
 *
 * One pass over the rows serves every pair of columns: the columns of a block
 * stay in the cache while all their products are taken, where one pass over
 * the whole arrays per pair reads the data K1 * K2 times from memory.
 */
#define CTOOLS_CROSS_BLOCK 512
/* A chunk is at least this many blocks, and there are at most
 * CTOOLS_CROSS_CHUNKS of them */
#define CTOOLS_CROSS_CHUNK_BLOCKS 8
#define CTOOLS_CROSS_CHUNKS 64
/* Chunk sums kept at once by the chunk arrangement (doubles): beyond this
 * the column arrangement runs, which keeps none */
#define CTOOLS_CROSS_SUMS ((size_t)8 << 20)
/* Columns of A per pass of the column arrangement */
#define CTOOLS_CROSS_STRIP 256
/* Element products below which the work stays on the calling thread */
#define CTOOLS_CROSS_PARALLEL 400000

/* Eight sums in the lane order ctools_simd_dot() uses, so that the compiler
 * can vectorize the loop without reassociating anything */
static inline ST_double cross_dot(const ST_double * restrict x, const ST_double * restrict y,
                                  ST_int n)
{
    ST_double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0;
    ST_int i = 0;
    for (; i + 8 <= n; i += 8) {
        s0 += x[i] * y[i];         s1 += x[i + 1] * y[i + 1];
        s2 += x[i + 2] * y[i + 2]; s3 += x[i + 3] * y[i + 3];
        s4 += x[i + 4] * y[i + 4]; s5 += x[i + 5] * y[i + 5];
        s6 += x[i + 6] * y[i + 6]; s7 += x[i + 7] * y[i + 7];
    }
    ST_double s = ((s0 + s4) + (s2 + s6)) + ((s1 + s5) + (s3 + s7));
    for (; i < n; i++) s += x[i] * y[i];
    return s;
}

/* cross_dot() of (w[i] * x[i]) and y: the value cross_dot(wx, y, n) has when
 * wx[i] = w[i] * x[i] was stored first */
static inline ST_double cross_dot_weighted(const ST_double * restrict w,
                                           const ST_double * restrict x,
                                           const ST_double * restrict y, ST_int n)
{
    ST_double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0;
    ST_int i = 0;
    for (; i + 8 <= n; i += 8) {
        s0 += (w[i] * x[i]) * y[i];             s1 += (w[i + 1] * x[i + 1]) * y[i + 1];
        s2 += (w[i + 2] * x[i + 2]) * y[i + 2]; s3 += (w[i + 3] * x[i + 3]) * y[i + 3];
        s4 += (w[i + 4] * x[i + 4]) * y[i + 4]; s5 += (w[i + 5] * x[i + 5]) * y[i + 5];
        s6 += (w[i + 6] * x[i + 6]) * y[i + 6]; s7 += (w[i + 7] * x[i + 7]) * y[i + 7];
    }
    ST_double s = ((s0 + s4) + (s2 + s6)) + ((s1 + s5) + (s3 + s7));
    for (; i < n; i++) s += (w[i] * x[i]) * y[i];
    return s;
}

/* Blocks [first, last) of chunk c out of nch, for nblk blocks in all */
static inline void cross_chunk_span(ST_int nblk, ST_int nch, ST_int c, ST_int *first, ST_int *last)
{
    *first = (ST_int)((long long)nblk * c / nch);
    *last = (ST_int)((long long)nblk * (c + 1) / nch);
}

/* The sums of blocks [first, last) for every pair of columns, into acc
 * (K1 x K2; with sym only the elements i <= j) */
static void cross_chunk(const ST_double *A, const ST_double *B, const ST_double *w, ST_int N,
                        ST_int K1, ST_int K2, int sym, ST_int first, ST_int last,
                        ST_double * restrict acc)
{
    ST_double wa[CTOOLS_CROSS_BLOCK];

    memset(acc, 0, (size_t)K1 * (size_t)K2 * sizeof(ST_double));
    for (ST_int blk = first; blk < last; blk++) {
        const ST_int i0 = blk * CTOOLS_CROSS_BLOCK;
        const ST_int len = (N - i0 < CTOOLS_CROSS_BLOCK) ? N - i0 : CTOOLS_CROSS_BLOCK;
        for (ST_int i = 0; i < K1; i++) {
            const ST_double *a = A + (size_t)i * N + i0;
            if (w) {
                const ST_double *wi = w + i0;
                for (ST_int t = 0; t < len; t++) wa[t] = wi[t] * a[t];
                a = wa;
            }
            for (ST_int j = sym ? i : 0; j < K2; j++) {
                acc[(size_t)j * K1 + i] += cross_dot(a, B + (size_t)j * N + i0, len);
            }
        }
    }
}

/* Column j of C by the column arrangement: no scratch beyond the stack, the
 * columns of A taken CTOOLS_CROSS_STRIP at a time */
static void cross_column(const ST_double *A, const ST_double *B, const ST_double *w, ST_int N,
                         ST_int K1, int sym, ST_int nblk, ST_int nch, ST_int j, ST_double *C)
{
    ST_double part[CTOOLS_CROSS_STRIP];
    const ST_int ilim = sym ? j + 1 : K1;
    const ST_double *b = B + (size_t)j * N;
    ST_double *c_col = C + (size_t)j * K1;

    for (ST_int s0 = 0; s0 < ilim; s0 += CTOOLS_CROSS_STRIP) {
        const ST_int sn = (ilim - s0 < CTOOLS_CROSS_STRIP) ? ilim - s0 : CTOOLS_CROSS_STRIP;
        for (ST_int i = 0; i < sn; i++) c_col[s0 + i] = 0.0;
        for (ST_int c = 0; c < nch; c++) {
            ST_int first, last;
            cross_chunk_span(nblk, nch, c, &first, &last);
            for (ST_int i = 0; i < sn; i++) part[i] = 0.0;
            for (ST_int blk = first; blk < last; blk++) {
                const ST_int i0 = blk * CTOOLS_CROSS_BLOCK;
                const ST_int len = (N - i0 < CTOOLS_CROSS_BLOCK) ? N - i0 : CTOOLS_CROSS_BLOCK;
                for (ST_int i = 0; i < sn; i++) {
                    const ST_double *a = A + (size_t)(s0 + i) * N + i0;
                    part[i] += w ? cross_dot_weighted(w + i0, a, b + i0, len)
                                 : cross_dot(a, b + i0, len);
                }
            }
            for (ST_int i = 0; i < sn; i++) c_col[s0 + i] += part[i];
        }
    }
}

/* C = A' diag(w) B. With sym (A and B the same K1 = K2 columns) the elements
 * i <= j are computed and mirrored. */
static void cross_product(const ST_double *A, const ST_double *B, const ST_double *w, ST_int N,
                          ST_int K1, ST_int K2, int sym, ST_double *C)
{
    const size_t cells = (size_t)K1 * (size_t)K2;
    if (K1 <= 0 || K2 <= 0) return;
    if (N < 0) N = 0;

    /* N + CTOOLS_CROSS_BLOCK - 1 would overflow int for N near INT_MAX */
    const ST_int nblk = N / CTOOLS_CROSS_BLOCK + (N % CTOOLS_CROSS_BLOCK != 0);
    ST_int nch = nblk / CTOOLS_CROSS_CHUNK_BLOCKS;
    if (nch > CTOOLS_CROSS_CHUNKS) nch = CTOOLS_CROSS_CHUNKS;
    if (nch < 1) nch = 1;
    const int threaded = (double)N * (double)cells >= (double)CTOOLS_CROSS_PARALLEL;
    ST_double *sums = NULL;
    (void)threaded; /* read only by the OpenMP clauses */

    if (nch == 1) {
        cross_chunk(A, B, w, N, K1, K2, sym, 0, nblk, C);
    } else if (cells <= CTOOLS_CROSS_SUMS / (size_t)nch &&
               (sums = (ST_double *)malloc((size_t)nch * cells * sizeof(ST_double))) != NULL) {
        /* Chunks in parallel, each into its own sums */
        ST_int c;
        #pragma omp parallel for schedule(dynamic, 1) if(threaded)
        for (c = 0; c < nch; c++) {
            ST_int first, last;
            cross_chunk_span(nblk, nch, c, &first, &last);
            cross_chunk(A, B, w, N, K1, K2, sym, first, last, sums + (size_t)c * cells);
        }
        for (size_t e = 0; e < cells; e++) {
            ST_double total = 0.0;
            for (c = 0; c < nch; c++) total += sums[(size_t)c * cells + e];
            C[e] = total;
        }
        free(sums);
    } else {
        /* Many elements, or no memory for the chunk sums: the columns of C in
         * parallel, each adding its chunks in order (the same sums) */
        ST_int j;
        #pragma omp parallel for schedule(dynamic, 1) if(threaded)
        for (j = 0; j < K2; j++) {
            cross_column(A, B, w, N, K1, sym, nblk, nch, j, C);
        }
    }

    if (sym) {
        for (ST_int j = 0; j < K2; j++) {
            for (ST_int i = 0; i < j; i++) {
                C[(size_t)i * K1 + j] = C[(size_t)j * K1 + i];
            }
        }
    }
}

/*
 * Matrix multiply C = A' * B
 * A is N x K1, B is N x K2, Result C is K1 x K2
 */
void ctools_matmul_atb(const ST_double * restrict A, const ST_double * restrict B,
                       ST_int N, ST_int K1, ST_int K2,
                       ST_double * restrict C)
{
    cross_product(A, B, NULL, N, K1, K2, A == B && K1 == K2, C);
}

/*
 * Gram matrix C = A' * A (upper triangle computed, then mirrored)
 */
void ctools_matmul_ata(const ST_double *A, ST_int N, ST_int K, ST_double *C)
{
    cross_product(A, A, NULL, N, K, K, 1, C);
}

/* Rows per block of the row-blocked kernels: a block of a few columns
 * stays in the L2 cache while the K passes over it run */
#define CTOOLS_ROW_BLOCK 2048

/*
 * C = X * B (and y += X * b): row blocks in parallel, columns in order
 */
static void matmul_xb_blocks(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                             ST_double *C, int zero_first, const ST_double *w)
{
    /* N + CTOOLS_ROW_BLOCK - 1 would overflow int for N within a block of
     * INT_MAX, leaving every row uncomputed */
    const ST_int num_blocks = N / CTOOLS_ROW_BLOCK + (N % CTOOLS_ROW_BLOCK != 0);
    ST_int blk;

    #pragma omp parallel for schedule(static) if((size_t)N * (size_t)K * (size_t)M > 400000)
    for (blk = 0; blk < num_blocks; blk++) {
        const ST_int i0 = blk * CTOOLS_ROW_BLOCK;
        const ST_int i1 = (N - i0 < CTOOLS_ROW_BLOCK) ? N : i0 + CTOOLS_ROW_BLOCK;
        /* with w, each row of the block is multiplied by sqrt(w) once its
         * sums are complete, while the block is still in the cache */
        ST_double sw[CTOOLS_ROW_BLOCK];
        if (w) {
            for (ST_int i = i0; i < i1; i++) sw[i - i0] = sqrt(w[i]);
        }
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
            if (w) {
                for (ST_int i = i0; i < i1; i++) c[i] = sw[i - i0] * c[i];
            }
        }
    }
}

void ctools_matvec_add(const ST_double *X, ST_int N, ST_int K, const ST_double *b, ST_double *y)
{
    matmul_xb_blocks(X, N, K, b, 1, y, 0, NULL);
}

void ctools_matmul_xb(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                      ST_double *C)
{
    matmul_xb_blocks(X, N, K, B, M, C, 1, NULL);
}

void ctools_matmul_xb_sqrtw(const ST_double *X, ST_int N, ST_int K, const ST_double *B, ST_int M,
                            const ST_double *w, ST_double *C)
{
    matmul_xb_blocks(X, N, K, B, M, C, 1, w);
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
 * Weighted matrix multiply C = A' * diag(w) * B (w == NULL: unit weights)
 */
void ctools_matmul_atdb(const ST_double * restrict A, const ST_double * restrict B,
                        const ST_double * restrict w, ST_int N, ST_int K1, ST_int K2,
                        ST_double * restrict C)
{
    cross_product(A, B, w, N, K1, K2, A == B && K1 == K2, C);
}

void ctools_matmul_atdb_full(const ST_double *A, const ST_double *B, const ST_double *w,
                             ST_int N, ST_int K1, ST_int K2, ST_double *C)
{
    cross_product(A, B, w, N, K1, K2, 0, C);
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

    /* Normalize weights for aweight/pweight (OLS convention): w * scale,
     * rounded to a double, formed where it is used */
    const int w_norm = d->normalize_weights && d->weights != NULL &&
                       (d->weight_type == 1 || d->weight_type == 3);
    ST_double scale = 1.0;
    if (w_norm) {
        ST_double sum_w = 0.0;
        for (idx = 0; idx < N; idx++) sum_w += d->weights[idx];
        if (!isfinite(sum_w) || sum_w <= 0) return 498;
        scale = (ST_double)N / sum_w;
    }

    /* Allocate high-precision meat and temp matrices */
    long double *meat = (long double *)calloc((size_t)K * K, sizeof(long double));
    long double *temp = (long double *)malloc((size_t)K * K * sizeof(long double));

    if (!meat || !temp) {
        if (meat) free(meat);
        if (temp) free(temp);
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
        return 920;
    }
    /* one value per row, rows in parallel */
    #pragma omp parallel for schedule(static) if(N > 200000)
    for (idx = 0; idx < N; idx++) {
        ST_double w_i;
        ST_double e = d->resid[idx];

        if (d->weights == NULL) {
            w_i = e * e;
        } else if (d->weight_type == 2) {
            w_i = e * e * d->weights[idx];
        } else if (w_norm) {
            const ST_double wn = d->weights[idx] * scale;
            ST_double ew = e * wn;
            w_i = ew * ew;
        } else {
            /* Raw weights (IV path): w * e^2 */
            w_i = d->weights[idx] * e * e;
        }
        obs_w[idx] = w_i;
    }

    /* meat[j,k] = sum_i (w_i * x_ij) * x_ik with compensated
     * (twice-working-precision) products and sums, as Mata's quadcross()
     * accumulates reghdfe's meat: with nearly collinear regressors D*meat*D
     * cancels most of the digits, and long double is no wider than double on
     * some platforms. Every element of the upper triangle is an independent
     * sum; ctools_kahan_cross_weighted() takes them all in one blocked pass
     * over the rows, in parallel. */
    ctools_kahan_cross_weighted(d->X_eff, obs_w, N, K, meat);
    free(obs_w);
    /* Mirror upper triangle to lower */
    for (j = 0; j < K; j++) {
        for (k = j + 1; k < K; k++) {
            meat[(size_t)k * K + j] = meat[(size_t)j * K + k];
        }
    }

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

/* Cluster sums (long doubles) the clustered VCE may hold at once: within
 * this, and with at most a quarter as many clusters as rows, every cluster's
 * e'X is formed in one pass per regressor over the rows; otherwise the
 * observations are sorted by cluster and the clusters taken in blocks */
#ifndef CTOOLS_VCE_CLUSTER_SUMS
#define CTOOLS_VCE_CLUSTER_SUMS ((size_t)16 << 20)
#endif

/* meat (upper triangle) += sums * sums' for one cluster's e'X */
static inline void vce_meat_add(long double *meat, const long double *sums, ST_int K)
{
    for (ST_int jj = 0; jj < K; jj++) {
        long double ecXj = sums[jj];
        meat[(size_t)jj * K + jj] += ecXj * ecXj;
        for (ST_int kk = jj + 1; kk < K; kk++) {
            meat[(size_t)jj * K + kk] += ecXj * sums[kk];
        }
    }
}

/* The clustered meat with the observations sorted by cluster: bounded scratch
 * whatever the number of clusters. Returns 0 or 920. */
static ST_retcode vce_cluster_meat_sorted(const ctools_vce_data *d, const ST_int *cluster_ids,
                                          ST_int num_clusters, const ST_double *e_w,
                                          long double *meat)
{
    const ST_int N = d->N, K = d->K;
    ST_int k;
    ST_int *sort_perm = (ST_int *)malloc((size_t)N * sizeof(ST_int));
    ST_int *boundaries = (ST_int *)malloc(((size_t)num_clusters + 1) * sizeof(ST_int));
    /* Clusters are processed in blocks: the e'X sums of a block are formed
     * one regressor per thread (each sum over its cluster's observations in
     * sorted order), then the outer products are added in cluster order. */
    ST_int block_clusters = 65536 / (K > 0 ? K : 1);
    if (block_clusters < 1) block_clusters = 1;
    if (block_clusters > num_clusters) block_clusters = num_clusters;
    long double *block_sums =
        (long double *)malloc((size_t)block_clusters * (size_t)K * sizeof(long double));

    /* Sort observations by cluster using counting sort */
    if (!sort_perm || !boundaries || !block_sums ||
        vce_sort_by_cluster(cluster_ids, N, num_clusters, sort_perm, boundaries) != 0) {
        free(sort_perm);
        free(boundaries);
        free(block_sums);
        return 920;
    }

    /* c0 advances to c1, which never passes num_clusters (c0 += block_clusters
     * could overflow int on the last block) */
    for (ST_int c0 = 0, c1; c0 < num_clusters; c0 = c1) {
        c1 = (num_clusters - c0 < block_clusters) ? num_clusters : c0 + block_clusters;
        #ifdef _OPENMP
        const size_t block_obs = (size_t)(boundaries[c1] - boundaries[c0]);
        #pragma omp parallel for schedule(static) if(K > 1 && block_obs * (size_t)K > 100000)
        #endif
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
            vce_meat_add(meat, block_sums + (size_t)(c - c0) * K, K);
        }
    }
    free(sort_perm);
    free(boundaries);
    free(block_sums);
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

    /* Normalize weights for aweight/pweight (OLS convention): w * scale,
     * rounded to a double, formed where it is used */
    const int w_norm = d->normalize_weights && d->weights != NULL &&
                       (d->weight_type == 1 || d->weight_type == 3);
    ST_double scale = 1.0;
    if (w_norm) {
        ST_double sum_w = 0.0;
        for (idx = 0; idx < N; idx++) sum_w += d->weights[idx];
        if (!isfinite(sum_w) || sum_w <= 0) return 498;
        scale = (ST_double)N / sum_w;
    }

    /* Allocate high-precision buffers */
    long double *meat = (long double *)calloc((size_t)K * K, sizeof(long double));
    long double *temp = (long double *)malloc((size_t)K * K * sizeof(long double));
    long double *ecX = (long double *)malloc(K * sizeof(long double));
    /* Weighted residual of every observation */
    ST_double *e_w = (ST_double *)malloc((size_t)N * sizeof(ST_double));

    if (!meat || !temp || !ecX || !e_w) {
        free(meat);
        free(temp);
        free(ecX);
        free(e_w);
        return 920;
    }
    /* one value per row, rows in parallel */
    #pragma omp parallel for schedule(static) if(N > 200000)
    for (idx = 0; idx < N; idx++) {
        if (d->weights == NULL) {
            e_w[idx] = d->resid[idx];
        } else if (d->weight_type == 2) {
            e_w[idx] = d->resid[idx] * d->weights[idx];
        } else if (w_norm) {
            const ST_double wn = d->weights[idx] * scale;
            e_w[idx] = d->resid[idx] * wn;
        } else {
            e_w[idx] = d->resid[idx] * d->weights[idx];
        }
    }

    /* e'X of every cluster: each element sums its cluster's observations in
     * row order. When all the sums fit, a regressor's pass over the rows adds
     * every observation to its cluster's sum (the regressors in parallel):
     * the rows are read in order and nothing is sorted. Otherwise the
     * observations are sorted by cluster and the clusters taken in blocks;
     * the sums, and the order in which their outer products enter the meat
     * (cluster order), are the same. */
    ST_retcode meat_rc = 0;
    long double *all_sums = NULL;
    /* The sums fit when there are few of them next to the data: at most a
     * quarter as many clusters as rows (the sums then take no more than half
     * the memory of X_eff), and CTOOLS_VCE_CLUSTER_SUMS in all */
    const int sums_fit = num_clusters <= N / 4 &&
                         (size_t)num_clusters <= CTOOLS_VCE_CLUSTER_SUMS / (size_t)K;
    if (sums_fit) {
        all_sums = (long double *)malloc((size_t)num_clusters * (size_t)K * sizeof(long double));
        if (!all_sums) meat_rc = 920;
    }
    if (all_sums) {
        #pragma omp parallel for schedule(static) if(K > 1 && (size_t)N * (size_t)K > 100000)
        for (k = 0; k < K; k++) {
            const ST_double *xk = d->X_eff + (size_t)k * N;
            long double *sums_k = all_sums + (size_t)k * num_clusters;
            for (ST_int c = 0; c < num_clusters; c++) sums_k[c] = 0.0L;
            /* A run of rows of one cluster (data sorted by it) is added in a
             * register: the same additions, without a store and a reload of
             * the sum per row */
            for (ST_int obs = 0; obs < N;) {
                const ST_int c = cluster_ids[obs];
                long double acc = sums_k[c];
                do {
                    acc += (long double)e_w[obs] * (long double)xk[obs];
                    obs++;
                } while (obs < N && cluster_ids[obs] == c);
                sums_k[c] = acc;
            }
        }
        for (ST_int c = 0; c < num_clusters; c++) {
            for (k = 0; k < K; k++) ecX[k] = all_sums[(size_t)k * num_clusters + c];
            vce_meat_add(meat, ecX, K);
        }
        free(all_sums);
    } else if (!sums_fit) {
        meat_rc = vce_cluster_meat_sorted(d, cluster_ids, num_clusters, e_w, meat);
    }
    free(e_w);
    free(ecX);
    if (meat_rc) {
        free(meat);
        free(temp);
        return meat_rc;
    }

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
