/*
 * creghdfe_vce.c
 *
 * Variance-Covariance Estimation: unadjusted, robust (HC1), clustered
 * Implements reghdfe-compatible VCE formulas
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_vce.h"
#include "creghdfe_utils.h"  /* For sort_by_cluster */

/* ========================================================================
 * Compute unadjusted VCE: V = sigma^2 * (X'X)^(-1)
 * where sigma^2 = RSS / df_r
 * ======================================================================== */

void compute_vce_unadjusted(
    const ST_double *inv_xx,  /* (X'X)^(-1), K x K */
    ST_double rss,
    ST_int df_r,
    ST_int K,
    ST_double *V              /* Output: K x K VCE matrix */
)
{
    ST_int i, j;
    ST_double sigma2 = rss / df_r;

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            V[i * K + j] = sigma2 * inv_xx[i * K + j];
        }
    }
}

/* ========================================================================
 * Compute robust VCE (HC1): V = D * M * D * dof_adj
 * where D = (X'X)^(-1), M = X'WX with W = diag(e^2), and X includes constant
 * This matches reghdfe.mata lines 3929-3930
 * ======================================================================== */

void compute_vce_robust(
    const ST_double *data,    /* N x (K_keep+2) matrix: y, X1...X_K_keep, constant */
    const ST_double *resid,   /* N x 1 pre-computed residuals from partialled X */
    const ST_double *inv_xx,  /* K_with_cons x K_with_cons inverse (X vars + constant) */
    const ST_double *weights, /* N x 1 weights (NULL for unweighted) */
    ST_int weight_type,       /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_int N,
    ST_int N_eff,             /* Effective N (= sum(weights) for fweight, else N) */
    ST_int K_with_cons,       /* Number of X vars + constant (excludes y) */
    ST_int df_a,              /* Degrees of freedom absorbed by FEs */
    ST_double *V              /* Output: K_with_cons x K_with_cons */
)
{
    ST_int i, j, k, idx;
    ST_double *XeeX = NULL;  /* X'WX = sum_i (e_i^2 * x_i * x_i') including constant */
    ST_double *temp = NULL;
    ST_double xi_j, xi_k;

    /* df_m is K_with_cons - 1 (X vars only, not constant) for consistency with reghdfe */
    ST_int df_m = K_with_cons - 1;
    /* HC1 adjustment: N_eff / (N_eff - df_m - df_a) - matches reghdfe.mata line 3924 */
    ST_double dof_adj = (ST_double)N_eff / (N_eff - df_m - df_a);

    /* For aweight/pweight: normalize weights to sum to N (reghdfe.mata line 3598) */
    ST_double *w_norm = NULL;
    if (weights != NULL && (weight_type == 1 || weight_type == 3)) {
        ST_double sum_w = 0.0;
        w_norm = (ST_double *)malloc(N * sizeof(ST_double));
        if (!w_norm) return;
        for (idx = 0; idx < N; idx++) sum_w += weights[idx];
        ST_double scale = (ST_double)N / sum_w;
        for (idx = 0; idx < N; idx++) w_norm[idx] = weights[idx] * scale;
    }

    XeeX = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    temp = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));

    if (!XeeX || !temp) {
        if (XeeX) free(XeeX);
        if (temp) free(temp);
        if (w_norm) free(w_norm);
        return;
    }

    /* Compute X'WX where W depends on weight type (reghdfe.mata lines 3914-3922):
     * - No weights:     w_i = e_i^2
     * - fweight:        w_i = e_i^2 * fw_i
     * - aweight/pweight: w_i = (e_i * w_norm_i)^2 */
    for (idx = 0; idx < N; idx++) {
        ST_double w_i;
        ST_double e = resid[idx];

        if (weights == NULL) {
            w_i = e * e;
        } else if (weight_type == 2) {
            /* fweight: e^2 * w */
            w_i = e * e * weights[idx];
        } else {
            /* aweight or pweight: (e * w_norm)^2 */
            ST_double ew = e * w_norm[idx];
            w_i = ew * ew;
        }

        for (j = 0; j < K_with_cons; j++) {
            xi_j = data[(j + 1) * N + idx];
            for (k = 0; k < K_with_cons; k++) {
                xi_k = data[(k + 1) * N + idx];
                XeeX[j * K_with_cons + k] += w_i * xi_j * xi_k;
            }
        }
    }

    if (w_norm) free(w_norm);

    /* Compute temp = D * M = (X'X)^(-1) * X'WX */
    for (i = 0; i < K_with_cons; i++) {
        for (j = 0; j < K_with_cons; j++) {
            temp[i * K_with_cons + j] = 0.0;
            for (k = 0; k < K_with_cons; k++) {
                temp[i * K_with_cons + j] += inv_xx[i * K_with_cons + k] * XeeX[k * K_with_cons + j];
            }
        }
    }

    /* Compute V = temp * D * dof_adj = D * M * D * dof_adj */
    for (i = 0; i < K_with_cons; i++) {
        for (j = 0; j < K_with_cons; j++) {
            V[i * K_with_cons + j] = 0.0;
            for (k = 0; k < K_with_cons; k++) {
                V[i * K_with_cons + j] += temp[i * K_with_cons + k] * inv_xx[k * K_with_cons + j];
            }
            V[i * K_with_cons + j] *= dof_adj;
        }
    }

    free(XeeX);
    free(temp);
}

/* ========================================================================
 * Compute clustered VCE using reghdfe's formula:
 * V = D * M * D * dof_adj
 * where D = (X'X)^(-1), M = sum_c (X_c'e_c)(e_c'X_c), X includes constant
 * dof_adj = (N-1)/(N - nested_adj - df_m - df_a) * M/(M-1)
 * where nested_adj = 1 if df_a_nested > 0, else 0
 * This matches reghdfe.mata lines 4021, 4026, 4031
 * ======================================================================== */

void compute_vce_cluster(
    const ST_double *data,      /* N x (K_keep+2) matrix: y, X1...X_K_keep, constant */
    const ST_double *resid,     /* N x 1 pre-computed residuals from partialled X */
    const ST_double *inv_xx,    /* K_with_cons x K_with_cons inverse (X vars + constant) */
    const ST_double *weights,   /* N x 1 weights (NULL for unweighted) */
    ST_int weight_type,         /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    const ST_int *cluster_ids,  /* N x 1 cluster IDs (0-indexed) */
    ST_int N,
    ST_int N_eff,               /* Effective N (= sum(weights) for fweight, else N) */
    ST_int K_with_cons,         /* Number of X vars + constant (excludes y) */
    ST_int num_clusters,
    ST_double *V,               /* Output: K_with_cons x K_with_cons */
    ST_int df_m,                /* Degrees of freedom model (number of X vars, excluding constant) */
    ST_int df_a,                /* Degrees of freedom absorbed by FEs */
    ST_int df_a_nested          /* Degrees of freedom nested within cluster (for adjustment) */
)
{
    ST_int i, j, k, idx;
    ST_double *XeeX = NULL;      /* Accumulates sum_c X_c'e_c e_c'X_c */
    ST_double *temp = NULL;
    ST_double *ecX = NULL;       /* Single cluster e'X buffer */
    ST_int *sort_perm = NULL;    /* Permutation to sort by cluster */
    ST_int *boundaries = NULL;   /* Cluster boundaries in sorted order */

    /* reghdfe's DOF adjustment: (N-1)/(N - nested_adj - df_m - S.df_a) * M/(M-1)
     * Use N_eff for fweight (sum of weights) */
    ST_int df_a_adjusted = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (N_eff - nested_adj - df_m - df_a_adjusted)) * ((ST_double)num_clusters / (num_clusters - 1));

    /* For aweight/pweight: normalize weights to sum to N (reghdfe.mata line 3598) */
    ST_double *w_norm = NULL;
    if (weights != NULL && (weight_type == 1 || weight_type == 3)) {
        ST_double sum_w = 0.0;
        w_norm = (ST_double *)malloc(N * sizeof(ST_double));
        if (!w_norm) return;
        for (idx = 0; idx < N; idx++) sum_w += weights[idx];
        ST_double scale = (ST_double)N / sum_w;
        for (idx = 0; idx < N; idx++) w_norm[idx] = weights[idx] * scale;
    }

    /* Allocate buffers - note: no all_ecX[num_clusters * K] needed!
     * We only need a single ecX[K] buffer and sort/stream by cluster. */
    XeeX = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    temp = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));
    ecX = (ST_double *)malloc(K_with_cons * sizeof(ST_double));
    sort_perm = (ST_int *)malloc(N * sizeof(ST_int));
    boundaries = (ST_int *)malloc((num_clusters + 1) * sizeof(ST_int));

    if (!XeeX || !temp || !ecX || !sort_perm || !boundaries) {
        if (XeeX) free(XeeX);
        if (temp) free(temp);
        if (ecX) free(ecX);
        if (sort_perm) free(sort_perm);
        if (boundaries) free(boundaries);
        if (w_norm) free(w_norm);
        return;
    }

    /* Sort observations by cluster using counting sort - O(N + num_clusters) */
    if (sort_by_cluster(cluster_ids, N, num_clusters, sort_perm, boundaries) != 0) {
        free(XeeX); free(temp); free(ecX);
        free(sort_perm); free(boundaries);
        if (w_norm) free(w_norm);
        return;
    }

    /* Stream through clusters, accumulating ecX for each cluster then updating XeeX.
     * This avoids allocating O(num_clusters * K) memory. */
    for (ST_int c = 0; c < num_clusters; c++) {
        ST_int start = boundaries[c];
        ST_int end = boundaries[c + 1];

        /* Reset ecX buffer for this cluster */
        for (k = 0; k < K_with_cons; k++) {
            ecX[k] = 0.0;
        }

        /* Accumulate e'X for this cluster */
        for (ST_int pos = start; pos < end; pos++) {
            idx = sort_perm[pos];  /* Original observation index */
            ST_double e_w;

            if (weights == NULL) {
                e_w = resid[idx];
            } else if (weight_type == 2) {
                /* fweight: e * w */
                e_w = resid[idx] * weights[idx];
            } else {
                /* aweight or pweight: e * w_norm */
                e_w = resid[idx] * w_norm[idx];
            }

            for (k = 0; k < K_with_cons; k++) {
                ecX[k] += e_w * data[(k + 1) * N + idx];
            }
        }

        /* Accumulate outer product ecX * ecX' into XeeX */
        for (ST_int jj = 0; jj < K_with_cons; jj++) {
            for (ST_int kk = 0; kk < K_with_cons; kk++) {
                XeeX[jj * K_with_cons + kk] += ecX[jj] * ecX[kk];
            }
        }
    }

    free(ecX);
    free(sort_perm);
    free(boundaries);
    if (w_norm) free(w_norm);

    /* Compute temp = D * M = (X'X)^(-1) * XeeX */
    for (i = 0; i < K_with_cons; i++) {
        for (j = 0; j < K_with_cons; j++) {
            temp[i * K_with_cons + j] = 0.0;
            for (k = 0; k < K_with_cons; k++) {
                temp[i * K_with_cons + j] += inv_xx[i * K_with_cons + k] * XeeX[k * K_with_cons + j];
            }
        }
    }

    /* Compute V = temp * D * dof_adj = D * M * D * dof_adj */
    for (i = 0; i < K_with_cons; i++) {
        for (j = 0; j < K_with_cons; j++) {
            V[i * K_with_cons + j] = 0.0;
            for (k = 0; k < K_with_cons; k++) {
                V[i * K_with_cons + j] += temp[i * K_with_cons + k] * inv_xx[k * K_with_cons + j];
            }
            V[i * K_with_cons + j] *= dof_adj;
        }
    }

    free(XeeX);
    free(temp);
}
