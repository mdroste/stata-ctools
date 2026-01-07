/*
 * creghdfe_vce.c
 *
 * Variance-Covariance Estimation: unadjusted, robust (HC1), clustered
 * Implements reghdfe-compatible VCE formulas
 * Part of the ctools Stata plugin suite
 */

#include "creghdfe_vce.h"

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
    ST_int N,
    ST_int K_with_cons,       /* Number of X vars + constant (excludes y) */
    ST_int df_a,              /* Degrees of freedom absorbed by FEs */
    ST_double *V              /* Output: K_with_cons x K_with_cons */
)
{
    ST_int i, j, k, idx;
    ST_double *XeeX = NULL;  /* X'WX = sum_i (e_i^2 * x_i * x_i') including constant */
    ST_double *temp = NULL;
    ST_double e_i, xi_j, xi_k;

    /* df_m is K_with_cons - 1 (X vars only, not constant) for consistency with reghdfe */
    ST_int df_m = K_with_cons - 1;
    /* HC1 adjustment: N / (N - df_m - df_a) - matches reghdfe.mata line 3924 */
    ST_double dof_adj = (ST_double)N / (N - df_m - df_a);

    XeeX = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    temp = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));

    if (!XeeX || !temp) {
        if (XeeX) free(XeeX);
        if (temp) free(temp);
        return;
    }

    /* Compute X'WX = sum_i (e_i^2 * x_i * x_i') where x_i includes constant
     * data layout: column 0 = y, columns 1..K_keep = X vars, column K_with_cons = constant
     * So X columns in data are at indices 1 to K_with_cons (inclusive) */
    for (idx = 0; idx < N; idx++) {
        e_i = resid[idx] * resid[idx];
        for (j = 0; j < K_with_cons; j++) {
            /* X vars are at data columns 1 to K_keep, constant at column K_with_cons */
            xi_j = data[(j + 1) * N + idx];
            for (k = 0; k < K_with_cons; k++) {
                xi_k = data[(k + 1) * N + idx];
                XeeX[j * K_with_cons + k] += e_i * xi_j * xi_k;
            }
        }
    }

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
    const ST_int *cluster_ids,  /* N x 1 cluster IDs (0-indexed) */
    ST_int N,
    ST_int K_with_cons,         /* Number of X vars + constant (excludes y) */
    ST_int num_clusters,
    ST_double *V,               /* Output: K_with_cons x K_with_cons */
    ST_int df_m,                /* Degrees of freedom model (number of X vars, excluding constant) */
    ST_int df_a,                /* Degrees of freedom absorbed by FEs */
    ST_int df_a_nested          /* Degrees of freedom nested within cluster (for adjustment) */
)
{
    ST_int i, j, k, idx, c;
    ST_double *XeeX = NULL;      /* Accumulates sum_c X_c'e_c e_c'X_c */
    ST_double *temp = NULL;
    ST_double *all_ecX = NULL;   /* Per-cluster e'X sums (includes constant) */

    /* reghdfe's DOF adjustment: (N-1)/(N - nested_adj - df_m - S.df_a) * M/(M-1)
     * where:
     *   - nested_adj is 0 or 1 (a flag to handle edge case, not the count)
     *   - S.df_a in reghdfe is ALREADY adjusted (= df_a_initial - df_a_redundant)
     *     where df_a_redundant includes the nested levels
     *   - Our df_a parameter is the UNADJUSTED value from the C code
     * So we need: df_a_adjusted = df_a - df_a_nested (to match reghdfe's S.df_a)
     * See reghdfe.mata lines 4021-4026 and the dof computation around 4265 */
    ST_int df_a_adjusted = df_a - df_a_nested;  /* Match reghdfe's S.df_a */
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;  /* Binary flag, not count */
    ST_double dof_adj = ((ST_double)(N - 1) / (N - nested_adj - df_m - df_a_adjusted)) * ((ST_double)num_clusters / (num_clusters - 1));

    XeeX = (ST_double *)calloc(K_with_cons * K_with_cons, sizeof(ST_double));
    temp = (ST_double *)malloc(K_with_cons * K_with_cons * sizeof(ST_double));
    all_ecX = (ST_double *)calloc(num_clusters * K_with_cons, sizeof(ST_double));

    if (!XeeX || !temp || !all_ecX) {
        if (XeeX) free(XeeX);
        if (temp) free(temp);
        if (all_ecX) free(all_ecX);
        return;
    }

    /* Accumulate e'X for each cluster (including constant column)
     * data layout: column 0 = y, columns 1..K_keep = X vars, column K_with_cons = constant
     * Matches reghdfe.mata line 3956: w = sol.resid :* w (unweighted: w = resid)
     * and meat computation at lines 4000-4014 */
    for (idx = 0; idx < N; idx++) {
        c = cluster_ids[idx];
        for (k = 0; k < K_with_cons; k++) {
            all_ecX[c * K_with_cons + k] += resid[idx] * data[(k + 1) * N + idx];
        }
    }

    /* Compute meat: M = sum_c (e'X)_c (e'X)_c' - parallelized with reduction */
    #pragma omp parallel for schedule(static) reduction(+:XeeX[:K_with_cons*K_with_cons])
    for (c = 0; c < num_clusters; c++) {
        ST_int jj, kk;
        for (jj = 0; jj < K_with_cons; jj++) {
            for (kk = 0; kk < K_with_cons; kk++) {
                XeeX[jj * K_with_cons + kk] += all_ecX[c * K_with_cons + jj] * all_ecX[c * K_with_cons + kk];
            }
        }
    }

    free(all_ecX);

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
