/*
 * creghdfe_vce.c
 *
 * Variance-Covariance Estimation: unadjusted, robust (HC1), clustered
 * Implements reghdfe-compatible VCE formulas using shared ctools_vce engine.
 * Part of the ctools Stata plugin suite
 */

#include <string.h>
#include "creghdfe_vce.h"
#include "../ctools_matrix.h"

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
    ST_int i;
    ST_double sigma2 = rss / df_r;

    for (i = 0; i < K * K; i++) {
        V[i] = sigma2 * inv_xx[i];
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
    /* df_m is K_with_cons - 1 (X vars only, not constant) for consistency with reghdfe */
    ST_int df_m = K_with_cons - 1;
    /* HC1 adjustment: N_eff / (N_eff - df_m - df_a) - matches reghdfe.mata line 3924 */
    ST_double dof_adj = (ST_double)N_eff / (N_eff - df_m - df_a);

    /* X_eff = data + N (skip y column, point to X1...X_K_keep, constant) */
    ctools_vce_data d;
    d.X_eff = data + N;
    d.D = inv_xx;
    d.resid = resid;
    d.weights = weights;
    d.weight_type = weight_type;
    d.N = N;
    d.K = K_with_cons;
    d.normalize_weights = 1;

    ctools_vce_robust(&d, dof_adj, V);
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
    /* Handle degenerate case: single cluster (VCE undefined, like reghdfe) */
    if (num_clusters <= 1) {
        memset(V, 0, K_with_cons * K_with_cons * sizeof(ST_double));
        return;
    }

    /* reghdfe's DOF adjustment: (N-1)/(N - nested_adj - df_m - S.df_a) * M/(M-1)
     * Use N_eff for fweight (sum of weights) */
    ST_int df_a_adjusted = df_a - df_a_nested;
    ST_int nested_adj = (df_a_nested > 0) ? 1 : 0;
    ST_double dof_adj = ((ST_double)(N_eff - 1) / (N_eff - nested_adj - df_m - df_a_adjusted)) * ((ST_double)num_clusters / (num_clusters - 1));

    /* X_eff = data + N (skip y column, point to X1...X_K_keep, constant) */
    ctools_vce_data d;
    d.X_eff = data + N;
    d.D = inv_xx;
    d.resid = resid;
    d.weights = weights;
    d.weight_type = weight_type;
    d.N = N;
    d.K = K_with_cons;
    d.normalize_weights = 1;

    ctools_vce_cluster(&d, cluster_ids, num_clusters, dof_adj, V);
}
