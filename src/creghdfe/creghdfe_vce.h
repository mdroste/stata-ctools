/*
 * creghdfe_vce.h
 *
 * Variance-Covariance Estimation: unadjusted, robust (HC1), clustered
 * Part of the ctools Stata plugin suite
 */

#ifndef CREGHDFE_VCE_H
#define CREGHDFE_VCE_H

#include "creghdfe_types.h"

/* ========================================================================
 * VCE Computation Functions
 * ======================================================================== */

/*
 * Compute unadjusted VCE: V = sigma^2 * (X'X)^(-1)
 * where sigma^2 = RSS / df_r
 */
ST_retcode compute_vce_unadjusted(
    const ST_double *inv_xx,  /* (X'X)^(-1), K x K */
    ST_double rss,
    ST_double df_r,
    ST_int K,
    ST_double *V              /* Output: K x K VCE matrix */
);

/*
 * Compute robust VCE (HC1): V = D * M * D * dof_adj
 * where D = (X'X)^(-1), M = X'WX, and X includes constant
 * Weight handling (matches reghdfe.mata lines 3914-3922):
 * - No weights:     W = diag(e^2)
 * - fweight:        W = diag(e^2 * w)
 * - aweight/pweight: W = diag((e * w_norm)^2) where w_norm is normalized to sum to N
 */
ST_retcode compute_vce_robust(
    const ST_double *data,    /* N x (K_keep+2) matrix: y, X1...X_K_keep, constant */
    const ST_double *resid,   /* N x 1 pre-computed residuals from partialled X */
    const ST_double *inv_xx,  /* K_with_cons x K_with_cons inverse (X vars + constant) */
    const ST_double *weights, /* N x 1 weights (NULL for unweighted) */
    ST_int weight_type,       /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    ST_int N,
    ST_double N_eff,             /* Effective N (= sum(weights) for fweight, else N) */
    ST_int K_with_cons,       /* Number of X vars + constant (excludes y) */
    ST_double df_a,              /* Degrees of freedom absorbed by FEs */
    ST_double *V              /* Output: K_with_cons x K_with_cons */
);

/*
 * Compute clustered VCE using reghdfe's formula:
 * V = D * M * D * dof_adj
 * where D = (X'X)^(-1), M = sum_c (X_c'e_c)(e_c'X_c), X includes constant
 * dof_adj = (N-1)/(N - nested_adj - df_m - df_a) * M/(M-1)
 * Weight handling (matches reghdfe.mata line 3956):
 * - Cluster sums use: e * w_norm (where w_norm is normalized to sum to N for aw/pw)
 */
ST_retcode compute_vce_cluster(
    const ST_double *data,      /* N x (K_keep+2) matrix: y, X1...X_K_keep, constant */
    const ST_double *resid,     /* N x 1 pre-computed residuals from partialled X */
    const ST_double *inv_xx,    /* K_with_cons x K_with_cons inverse (X vars + constant) */
    const ST_double *weights,   /* N x 1 weights (NULL for unweighted) */
    ST_int weight_type,         /* 0=none, 1=aweight, 2=fweight, 3=pweight */
    const ST_int *cluster_ids,  /* N x 1 cluster IDs (0-indexed) */
    ST_int N,
    ST_double N_eff,               /* Effective N (= sum(weights) for fweight, else N) */
    ST_int K_with_cons,         /* Number of X vars + constant (excludes y) */
    ST_int num_clusters,
    ST_double *V,               /* Output: K_with_cons x K_with_cons */
    ST_int df_m,                /* Degrees of freedom model (number of X vars, excluding constant) */
    ST_double df_a,                /* Degrees of freedom absorbed by FEs */
    ST_double df_a_nested          /* Degrees of freedom nested within cluster (for adjustment) */
);

/*
 * Multi-way clustered VCE (reghdfe_vce_cluster with Q > 1 cluster
 * variables): the sum over every intersection of the cluster variables of
 * its clustered sandwich, with sign + for an odd and - for an even number of
 * variables, one small-sample factor with the fewest clusters, and
 * reghdfe_fix_psd() (*psd_fixed tells whether it changed V). cluster_ids
 * holds Q columns of N dense 0-based IDs, num_clusters their counts. When
 * reghdfe's corrected V is indefinite beyond rounding, V_psd (K_with_cons
 * square, may be NULL) gets a positive semidefinite version with the same
 * coefficient block and *psd_made is 1 (see fix_psd in creghdfe_vce.c).
 */
ST_retcode compute_vce_multicluster(
    const ST_double *data, const ST_double *resid, const ST_double *inv_xx,
    const ST_double *weights, ST_int weight_type, const ST_int *cluster_ids,
    const ST_int *num_clusters, ST_int Q, ST_int N, ST_double N_eff, ST_int K_with_cons,
    ST_double *V, ST_int df_m, ST_double df_a, ST_double df_a_nested, ST_int *psd_fixed,
    ST_double *V_psd, ST_int *psd_made);

/*
 * Driscoll-Kraay VCE (reghdfe_vce_dkraay): the moments summed within each
 * time period (time_ids: dense 0-based in time order, T periods), a Bartlett
 * kernel over bandwidth - 1 lags (bandwidth 0: reghdfe's default; at most T
 * - 1 lags, *bandwidth_used = lags + 1), the factor (N-1)/(N-df_m-df_a) *
 * T/(T-1) with df_a net of nested effects, and reghdfe_fix_psd() (V_psd and
 * *psd_made as for compute_vce_multicluster).
 */
ST_retcode compute_vce_dkraay(
    const ST_double *data, const ST_double *resid, const ST_double *inv_xx,
    const ST_double *weights, ST_int weight_type, const ST_int *time_ids, ST_int T,
    ST_int bandwidth, ST_int N, ST_double N_eff, ST_int K_with_cons, ST_double *V,
    ST_int df_m, ST_double df_a, ST_int *bandwidth_used, ST_double *V_psd, ST_int *psd_made);

#endif /* CREGHDFE_VCE_H */
