/*
    civreghdfe_tests.h
    Diagnostic test statistics for IV regression

    Implements test statistics for instrument validity and endogeneity:
    - First-stage F statistics
    - Anderson canonical correlation / Kleibergen-Paap rk LM (underidentification)
    - Cragg-Donald / Kleibergen-Paap rk Wald F (weak instruments)
    - Sargan / Hansen J (overidentification)
    - Durbin-Wu-Hausman (endogeneity)
*/

#ifndef CIVREGHDFE_TESTS_H
#define CIVREGHDFE_TESTS_H

#include "../stplugin.h"

/*
    Structure to hold all IV diagnostic test results.
    Populated by civreghdfe_compute_diagnostics().
*/
typedef struct {
    /* First-stage statistics */
    ST_double *first_stage_F;   /* Array of K_endog F-statistics */

    /* Underidentification test (rank condition) */
    ST_double underid_stat;     /* Anderson LM or KP rk LM statistic (chi-sq) */
    ST_int underid_df;          /* Degrees of freedom */

    /* Weak instruments test */
    ST_double cd_f;             /* Cragg-Donald Wald F (homoskedastic) */
    ST_double kp_f;             /* Kleibergen-Paap rk Wald F (robust) */

    /* Overidentification test */
    ST_double sargan_stat;      /* Sargan/Hansen J statistic (chi-sq) */
    ST_int overid_df;           /* Degrees of freedom */

    /* Endogeneity test (Durbin-Wu-Hausman) */
    ST_double endog_chi2;       /* DWH chi-squared statistic */
    ST_double endog_f;          /* DWH F statistic */
    ST_int endog_df;            /* Degrees of freedom */
} civreghdfe_diagnostics_t;

/*
    Compute first-stage F statistics for each endogenous variable.

    Parameters:
    - X_endog: Endogenous regressors (N x K_endog, column-major)
    - Z: All instruments including exogenous (N x K_iv, column-major)
    - ZtZ_inv: Precomputed (Z'Z)^-1 (K_iv x K_iv)
    - ZtX: Precomputed Z'X for all X (K_iv x K_total)
    - weights: Optional weights (may be NULL)
    - weight_type: 0=none, 1=aweight, 2=fweight, 3=pweight
    - N, K_exog, K_endog, K_iv: Dimensions
    - df_a: Absorbed degrees of freedom
    - first_stage_F: Output array (K_endog elements, must be preallocated)
*/
void civreghdfe_compute_first_stage_F(
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *ZtZ_inv,
    const ST_double *ZtX,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_double *first_stage_F
);

/*
    Compute underidentification test (Anderson LM / Kleibergen-Paap rk LM).

    For homoskedastic VCE: Anderson canonical correlation LM
    For robust/cluster VCE: Kleibergen-Paap rk LM

    Parameters:
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv)
    - ZtZ: Precomputed Z'Z (K_iv x K_iv)
    - ZtZ_inv: Precomputed (Z'Z)^-1 (K_iv x K_iv)
    - temp1: First-stage coefficients (K_iv x K_total)
    - first_stage_F: First-stage F statistics (K_endog)
    - weights, weight_type, N, K_exog, K_endog, K_iv, df_a: Parameters
    - vce_type: 0=unadjusted, 1=robust, 2=cluster, 3=dkraay
    - cluster_ids: Cluster identifiers (may be NULL)
    - num_clusters: Number of clusters
    - kernel_type: HAC kernel type (0=none, 1=Bartlett, etc.)
    - bw: HAC bandwidth
    - kiefer: 1 if Kiefer VCE, 0 otherwise
    - hac_panel_ids: Panel IDs for panel-aware HAC (NULL if not used)
    - num_hac_panels: Number of panels for HAC
    - underid_stat: Output - underidentification test statistic
    - underid_df: Output - degrees of freedom
    - cd_f: Output - Cragg-Donald F statistic
    - kp_f: Output - Kleibergen-Paap rk Wald F statistic
*/
void civreghdfe_compute_underid_test(
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *ZtZ,
    const ST_double *ZtZ_inv,
    const ST_double *temp1,
    const ST_double *first_stage_F,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    const ST_int *cluster2_ids,
    ST_int num_clusters2,
    ST_int kernel_type,
    ST_int bw,
    ST_int kiefer,
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels,
    ST_double *underid_stat,
    ST_int *underid_df,
    ST_double *cd_f,
    ST_double *kp_f
);

/*
    Compute Sargan/Hansen J overidentification test.

    For homoskedastic VCE: Sargan statistic
    For robust/cluster VCE: Hansen J statistic

    Parameters:
    - resid: 2SLS residuals (N x 1)
    - Z: All instruments (N x K_iv)
    - ZtZ_inv: Precomputed (Z'Z)^-1 (K_iv x K_iv)
    - sigma2: Residual variance
    - weights, weight_type, N, K_exog, K_endog, K_iv: Parameters
    - vce_type: 0=unadjusted, 1=robust, 2=cluster, 4=dkraay
    - cluster_ids: Cluster identifiers (may be NULL)
    - num_clusters: Number of clusters
    - kernel_type: HAC kernel type (0=none, 1=Bartlett, etc.)
    - bw: HAC bandwidth
    - kiefer: 1 if Kiefer VCE, 0 otherwise
    - hac_panel_ids: Panel IDs for panel-aware HAC (NULL if not used)
    - num_hac_panels: Number of panels for HAC
    - sargan_stat: Output - overidentification test statistic
    - overid_df: Output - degrees of freedom (L - K_endog)
*/
void civreghdfe_compute_sargan_j(
    const ST_double *resid,
    const ST_double *Z,
    const ST_double *ZtZ_inv,
    ST_double sigma2,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_int vce_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_int kernel_type,
    ST_int bw,
    ST_int kiefer,
    const ST_int *hac_panel_ids,
    ST_int num_hac_panels,
    ST_double *sargan_stat,
    ST_int *overid_df
);

/*
    Compute Durbin-Wu-Hausman endogeneity test.

    Tests H0: X_endog are exogenous (no endogeneity)
    Uses augmented regression approach with first-stage residuals.

    Parameters:
    - y: Dependent variable (N x 1)
    - X_exog: Exogenous regressors (N x K_exog, may be NULL)
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv)
    - temp1: First-stage coefficients (K_iv x K_total)
    - N, K_exog, K_endog, K_iv, df_a: Parameters
    - endog_chi2: Output - DWH chi-squared statistic
    - endog_f: Output - DWH F statistic
    - endog_df: Output - degrees of freedom
*/
void civreghdfe_compute_dwh_test(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *temp1,
    ST_int N,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    ST_double *endog_chi2,
    ST_double *endog_f,
    ST_int *endog_df
);

/*
    Compute C-statistic (orthogonality test for specified instruments).

    Tests H0: specified instruments are exogenous.
    C = J_full - J_restricted (difference in overidentification statistics).
    df = number of tested instruments.

    Parameters:
    - y: Dependent variable (N x 1)
    - X_exog: Exogenous regressors (N x K_exog)
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv)
    - K_exog: Number of exogenous regressors (included in instruments)
    - K_endog: Number of endogenous regressors
    - K_iv: Total instruments
    - orthog_indices: 1-based indices of excluded instruments to test
    - n_orthog: Number of instruments to test
    - vce_type, weights, weight_type, cluster_ids, num_clusters: VCE parameters
    - sargan_full: Full model Sargan statistic (already computed)
    - cstat: Output - C statistic
    - cstat_df: Output - degrees of freedom
*/
void civreghdfe_compute_cstat(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const ST_int *orthog_indices,
    ST_int n_orthog,
    ST_int vce_type,
    const ST_double *weights,
    ST_int weight_type,
    const ST_int *cluster_ids,
    ST_int num_clusters,
    ST_double sargan_full,
    ST_double rss_full,
    ST_double *cstat,
    ST_int *cstat_df
);

/*
    Compute endogeneity test for specified subset of endogenous regressors.

    Tests H0: specified regressors are exogenous (can be treated as exogenous).
    Uses difference-in-Sargan C-statistic approach (matching ivreg2):
    C = J_exog - J_orig, where:
    - J_exog: Sargan from model where tested vars are treated as exogenous
      (added to instruments), using sigma^2 from that model
    - J_orig: Sargan from original model, using sigma^2 from the exogenous model

    Parameters:
    - y: Dependent variable (N x 1)
    - X_exog: Exogenous regressors (N x K_exog, may be NULL)
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv)
    - resid: Original model residuals (N x 1)
    - ZtZ_inv: Original model (Z'Z)^-1 (K_iv x K_iv)
    - N, K_exog, K_endog, K_iv: Parameters
    - endogtest_indices: 1-based indices of endogenous regressors to test
    - n_endogtest: Number of regressors to test
    - endogtest_stat: Output - C statistic (chi-sq)
    - endogtest_df: Output - degrees of freedom
*/
void civreghdfe_compute_endogtest_subset(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *resid,
    const ST_double *ZtZ_inv,
    ST_int N,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const ST_int *endogtest_indices,
    ST_int n_endogtest,
    ST_double *endogtest_stat,
    ST_int *endogtest_df
);

/*
    Compute instrument redundancy test.

    Tests H0: specified instruments add no information beyond other instruments.
    LM test based on partial correlation.

    Parameters:
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv)
    - K_exog, K_endog, K_iv: Dimensions
    - redund_indices: 1-based indices of excluded instruments to test
    - n_redund: Number of instruments to test
    - redund_stat: Output - LM test statistic (chi-sq)
    - redund_df: Output - degrees of freedom (K_endog * n_redund)
*/
void civreghdfe_compute_redundant(
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const ST_int *redund_indices,
    ST_int n_redund,
    ST_double *redund_stat,
    ST_int *redund_df
);

#endif /* CIVREGHDFE_TESTS_H */
