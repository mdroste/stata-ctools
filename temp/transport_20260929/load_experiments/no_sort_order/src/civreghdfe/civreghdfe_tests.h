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
#include "civreghdfe_omega.h"

/*
    Compute underidentification test (Anderson LM / Kleibergen-Paap rk LM)
    and weak identification statistics (Cragg-Donald / KP rk Wald F).

    As ivreghdfe calls ranktest: with an iid rk_om (no robust, cluster or
    bw) the Anderson canonical-correlation LM and the Cragg-Donald F (which
    is then also the reported rk Wald F); otherwise the Kleibergen-Paap rk LM
    and rk Wald statistics for any number of endogenous regressors (ranktest's
    SVD algorithm) with rk_om's moment covariance.

    Parameters:
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv), included exogenous first
    - weights, weight_type, N, N_eff, K_exog, K_endog, K_iv, df_a: Parameters
    - rk_om: ranktest's S structure (a kernel becomes Bartlett, Kiefer iid)
    - rk_cluster, rk_nclust: cluster-type F conversion and its cluster count
    - underid_stat: Output - underidentification test statistic
    - underid_df: Output - degrees of freedom
    - cd_f: Output - Cragg-Donald F statistic
    - kp_f: Output - Kleibergen-Paap rk Wald F statistic
*/
void civreghdfe_compute_underid_test(
    const ST_double *X_endog,
    const ST_double *Z,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    ST_int N_eff,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int df_a,
    const civ_omega *rk_om,
    ST_int rk_cluster,
    ST_int rk_nclust,
    ST_double *underid_stat,
    ST_int *underid_df,
    ST_double *cd_f,
    ST_double *kp_f
);

/*
    Compute Sargan/Hansen J overidentification test.

    J = g(e)' S^-1 g(e), g = Z'We, where S is the model's moment covariance
    (sigma^2 Z'WZ when iid: the Sargan statistic).

    Parameters:
    - y, X_all (N x K_total), Z (N x K_iv): Model data
    - weights, weight_type, N: Weighting and sample size
    - S: Moment covariance (K_iv x K_iv, raw)
    - resid: residuals for J (the model's own when it is efficient for S:
      iid, GMM2S, CUE); NULL = the residuals of efficient GMM for S (Hansen J
      of an inefficient 2SLS / k-class fit)
    - J: Output - overidentification test statistic
    - overid_df: Output - degrees of freedom (K_iv - K_total)
*/
void civreghdfe_compute_hansen_j(
    const ST_double *y,
    const ST_double *X_all,
    ST_int K_total,
    const ST_double *Z,
    ST_int K_iv,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    const ST_double *S,
    const ST_double *resid,
    ST_double *J,
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
    C = J_full - J_r, where J_r is the efficient GMM J without the tested
    instruments, weighted by their complement's block of the full model's S
    (ivreg2's smatrix approach). df = number of tested instruments.

    Parameters:
    - y, X_all (N x K_total), Z (N x K_iv, K_exog included exogenous first)
    - weights, weight_type, N: Weighting and sample size
    - S: Full model's moment covariance (K_iv x K_iv, raw)
    - J_full: Full model's J (already computed with S)
    - orthog_indices: distinct 1-based indices of excluded instruments to test
    - n_orthog: Number of instruments to test
    - cstat: Output - C statistic
    - cstat_df: Output - degrees of freedom
*/
void civreghdfe_compute_cstat(
    const ST_double *y,
    const ST_double *X_all,
    ST_int K_total,
    const ST_double *Z,
    ST_int K_iv,
    ST_int K_exog,
    const ST_double *weights,
    ST_int weight_type,
    ST_int N,
    const ST_double *S,
    ST_double J_full,
    const ST_int *orthog_indices,
    ST_int n_orthog,
    ST_double *cstat,
    ST_int *cstat_df
);

/*
    Compute endogeneity test for specified subset of endogenous regressors.

    Tests H0: specified regressors are exogenous (can be treated as exogenous).
    As ivreg2's endog(): the C statistic of the re-estimated model that
    treats the tested regressors as exogenous (they join the instruments),
    stat = J_exog - J_r, J_r = efficient GMM J of the original model weighted
    by the exogenous model's S block for the original instruments.

    Parameters:
    - y, X_exog (may be NULL), X_endog, Z (N x K_iv): Model data
    - N, weights, weight_type, K_exog, K_endog, K_iv: Parameters
    - est_type: estimator of the re-estimated model (0 2SLS, 1 LIML, 4 GMM2S)
    - om: S structure of the re-estimated model
    - endogtest_indices: distinct 1-based indices of endogenous regressors to test
    - n_endogtest: Number of regressors to test
    - endogtest_stat: Output - C statistic (chi-sq)
    - endogtest_df: Output - degrees of freedom
*/
void civreghdfe_compute_endogtest_subset(
    const ST_double *y,
    const ST_double *X_exog,
    const ST_double *X_endog,
    const ST_double *Z,
    ST_int N,
    const ST_double *weights,
    ST_int weight_type,
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    ST_int est_type,
    const civ_omega *om,
    const ST_int *endogtest_indices,
    ST_int n_endogtest,
    ST_double *endogtest_stat,
    ST_int *endogtest_df
);

/*
    Compute instrument redundancy test.

    Tests H0: specified instruments add no information beyond other instruments:
    ranktest's LM test of rank 0 for the endogenous regressors and the tested
    instruments, with om's moment covariance.

    Parameters:
    - X_endog: Endogenous regressors (N x K_endog)
    - Z: All instruments (N x K_iv)
    - K_exog, K_endog, K_iv: Dimensions
    - om: ranktest's S structure
    - redund_indices: distinct 1-based indices of excluded instruments to test
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
    ST_int K_exog,
    ST_int K_endog,
    ST_int K_iv,
    const civ_omega *om,
    const ST_int *redund_indices,
    ST_int n_redund,
    ST_double *redund_stat,
    ST_int *redund_df
);

#endif /* CIVREGHDFE_TESTS_H */
