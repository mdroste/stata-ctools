#ifndef CIVREGHDFE_RESULT_H
#define CIVREGHDFE_RESULT_H
#include "civreghdfe_impl.h"
/* Borrowed result view; arrays stay owned by the estimator through publication. */
typedef struct
{
    ST_double df_r_val;   /* N_eff-based: may exceed INT_MAX with fweights */
    ST_int df_a;
    ST_int df_a_for_vce;
    ST_int K_total;
    ST_int n_abs_x;
    ST_int has_cluster2;
    ST_int num_clusters;
    ST_int num_clusters2;
    ST_int has_cluster;
    ST_int est_method;
    ST_int K_exog_orig;
    ST_int K_endog_orig;
    ST_int K_endog;
    ST_double N_eff;
    ST_double rss;
    ST_double tss;
    ST_double r2;
    ST_double rmse;
    ST_double lambda;
    const ST_double *beta;
    const ST_double *V;
    const ST_double *first_stage_F;
    const ST_int *is_collinear_x;
    ST_retcode residual_rc;
} civreghdfe_result;
ST_retcode civreghdfe_post_results(const civreghdfe_result *result);
#endif
