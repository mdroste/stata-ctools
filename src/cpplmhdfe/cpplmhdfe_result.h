#ifndef CPPLMHDFE_RESULT_H
#define CPPLMHDFE_RESULT_H
#include "cpplmhdfe_types.h"
/* Borrowed result view; arrays stay owned by the estimator through publication. */
typedef struct
{
    ST_int N;
    ST_int weight_type;
    ST_int has_weights;
    ST_int num_singletons;
    ST_int num_separated;
    ST_int K_keep;
    ST_int df_a;
    ST_int mobility_groups;
    ST_int irls_iter;
    ST_int irls_converged;
    ST_int G;
    ST_int cluster_terms;
    ST_int vcetype;
    ST_int num_clusters;
    ST_double deviance;
    ST_double stdev_y;
    ST_double user_weight_scale;
    ST_double ll;
    ST_double ll_0;
    const ST_double *w_user_compact;
    const ST_double *beta_final;
    const ST_double *V_keep;
    const PPMLFactorData *factors;
} cpplmhdfe_result;
ST_retcode cpplmhdfe_post_results(const cpplmhdfe_result *result);
#endif
