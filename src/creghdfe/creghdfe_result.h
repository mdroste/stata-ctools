#ifndef CREGHDFE_RESULT_H
#define CREGHDFE_RESULT_H
#include "creghdfe_types.h"
/* Borrowed result view; arrays stay owned by the estimator through publication. */
typedef struct
{
    ST_int weight_type;
    ST_int has_weights;
    ST_int N;
    ST_int num_singletons;
    ST_int G;
    ST_int mobility_groups;
    ST_int K_keep;
    ST_int K_with_cons;
    ST_int vcetype;
    ST_int num_clusters;
    ST_double df_a;
    ST_double df_a_nested;
    ST_double rss;
    ST_double tss_within;
    const ST_double *tss;
    const ST_double *beta_keep;
    const ST_double *V_keep;
    const FactorData *factors;
    const HDFE_State *g_state;
} creghdfe_result;
ST_retcode creghdfe_post_results(const creghdfe_result *result);
#endif
