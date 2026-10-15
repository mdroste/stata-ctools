#ifndef CPPLMHDFE_INPUT_H
#define CPPLMHDFE_INPUT_H
#include "stplugin.h"
/* Wire names/defaults belong here; numerical phases consume a typed snapshot. */
typedef struct
{
    ST_int K, G, verbose, maxiter, vcetype, cluster_terms, has_weights, weight_type, has_offset,
        irls_maxiter, dof_adjust_type, dof_flags, keep_singletons, check_fe, check_relu;
    ST_double tolerance, irls_tol;
    ST_int slope_only[10], slope_idx[10];
} cpplmhdfe_input;
ST_retcode cpplmhdfe_read_input(cpplmhdfe_input *input);
#endif
