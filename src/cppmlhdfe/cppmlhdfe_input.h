#ifndef CPPMLHDFE_INPUT_H
#define CPPMLHDFE_INPUT_H
#include "stplugin.h"
/* Wire names/defaults belong here; numerical phases consume a typed snapshot. */
typedef struct
{
    ST_int K, G, verbose, maxiter, vcetype, cluster_terms, has_weights, weight_type, has_offset,
        irls_maxiter, dof_adjust_type, dof_flags, keep_singletons, check_fe, check_relu, touse;
    ST_double tolerance, irls_tol;
    ST_int slope_only[10], slope_idx[10];
} cppmlhdfe_input;
ST_retcode cppmlhdfe_read_input(cppmlhdfe_input *input);
#endif
