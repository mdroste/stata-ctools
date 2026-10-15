#ifndef CREGHDFE_INPUT_H
#define CREGHDFE_INPUT_H
#include "stplugin.h"
/* Wire names/defaults belong here; numerical phases consume a typed snapshot. */
typedef struct
{
    ST_int K, G, drop_singletons, verbose, maxiter, standardize, vcetype, has_weights, weight_type,
        compute_dof, dof_adjust_type, dof_clusters, savefe, savefe_var_idx, compute_groupvar,
        groupvar_var_idx, compute_resid, resid_var_idx, use_quad;
    ST_double tolerance, df_a_nested;
} creghdfe_input;
ST_retcode creghdfe_read_input(creghdfe_input *input);
#endif
