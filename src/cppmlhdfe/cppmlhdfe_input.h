#ifndef CPPMLHDFE_INPUT_H
#define CPPMLHDFE_INPUT_H
#include "stplugin.h"
/* Largest regressor count whose square matrices, with the constant column
 * added, fit int products: (46339 + 1)^2 < 2^31 <= 46341^2 */
#define CPPMLHDFE_MAX_REGRESSORS 46339
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
