#ifndef CIVREGHDFE_INPUT_H
#define CIVREGHDFE_INPUT_H
#include "stplugin.h"
/* Wire names/defaults belong here; numerical phases consume a typed snapshot. */
typedef struct
{
    ST_int K_endog, K_exog, K_iv, G, has_cluster, has_cluster2, has_weights, weight_type, vce_type,
        maxiter, verbose, est_method, kernel_type, bw, kiefer, ts_mode, ts_cols, dofminus,
        sdofminus_opt, nopartialsmall, center;
    ST_double tolerance, kclass_user, fuller_alpha, tdelta;
} civreghdfe_input;
ST_retcode civreghdfe_read_input(civreghdfe_input *input);
#endif
