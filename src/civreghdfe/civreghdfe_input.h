#ifndef CIVREGHDFE_INPUT_H
#define CIVREGHDFE_INPUT_H
#include "stplugin.h"
/* Wire names/defaults belong here; numerical phases consume a typed snapshot. */
typedef struct
{
    ST_int K_endog, K_exog, K_iv, G, has_cluster, has_cluster2, has_weights, weight_type, vce_type,
        maxiter, verbose, est_method, kernel_type, kiefer, ts_mode, ts_cols, dofminus,
        sdofminus_opt, nopartialsmall, center, touse;
    ST_double tolerance, kclass_user, fuller_alpha, tdelta, bw;
    /* Heterogeneous slopes: factor g is a slope (slope_flag[g]) of the
     * absorbed effect fe_absvar[g] (1-based); the num_slopes slope variables
     * follow the tsset columns. dof_continuous: dofadjustments(continuous) */
    ST_int num_slopes, slope_flag[10], fe_absvar[10], dof_continuous;
} civreghdfe_input;
ST_retcode civreghdfe_read_input(civreghdfe_input *input);
#endif
