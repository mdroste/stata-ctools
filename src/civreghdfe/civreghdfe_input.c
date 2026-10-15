#include "civreghdfe_input.h"
#include <string.h>
#include <stdio.h>

ST_retcode civreghdfe_read_input(civreghdfe_input *input)
{
    ST_double dval = 0;
    /* Read configuration scalars */
    ST_int K_endog, K_exog, K_iv, G;
    ST_int has_cluster, has_weights, weight_type;
    ST_int vce_type, maxiter, verbose;
    ST_double tolerance;

    SF_scal_use("__civreghdfe_K_endog", &dval);
    K_endog = (ST_int)dval;
    SF_scal_use("__civreghdfe_K_exog", &dval);
    K_exog = (ST_int)dval;
    SF_scal_use("__civreghdfe_K_iv", &dval);
    K_iv = (ST_int)dval;
    SF_scal_use("__civreghdfe_G", &dval);
    G = (ST_int)dval;
    SF_scal_use("__civreghdfe_has_cluster", &dval);
    has_cluster = (ST_int)dval;
    ST_int has_cluster2 = 0;
    SF_scal_use("__civreghdfe_has_cluster2", &dval);
    has_cluster2 = (ST_int)dval;
    SF_scal_use("__civreghdfe_has_weights", &dval);
    has_weights = (ST_int)dval;
    SF_scal_use("__civreghdfe_weight_type", &dval);
    weight_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_vce_type", &dval);
    vce_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_maxiter", &dval);
    maxiter = (ST_int)dval;
    SF_scal_use("__civreghdfe_tolerance", &dval);
    tolerance = dval;
    SF_scal_use("__civreghdfe_verbose", &dval);
    verbose = (ST_int)dval;
    /* nested_fe_index from .ado is no longer used — data-based detection in C handles it */

    /* New scalars for estimation method */
    ST_int est_method = 0;
    ST_double kclass_user = 0, fuller_alpha = 0;
    SF_scal_use("__civreghdfe_est_method", &dval);
    est_method = (ST_int)dval;
    SF_scal_use("__civreghdfe_kclass", &dval);
    kclass_user = dval;
    SF_scal_use("__civreghdfe_fuller", &dval);
    fuller_alpha = dval;

    /* HAC parameters */
    ST_int kernel_type = 0, bw = 0, kiefer = 0;
    SF_scal_use("__civreghdfe_kernel", &dval);
    kernel_type = (ST_int)dval;
    SF_scal_use("__civreghdfe_bw", &dval);
    bw = (ST_int)dval;
    SF_scal_use("__civreghdfe_kiefer", &dval);
    kiefer = (ST_int)dval;

    /* tsset structure for the kernel estimators (HAC, AC, Driscoll-Kraay,
       Kiefer): ts_mode 1 = the time variable, 2 = the panel and time
       variables, passed after the weight; tdelta = tsset delta */
    ST_int ts_mode = 0;
    ST_double tdelta = 1.0;
    if (SF_scal_use("__civreghdfe_ts_mode", &dval) == 0 && dval > 0)
        ts_mode = (dval >= 2) ? 2 : 1;
    if (SF_scal_use("__civreghdfe_tdelta", &dval) == 0 && dval > 0 && !SF_is_missing(dval))
        tdelta = dval;
    ST_int ts_cols = (ts_mode == 2) ? 2 : ts_mode;

    /* DOF adjustment parameters */
    ST_int dofminus = 0, sdofminus_opt = 0, nopartialsmall = 0, center = 0;
    SF_scal_use("__civreghdfe_dofminus", &dval);
    dofminus = (ST_int)dval;
    SF_scal_use("__civreghdfe_sdofminus", &dval);
    sdofminus_opt = (ST_int)dval;
    SF_scal_use("__civreghdfe_nopartialsmall", &dval);
    nopartialsmall = (ST_int)dval;
    SF_scal_use("__civreghdfe_center", &dval);
    center = (ST_int)dval;

    input->K_endog = K_endog;
    input->K_exog = K_exog;
    input->K_iv = K_iv;
    input->G = G;
    input->has_cluster = has_cluster;
    input->has_cluster2 = has_cluster2;
    input->has_weights = has_weights;
    input->weight_type = weight_type;
    input->vce_type = vce_type;
    input->maxiter = maxiter;
    input->verbose = verbose;
    input->est_method = est_method;
    input->kernel_type = kernel_type;
    input->bw = bw;
    input->kiefer = kiefer;
    input->ts_mode = ts_mode;
    input->ts_cols = ts_cols;
    input->dofminus = dofminus;
    input->sdofminus_opt = sdofminus_opt;
    input->nopartialsmall = nopartialsmall;
    input->center = center;
    input->tolerance = tolerance;
    input->kclass_user = kclass_user;
    input->fuller_alpha = fuller_alpha;
    input->tdelta = tdelta;
    return 0;
}
