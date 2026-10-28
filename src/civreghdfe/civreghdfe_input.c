#include "civreghdfe_input.h"
#include <string.h>
#include <stdio.h>

/* A scalar the wrapper always sets: a missing one is a protocol error, not a
 * value left over from the previous read. */
static ST_retcode read_required(const char *name, ST_double *value)
{
    if (SF_scal_use((char *)name, value) == 0)
        return 0;
    char msg[96];
    snprintf(msg, sizeof(msg), "civreghdfe: plugin scalar %s is not set\n", name);
    SF_error(msg);
    return 198;
}

ST_retcode civreghdfe_read_input(civreghdfe_input *input)
{
    ST_double dval = 0;
    ST_retcode rc;
    /* Read configuration scalars */
    ST_int K_endog, K_exog, K_iv, G;
    ST_int has_cluster, has_weights, weight_type;
    ST_int vce_type, maxiter, verbose;
    ST_double tolerance;

    if ((rc = read_required("__civreghdfe_K_endog", &dval)))
        return rc;
    K_endog = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_K_exog", &dval)))
        return rc;
    K_exog = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_K_iv", &dval)))
        return rc;
    K_iv = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_G", &dval)))
        return rc;
    G = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_has_cluster", &dval)))
        return rc;
    has_cluster = (ST_int)dval;
    ST_int has_cluster2 = 0;
    if ((rc = read_required("__civreghdfe_has_cluster2", &dval)))
        return rc;
    has_cluster2 = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_has_weights", &dval)))
        return rc;
    has_weights = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_weight_type", &dval)))
        return rc;
    weight_type = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_vce_type", &dval)))
        return rc;
    vce_type = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_maxiter", &dval)))
        return rc;
    maxiter = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_tolerance", &dval)))
        return rc;
    tolerance = dval;
    if ((rc = read_required("__civreghdfe_verbose", &dval)))
        return rc;
    verbose = (ST_int)dval;
    /* nested_fe_index from .ado is no longer used — data-based detection in C handles it */

    /* New scalars for estimation method */
    ST_int est_method = 0;
    ST_double kclass_user = 0, fuller_alpha = 0;
    if ((rc = read_required("__civreghdfe_est_method", &dval)))
        return rc;
    est_method = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_kclass", &dval)))
        return rc;
    kclass_user = dval;
    if ((rc = read_required("__civreghdfe_fuller", &dval)))
        return rc;
    fuller_alpha = dval;

    /* HAC parameters */
    ST_int kernel_type = 0, bw = 0, kiefer = 0;
    if ((rc = read_required("__civreghdfe_kernel", &dval)))
        return rc;
    kernel_type = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_bw", &dval)))
        return rc;
    bw = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_kiefer", &dval)))
        return rc;
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
    if ((rc = read_required("__civreghdfe_dofminus", &dval)))
        return rc;
    dofminus = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_sdofminus", &dval)))
        return rc;
    sdofminus_opt = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_nopartialsmall", &dval)))
        return rc;
    nopartialsmall = (ST_int)dval;
    if ((rc = read_required("__civreghdfe_center", &dval)))
        return rc;
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
