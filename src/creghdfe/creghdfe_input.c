#include "creghdfe_input.h"
#include <string.h>
#include <stdio.h>

/* A scalar the wrapper always sets: a missing one is a protocol error, not a
 * value left over from the previous read. */
static ST_retcode read_required(const char *name, ST_double *value)
{
    if (SF_scal_use((char *)name, value) == 0)
        return 0;
    char msg[96];
    snprintf(msg, sizeof(msg), "creghdfe: plugin scalar %s is not set\n", name);
    SF_error(msg);
    return 198;
}

ST_retcode creghdfe_read_input(creghdfe_input *input)
{
    ST_int K = 0;
    ST_int G = 0;
    ST_int drop_singletons = 0;
    ST_int verbose = 0;
    ST_int maxiter = 0;
    ST_int standardize = 0;
    ST_int vcetype = 0;
    ST_int compute_dof = 0;
    ST_int compute_resid = 0;
    ST_int resid_var_idx = 0;
    ST_double df_a_nested = 0;
    ST_double val = 0;
    ST_retcode rc;
    /* Read parameters */
    if ((rc = read_required("__creghdfe_K", &val)))
        return rc;
    K = (ST_int)val;
    if ((rc = read_required("__creghdfe_G", &val)))
        return rc;
    G = (ST_int)val;
    if ((rc = read_required("__creghdfe_drop_singletons", &val)))
        return rc;
    drop_singletons = (ST_int)val;
    if ((rc = read_required("__creghdfe_verbose", &val)))
        return rc;
    verbose = (ST_int)val;
    if ((rc = read_required("__creghdfe_maxiter", &val)))
        return rc;
    maxiter = (ST_int)val;
    if ((rc = read_required("__creghdfe_tolerance", &val)))
        return rc;
    ST_double tolerance = val;
    if ((rc = read_required("__creghdfe_standardize", &val)))
        return rc;
    standardize = (ST_int)val;
    if ((rc = read_required("__creghdfe_vce_type", &val)))
        return rc;
    vcetype = (ST_int)val;
    if ((rc = read_required("__creghdfe_df_a_nested", &val)))
        return rc;
    df_a_nested = val;

    /* Read weight parameters */
    ST_int has_weights = 0, weight_type = 0;
    if (SF_scal_use("__creghdfe_has_weights", &val) == 0)
    {
        has_weights = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_weight_type", &val) == 0)
    {
        weight_type = (ST_int)val;
    }

    compute_dof = 0;
    if (SF_scal_use("__creghdfe_compute_dof", &val) == 0)
    {
        compute_dof = (ST_int)val;
    }

    /* Read DOF adjustment type: 0=all, 1=none, 2=firstpair, 3=pairwise.
     * "all" includes reghdfe's cluster-nesting adjustment; the others apply
     * it only when the wrapper also requests it. */
    ST_int dof_adjust_type = 0;
    if (SF_scal_use("__creghdfe_dof_adjust_type", &val) == 0)
    {
        dof_adjust_type = (ST_int)val;
    }
    ST_int dof_clusters = dof_adjust_type == 0;
    if (dof_adjust_type != 0 && SF_scal_use("__creghdfe_dof_clusters", &val) == 0)
    {
        dof_clusters = val != 0;
    }

    /* Read savefe flag */
    ST_int savefe = 0;
    ST_int savefe_var_idx = 0;
    if (SF_scal_use("__creghdfe_savefe", &val) == 0)
    {
        savefe = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_savefe_idx", &val) == 0)
    {
        savefe_var_idx = (ST_int)val;
    }

    /* reghdfe judges collinearity with the fixed effects by the tolerance()
     * the user gave, so it is taken before the savefe override below. */
    input->collinear_tol = (tolerance / 10.0 < 1e-6) ? tolerance / 10.0 : 1e-6;

    /* Individual FE estimates require a more precise projection than slopes. */
    if (savefe && tolerance > 1e-12)
        tolerance = 1e-12;

    /* Read groupvar flag */
    ST_int compute_groupvar = 0;
    ST_int groupvar_var_idx = 0;
    if (SF_scal_use("__creghdfe_compute_groupvar", &val) == 0)
    {
        compute_groupvar = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_groupvar_idx", &val) == 0)
    {
        groupvar_var_idx = (ST_int)val;
    }

    /* Check if we should compute and store residuals */
    if (SF_scal_use("__creghdfe_compute_resid", &val) == 0)
    {
        compute_resid = (ST_int)val;
    }
    if (SF_scal_use("__creghdfe_resid_var_idx", &val) == 0)
    {
        resid_var_idx = (ST_int)val;
    }

    /* Check if quad-precision accumulation is requested */
    ST_int use_quad = 0;
    if (SF_scal_use("__creghdfe_use_quad", &val) == 0)
    {
        use_quad = (ST_int)val;
    }

    /* The sample marker follows every other variable in place of an if */
    ST_int touse = 0;
    if (SF_scal_use("__creghdfe_touse", &val) == 0)
    {
        touse = val != 0;
    }

    ST_int bs_unit_idx = 0, bs_strata_idx = 0;
    if (SF_scal_use("__creghdfe_bs_unit_idx", &val) == 0)
        bs_unit_idx = (ST_int)val;
    if (SF_scal_use("__creghdfe_bs_strata_idx", &val) == 0)
        bs_strata_idx = (ST_int)val;
    input->bs_retain = 0;
    input->bs_unit_idx = bs_unit_idx;
    input->bs_strata_idx = bs_strata_idx;

    /* One cluster variable unless the wrapper says otherwise */
    input->num_cluster = 1;
    if (SF_scal_use("__creghdfe_num_cluster", &val) == 0 && val >= 1)
        input->num_cluster = (ST_int)val;
    input->dkraay_bw = 0;
    if (SF_scal_use("__creghdfe_dkraay_bw", &val) == 0 && val >= 1 && val < 2e9)
        input->dkraay_bw = (ST_int)val;
    input->dof_continuous = dof_adjust_type == 0;
    if (SF_scal_use("__creghdfe_dof_continuous", &val) == 0)
        input->dof_continuous = val != 0;
    input->num_slopes = 0;
    for (ST_int g = 0; g < 10; g++)
    {
        char name[48];
        input->slope_idx[g] = 0;
        input->fe_absvar[g] = g + 1;
        if (g >= G)
            continue;
        snprintf(name, sizeof(name), "__creghdfe_slope_%d", (int)g + 1);
        if (SF_scal_use(name, &val) == 0 && val > 0)
            input->slope_idx[g] = (ST_int)val;
        snprintf(name, sizeof(name), "__creghdfe_absvar_%d", (int)g + 1);
        if (SF_scal_use(name, &val) == 0 && val >= 1)
            input->fe_absvar[g] = (ST_int)val;
        input->num_slopes += input->slope_idx[g] > 0;
    }

    input->indiv = 0;
    input->indiv_mean = 1;
    input->indiv_gidx = input->indiv_iidx = input->indiv_tidx = 0;
    if (SF_scal_use("__creghdfe_indiv", &val) == 0 && val != 0)
    {
        ST_double gi = 0, ii = 0, ti = 0, mv = 1;
        if (SF_scal_use("__creghdfe_indiv_gidx", &gi) || SF_scal_use("__creghdfe_indiv_iidx", &ii) ||
            SF_scal_use("__creghdfe_indiv_tidx", &ti) || SF_scal_use("__creghdfe_indiv_mean", &mv) ||
            gi < 1 || ii < 1 || ti < 1)
        {
            SF_error("creghdfe: the individual-effect variables are not set\n");
            return 198;
        }
        input->indiv = 1;
        input->indiv_mean = mv != 0;
        input->indiv_gidx = (ST_int)gi;
        input->indiv_iidx = (ST_int)ii;
        input->indiv_tidx = (ST_int)ti;
    }

    input->transform = 0;
    input->acceleration = 1;
    input->abort_nonconv = 1;
    if (SF_scal_use("__creghdfe_transform", &val) == 0 && val >= 0 && val <= 2)
        input->transform = (ST_int)val;
    if (SF_scal_use("__creghdfe_acceleration", &val) == 0 && val >= 0 && val <= 3)
        input->acceleration = (ST_int)val;
    if (SF_scal_use("__creghdfe_abort", &val) == 0)
        input->abort_nonconv = val != 0;

    input->K = K;
    input->G = G;
    input->drop_singletons = drop_singletons;
    input->verbose = verbose;
    input->maxiter = maxiter;
    input->standardize = standardize;
    input->vcetype = vcetype;
    input->has_weights = has_weights;
    input->weight_type = weight_type;
    input->compute_dof = compute_dof;
    input->dof_adjust_type = dof_adjust_type;
    input->dof_clusters = dof_clusters;
    input->savefe = savefe;
    input->savefe_var_idx = savefe_var_idx;
    input->compute_groupvar = compute_groupvar;
    input->groupvar_var_idx = groupvar_var_idx;
    input->compute_resid = compute_resid;
    input->resid_var_idx = resid_var_idx;
    input->use_quad = use_quad;
    input->touse = touse;
    input->tolerance = tolerance;
    input->df_a_nested = df_a_nested;
    return 0;
}
