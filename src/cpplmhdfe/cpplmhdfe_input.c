#include "cpplmhdfe_input.h"
#include <string.h>
#include <stdio.h>

ST_retcode cpplmhdfe_read_input(cpplmhdfe_input *input)
{
    ST_int K = 0;
    ST_int G = 0;
    ST_int verbose = 0;
    ST_double val = 0;
    ST_int g;
    char scalar_name[64];
    /* Read parameters */
    SF_scal_use("__cpplmhdfe_K", &val);
    K = (ST_int)val;
    SF_scal_use("__cpplmhdfe_G", &val);
    G = (ST_int)val;
    SF_scal_use("__cpplmhdfe_verbose", &val);
    verbose = (ST_int)val;

    ST_int maxiter;
    ST_double tolerance;
    SF_scal_use("__cpplmhdfe_maxiter", &val);
    maxiter = (ST_int)val;
    SF_scal_use("__cpplmhdfe_tolerance", &val);
    tolerance = val;

    ST_int vcetype;
    SF_scal_use("__cpplmhdfe_vce_type", &val);
    vcetype = (ST_int)val;

    ST_int cluster_terms = vcetype == 2 ? 1 : 0;
    if (SF_scal_use("__cpplmhdfe_cluster_terms", &val) == 0)
        cluster_terms = (ST_int)val;
    if (vcetype == 2 && (cluster_terms < 1 || cluster_terms > 1023))
        return 198;

    ST_int has_weights = 0, weight_type = 0;
    if (SF_scal_use("__cpplmhdfe_has_weights", &val) == 0)
        has_weights = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_weight_type", &val) == 0)
        weight_type = (ST_int)val;

    ST_int has_offset = 0;
    if (SF_scal_use("__cpplmhdfe_has_offset", &val) == 0)
        has_offset = (ST_int)val;

    /* IRLS parameters */
    ST_int irls_maxiter = 10000;
    ST_double irls_tol = 1e-8;
    if (SF_scal_use("__cpplmhdfe_irls_maxiter", &val) == 0)
        irls_maxiter = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_irls_tol", &val) == 0)
        irls_tol = val;

    /* DOF adjustment */
    ST_int dof_adjust_type = 0;
    if (SF_scal_use("__cpplmhdfe_dof_adjust_type", &val) == 0)
        dof_adjust_type = (ST_int)val;

    ST_int dof_flags = 13; /* pairwise=1, firstpair=2, clusters=4, continuous=8 */
    if (SF_scal_use("__cpplmhdfe_dof_flags", &val) == 0)
        dof_flags = (ST_int)val;

    ST_int slope_only[10] = {0}, slope_idx[10] = {0};
    ST_int keep_singletons = 0;
    if (SF_scal_use("__cpplmhdfe_keepsingletons", &val) == 0)
        keep_singletons = (ST_int)val;
    for (g = 0; g < G && g < 10; g++)
    {
        snprintf(scalar_name, sizeof(scalar_name), "__cpplmhdfe_slope_%d", g + 1);
        if (SF_scal_use(scalar_name, &val) == 0)
            slope_idx[g] = (ST_int)val;
        slope_only[g] = slope_idx[g] > 0;
    }

    ST_int check_fe = 1, check_relu = 1;
    if (SF_scal_use("__cpplmhdfe_check_fe", &val) == 0)
        check_fe = (ST_int)val;
    if (SF_scal_use("__cpplmhdfe_check_relu", &val) == 0)
        check_relu = (ST_int)val;

    input->K = K;
    input->G = G;
    input->verbose = verbose;
    input->maxiter = maxiter;
    input->vcetype = vcetype;
    input->cluster_terms = cluster_terms;
    input->has_weights = has_weights;
    input->weight_type = weight_type;
    input->has_offset = has_offset;
    input->irls_maxiter = irls_maxiter;
    input->dof_adjust_type = dof_adjust_type;
    input->dof_flags = dof_flags;
    input->keep_singletons = keep_singletons;
    input->check_fe = check_fe;
    input->check_relu = check_relu;
    input->tolerance = tolerance;
    input->irls_tol = irls_tol;
    memcpy(input->slope_only, slope_only, sizeof(slope_only));
    memcpy(input->slope_idx, slope_idx, sizeof(slope_idx));
    return 0;
}
