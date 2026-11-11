#ifndef CREGHDFE_INPUT_H
#define CREGHDFE_INPUT_H
#include "stplugin.h"
/* Wire names/defaults belong here; numerical phases consume a typed snapshot. */
typedef struct
{
    ST_int K, G, drop_singletons, verbose, maxiter, standardize, vcetype, has_weights, weight_type,
        compute_dof, dof_adjust_type, dof_clusters, savefe, savefe_var_idx, compute_groupvar,
        groupvar_var_idx, compute_resid, resid_var_idx, use_quad, touse;
    ST_double tolerance, df_a_nested;
    /* reghdfe's collinear-with-FE tolerance, min(1e-6, tolerance()/10) of the
     * user's tolerance (before the savefe override of tolerance) */
    ST_double collinear_tol;
    /* vce(bootstrap): keep the model for replicates (full_regression retain)
     * and the plugin variable indices of the resampling units and strata
     * (0: observations, one stratum) */
    ST_int bs_retain, bs_unit_idx, bs_strata_idx;
    /* Cluster variables after the FE variables: the dimensions of
     * vce(cluster) (vce type 2), or the time variable of vce(dkraay) (vce
     * type 3, with its bandwidth, 0 for the default) */
    ST_int num_cluster, dkraay_bw;
    /* dofadjustments(continuous): slopes constant (or zero) within a level
     * are redundant */
    ST_int dof_continuous;
    /* Heterogeneous slopes: factor g's slope variable in the plugin
     * varlist (0 for an intercept) and the absorbed effect (1-based) it
     * belongs to; num_slopes counts the slope factors */
    ST_int slope_idx[10], fe_absvar[10], num_slopes;
    /* group() and individual(): an individual effect aggregated within the
     * groups (indiv_mean: aggregation(mean), else sum), from the plugin
     * variables of the group number, the individual number and the
     * individual-level sample marker (1-based plugin variable indices) */
    ST_int indiv, indiv_mean, indiv_gidx, indiv_iidx, indiv_tidx;
    /* reghdfe's transform() and acceleration() of the main projection
     * (partial_out_columns_options codes; 0 and 1, symmetric Kaczmarz with
     * CG, by default), and its abort(): 0 keeps a projection that reached
     * iterate() without converging, reported as e(converged) = 0 */
    ST_int transform, acceleration, abort_nonconv;
} creghdfe_input;
ST_retcode creghdfe_read_input(creghdfe_input *input);
#endif
