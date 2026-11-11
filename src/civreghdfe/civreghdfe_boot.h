/*
 * civreghdfe_boot.h
 *
 * vce(bootstrap) and fevce(bootstrap) for civreghdfe: the model retained
 * after the full-sample fit and the plugin calls that run replicates on it
 * (docs/development/BOOTSTRAP_PLAN.md, section 9.3).
 * Part of the ctools Stata plugin suite
 */

#ifndef CIVREGHDFE_BOOT_H
#define CIVREGHDFE_BOOT_H

#include "../stplugin.h"
#include "../ctools_types.h"
#include "../ctools_hdfe_utils.h"

/* What the full-sample fit hands over (borrowed; retain copies it). The
 * columns are the compact estimation rows, in the plugin's layout. */
typedef struct civreghdfe_boot_source_s
{
    ST_int N, G;
    ST_int K_endog, K_exog, K_excl; /* as sent: K_exog counts partial() columns */
    const ST_double *y, *X_endog, *X_exog, *Z_excl; /* raw columns, N x K each */
    /* The columns with the fixed effects partialled out (before partial()),
     * N x (K_excl + 1 + K_endog + K_exog): excluded instruments, y,
     * endogenous, exogenous */
    const ST_double *tilde;
    const ST_double *weights; /* N weights as the fit used them, or NULL */
    ST_int weight_type;       /* 1 aw, 2 fw, 3 pw */
    const FE_Factor *factors; /* G factors on the compact rows */
    const ST_int *partial;    /* K_exog partial() flags as the ado sent them, or NULL */
    /* The fit's omitted columns: K_x = K_exog less partial() columns, then
     * K_endog ([exogenous, endogenous]); and the excluded instruments */
    const ST_int *collin_x, *collin_z;
    ST_int endo_in_z; /* endogenous regressors collinear with the instruments */
    ST_int K_b;       /* coefficients posted: [endogenous, exogenous], zeros where omitted */
    const ST_double *b_hat;
    ST_int est_method;
    ST_double kclass_user, fuller_alpha;
    ST_int center, dofminus, maxiter, keep_singletons;
    ST_double tolerance;
    ST_int vce_type, kernel_type, kiefer; /* the fit's VCE: the GMM moment covariance */
    const ST_int *cluster_ids;            /* N 1-based clusters of vce_type 2, or NULL */
    const perm_idx_t *obs_map;            /* N 1-based Stata observations */
    ST_double N_eff;                      /* observations (sum of fweights) */
    ST_int K_iv;                          /* instruments kept (Fuller's N_eff - K_iv) */
} civreghdfe_boot_source;

/* Hooks the fit calls (NULL unless civreghdfe_main installed them), so that
 * the native harness links the fit without the bootstrap */
extern ST_retcode (*civreghdfe_boot_retain_hook)(const civreghdfe_boot_source *src);
extern void (*civreghdfe_boot_release_hook)(void);
extern int civreghdfe_boot_retain_requested;

/* Keep the model for boot_chunk. Reads scalars __civreghdfe_bs_unit_idx and
 * __civreghdfe_bs_strata_idx (plugin variables of the resampling unit and
 * strata, 0 = none), __civreghdfe_bs_s_units (1: the GMM clusters are the
 * resampled units, each drawn copy its own cluster). Sets scalars
 * __civreghdfe_bs_U, __civreghdfe_bs_S and __civreghdfe_bs_par. */
ST_retcode civreghdfe_boot_retain(const civreghdfe_boot_source *src);
void civreghdfe_boot_release(void);

/* boot_chunk first=# count=# total=# (seedhi=# seedlo=# | stata | supplied)
 * [fe] [dse] [ci] [rows]: replicates first .. first+count-1, as creghdfe's
 * (creghdfe_boot.h).
 * fe: matrix __civreghdfe_bs_fe_targets (1 x G) flags the saved effects.
 * Statuses in matrix __civreghdfe_bs_status (1 x count). */
ST_retcode civreghdfe_boot_chunk(const char *args);

/* boot_finish level=# [coef] [mse] [nreps]: the coefficient summary into
 * __civreghdfe_bs_V, _b, _se, _bias, _z0, _ci_normal, _ci_percentile,
 * _ci_bc, _counts and scalars _N_reps, _N_misreps (coef); with fe, the
 * saved effects' standard errors into the call's variables (one per saved
 * effect, the same again with valid-replicate counts with nreps, then d's
 * with dse), matrix __civreghdfe_bs_fe_table (saved effects x 5) and
 * scalars __civreghdfe_bs_fe_split_reps, __civreghdfe_bs_fe_N_reps. */
ST_retcode civreghdfe_boot_finish(const char *args);

/*
 * wildbootstrap(): the WRE bootstrap on the retained fit (civreghdfe_wre.h).
 *
 * wild_setup vcetype=# [vce] [boot] [nonull] [partialled] [projected]:
 * plugin variables the vce() cluster (vce; vcetype 2) and the bootstrap
 * cluster (boot; the vce() cluster without it, each row without either),
 * numbered in ascending order of their values. Tested coefficients: their
 * positions in e(b) in matrix __civreghdfe_wb_test. Sets scalars
 * __civreghdfe_wb_C and __civreghdfe_wb_path (0 nested, 1 projected,
 * 2 partialled, 3 projected although nested).
 *
 * wild_chunk first=# count=# total=# [enum]: the weights of replicates
 * first .. first+count-1, packed in the call's variables (or the
 * enumerated sign vectors).
 *
 * wild_finish level=# ptype=# [ci] [keep]: matrices __civreghdfe_wb_p,
 * _ci, _cistat (the interval search scales by __civreghdfe_wb_scale, the
 * analytic standard errors), scalars __civreghdfe_wb_reps, _misreps and
 * _cached. keep leaves the fit for fevce(bootstrap).
 */
ST_retcode civreghdfe_wild_setup(const char *args);
ST_retcode civreghdfe_wild_chunk(const char *args);
ST_retcode civreghdfe_wild_finish(const char *args);

#endif /* CIVREGHDFE_BOOT_H */
