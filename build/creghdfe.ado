*! version 1.2.4 20261103 github.com/mdroste/stata-ctools

program define creghdfe, eclass
    version 14.1
    * Replay never enters the estimation transaction or clears existing e().
    if replay() {
        if "`e(cmd)'" != "creghdfe" error 301
        _creghdfe_display `0'
        exit
    }
    * A failed or interrupted regression leaves the data as they were without
    * a snapshot of the dataset: creghdfe only adds variables at the end
    * (residuals, groupvar(), savefe and temporary variables) and writes into
    * nothing else, so _ctools_undo drops them and resets c(changed). Option
    * snapshot preserves the whole dataset instead (on disk above
    * c(max_preservemem)); option fast undoes nothing.
    _ctools_mode `0'
    * Break is recognized only inside _creghdfe_impl.
    nobreak {
        if "`mode'" == "snapshot" preserve
        else if "`mode'" == "" _ctools_undo record CTOOLS_u_creghdfe_
        capture noisily break _creghdfe_impl `0'
        local rc = _rc
        if `rc' {
            if "`mode'" == "snapshot" restore
            else if "`mode'" == "" _ctools_undo rollback CTOOLS_u_creghdfe_
            ereturn clear
            exit `rc'
        }
        if "`mode'" == "snapshot" restore, not
        else if "`mode'" == "" _ctools_undo clear CTOOLS_u_creghdfe_
    }
end

program define _creghdfe_impl, eclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start wall clock timer immediately
    local __do_timing = 0
    quietly _ctools_clock __ct_lap98

    * Drop whatever an interrupted run left behind
    _ctools_cleanup creghdfe

    * A refit on the singleton-free sample (see after the plugin call) passes
    * that sample and the first run's singleton count here
    local __refit_sample : word 1 of $CTOOLS_creghdfe_refit_
    local __refit_singletons : word 2 of $CTOOLS_creghdfe_refit_
    global CTOOLS_creghdfe_refit_
    local __creghdfe_args `"`0'"'

    syntax varlist(min=1 fv ts) [aw fw pw] [if] [in] [, Absorb(string) NOABsorb VCE(string) ///
        CLuster(string) Verbose ///
        TOLerance(real 1e-8) ITERATE(integer 16000) NOSTANDardize RESID RESID2(name) RESIDuals(name) ///
        DOFadjustments(string) GROUPvar(name) THReads(integer 0) QUAD SNAPshot FAST ///
        noCONstant KEEPSINgletons noHEADer noTABle noFOOTnote *]

    * As in reghdfe: verbose(#) (a positive level gives the timing
    * breakdown), and the display options of estimation commands
    _creghdfe_verbose_opt, `options'
    local verbose = cond("`verbose'" != "" | `__verbose_level' > 0, "verbose", "")
    _get_diopts diopts, `__options'

    * cluster(varlist) is vce(cluster varlist); vce(conventional) is the
    * unadjusted VCE; noabsorb is the default without absorb()
    if `"`cluster'"' != "" {
        if `"`vce'"' != "" {
            di as error "only one of cluster() and vce() can be specified"
            exit 198
        }
        local vce `"cluster `cluster'"'
    }
    if inlist(lower(strtrim(`"`vce'"')), "conventional", "conv") local vce unadjusted
    if "`keepsingletons'" != "" {
        di as error "WARNING: Singleton observations not dropped; statistical significance is biased"
    }

    if `iterate' < 1 | missing(`tolerance') | `tolerance' <= 0 {
        di as error "creghdfe: iterate() and tolerance() must be positive"
        exit 198
    }
    local __do_timing = ("`verbose'" != "")

    * Handle residuals() as alias for resid2()
    if "`residuals'" != "" & "`resid2'" == "" {
        local resid2 "`residuals'"
    }

    * If absorb is not specified, create a constant FE (absorb the intercept only)
    local __noabsorb = 0
    if `"`absorb'"' == "" {
        local __noabsorb = 1
        tempvar _cons_fe
        quietly gen byte `_cons_fe' = 1
        local absorb "`_cons_fe'"
    }

    * Parse absorb option for suboptions (savefe)
    local savefe = 0
    local absorb_vars "`absorb'"
    if strpos(`"`absorb'"', ",") > 0 {
        gettoken absorb_vars absorb_opts : absorb, parse(",")
        _creghdfe_absorb_opts `absorb_opts'
    }

    * Parse dofadjustments option as reghdfe does
    * 0=all (default), 1=none, 2=firstpair, 3=pairwise; all implies clusters
    local dof_lower = lower(`"`dofadjustments'"')
    _ctools_hdfe_dof, `dof_lower'
    local dof_adjust_type = `s(type)'
    local dof_clusters = `s(clusters)'
    local dofmethods `s(methods)'

    * Mark sample - include absorb and cluster variables
    marksample touse
    if "`__refit_sample'" != "" quietly replace `touse' = 0 if !`__refit_sample'

    * One variable per absorbed effect, in order: the FE count and every
    * plugin variable position are computed from this list. A numeric
    * variable goes to the plugin as it is; string variables and
    * interactions are numbered first (the plugin renumbers every effect on
    * the estimation sample). absorb_names holds the effects as reghdfe
    * names them, for e(absorb) and the output.
    _ctools_hdfe_absorb `absorb_vars', touse(`touse') cmd(creghdfe)
    local nfe = r(n)
    forvalues g = 1/`nfe' {
        local absorb_name_`g' `r(name_`g')'
        local absorb_vars_`g' `r(vars_`g')'
        local absorb_group_`g' = r(group_`g')
    }
    local absorb
    local absorb_names
    forvalues g = 1/`nfe' {
        local absorb_names `absorb_names' `absorb_name_`g''
        if `absorb_group_`g'' {
            tempvar absorb_id_`g'
            quietly _ctools_group `absorb_vars_`g'' if `touse', generate(`absorb_id_`g'')
            local absorb `absorb' `absorb_id_`g''
        }
        else local absorb `absorb' `absorb_vars_`g''
    }

    * Parse weights
    local weight_var ""
    local weight_type = 0   // 0=none, 1=aweight, 2=fweight, 3=pweight
    if "`weight'" != "" {
        * Evaluate the weight expression; marks out and validates the weights
        tempvar evaluated_weight
        _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"`exp'"')
        local weight_var "`evaluated_weight'"
        local weight_type = cond("`weight'" == "aweight", 1, cond("`weight'" == "fweight", 2, 3))
    }

    * Parse VCE option. Its cluster variables are marked out before the
    * factor variables are expanded, so that their levels come from the
    * estimation sample.
    local vcetype = 0   // 0=unadjusted, 1=robust, 2=cluster
    local clustervar ""
    local clustervar_orig ""

    * pweight forces robust VCE (unless clustering)
    if `weight_type' == 3 {
        local vcetype = 1
    }

    if `"`vce'"' != "" {
        _ctools_hdfe_vce, spec(`"`vce'"') touse(`touse') cmd(creghdfe)
        local vce_kind `r(kind)'
        if "`vce_kind'" == "cluster" {
            local vcetype = 2
            if r(ndims) == 0 {
                di as error "creghdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            if r(ndims) > 1 {
                di as error "creghdfe: multi-way clustering is not supported; specify one cluster variable in vce(cluster)"
                exit 198
            }
            * Reported as reghdfe names it; a string variable or an
            * interaction is numbered first, as egen group() would
            local clustervar_orig `r(name_1)'
            local clustervar `r(vars_1)'
            if r(group_1) {
                tempvar clustervar_numeric
                quietly _ctools_group `clustervar' if `touse', generate(`clustervar_numeric')
                local clustervar `clustervar_numeric'
            }
        }
        else if "`vce_kind'" == "robust" {
            local vcetype = 1
        }
        else if inlist("`vce_kind'", "unadjusted", "ols") {
            * As in reghdfe: pweights need robust or clustered standard errors
            if `weight_type' == 3 {
                di as error "pweights do not work with vce(unadjusted), use a different vce()"
                exit 9
            }
            local vcetype = 0
        }
        else {
            di as error "creghdfe: unrecognized vce() option"
            exit 9
        }
    }

    * Parse variable list
    gettoken depvar indepvars : varlist

    * Save original depvar name for labeling (before ts expansion)
    local depvar_orig "`depvar'"

    * Expand time-series operators on depvar (creates temp var if needed)
    tsrevar `depvar'
    local depvar "`r(varlist)'"

    * Expand factor variables to temporary numeric variables
    if "`indepvars'" != "" {
        * Expand factor variables to temp numeric vars
        fvrevar `indepvars' if `touse'
        local __indepvars_all "`r(varlist)'"

        * Get coefficient names (same count as fvrevar, one per level)
        fvexpand `indepvars' if `touse'
        local coef_names_all "`r(varlist)'"

        * Base/omitted placeholders are not sent to the plugin, which also
        * detects collinearity (detect_collinearity)
        _ctools_fvsplit, names(`coef_names_all') vars(`__indepvars_all')
        local indepvars_expanded "`r(vars)'"
        local coef_names "`r(names)'"
        local coef_is_base "`r(base)'"

        * Build expanded varlist for plugin
        local varlist_expanded "`depvar' `indepvars_expanded'"
    }
    else {
        local indepvars_expanded ""
        local varlist_expanded "`depvar'"
        local coef_names ""
        local coef_names_all ""
        local coef_is_base ""
    }

    local nvars : word count `varlist_expanded'

    * Count observations
    qui count if `touse'
    local nobs = r(N)

    if `nobs' == 0 {
        di as error "creghdfe: no observations"
        exit 2000
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Handle resid option - determine residual variable name
    local resid_varname ""
    local compute_resid = 0
    if "`resid'" != "" | "`resid2'" != "" {
        local compute_resid = 1
        if "`resid2'" != "" {
            local resid_varname "`resid2'"
        }
        else {
            local resid_varname "_reghdfe_resid"
        }
    }

    * Validate the entire output list before creating any destination.
    local output_names `resid_varname' `groupvar'
    if `savefe' {
        forvalues g=1/`nfe' {
            local output_names `output_names' __hdfe`g'__
        }
    }
    if "`output_names'" != "" _ctools_newvars `output_names'
    if `compute_resid' quietly gen double `resid_varname' = .

    * Set up parameters via Stata scalars (expected by C plugin)
    scalar __creghdfe_K = `nvars'           // depvar + indepvars
    scalar __creghdfe_G = `nfe'             // number of FE groups
    scalar __creghdfe_drop_singletons = ("`keepsingletons'" == "")
    scalar __creghdfe_verbose = ("`verbose'" != "")
    scalar __creghdfe_maxiter = `iterate'
    scalar __creghdfe_tolerance = `tolerance'
    scalar __creghdfe_standardize = ("`nostandardize'" == "")
    scalar __creghdfe_vce_type = `vcetype'
    scalar __creghdfe_compute_dof = 1
    scalar __creghdfe_df_a_nested = 0
    scalar __creghdfe_compute_resid = `compute_resid'
    scalar __creghdfe_has_weights = (`weight_type' > 0)
    scalar __creghdfe_weight_type = `weight_type'
    scalar __creghdfe_dof_adjust_type = `dof_adjust_type'
    scalar __creghdfe_dof_clusters = `dof_clusters'
    scalar __creghdfe_savefe = `savefe'
    scalar __creghdfe_compute_groupvar = ("`groupvar'" != "")
    scalar __creghdfe_use_quad = ("`quad'" != "")

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record setup time
    _ctools_lap __time_setup __ct_lap98 `__ct_lap98'

    * Build varlist for plugin: depvar indepvars fe_vars [cluster_var] [weight_var] [resid_var]
    * Use expanded varlist (factor variables converted to temp numeric vars)
    local plugin_varlist `varlist_expanded' `absorb'
    if `vcetype' == 2 {
        local plugin_varlist `plugin_varlist' `clustervar'
    }
    * Add weight variable to plugin varlist if using weights
    if `weight_type' > 0 {
        local plugin_varlist `plugin_varlist' `weight_var'
    }
    * Add residual variable to the plugin varlist if computing residuals
    if `compute_resid' {
        * Position in varlist = K + G + (1 if cluster) + (1 if weights) + 1
        local resid_var_pos = `nvars' + `nfe' + (`vcetype' == 2) + (`weight_type' > 0) + 1
        local plugin_varlist `plugin_varlist' `resid_varname'
        scalar __creghdfe_resid_var_idx = `resid_var_pos'
    }
    else {
        scalar __creghdfe_resid_var_idx = 0
    }

    * Handle groupvar - create variable and add to plugin varlist
    local groupvar_var_pos = 0
    if "`groupvar'" != "" {
        quietly gen long `groupvar' = .
        local groupvar_var_pos = `nvars' + `nfe' + (`vcetype' == 2) + (`weight_type' > 0) + `compute_resid' + 1
        local plugin_varlist `plugin_varlist' `groupvar'
        scalar __creghdfe_groupvar_idx = `groupvar_var_pos'
    }
    else {
        scalar __creghdfe_groupvar_idx = 0
    }

    * Handle savefe - create FE variables and add to plugin varlist
    if `savefe' {
        forval g = 1/`nfe' {
            quietly gen double __hdfe`g'__ = .
            local plugin_varlist `plugin_varlist' __hdfe`g'__
        }
        * First savefe var position
        local savefe_var_pos = `nvars' + `nfe' + (`vcetype' == 2) + (`weight_type' > 0) + `compute_resid' + ("`groupvar'" != "") + 1
        scalar __creghdfe_savefe_idx = `savefe_var_pos'
    }
    else {
        scalar __creghdfe_savefe_idx = 0
    }

    * The plugin returns the final sample after internal row removals.
    tempvar final_sample
    quietly generate byte `final_sample' = 0
    local plugin_varlist `plugin_varlist' `final_sample'
    local sample_idx : word count `plugin_varlist'
    scalar __creghdfe_sample_idx = `sample_idx'

    * The sample marker goes last, in place of if `touse': the plugin reads
    * it in parallel, where Stata answers an if one observation at a time
    local plugin_varlist `plugin_varlist' `touse'
    scalar __creghdfe_touse = 1

    * Create matrix to store results (plugin will fill __creghdfe_V)
    local K_x = `nvars' - 1  // number of indepvars (excluding depvar)
    matrix __creghdfe_V = J(`K_x' + 1, `K_x' + 1, 0)

    * Validation hook: fail as Break would once the outputs exist
    if "$CTOOLS_creghdfe_fault" == "outputs" error 1

    * Call the C plugin with full_regression subcommand
    capture noisily plugin call ctools_plugin `plugin_varlist', ///
        "creghdfe `threads_code' full_regression"

    local reg_rc = _rc
    if `reg_rc' {
        _ctools_cleanup creghdfe
        di as error "Error in regression (rc=`reg_rc')"
        exit `reg_rc'
    }

    * Record plugin call time
    _ctools_lap __time_plugin_call __ct_lap98 `__ct_lap98'

    * Retrieve results from scalars set by C plugin
    local N_final = __creghdfe_N
    local num_singletons = __creghdfe_num_singletons

    * Factor variables were expanded on the sample before singleton removal.
    * As reghdfe (which expands after its singleton pass), refit on the final
    * sample when that changes the terms, e.g. a level or the base level that
    * occurs only on singleton rows. Singletons are a fixed point, so the
    * refit drops none; it reports the first run's count.
    if `num_singletons' > 0 & "`indepvars'" != "" & "`__refit_sample'" == "" {
        fvexpand `indepvars' if `final_sample'
        if "`r(varlist)'" != "`coef_names_all'" {
            if "`resid_varname'" != "" drop `resid_varname'
            if "`groupvar'" != "" drop `groupvar'
            if `savefe' {
                forvalues g = 1/`nfe' {
                    drop __hdfe`g'__
                }
            }
            _ctools_cleanup creghdfe
            global CTOOLS_creghdfe_refit_ "`final_sample' `num_singletons'"
            _creghdfe_impl `__creghdfe_args'
            exit
        }
    }
    if "`__refit_singletons'" != "" local num_singletons = `__refit_singletons'

    * reghdfe's notes on the regressors it omits, and its warning when the
    * fixed effects explain the dependent variable (which is then zeroed)
    local __coltol = min(1e-6, `tolerance'/10)
    if scalar(__creghdfe_depvar_explained) == 1 {
        di as text "warning: dependent variable {res}`depvar_orig'{txt} is likely perfectly explained by the fixed effects (tol =" %3.1e `__coltol' ")"
    }
    local __k = 0
    foreach __v of local coef_names {
        local ++__k
        if scalar(__creghdfe_fe_collinear_`__k') == 1 {
            di as text "note: {res}`__v'{txt} is probably collinear with the fixed effects (all partialled-out values are close to zero; tol =" %3.1e `__coltol' ")"
        }
        else if scalar(__creghdfe_collinear_`__k') == 1 {
            di as text "note: {res}`__v'{txt} omitted because of collinearity"
        }
    }
    local K_keep = __creghdfe_K_keep
    local df_a = __creghdfe_df_a
    capture local num_iters = __creghdfe_iterations
    if _rc != 0 local num_iters = 1

    * Get number of levels and nested status for each FE (computed by C plugin)
    local df_a_nested = 0
    forval g = 1/`nfe' {
        capture local fe_levels_`g' = __creghdfe_num_levels_`g'
        if _rc != 0 local fe_levels_`g' = .

        * Get nested status from C plugin (replaces slow egen-based check)
        capture local fe_nested_`g' = __creghdfe_fe_nested_`g'
        if _rc != 0 local fe_nested_`g' = 0

        * Accumulate df_a_nested
        if `fe_nested_`g'' == 1 {
            local df_a_nested = `df_a_nested' + `fe_levels_`g''
        }
    }

    * Record post-processing time (nested FE now computed in C plugin)
    _ctools_lap __time_postproc __ct_lap98 `__ct_lap98'

    * Adjust df_a for nested FEs
    local df_a_adjusted = `df_a' - `df_a_nested'
    if `df_a_adjusted' < 0 {
        local df_a_adjusted = 0
    }

    * Calculate R-squared and adjusted R-squared
    * Use scalars for precision-critical computations to avoid local precision loss
    tempname cr_rss cr_tss cr_tss_w cr_mss cr_mss_w
    scalar `cr_rss' = scalar(__creghdfe_rss)
    scalar `cr_tss' = scalar(__creghdfe_tss)
    scalar `cr_tss_w' = scalar(__creghdfe_tss_within)
    local r2 = 1 - scalar(`cr_rss') / scalar(`cr_tss')
    local r2_within = 1 - scalar(`cr_rss') / scalar(`cr_tss_w')
    local df_m = `K_keep'

    * df_r depends on VCE type:
    * - For unadjusted/robust: N - K - df_a
    * - For cluster: min(num_clusters - 1, N - K - df_a) (like reghdfe)
    local df_r_ols_base = `N_final' - `K_keep' - `df_a_adjusted'
    if `vcetype' == 2 {
        * Cluster VCE: df_r = min(num_clusters - 1, df_r_ols)
        * reghdfe uses the minimum to ensure internal consistency
        local num_clusters = __creghdfe_N_clust
        local df_r_cluster = `num_clusters' - 1
        if `df_r_ols_base' < `df_r_cluster' {
            local df_r = `df_r_ols_base'
        }
        else {
            local df_r = `df_r_cluster'
        }
    }
    else {
        * Unadjusted or robust: df_r = N - K - df_a
        local df_r = `df_r_ols_base'
    }
    * A saturated model has no residual degrees of freedom, as in reghdfe
    if `df_r' < 0 {
        local df_r = 0
    }

    * R2_a and RMSE use df_r adjusted for nested FEs (like reghdfe)
    * used_df_r = N - df_a - df_m - df_a_nested
    * (both missing or rss itself without residual degrees of freedom)
    local df_r_ols = `N_final' - `df_a_adjusted' - `K_keep' - `df_a_nested'
    local r2_a = 1 - (scalar(`cr_rss')/`df_r_ols') / (scalar(`cr_tss')/(`N_final'-1))
    local rmse = cond(`df_r_ols' == 0, sqrt(scalar(`cr_rss')), sqrt(scalar(`cr_rss') / `df_r_ols'))
    * reghdfe's r2_a_within, ll and ll_0
    local r2_a_within = 1 - (scalar(`cr_rss')/`df_r_ols') / (scalar(`cr_tss_w')/(`df_r_ols' + `K_keep'))
    local ll = -0.5 * `N_final' * (1 + ln(2*_pi) + ln(scalar(`cr_rss') / `N_final'))
    local ll_0 = -0.5 * `N_final' * (1 + ln(2*_pi) + ln(scalar(`cr_tss_w') / `N_final'))

    * Calculate F-statistic
    * For unadjusted VCE: use traditional MSS/MSE formula with within-group variation
    * For robust/cluster VCE: use Wald test F = b'*V^-1*b / k
    * e(mss) = TSS - RSS (total model sum of squares for reporting)
    * But F uses within-group MSS = tss_within - rss
    scalar `cr_mss' = scalar(`cr_tss') - scalar(`cr_rss')
    scalar `cr_mss_w' = scalar(`cr_tss_w') - scalar(`cr_rss')

    * e(b) and e(V) as reghdfe posts them: every term in order, with zeros for
    * omitted regressors (o. prefix) and for factor-variable base levels
    tempname b V b_kept
    _ctools_ebv `b' `V', names(`coef_names_all') base(`coef_is_base') ///
        stub(__creghdfe) kkeep(`K_keep') vmatrix(__creghdfe_V) cons(__creghdfe_cons) kept(`b_kept')

    * Label the residual variable (already filled by C plugin)
    if "`resid_varname'" != "" {
        label var `resid_varname' "Residuals"
    }

    * Label the groupvar variable (mobility group identifier), as reghdfe
    * labels the groups of the first pair of absorbed effects not nested in
    * the clusters. Where reghdfe forms no such pair (effects nested in the
    * clusters, dofadjustments() without mobility groups) it creates no
    * variable, and neither does creghdfe.
    if "`groupvar'" != "" {
        if scalar(__creghdfe_groupvar_created) == 0 {
            drop `groupvar'
        }
        else if scalar(__creghdfe_groupvar_fe2) > 0 {
            local __fe1 : word `=scalar(__creghdfe_groupvar_fe1)' of `absorb_names'
            local __fe2 : word `=scalar(__creghdfe_groupvar_fe2)' of `absorb_names'
            label var `groupvar' "Mobility group between `__fe1' and `__fe2'"
        }
        else label var `groupvar' "Mobility group"
    }

    * Label the savefe variables (fixed effect estimates)
    if `savefe' {
        forval g = 1/`nfe' {
            local fevar : word `g' of `absorb_names'
            label var __hdfe`g'__ "FE for `fevar'"
        }
    }

    * The Wald F-statistic of all kept regressors, computed by the plugin as
    * reghdfe computes it for every VCE (in the units of the standardized
    * regression; missing when V is rank deficient: too few clusters, a
    * perfect fit, a saturated model)
    local F = scalar(__creghdfe_F)

    * noconstant: as reghdfe, _cons is not reported (the constant is still
    * absorbed with the fixed effects)
    if "`constant'" == "noconstant" & colsof(`b') > 1 {
        local __kb = colsof(`b') - 1
        matrix `b' = `b'[1, 1..`__kb']
        matrix `V' = `V'[1..`__kb', 1..`__kb']
    }

    * Record results building time
    _ctools_lap __time_build __ct_lap98 `__ct_lap98'

    * Post results. With factor variables, buildfvinfo lets margins check
    * estimability, as after reghdfe.
    * buildfvinfo costs a pass over the data and is only
    * needed when a coefficient other than a base level is omitted; with a
    * full-rank design every margin is estimable.
    local bfv
    _ms_op_info `b'
    if r(fvops) {
        local __bfv_names : colnames `b'
        foreach __bfv_name of local __bfv_names {
            if strpos("`__bfv_name'", "o.") & !strpos("`__bfv_name'", "b.") local bfv buildfvinfo
        }
    }
    * As reghdfe: a coefficient with zero b and zero variance is omitted
    _ms_findomitted `b' `V'
    ereturn post `b' `V', esample(`final_sample') depname(`depvar_orig') obs(`N_final') `bfv'
    * The weights, which margins and other postestimation commands reuse
    if "`weight'" != "" {
        ereturn local wtype "`weight'"
        ereturn local wexp `"`exp'"'
    }

    * Store additional e() results
    ereturn scalar N = `N_final'
    ereturn scalar df_m = `df_m'
    ereturn scalar df_r = `df_r'
    ereturn scalar r2 = `r2'
    ereturn scalar r2_within = `r2_within'
    ereturn scalar r2_a = `r2_a'
    ereturn scalar rmse = `rmse'
    ereturn scalar rss = scalar(`cr_rss')
    ereturn scalar mss = scalar(`cr_mss')
    ereturn scalar tss = scalar(`cr_tss')
    ereturn scalar tss_within = scalar(`cr_tss_w')
    ereturn scalar df_a = `df_a_adjusted'
    ereturn scalar df_a_nested = `df_a_nested'
    ereturn scalar F = `F'
    ereturn scalar rank = `df_m'
    ereturn scalar num_singletons = `num_singletons'
    ereturn scalar drop_singletons = ("`keepsingletons'" == "")
    * Likelihood, fit and solver results reghdfe posts
    ereturn scalar r2_a_within = `r2_a_within'
    ereturn scalar ll = `ll'
    ereturn scalar ll_0 = `ll_0'
    ereturn scalar sumweights = scalar(__creghdfe_sumweights)
    ereturn scalar ic = `num_iters'
    ereturn scalar converged = 1
    ereturn scalar report_constant = ("`constant'" != "noconstant")

    * The absorbed effects (creghdfe has no slopes to add), the absorbed
    * degrees of freedom table and the covariance labels, as reghdfe posts
    * them; the model without absorb() absorbs _cons
    if `__noabsorb' {
        ereturn local absorb ""
        local absvars _cons
    }
    else {
        ereturn local absorb "`absorb_names'"
        local absvars `absorb_names'
    }
    local vcename = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", "cluster"))
    local cluster_post
    if `vcetype' == 2 {
        local cluster_post clustvars(`clustervar_orig') nclust(`num_clusters')
    }
    _ctools_hdfe_post, stub(__creghdfe) absvars(`absvars') nhdfe(`nfe') vce(`vcename') ///
        dofmethod(`dofmethods') `cluster_post'
    * reghdfe calls the unadjusted VCE ols
    if `vcetype' == 0 ereturn local vce "ols"
    ereturn local model "ols"
    ereturn local marginsnotok "Residuals SCore"
    ereturn local title "HDFE Linear regression"
    ereturn local title2 = "Absorbing `nfe' HDFE " + cond(`nfe' > 1, "groups", "group")
    ereturn local depvar "`depvar_orig'"
    ereturn local indepvars "`indepvars'"
    if "`resid_varname'" != "" {
        ereturn local resid "`resid_varname'"
    }
    ereturn local predict "creghdfe_p"
    ereturn local cmd "creghdfe"
    ereturn local cmdline "creghdfe `0'"

    * Record ereturn posting time
    _ctools_lap __time_ereturn __ct_lap98 `__ct_lap98'

    * Display singleton and convergence messages (before header, like reghdfe)
    if `num_singletons' > 0 {
        di as text "(dropped " as result `num_singletons' as text " singleton observations)"
    }
    di as text `"({browse "http://scorreia.com/research/hdfe.pdf":MWFE estimator} converged in "' as result `num_iters' as text " iterations)"

    * Header, coefficient table and absorbed degrees of freedom, as replay
    * shows them
    _creghdfe_display, `header' `table' `footnote' `diopts'

    * Record display time and compute total
    _ctools_lap __time_display __ct_lap98 `__ct_lap98'
    local __time_total = `__time_setup' + `__time_plugin_call' + `__time_postproc' + `__time_build' + `__time_ereturn' + `__time_display'

    * Display timing if requested
    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "creghdfe timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load (SPI):        " as result %8.4f _creghdfe_time_load " sec"
        di as text "    Data copy + sums:       " as result %8.4f _creghdfe_time_copy " sec"
        di as text "    FE remap + counts:      " as result %8.4f _creghdfe_time_remap " sec"
        di as text "    Singleton removal:      " as result %8.4f _creghdfe_time_singleton " sec"
        di as text "    DOF computation:        " as result %8.4f _creghdfe_time_dof " sec"
        di as text "    HDFE partial out:       " as result %8.4f _creghdfe_time_partial " sec"
        di as text "    OLS estimation:         " as result %8.4f _creghdfe_time_ols " sec"
        di as text "    VCE computation:        " as result %8.4f _creghdfe_time_vce " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _creghdfe_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Setup (parsing, etc):   " as result %8.4f `__time_setup' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f (`__time_plugin_call' - _creghdfe_time_total) " sec"
        di as text "    Post-processing:        " as result %8.4f `__time_postproc' " sec"
        di as text "    Results building:       " as result %8.4f `__time_build' " sec"
        di as text "    ereturn posting:        " as result %8.4f `__time_ereturn' " sec"
        di as text "    Display output:         " as result %8.4f `__time_display' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_setup' + (`__time_plugin_call' - _creghdfe_time_total) + `__time_postproc' + `__time_build' + `__time_ereturn' + `__time_display'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = _creghdfe_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _creghdfe_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    _ctools_cleanup creghdfe
    * Validation hook: fail after the plugin filled the outputs and e() was posted
    if "$CTOOLS_creghdfe_fault" == "post" error 1
end

* Results as reghdfe shows them, after estimation and on replay: the header,
* the coefficient table (display options as for any estimation command) and
* the absorbed degrees of freedom (e(dof_table)), each of which noheader,
* notable and nofootnote suppress
program define _creghdfe_display
    version 14.1
    syntax [, noHEADer noTABle noFOOTnote *]
    _get_diopts diopts, `options'
    if "`header'" == "" {
        di as text ""
        di as text "`e(title)'" _col(49) "Number of obs" _col(67) "= " as result %10.0fc e(N)
        di as text "`e(title2)'" _col(49) "F(" as result %3.0f e(df_m) as text "," as result %8.0f e(df_r) as text ")" _col(67) "= " as result %10.2f e(F)
        di as text _col(49) "Prob > F" _col(67) "= " as result %10.4f Ftail(e(df_m), e(df_r), e(F))
        di as text _col(49) "R-squared" _col(67) "= " as result %10.4f e(r2)
        di as text _col(49) "Adj R-squared" _col(67) "= " as result %10.4f e(r2_a)
        di as text _col(49) "Within R-sq." _col(67) "= " as result %10.4f e(r2_within)
        di as text _col(49) "Root MSE" _col(67) "= " as result %10.4f e(rmse)
    }
    if "`table'" == "" {
        di as text ""
        ereturn display, `diopts'
    }
    * The absorbed degrees of freedom (only when absorb() was specified)
    if "`footnote'" == "" & "`e(absorb)'" != "" {
        tempname dof
        matrix `dof' = e(dof_table)
        di as text ""
        di as text "Absorbed degrees of freedom:"
        di as text "{hline 13}{c TT}{hline 39}{c TRC}"
        di as text " Absorbed FE {c |} Categories  - Redundant  = Num. Coefs {c |}"
        di as text "{hline 13}{c +}{hline 39}{c RT}"
        local has_nested = 0
        forvalues g = 1/`=rowsof(`dof')' {
            local fevar : word `g' of `e(absvars)'
            * Redundant coefficients as counted by the plugin: every level of
            * an FE nested within the clusters, otherwise reghdfe's
            * mobility-group count
            local nested_marker " "
            if `dof'[`g', 5] == 1 {
                local has_nested = 1
                local nested_marker "*"
            }
            di as text %12s abbrev("`fevar'", 12) " {c |}" as result %10.0fc `dof'[`g', 1] as text "  " as result %10.0fc `dof'[`g', 2] as text "  " as result %10.0fc `dof'[`g', 3] as text "    `nested_marker'{c |}"
        }
        di as text "{hline 13}{c BT}{hline 39}{c BRC}"
        if `has_nested' {
            di as text "* = FE nested within cluster; treated as redundant for DoF computation"
        }
    }
end

* verbose(#) among the remaining options, returned in caller locals
program define _creghdfe_verbose_opt
    version 14.1
    syntax [, VERbose(integer 0) *]
    c_local __verbose_level `verbose'
    c_local __options `"`options'"'
end

* absorb(absvars, savefe): reghdfe's only suboption; anything else is an error
program define _creghdfe_absorb_opts
    version 14.1
    syntax [, SAVEfe]
    c_local savefe = ("`savefe'" != "")
end
