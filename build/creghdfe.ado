*! version 1.0.2 9feb2026 github.com/mdroste/stata-ctools

program define creghdfe, eclass
    version 14.1
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

    syntax varlist(min=2 fv ts) [aw fw pw] [if] [in] [, Absorb(string) VCE(string) Verbose ///
        TOLerance(real 1e-8) ITERATE(integer 10000) NOSTANDardize RESID RESID2(name) RESIDuals(name) ///
        DOFadjustments(string) GROUPvar(name) THReads(integer 0) QUAD SNAPshot FAST]

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
    local absorb_vars ""
    if strpos(`"`absorb'"', ",") > 0 {
        gettoken absorb_vars absorb_opts : absorb, parse(",")
        local absorb_opts = subinstr("`absorb_opts'", ",", "", 1)
        local absorb_opts = trim("`absorb_opts'")
        if strpos(lower("`absorb_opts'"), "savefe") > 0 {
            local savefe = 1
        }
        local absorb_vars = trim("`absorb_vars'")
    }
    else {
        local absorb_vars "`absorb'"
    }
    * Reassign absorb to the varlist portion only. Expand wildcards and
    * ranges (fe*, fe1-fe3): the FE count and every plugin variable position
    * are computed from this list, while markout and plugin call expand it.
    unab absorb : `absorb_vars'

    * Parse dofadjustments option as reghdfe does
    * 0=all (default), 1=none, 2=firstpair, 3=pairwise; all implies clusters
    local dof_lower = lower(`"`dofadjustments'"')
    _creghdfe_parse_dof, `dof_lower'
    local dof_adjust_type = `s(type)'
    local dof_clusters = `s(clusters)'

    * Mark sample - include absorb and cluster variables
    marksample touse
    markout `touse' `absorb'

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
    local nfe : word count `absorb'

    * Parse VCE option
    local vcetype = 0   // 0=unadjusted, 1=robust, 2=cluster
    local clustervar ""
    local clustervar_orig ""
    local clustervar_is_temp = 0

    * pweight forces robust VCE (unless clustering)
    if `weight_type' == 3 {
        local vcetype = 1
    }

    if `"`vce'"' != "" {
        _ctools_vce vce_kind clustervar : `vce'
        if "`vce_kind'" == "cluster" {
            local vcetype = 2
            if "`clustervar'" == "" {
                di as error "creghdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            unab clustervar : `clustervar'
            if `: word count `clustervar'' > 1 {
                di as error "creghdfe: multi-way clustering is not supported; specify one cluster variable in vce(cluster)"
                exit 198
            }
            * Save original cluster variable name for reporting
            local clustervar_orig "`clustervar'"
            * Check if cluster variable is string - if so, convert to numeric
            capture confirm numeric variable `clustervar'
            if _rc != 0 {
                * String variable - number its values in C, as egen group() would
                tempvar clustervar_numeric
                quietly _ctools_group `clustervar', generate(`clustervar_numeric')
                local clustervar "`clustervar_numeric'"
                local clustervar_is_temp = 1
            }
            * Mark out cluster variable from sample (only for numeric vars)
            * String variables are handled by marking missing string values
            capture confirm numeric variable `clustervar_orig'
            if _rc == 0 {
                markout `touse' `clustervar_orig'
            }
            else {
                * For string variables, mark out empty strings
                quietly replace `touse' = 0 if `clustervar_orig' == "" & `touse' == 1
            }
        }
        else if "`vce_kind'" == "robust" {
            local vcetype = 1
        }
        else if inlist("`vce_kind'", "unadjusted", "ols") {
            local vcetype = 0
        }
        else {
            di as error "creghdfe: unrecognized vce() option"
            exit 9
        }
    }

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
    scalar __creghdfe_drop_singletons = 1   // always drop singletons
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

    * Create matrix to store results (plugin will fill __creghdfe_V)
    local K_x = `nvars' - 1  // number of indepvars (excluding depvar)
    matrix __creghdfe_V = J(`K_x' + 1, `K_x' + 1, 0)

    * Validation hook: fail as Break would once the outputs exist
    if "$CTOOLS_creghdfe_fault" == "outputs" error 1

    * Call the C plugin with full_regression subcommand
    capture noisily plugin call ctools_plugin `plugin_varlist' if `touse', ///
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
    local K_keep = __creghdfe_K_keep
    local df_a = __creghdfe_df_a
    local mobility_groups = __creghdfe_mobility_groups
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
    if `df_r' < 1 {
        local df_r = 1
    }

    * R2_a and RMSE use df_r adjusted for nested FEs (like reghdfe)
    * used_df_r = N - df_a - df_m - df_a_nested
    local df_r_ols = `N_final' - `df_a_adjusted' - `K_keep' - `df_a_nested'
    if `df_r_ols' < 1 {
        local df_r_ols = 1
    }
    local r2_a = 1 - (scalar(`cr_rss')/`df_r_ols') / (scalar(`cr_tss')/(`N_final'-1))
    local rmse = sqrt(scalar(`cr_rss') / `df_r_ols')

    * Calculate F-statistic
    * For unadjusted VCE: use traditional MSS/MSE formula with within-group variation
    * For robust/cluster VCE: use Wald test F = b'*V^-1*b / k
    * e(mss) = TSS - RSS (total model sum of squares for reporting)
    * But F uses within-group MSS = tss_within - rss
    scalar `cr_mss' = scalar(`cr_tss') - scalar(`cr_rss')
    scalar `cr_mss_w' = scalar(`cr_tss_w') - scalar(`cr_rss')
    if `vcetype' == 0 {
        * Traditional F-statistic using within-group variation
        local F = (scalar(`cr_mss_w') / `df_m') / (scalar(`cr_rss') / `df_r')
    }
    else {
        * Wald F-statistic for robust/cluster VCE - computed after V matrix is built
        local F = .
    }

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
    * labels the groups of its first pair of absorbed effects
    if "`groupvar'" != "" {
        local __fe1 : word 1 of `absorb'
        local __fe2 : word 2 of `absorb'
        if "`__fe2'" != "" label var `groupvar' "Mobility group between `__fe1' and `__fe2'"
        else label var `groupvar' "Mobility group"
    }

    * Label the savefe variables (fixed effect estimates)
    if `savefe' {
        forval g = 1/`nfe' {
            local fevar : word `g' of `absorb'
            label var __hdfe`g'__ "FE for `fevar'"
        }
    }

    * Compute Wald F-statistic for robust/cluster VCE
    * Test all kept regressors jointly, as reghdfe does. Omitted regressors may
    * precede kept ones in b/V, so use the plugin's compact results: the kept
    * coefficients in order and the top-left K_keep x K_keep block of its V.
    if `vcetype' != 0 & `K_keep' > 0 {
        tempname V_noconstant Vinv
        matrix `V_noconstant' = __creghdfe_V[1..`K_keep', 1..`K_keep']
        matrix `Vinv' = syminv(`V_noconstant')
        * A rank-deficient V (e.g. too few clusters) gives a missing F, as in reghdfe
        if diag0cnt(`Vinv') > 0 {
            local F = .
        }
        else {
            * Wald statistic = b' * V^-1 * b
            tempname Wald
            matrix `Wald' = `b_kept' * `Vinv' * `b_kept''
            local F = `Wald'[1,1] / `df_m'
        }
    }

    * Record results building time
    _ctools_lap __time_build __ct_lap98 `__ct_lap98'

    * Post results
    ereturn post `b' `V', esample(`final_sample') depname(`depvar_orig') obs(`N_final')

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
    ereturn scalar df_a_initial = `df_a'
    ereturn scalar df_a_nested = `df_a_nested'
    ereturn scalar F = `F'
    ereturn scalar rank = `df_m'
    ereturn scalar N_hdfe = `nfe'
    ereturn scalar N_hdfe_extended = `mobility_groups'
    ereturn scalar num_singletons = `num_singletons'

    if `vcetype' == 2 {
        ereturn scalar N_clust = __creghdfe_N_clust
        ereturn local clustvar "`clustervar_orig'"
    }

    if `__noabsorb' {
        ereturn local absorb ""
    }
    else {
        ereturn local absorb "`absorb'"
    }
    ereturn local depvar "`depvar_orig'"
    ereturn local indepvars "`indepvars'"
    ereturn local vce = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", "cluster"))
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

    * Display results header
    di as text ""
    di as text "HDFE Linear regression" _col(49) "Number of obs" _col(67) "= " as result %10.0fc `N_final'
    di as text "Absorbing " as result `nfe' as text " HDFE " cond(`nfe'>1, "groups", "group") _col(49) "F(" as result %3.0f `df_m' as text "," as result %8.0f `df_r' as text ")" _col(67) "= " as result %10.2f e(F)
    di as text _col(49) "Prob > F" _col(67) "= " as result %10.4f Ftail(`df_m', `df_r', e(F))
    di as text _col(49) "R-squared" _col(67) "= " as result %10.4f `r2'
    di as text _col(49) "Adj R-squared" _col(67) "= " as result %10.4f `r2_a'
    di as text _col(49) "Within R-sq." _col(67) "= " as result %10.4f `r2_within'
    di as text _col(49) "Root MSE" _col(67) "= " as result %10.4f `rmse'

    * Display coefficient table
    di as text ""
    ereturn display

    * Display absorbed degrees of freedom table (only when absorb was specified)
    if !`__noabsorb' {
    di as text ""
    di as text "Absorbed degrees of freedom:"
    di as text "{hline 13}{c TT}{hline 39}{c TRC}"
    di as text " Absorbed FE {c |} Categories  - Redundant  = Num. Coefs {c |}"
    di as text "{hline 13}{c +}{hline 39}{c RT}"
    local total_coefs = 0
    local has_nested = 0
    forval g = 1/`nfe' {
        local fevar : word `g' of `absorb'
        local cats = `fe_levels_`g''
        * Redundant coefficients as counted by the plugin: every level of an
        * FE nested within the clusters, otherwise reghdfe's mobility-group count
        local redundant = __creghdfe_fe_redundant_`g'
        local nested_marker " "
        if `fe_nested_`g'' == 1 {
            local has_nested = 1
            local nested_marker "*"
        }
        local coefs = `cats' - `redundant'
        local total_coefs = `total_coefs' + `coefs'
        di as text %12s abbrev("`fevar'", 12) " {c |}" as result %10.0fc `cats' as text "  " as result %10.0fc `redundant' as text "  " as result %10.0fc `coefs' as text "    `nested_marker'{c |}"
    }
    di as text "{hline 13}{c BT}{hline 39}{c BRC}"
    if `has_nested' {
        di as text "* = FE nested within cluster; treated as redundant for DoF computation"
    }
    } /* end if !__noabsorb */

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

* dofadjustments(): all (the default), none, or any of firstpair|pairwise,
* clusters and continuous, with reghdfe's abbreviations and exclusions.
* Returns s(type): 0=all, 1=none, 2=firstpair, 3=pairwise, and s(clusters).
* continuous concerns slope coefficients, which creghdfe does not absorb.
program define _creghdfe_parse_dof, sclass
    version 14.1
    syntax [, ALL NONE FIRSTpair PAIRwise CLusters CONTinuous]
    if "`none'`firstpair'`pairwise'`clusters'`continuous'" == "" local all all
    opts_exclusive "`all' `none' `firstpair' `pairwise'" dofadjustments
    opts_exclusive "`all' `none' `clusters'" dofadjustments
    opts_exclusive "`all' `none' `continuous'" dofadjustments
    sreturn local type = cond("`all'" != "", 0, cond("`firstpair'" != "", 2, cond("`pairwise'" != "", 3, 1)))
    sreturn local clusters = "`all'`clusters'" != ""
end
