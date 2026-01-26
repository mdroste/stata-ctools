*! version 0.9.0 26Jan2026
*! creghdfe: C-accelerated high-dimensional fixed effects regression
*! Part of the ctools suite
*!
*! Description:
*!   High-performance replacement for reghdfe using a C plugin
*!   with optimized fixed effects absorption.
*!
*! Syntax:
*!   creghdfe depvar indepvars [if] [in], Absorb(varlist) [options]
*!
*! Options:
*!   absorb(varlist)     - Fixed effects to absorb (required)
*!   vce(cluster varlist) - Clustered standard errors
*!   verbose             - Display progress information

program define creghdfe, eclass
    version 14.0

    * Start wall clock timer immediately
    local __do_timing = 0
    timer clear 98
    timer on 98

    syntax varlist(min=2 fv) [aw fw pw] [if] [in], Absorb(string) [VCE(string) Verbose TIMEit ///
        TOLerance(real 1e-8) MAXiter(integer 10000) NOSTANDardize RESID RESID2(name) RESIDuals(name) ///
        DOFadjustments(string) GROUPvar(name) THReads(integer 0)]

    local __do_timing = ("`verbose'" != "" | "`timeit'" != "")

    * Handle residuals() as alias for resid2()
    if "`residuals'" != "" & "`resid2'" == "" {
        local resid2 "`residuals'"
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
    * Reassign absorb to the varlist portion only
    local absorb "`absorb_vars'"

    * Parse dofadjustments option
    * 0=all (default), 1=none, 2=firstpair, 3=pairwise
    local dof_adjust_type = 0
    if "`dofadjustments'" != "" {
        local dof_lower = lower("`dofadjustments'")
        if "`dof_lower'" == "none" {
            local dof_adjust_type = 1
        }
        else if strpos("`dof_lower'", "first") > 0 {
            local dof_adjust_type = 2
        }
        else if strpos("`dof_lower'", "pair") > 0 {
            local dof_adjust_type = 3
        }
    }

    * Mark sample - include absorb and cluster variables
    marksample touse
    markout `touse' `absorb'

    * Parse weights
    local weight_var ""
    local weight_type = 0   // 0=none, 1=aweight, 2=fweight, 3=pweight
    if "`weight'" != "" {
        * Extract weight variable from expression (remove leading =)
        local weight_var "`exp'"
        local weight_var = subinstr("`weight_var'", "=", "", .)
        local weight_var = trim("`weight_var'")

        if "`weight'" == "aweight" {
            local weight_type = 1
        }
        else if "`weight'" == "fweight" {
            local weight_type = 2
            * fweight must be positive integer
            capture confirm numeric variable `weight_var'
            if _rc != 0 {
                di as error "creghdfe: fweight variable must be numeric"
                exit 198
            }
        }
        else if "`weight'" == "pweight" {
            local weight_type = 3
        }

        * Mark out missing weights
        markout `touse' `weight_var'

        * Check for non-positive weights
        qui count if `weight_var' <= 0 & `touse'
        if r(N) > 0 {
            di as error "creghdfe: weights must be positive"
            exit 198
        }
    }

    * Parse variable list
    gettoken depvar indepvars : varlist

    * Expand factor variables to temporary numeric variables
    if "`indepvars'" != "" {
        * Step 1: Expand factor variables to temp numeric vars
        fvrevar `indepvars' if `touse'
        local indepvars_fvrevar "`r(varlist)'"

        * Step 2: Remove collinear variables (including base levels of factors)
        _rmcoll `indepvars_fvrevar' if `touse', forcedrop
        local indepvars_expanded "`r(varlist)'"

        * Step 3: Get proper coefficient names using fvexpand
        fvexpand `indepvars' if `touse'
        local coef_names_all "`r(varlist)'"

        * Filter out base levels (those with 'b' or 'o' notation)
        local coef_names ""
        foreach v of local coef_names_all {
            if !regexm("`v'", "^[0-9]+b\.") & !regexm("`v'", "^o\.") & "`v'" != "" {
                local coef_names `coef_names' `v'
            }
        }

        * Build expanded varlist for plugin
        local varlist_expanded "`depvar' `indepvars_expanded'"
    }
    else {
        local indepvars_expanded ""
        local varlist_expanded "`depvar'"
        local coef_names ""
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
        local vce_lower = lower(trim(`"`vce'"'))
        if substr("`vce_lower'", 1, 2) == "cl" {
            * cluster option - extract variable name
            local vcetype = 2
            * Parse "cluster varname" or "cl varname"
            gettoken vce_type clustervar : vce
            local clustervar = trim("`clustervar'")
            if "`clustervar'" == "" {
                di as error "creghdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            * Save original cluster variable name for reporting
            local clustervar_orig "`clustervar'"
            * Check if cluster variable is string - if so, convert to numeric
            capture confirm numeric variable `clustervar'
            if _rc != 0 {
                * String variable - create temporary numeric grouping
                * Use gegen if gtools is available (faster), otherwise fall back to egen
                tempvar clustervar_numeric
                capture which gegen
                if _rc == 0 {
                    quietly gegen `clustervar_numeric' = group(`clustervar')
                }
                else {
                    quietly egen `clustervar_numeric' = group(`clustervar')
                }
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
        else if "`vce_lower'" == "robust" | "`vce_lower'" == "r" {
            local vcetype = 1
        }
        else if "`vce_lower'" == "unadjusted" | "`vce_lower'" == "ols" {
            local vcetype = 0
        }
        else {
            di as error "creghdfe: unrecognized vce() option"
            exit 198
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
    capture program list ctools_plugin
    if _rc != 0 {
        local __os = c(os)
        local __machine = c(machine_type)
        local __is_mac = 0
        if "`__os'" == "MacOSX" {
            local __is_mac = 1
        }
        else if strpos(lower("`__machine'"), "mac") > 0 {
            local __is_mac = 1
        }
        local __plugin = ""
        if "`__os'" == "Windows" {
            local __plugin "ctools_windows.plugin"
        }
        else if `__is_mac' {
            local __is_arm = 0
            if strpos(lower("`__machine'"), "apple") > 0 | strpos(lower("`__machine'"), "arm") > 0 | strpos(lower("`__machine'"), "silicon") > 0 {
                local __is_arm = 1
            }
            if `__is_arm' == 0 {
                tempfile __archfile
                quietly shell uname -m > "`__archfile'" 2>&1
                tempname __fh
                file open `__fh' using "`__archfile'", read text
                file read `__fh' __archline
                file close `__fh'
                capture erase "`__archfile'"
                if strpos("`__archline'", "arm64") > 0 {
                    local __is_arm = 1
                }
            }
            if `__is_arm' {
                local __plugin "ctools_mac_arm.plugin"
            }
            else {
                local __plugin "ctools_mac_x86.plugin"
            }
        }
        else if "`__os'" == "Unix" {
            local __plugin "ctools_linux.plugin"
        }
        else {
            local __plugin "ctools.plugin"
        }
        capture program ctools_plugin, plugin using("`__plugin'")
        if _rc != 0 & _rc != 110 & "`__plugin'" != "ctools.plugin" {
            capture program ctools_plugin, plugin using("ctools.plugin")
        }
        if _rc != 0 & _rc != 110 {
            di as error "creghdfe: Could not load ctools plugin"
            exit 601
        }
    }

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
        * Drop existing variable if it exists
        capture drop `resid_varname'
        * Create the residual variable (will be filled by C plugin)
        quietly gen double `resid_varname' = .
    }

    * Set up parameters via Stata scalars (expected by C plugin)
    scalar __creghdfe_K = `nvars'           // depvar + indepvars
    scalar __creghdfe_G = `nfe'             // number of FE groups
    scalar __creghdfe_drop_singletons = 1   // always drop singletons
    scalar __creghdfe_verbose = ("`verbose'" != "")
    scalar __creghdfe_maxiter = `maxiter'
    scalar __creghdfe_tolerance = `tolerance'
    scalar __creghdfe_standardize = ("`nostandardize'" == "")
    scalar __creghdfe_vce_type = `vcetype'
    scalar __creghdfe_compute_dof = 1
    scalar __creghdfe_df_a_nested = 0
    scalar __creghdfe_compute_resid = `compute_resid'
    scalar __creghdfe_has_weights = (`weight_type' > 0)
    scalar __creghdfe_weight_type = `weight_type'
    scalar __creghdfe_dof_adjust_type = `dof_adjust_type'
    scalar __creghdfe_savefe = `savefe'
    scalar __creghdfe_compute_groupvar = ("`groupvar'" != "")

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record setup time
    timer off 98
    quietly timer list 98
    local __time_setup = r(t98)
    timer clear 98
    timer on 98

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
        capture drop `groupvar'
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
            capture drop __hdfe`g'__
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

    * Create matrix to store results (plugin will fill __creghdfe_V)
    local K_x = `nvars' - 1  // number of indepvars (excluding depvar)
    matrix __creghdfe_V = J(`K_x' + 1, `K_x' + 1, 0)

    * Call the C plugin with full_regression subcommand
    capture noisily plugin call ctools_plugin `plugin_varlist' if `touse', ///
        "creghdfe `threads_code' full_regression"

    local reg_rc = _rc
    if `reg_rc' {
        * Clean up scalars on error
        capture scalar drop __creghdfe_K __creghdfe_G __creghdfe_drop_singletons
        capture scalar drop __creghdfe_verbose __creghdfe_maxiter __creghdfe_tolerance
        capture scalar drop __creghdfe_standardize __creghdfe_vce_type __creghdfe_compute_dof
        capture scalar drop __creghdfe_df_a_nested __creghdfe_compute_resid
        capture scalar drop __creghdfe_resid_var_idx
        capture scalar drop __creghdfe_has_weights __creghdfe_weight_type
        capture scalar drop __creghdfe_dof_adjust_type __creghdfe_savefe
        capture scalar drop __creghdfe_compute_groupvar __creghdfe_groupvar_idx __creghdfe_savefe_idx
        capture matrix drop __creghdfe_V
        di as error "Error in regression (rc=`reg_rc')"
        exit `reg_rc'
    }

    * Record plugin call time
    timer off 98
    quietly timer list 98
    local __time_plugin_call = r(t98)
    timer clear 98
    timer on 98

    * Retrieve results from scalars set by C plugin
    local N_final = __creghdfe_N
    local num_singletons = __creghdfe_num_singletons
    local K_keep = __creghdfe_K_keep
    local df_a = __creghdfe_df_a
    local mobility_groups = __creghdfe_mobility_groups
    local rss = __creghdfe_rss
    local tss = __creghdfe_tss
    local tss_within = __creghdfe_tss_within
    local has_cons = __creghdfe_has_cons
    local cons = __creghdfe_cons
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
    timer off 98
    quietly timer list 98
    local __time_postproc = r(t98)
    timer clear 98
    timer on 98

    * Adjust df_a for nested FEs
    local df_a_adjusted = `df_a' - `df_a_nested'
    if `df_a_adjusted' < 0 {
        local df_a_adjusted = 0
    }

    * Calculate R-squared and adjusted R-squared
    local r2 = 1 - `rss' / `tss'
    local r2_within = 1 - `rss' / `tss_within'
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
    local r2_a = 1 - (`rss'/`df_r_ols') / (`tss'/(`N_final'-1))
    local rmse = sqrt(`rss' / `df_r_ols')

    * Calculate F-statistic
    * For unadjusted VCE: use traditional MSS/MSE formula with within-group variation
    * For robust/cluster VCE: use Wald test F = b'*V^-1*b / k
    * e(mss) = TSS - RSS (total model sum of squares for reporting)
    * But F uses within-group MSS = tss_within - rss
    local mss = `tss' - `rss'
    local mss_within = `tss_within' - `rss'
    if `vcetype' == 0 {
        * Traditional F-statistic using within-group variation
        local F = (`mss_within' / `df_m') / (`rss' / `df_r')
    }
    else {
        * Wald F-statistic for robust/cluster VCE - computed after V matrix is built
        local F = .
    }

    * Build coefficient vector and names
    tempname b V

    * Always build full matrix with all X vars (including omitted) to match reghdfe
    local K_full = `K_x' + 1  /* All X vars + constant */
    matrix `b' = J(1, `K_full', 0)
    matrix `V' = J(`K_full', `K_full', 0)

    * Fill in coefficients and VCE for non-collinear variables
    if `K_keep' > 0 {
        * Build mapping from kept index to full index
        local kept_idx = 1
        forval k = 1/`K_x' {
            local is_collin = __creghdfe_collinear_`k'
            if `is_collin' == 0 {
                * Fill coefficient
                matrix `b'[1, `k'] = __creghdfe_beta_`kept_idx'
                * Fill VCE row/column
                local kept_idx2 = 1
                forval j = 1/`K_x' {
                    local is_collin_j = __creghdfe_collinear_`j'
                    if `is_collin_j' == 0 {
                        matrix `V'[`k', `j'] = __creghdfe_V[`kept_idx', `kept_idx2']
                        local ++kept_idx2
                    }
                }
                * Cross terms with constant
                matrix `V'[`k', `K_full'] = __creghdfe_V[`kept_idx', `K_keep'+1]
                matrix `V'[`K_full', `k'] = __creghdfe_V[`K_keep'+1, `kept_idx']
                local ++kept_idx
            }
        }
        * Constant variance
        matrix `V'[`K_full', `K_full'] = __creghdfe_V[`K_keep'+1, `K_keep'+1]
    }
    else {
        * K_keep == 0: set variance for constant only
        local var_cons = `rss' / (`N_final' * `df_r')
        matrix `V'[`K_full', `K_full'] = `var_cons'
    }
    * Constant coefficient
    matrix `b'[1, `K_full'] = `cons'

    * Label the residual variable (already filled by C plugin)
    if "`resid_varname'" != "" {
        label var `resid_varname' "Residuals"
    }

    * Label the groupvar variable (mobility group identifier)
    if "`groupvar'" != "" {
        label var `groupvar' "Mobility group"
    }

    * Label the savefe variables (fixed effect estimates)
    if `savefe' {
        forval g = 1/`nfe' {
            local fevar : word `g' of `absorb'
            label var __hdfe`g'__ "FE for `fevar'"
        }
    }

    * Build column names - include omitted vars with o. prefix
    * Use coef_names which has proper factor variable notation
    local colnames ""
    forval k = 1/`K_x' {
        local vname : word `k' of `coef_names'
        if `K_keep' == 0 {
            local colnames `colnames' o.`vname'
        }
        else {
            local is_collin = __creghdfe_collinear_`k'
            if `is_collin' == 1 {
                local colnames `colnames' o.`vname'
            }
            else {
                local colnames `colnames' `vname'
            }
        }
    }
    local colnames `colnames' _cons

    matrix colnames `b' = `colnames'
    matrix rownames `V' = `colnames'
    matrix colnames `V' = `colnames'

    * Compute Wald F-statistic for robust/cluster VCE
    if `vcetype' != 0 & `K_keep' > 0 {
        * Extract coefficients (excluding constant) and corresponding V submatrix
        tempname b_noconstant V_noconstant Vinv
        matrix `b_noconstant' = `b'[1, 1..`K_keep']
        matrix `V_noconstant' = `V'[1..`K_keep', 1..`K_keep']
        matrix `Vinv' = syminv(`V_noconstant')
        * Wald statistic = b' * V^-1 * b
        tempname Wald
        matrix `Wald' = `b_noconstant' * `Vinv' * `b_noconstant''
        local F = `Wald'[1,1] / `df_m'
    }

    * Record results building time
    timer off 98
    quietly timer list 98
    local __time_build = r(t98)
    timer clear 98
    timer on 98

    * Post results
    ereturn post `b' `V', esample(`touse') depname(`depvar') obs(`N_final')

    * Store additional e() results
    ereturn scalar N = `N_final'
    ereturn scalar df_m = `df_m'
    ereturn scalar df_r = `df_r'
    ereturn scalar r2 = `r2'
    ereturn scalar r2_within = `r2_within'
    ereturn scalar r2_a = `r2_a'
    ereturn scalar rmse = `rmse'
    ereturn scalar rss = `rss'
    ereturn scalar mss = `mss'
    ereturn scalar tss = `tss'
    ereturn scalar tss_within = `tss_within'
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

    ereturn local absorb "`absorb'"
    ereturn local depvar "`depvar'"
    ereturn local indepvars "`indepvars'"
    ereturn local vce = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", "cluster"))
    if "`resid_varname'" != "" {
        ereturn local resid "`resid_varname'"
    }
    ereturn local predict "creghdfe_p"
    ereturn local cmd "creghdfe"
    ereturn local cmdline "creghdfe `0'"

    * Record ereturn posting time
    timer off 98
    quietly timer list 98
    local __time_ereturn = r(t98)
    timer clear 98
    timer on 98

    * Display singleton and convergence messages (before header, like reghdfe)
    if `num_singletons' > 0 {
        di as text "(dropped " as result `num_singletons' as text " singleton observations)"
    }
    di as text "(HDFE estimator converged in " as result `num_iters' as text " iterations)"

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

    * Display absorbed degrees of freedom table
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
        * For nested FEs: all levels are redundant
        * For first non-nested FE: no redundant
        * For subsequent non-nested FEs: 1 redundant (identification)
        if `fe_nested_`g'' == 1 {
            local redundant = `cats'
            local has_nested = 1
            local nested_marker "*"
        }
        else {
            if `g' == 1 | (`has_nested' == 0 & `g' == 1) {
                * First non-nested FE: check if any prior FE was nested
                local first_nonnested = 1
                forval prev = 1/`=`g'-1' {
                    if `fe_nested_`prev'' == 0 {
                        local first_nonnested = 0
                    }
                }
                if `first_nonnested' {
                    local redundant = 0
                }
                else {
                    local redundant = 1
                }
            }
            else {
                local redundant = 1
            }
            local nested_marker " "
        }
        local coefs = `cats' - `redundant'
        local total_coefs = `total_coefs' + `coefs'
        di as text %12s abbrev("`fevar'", 12) " {c |}" as result %10.0fc `cats' as text "  " as result %10.0fc `redundant' as text "  " as result %10.0fc `coefs' as text "   `nested_marker'{c |}"
    }
    di as text "{hline 13}{c BT}{hline 39}{c BRC}"
    if `has_nested' {
        di as text "* = FE nested within cluster; treated as redundant for DoF computation"
    }

    * Clean up scalars
    capture scalar drop __creghdfe_K __creghdfe_G __creghdfe_drop_singletons
    capture scalar drop __creghdfe_verbose __creghdfe_maxiter __creghdfe_tolerance
    capture scalar drop __creghdfe_standardize __creghdfe_vce_type __creghdfe_compute_dof
    capture scalar drop __creghdfe_df_a_nested __creghdfe_compute_resid
    capture scalar drop __creghdfe_resid_var_idx
    capture scalar drop __creghdfe_has_weights __creghdfe_weight_type
    capture scalar drop __creghdfe_dof_adjust_type __creghdfe_savefe
    capture scalar drop __creghdfe_compute_groupvar __creghdfe_groupvar_idx __creghdfe_savefe_idx
    capture scalar drop __creghdfe_N __creghdfe_num_singletons __creghdfe_K_keep
    capture scalar drop __creghdfe_df_a __creghdfe_mobility_groups
    capture scalar drop __creghdfe_rss __creghdfe_tss __creghdfe_tss_within
    capture scalar drop __creghdfe_has_cons __creghdfe_cons
    capture scalar drop __creghdfe_ols_N __creghdfe_df_a_nested_computed
    capture scalar drop __creghdfe_num_collinear __creghdfe_N_clust __creghdfe_iterations
    forval k = 1/20 {
        capture scalar drop __creghdfe_beta_`k'
        capture scalar drop __creghdfe_collinear_`k'
        capture scalar drop __creghdfe_num_levels_`k'
        capture scalar drop __creghdfe_collinear_varnum_`k'
    }
    capture matrix drop __creghdfe_V

    * Clean up nested FE scalars
    forval k = 1/20 {
        capture scalar drop __creghdfe_fe_nested_`k'
    }

    * Record display time and compute total
    timer off 98
    quietly timer list 98
    local __time_display = r(t98)
    local __time_total = `__time_setup' + `__time_plugin_call' + `__time_postproc' + `__time_build' + `__time_ereturn' + `__time_display'

    * Display timing if requested
    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "creghdfe timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f _creghdfe_time_read " sec"
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

        * Clean up timing scalars
        capture scalar drop _creghdfe_time_read _creghdfe_time_singleton
        capture scalar drop _creghdfe_time_dof _creghdfe_time_partial
        capture scalar drop _creghdfe_time_ols _creghdfe_time_vce _creghdfe_time_total
        capture scalar drop _creghdfe_threads_max _creghdfe_openmp_enabled
    }

    timer clear 98

end
