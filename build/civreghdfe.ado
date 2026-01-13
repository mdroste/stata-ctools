*! version 1.1.0 13jan2026
*! civreghdfe: C-accelerated instrumental variables regression with HDFE
*! Implements 2SLS/IV/LIML/GMM2S with high-dimensional fixed effects absorption

program define civreghdfe, eclass sortpreserve
    version 14.0

    * Parse syntax
    * Basic syntax: civreghdfe depvar (endogvars = instruments) [exogvars], absorb() [options]
    syntax anything(name=0 equalok) [if] [in] [aw fw pw/], ///
        Absorb(varlist) ///
        [vce(string) ///
        TOLerance(real 1e-8) ///
        MAXiter(integer 500) ///
        Verbose TIMEit ///
        FIRST ///
        SMALL ///
        LIML ///
        FULLER(real 0) ///
        Kclass(real 0) ///
        GMM2s ///
        CUE ///
        COVIV ///
        RESiduals(name) ///
        NOConstant]

    * Start timing
    if "`timeit'" != "" {
        local t0 = clock("", "hms")
    }

    * Parse the varlist with parentheses for endogenous vars and instruments
    * Format: depvar (endog1 endog2 = inst1 inst2 inst3) exog1 exog2
    local 0 `"`0'"'

    * Get dependent variable (first token)
    gettoken depvar rest : 0, parse("(")
    local depvar = strtrim("`depvar'")

    * Check for parentheses containing endogenous = instruments
    if strpos("`rest'", "(") == 0 {
        di as error "civreghdfe requires (endogvars = instruments) specification"
        exit 198
    }

    * Extract content within parentheses
    local rest = strtrim("`rest'")
    if substr("`rest'", 1, 1) == "(" {
        local rest = substr("`rest'", 2, .)
    }

    local paren_end = strpos("`rest'", ")")
    if `paren_end' == 0 {
        di as error "Unmatched parenthesis in variable specification"
        exit 198
    }

    local paren_content = substr("`rest'", 1, `paren_end' - 1)
    local exogvars = strtrim(substr("`rest'", `paren_end' + 1, .))

    * Parse endogenous = instruments
    local eq_pos = strpos("`paren_content'", "=")
    if `eq_pos' == 0 {
        di as error "Must specify instruments after '=' in parentheses"
        exit 198
    }

    local endogvars = strtrim(substr("`paren_content'", 1, `eq_pos' - 1))
    local instruments = strtrim(substr("`paren_content'", `eq_pos' + 1, .))

    * Validate variables exist
    foreach v of local depvar {
        confirm numeric variable `v'
    }
    foreach v of local endogvars {
        confirm numeric variable `v'
    }
    foreach v of local instruments {
        confirm numeric variable `v'
    }
    foreach v of local exogvars {
        confirm numeric variable `v'
    }
    foreach v of local absorb {
        confirm variable `v'
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
            di as error "civreghdfe: Could not load ctools plugin"
            exit 601
        }
    }

    * Count variables
    local K_endog : word count `endogvars'
    local K_exog : word count `exogvars'
    local K_inst_excl : word count `instruments'
    local G : word count `absorb'

    * Build full instrument list (exogenous + excluded instruments)
    local all_instruments `exogvars' `instruments'
    local K_iv : word count `all_instruments'

    * Check identification
    if `K_iv' < `K_endog' + `K_exog' {
        di as error "Model is underidentified"
        di as error "  Endogenous variables: `K_endog'"
        di as error "  Exogenous variables: `K_exog'"
        di as error "  Total instruments: `K_iv'"
        exit 481
    }

    * Parse VCE
    local vce_type = 0
    local cluster_var ""
    if `"`vce'"' != "" {
        local vce_lower = lower(`"`vce'"')
        if `"`vce_lower'"' == "robust" | `"`vce_lower'"' == "r" {
            local vce_type = 1
        }
        else if substr(`"`vce_lower'"', 1, 2) == "cl" {
            local vce_type = 2
            * Extract cluster variable
            local cluster_spec `"`vce'"'
            gettoken vce_cmd cluster_var : cluster_spec, parse(" (")
            local cluster_var = strtrim(subinstr("`cluster_var'", "(", "", .))
            local cluster_var = strtrim(subinstr("`cluster_var'", ")", "", .))
            confirm variable `cluster_var'
        }
        else if `"`vce_lower'"' == "unadjusted" | `"`vce_lower'"' == "ols" {
            local vce_type = 0
        }
        else {
            di as error "vce(`vce') not supported. Use vce(robust) or vce(cluster varname)"
            exit 198
        }
    }

    * Handle weights
    local has_weights = 0
    local weight_type = 0
    local weight_var ""
    if "`weight'" != "" {
        local has_weights = 1
        if "`weight'" == "aweight" local weight_type = 1
        else if "`weight'" == "fweight" local weight_type = 2
        else if "`weight'" == "pweight" local weight_type = 3
        else {
            di as error "Weight type `weight' not supported"
            exit 198
        }
        local weight_var `"`exp'"'
        confirm numeric variable `weight_var'
    }

    * Mark sample
    marksample touse
    markout `touse' `depvar' `endogvars' `exogvars' `instruments' `absorb'
    if "`cluster_var'" != "" {
        markout `touse' `cluster_var'
    }
    if `has_weights' {
        markout `touse' `weight_var'
    }

    * Count observations
    qui count if `touse'
    local N = r(N)

    if `N' == 0 {
        di as error "No observations"
        exit 2000
    }

    * Verbose output
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "civreghdfe: C-Accelerated IV Regression with HDFE"
        di as text "{hline 60}"
        di as text "Dependent variable:    `depvar'"
        di as text "Endogenous variables:  `endogvars' (`K_endog')"
        di as text "Exogenous variables:   `exogvars' (`K_exog')"
        di as text "Excluded instruments:  `instruments' (`K_inst_excl')"
        di as text "Absorbed FE:           `absorb' (`G')"
        di as text "Observations:          `N'"
        if `vce_type' == 1 {
            di as text "VCE:                   Robust"
        }
        else if `vce_type' == 2 {
            di as text "VCE:                   Cluster (`cluster_var')"
        }
        else {
            di as text "VCE:                   Unadjusted"
        }
        di as text "{hline 60}"
    }

    * Detect if cluster variable is nested within (same as) an absorb variable
    * When cluster var = FE var, that FE shouldn't reduce DOF for VCE
    local nested_fe_index = 0
    if `vce_type' == 2 & "`cluster_var'" != "" {
        local fe_idx = 1
        foreach fe_var of local absorb {
            if "`fe_var'" == "`cluster_var'" {
                local nested_fe_index = `fe_idx'
            }
            local fe_idx = `fe_idx' + 1
        }
    }

    * Determine estimation method
    * 0 = 2SLS (default)
    * 1 = LIML
    * 2 = Fuller's modified LIML
    * 3 = k-class
    * 4 = GMM2S
    * 5 = CUE
    local est_method = 0
    local kclass_val = 0
    local fuller_val = 0

    if "`liml'" != "" {
        local est_method = 1
        if "`verbose'" != "" {
            di as text "Estimation method:     LIML"
        }
    }
    else if `fuller' != 0 {
        local est_method = 2
        local fuller_val = `fuller'
        if "`verbose'" != "" {
            di as text "Estimation method:     Fuller LIML (alpha=`fuller')"
        }
    }
    else if `kclass' != 0 {
        local est_method = 3
        local kclass_val = `kclass'
        if "`verbose'" != "" {
            di as text "Estimation method:     k-class (k=`kclass')"
        }
    }
    else if "`gmm2s'" != "" {
        local est_method = 4
        if "`verbose'" != "" {
            di as text "Estimation method:     Two-step GMM"
        }
    }
    else if "`cue'" != "" {
        local est_method = 5
        if "`verbose'" != "" {
            di as text "Estimation method:     CUE (Continuously Updated GMM)"
        }
    }
    else {
        if "`verbose'" != "" {
            di as text "Estimation method:     2SLS"
        }
    }

    * Set up scalars for plugin
    scalar __civreghdfe_K_endog = `K_endog'
    scalar __civreghdfe_K_exog = `K_exog'
    scalar __civreghdfe_K_iv = `K_iv'
    scalar __civreghdfe_G = `G'
    scalar __civreghdfe_has_cluster = (`vce_type' == 2)
    scalar __civreghdfe_has_weights = `has_weights'
    scalar __civreghdfe_weight_type = `weight_type'
    scalar __civreghdfe_vce_type = `vce_type'
    scalar __civreghdfe_maxiter = `maxiter'
    scalar __civreghdfe_tolerance = `tolerance'
    scalar __civreghdfe_verbose = ("`verbose'" != "")
    scalar __civreghdfe_nested_fe_index = `nested_fe_index'
    scalar __civreghdfe_est_method = `est_method'
    scalar __civreghdfe_kclass = `kclass_val'
    scalar __civreghdfe_fuller = `fuller_val'
    scalar __civreghdfe_coviv = ("`coviv'" != "")

    * Create output matrices
    local K_total = `K_endog' + `K_exog'
    tempname b_temp V_temp

    matrix `b_temp' = J(1, `K_total', 0)
    matrix `V_temp' = J(`K_total', `K_total', 0)
    matrix __civreghdfe_b = `b_temp'
    matrix __civreghdfe_V = `V_temp'

    * Build variable list for plugin
    * Order: depvar, endogvars, exogvars, all_instruments, absorb, [cluster], [weight]
    local plugin_vars `depvar' `endogvars' `exogvars' `all_instruments' `absorb'
    if `vce_type' == 2 {
        local plugin_vars `plugin_vars' `cluster_var'
    }
    if `has_weights' {
        local plugin_vars `plugin_vars' `weight_var'
    }

    * Call the C plugin
    plugin call ctools_plugin `plugin_vars' if `touse', "civreghdfe iv_regression"

    * Retrieve results
    local N_used = __civreghdfe_N
    local df_r = __civreghdfe_df_r
    local df_a = __civreghdfe_df_a
    local K = __civreghdfe_K

    if `vce_type' == 2 {
        local N_clust = __civreghdfe_N_clust
    }

    * Get coefficient vector and VCE
    matrix `b_temp' = __civreghdfe_b
    matrix `V_temp' = __civreghdfe_V

    * Build variable names for matrices
    * ivreghdfe convention: [endog, exog] - we need to reorder from C's [exog, endog]
    local varnames_c ""
    foreach v of local exogvars {
        local varnames_c `varnames_c' `v'
    }
    foreach v of local endogvars {
        local varnames_c `varnames_c' `v'
    }

    * Reorder to match ivreghdfe: [endog, exog]
    local varnames ""
    foreach v of local endogvars {
        local varnames `varnames' `v'
    }
    foreach v of local exogvars {
        local varnames `varnames' `v'
    }

    * Reorder coefficient vector and VCE matrix
    * C code returns [exog, endog], we want [endog, exog]
    local K_total = `K_exog' + `K_endog'

    if `K_endog' > 0 & `K_exog' > 0 {
        * Need to reorder
        tempname b_reorder V_reorder
        matrix `b_reorder' = J(1, `K_total', 0)
        matrix `V_reorder' = J(`K_total', `K_total', 0)

        * Copy endogenous coefficients (from end of C output to start of reordered)
        forvalues e = 1/`K_endog' {
            local c_idx = `K_exog' + `e'
            matrix `b_reorder'[1, `e'] = `b_temp'[1, `c_idx']
        }

        * Copy exogenous coefficients (from start of C output to end of reordered)
        forvalues x = 1/`K_exog' {
            local new_idx = `K_endog' + `x'
            matrix `b_reorder'[1, `new_idx'] = `b_temp'[1, `x']
        }

        * Reorder VCE matrix
        * Build mapping: new_order[i] = old_index
        * For i = 1..K_endog: old index is K_exog + i
        * For i = K_endog+1..K_total: old index is i - K_endog
        forvalues i = 1/`K_total' {
            if `i' <= `K_endog' {
                local old_i = `K_exog' + `i'
            }
            else {
                local old_i = `i' - `K_endog'
            }
            forvalues j = 1/`K_total' {
                if `j' <= `K_endog' {
                    local old_j = `K_exog' + `j'
                }
                else {
                    local old_j = `j' - `K_endog'
                }
                matrix `V_reorder'[`i', `j'] = `V_temp'[`old_i', `old_j']
            }
        }

        matrix `b_temp' = `b_reorder'
        matrix `V_temp' = `V_reorder'
    }

    * Add column/row names to matrices
    matrix colnames `b_temp' = `varnames'
    matrix rownames `b_temp' = `depvar'
    matrix colnames `V_temp' = `varnames'
    matrix rownames `V_temp' = `varnames'

    * Compute additional statistics
    * RSS is not directly available, compute from residuals if needed

    * Post results
    ereturn clear
    ereturn post `b_temp' `V_temp', esample(`touse') depname(`depvar')

    * Store scalars
    ereturn scalar N = `N_used'
    ereturn scalar df_r = `df_r'
    ereturn scalar df_a = `df_a'
    ereturn scalar K = `K'
    ereturn scalar K_endog = `K_endog'
    ereturn scalar K_exog = `K_exog'
    ereturn scalar K_iv = `K_iv'
    ereturn scalar G = `G'

    if `vce_type' == 2 {
        ereturn scalar N_clust = `N_clust'
    }

    * Store first-stage F statistics
    forvalues e = 1/`K_endog' {
        capture scalar temp_F = __civreghdfe_F1_`e'
        if _rc == 0 {
            ereturn scalar F_first`e' = temp_F
        }
    }

    * Store macros
    ereturn local cmd "civreghdfe"
    ereturn local cmdline `"civreghdfe `0'"'
    ereturn local depvar "`depvar'"
    ereturn local endogvars "`endogvars'"
    ereturn local exogvars "`exogvars'"
    ereturn local instruments "`instruments'"
    ereturn local absorb "`absorb'"

    * Set title and estimation method based on est_method
    if `est_method' == 1 {
        ereturn local title "LIML regression with HDFE"
        ereturn local model "liml"
        capture ereturn scalar lambda = __civreghdfe_lambda
    }
    else if `est_method' == 2 {
        ereturn local title "Fuller LIML regression with HDFE"
        ereturn local model "fuller"
        ereturn scalar fuller = `fuller_val'
        capture ereturn scalar lambda = __civreghdfe_lambda
    }
    else if `est_method' == 3 {
        ereturn local title "k-class regression with HDFE"
        ereturn local model "kclass"
        ereturn scalar kclass = `kclass_val'
    }
    else if `est_method' == 4 {
        ereturn local title "2-step GMM with HDFE"
        ereturn local model "gmm2s"
    }
    else if `est_method' == 5 {
        ereturn local title "CUE regression with HDFE"
        ereturn local model "cue"
    }
    else {
        ereturn local title "IV (2SLS) regression with HDFE"
        ereturn local model "2sls"
    }

    if `vce_type' == 0 {
        ereturn local vcetype "Unadjusted"
    }
    else if `vce_type' == 1 {
        ereturn local vcetype "Robust"
    }
    else if `vce_type' == 2 {
        ereturn local vcetype "Cluster"
        ereturn local clustvar "`cluster_var'"
    }

    * Display results
    di as text ""
    local disp_title = e(title)
    if `vce_type' == 2 {
        di as text "`disp_title'" _col(49) "Number of obs" _col(68) "=" _col(70) %9.0fc `N_used'
        di as text "Absorbed " `G' " HDFE groups" _col(49) "Num. clusters" _col(68) "=" _col(70) %9.0fc `N_clust'
    }
    else {
        di as text "`disp_title'" _col(49) "Number of obs" _col(68) "=" _col(70) %9.0fc `N_used'
        di as text "Absorbed " `G' " HDFE groups"
    }
    di as text ""

    * Display coefficient table
    ereturn display, level(95)

    * Display first-stage F statistics if requested
    if "`first'" != "" | `K_endog' > 0 {
        di as text ""
        di as text "First-stage regression summary:"
        di as text "{hline 60}"
        forvalues e = 1/`K_endog' {
            local endogvar : word `e' of `endogvars'
            capture local F_val = e(F_first`e')
            if _rc == 0 {
                di as text "  `endogvar'" _col(30) "F(" %3.0f `K_inst_excl' ", " %6.0f `df_r' ") = " _col(50) as result %8.2f `F_val'
            }
        }
        di as text "{hline 60}"
        di as text "  Stock-Yogo weak ID test critical values (10% maximal IV size): 16.38"
    }

    * Timing
    if "`timeit'" != "" {
        local t1 = clock("", "hms")
        local elapsed = (`t1' - `t0') / 1000
        di as text ""
        di as text "Total time: " as result %6.3f `elapsed' " seconds"
    }

    * Clean up temp scalars and matrices
    capture scalar drop __civreghdfe_K_endog
    capture scalar drop __civreghdfe_K_exog
    capture scalar drop __civreghdfe_K_iv
    capture scalar drop __civreghdfe_G
    capture scalar drop __civreghdfe_has_cluster
    capture scalar drop __civreghdfe_has_weights
    capture scalar drop __civreghdfe_weight_type
    capture scalar drop __civreghdfe_vce_type
    capture scalar drop __civreghdfe_maxiter
    capture scalar drop __civreghdfe_tolerance
    capture scalar drop __civreghdfe_verbose
    capture scalar drop __civreghdfe_N
    capture scalar drop __civreghdfe_df_r
    capture scalar drop __civreghdfe_df_a
    capture scalar drop __civreghdfe_K
    capture scalar drop __civreghdfe_N_clust
    capture scalar drop __civreghdfe_nested_fe_index
    capture scalar drop __civreghdfe_est_method
    capture scalar drop __civreghdfe_kclass
    capture scalar drop __civreghdfe_fuller
    capture scalar drop __civreghdfe_coviv
    capture scalar drop __civreghdfe_lambda
    forvalues e = 1/`K_endog' {
        capture scalar drop __civreghdfe_F1_`e'
    }
    capture matrix drop __civreghdfe_b
    capture matrix drop __civreghdfe_V

end
