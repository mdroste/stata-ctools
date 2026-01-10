*! version 1.0.0 09jan2026
*! civreghdfe: C-accelerated instrumental variables regression with HDFE
*! Implements 2SLS/IV with high-dimensional fixed effects absorption

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
        SMALL]

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
    local varnames ""
    foreach v of local exogvars {
        local varnames `varnames' `v'
    }
    foreach v of local endogvars {
        local varnames `varnames' `v'
    }

    * Note: The C code returns coefficients in order [exog, endog]
    * But Stata convention is often [endog, exog] - let's keep C order for simplicity

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
    ereturn local title "IV regression with HDFE"

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
    if `vce_type' == 2 {
        di as text "IV regression with HDFE" _col(49) "Number of obs" _col(68) "=" _col(70) %9.0fc `N_used'
        di as text "Absorbed " `G' " HDFE groups" _col(49) "Num. clusters" _col(68) "=" _col(70) %9.0fc `N_clust'
    }
    else {
        di as text "IV regression with HDFE" _col(49) "Number of obs" _col(68) "=" _col(70) %9.0fc `N_used'
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
    forvalues e = 1/`K_endog' {
        capture scalar drop __civreghdfe_F1_`e'
    }
    capture matrix drop __civreghdfe_b
    capture matrix drop __civreghdfe_V

end
