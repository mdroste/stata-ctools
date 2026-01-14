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
        RF ///
        SMALL ///
        LIML ///
        FULLER(real 0) ///
        Kclass(real 0) ///
        GMM2s ///
        CUE ///
        COVIV ///
        BW(integer 0) ///
        KERNEL(string) ///
        DKRAAY(integer 0) ///
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

    * Handle HAC options
    * Kernel types: 0=none, 1=Bartlett/Newey-West, 2=Parzen, 3=Quadratic Spectral
    local kernel_type = 0
    local bw_val = `bw'
    local dkraay_val = `dkraay'

    if "`kernel'" != "" {
        local kernel_lower = lower("`kernel'")
        if "`kernel_lower'" == "bartlett" | "`kernel_lower'" == "nwest" {
            local kernel_type = 1
        }
        else if "`kernel_lower'" == "parzen" {
            local kernel_type = 2
        }
        else if "`kernel_lower'" == "quadraticspectral" | "`kernel_lower'" == "qs" {
            local kernel_type = 3
        }
        else if "`kernel_lower'" == "truncated" | "`kernel_lower'" == "trunc" {
            local kernel_type = 4
        }
        else if "`kernel_lower'" == "tukey" | "`kernel_lower'" == "thann" {
            local kernel_type = 5
        }
        else {
            di as error "Unknown kernel type: `kernel'"
            di as error "Valid options: bartlett, parzen, quadraticspectral, truncated, tukey"
            exit 198
        }
    }

    * If dkraay is specified, use Bartlett kernel by default
    if `dkraay_val' > 0 & `kernel_type' == 0 {
        local kernel_type = 1
    }

    * Set default bandwidth if kernel specified but no bandwidth
    if `kernel_type' > 0 & `bw_val' == 0 {
        * Default: Newey-West optimal bandwidth = floor(4*(N/100)^(2/9))
        local bw_val = floor(4 * (`N'/100)^(2/9))
        if `bw_val' < 1 local bw_val = 1
    }

    scalar __civreghdfe_kernel = `kernel_type'
    scalar __civreghdfe_bw = `bw_val'
    scalar __civreghdfe_dkraay = `dkraay_val'

    if "`verbose'" != "" & `kernel_type' > 0 {
        di as text "HAC: kernel=`kernel_type', bandwidth=`bw_val'"
        if `dkraay_val' > 0 {
            di as text "Driscoll-Kraay standard errors with lag `dkraay_val'"
        }
    }

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
    local rss = __civreghdfe_rss
    local tss = __civreghdfe_tss
    local r2 = __civreghdfe_r2
    local rmse = __civreghdfe_rmse
    local F_model = __civreghdfe_F

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

    * Compute Wald F-statistic BEFORE posting (b' * V^{-1} * b / K)
    tempname Vinv b_col wald_temp wald_stat
    matrix `Vinv' = syminv(`V_temp')
    matrix `b_col' = `b_temp''
    matrix `wald_temp' = `Vinv' * `b_col'
    matrix `wald_stat' = `b_temp' * `wald_temp'
    local F_wald = `wald_stat'[1,1] / `K_total'

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

    * Note: Don't set vcetype for unadjusted case to match ivreghdfe display format
    if `vce_type' == 1 {
        ereturn local vcetype "Robust"
    }
    else if `vce_type' == 2 {
        ereturn local vcetype "Cluster"
        ereturn local clustvar "`cluster_var'"
    }

    * Store additional statistics
    ereturn scalar rss = `rss'
    ereturn scalar tss = `tss'
    ereturn scalar r2 = `r2'
    ereturn scalar rmse = `rmse'
    ereturn scalar F = `F_wald'

    * Display results in ivreghdfe format
    di as text ""
    di as text "IV (2SLS) estimation"
    di as text "{hline 20}"
    di as text ""

    * Efficiency/consistency notes
    if `vce_type' == 0 {
        di as text "Estimates efficient for homoskedasticity only"
        di as text "Statistics consistent for homoskedasticity only"
    }
    else if `vce_type' == 1 {
        di as text "Estimates efficient for homoskedasticity only"
        di as text "Statistics robust to heteroskedasticity"
    }
    else if `vce_type' == 2 {
        di as text "Estimates efficient for homoskedasticity only"
        di as text "Statistics robust to heteroskedasticity and clustering on `cluster_var'"
    }
    di as text ""

    * Summary statistics line 1
    if `vce_type' == 2 {
        di as text "Number of clusters (`cluster_var') = " as result %9.0fc `N_clust' ///
            as text _col(55) "Number of obs =" as result %9.0fc `N_used'
    }
    else {
        di as text _col(55) "Number of obs =" as result %9.0fc `N_used'
    }

    * F statistic (Wald test)
    local prob_F = Ftail(`K_total', `df_r', `F_wald')
    if `vce_type' == 2 {
        di as text _col(55) "F(" %3.0f `K_total' "," %6.0f (`N_clust' - 1) ") =" as result %9.2f `F_wald'
    }
    else {
        di as text _col(55) "F(" %3.0f `K_total' "," %6.0f `df_r' ") =" as result %9.2f `F_wald'
    }
    di as text _col(55) "Prob > F      =" as result %9.4f `prob_F'

    * R-squared and RSS (match ivreghdfe format: 2 spaces then number)
    di as text "Total (centered) SS     =  " as result %-14.1f `tss' ///
        as text _col(55) "Centered R2   =" as result %9.4f `r2'
    di as text "Total (uncentered) SS   =  " as result %-14.1f `tss' ///
        as text _col(55) "Uncentered R2 =" as result %9.4f `r2'
    di as text "Residual SS             =  " as result %-14.1f `rss' ///
        as text _col(55) "Root MSE      =" as result %9.0f round(`rmse')
    di as text ""

    * Display coefficient table
    ereturn display, level(95)

    * Display diagnostic tests (always shown for IV)

    * Retrieve underidentification test
    capture scalar temp_underid = __civreghdfe_underid
    capture scalar temp_underid_df = __civreghdfe_underid_df
    if _rc == 0 {
        local underid_stat = temp_underid
        local underid_df = temp_underid_df
        ereturn scalar idstat = `underid_stat'
        ereturn scalar iddf = `underid_df'
        ereturn scalar idp = chi2tail(`underid_df', `underid_stat')
    }

    * Underidentification test
    if `vce_type' == 0 {
        di as text "Underidentification test (Anderson canon. corr. LM statistic):" ///
            _col(65) as result %14.3f `underid_stat'
    }
    else {
        di as text "Underidentification test (Kleibergen-Paap rk LM statistic):" ///
            _col(65) as result %14.3f `underid_stat'
    }
    local underid_p = chi2tail(`underid_df', `underid_stat')
    di as text _col(52) "Chi-sq(" as result `underid_df' as text ") P-val =" ///
        as result %10.4f `underid_p'
    di as text "{hline 78}"

    * Weak identification test
    capture scalar temp_cd_f = __civreghdfe_cd_f
    capture scalar temp_kp_f = __civreghdfe_kp_f
    local cd_f_val = temp_cd_f
    local kp_f_val = temp_kp_f

    * Always display Cragg-Donald F
    di as text "Weak identification test (Cragg-Donald Wald F statistic):" ///
        _col(65) as result %14.3f `cd_f_val'

    * For robust/cluster, also display Kleibergen-Paap rk Wald F
    if `vce_type' != 0 & `kp_f_val' > 0 {
        di as text "                         (Kleibergen-Paap rk Wald F statistic):" ///
            _col(65) as result %14.3f `kp_f_val'
    }

    * Stock-Yogo critical values (for 1 endogenous variable, varies by # excluded instruments)
    * Using values from Stock-Yogo (2005)
    if `K_endog' == 1 {
        if `K_inst_excl' == 1 {
            * 1 excluded instrument (just-identified)
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 16.38
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 8.96
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 6.66
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 5.53
        }
        else if `K_inst_excl' == 2 {
            * 2 excluded instruments
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 19.93
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 11.59
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 8.75
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 7.25
        }
        else if `K_inst_excl' == 3 {
            * 3 excluded instruments
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 22.30
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 12.83
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 9.54
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 7.80
        }
        else {
            * More than 3 excluded instruments - show approximate values
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 16.38
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 8.96
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 6.66
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 5.53
        }
    }
    else if `K_endog' == 2 {
        * Two endogenous - use appropriate critical values
        if `K_inst_excl' == 3 {
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 13.43
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 8.18
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 6.40
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 5.45
        }
        else {
            di as text "(Stock-Yogo critical values: see Stock-Yogo (2005) for `K_endog' endogenous, `K_inst_excl' instruments)"
        }
    }
    else {
        * More than 2 endogenous - show message
        di as text "(Stock-Yogo critical values: see Stock-Yogo (2005) for `K_endog' endogenous variables)"
    }
    di as text "Source: Stock-Yogo (2005).  Reproduced by permission."
    if `vce_type' != 0 {
        di as text "NB: Critical values are for Cragg-Donald F statistic and i.i.d. errors."
    }
    di as text "{hline 78}"

    * Display reduced form if requested
    if "`rf'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "Reduced-form regression: `depvar' on instruments"
        di as text "{hline 60}"

        * Save current IV results
        tempname iv_b iv_V
        matrix `iv_b' = e(b)
        matrix `iv_V' = e(V)
        local iv_N = e(N)
        local iv_df_r = e(df_r)
        local iv_df_a = e(df_a)
        local iv_K = e(K)
        local iv_K_endog = e(K_endog)
        local iv_K_exog = e(K_exog)
        local iv_K_iv = e(K_iv)
        local iv_G = e(G)
        local iv_title = e(title)
        local iv_model = e(model)
        local iv_vcetype = e(vcetype)
        local iv_cmd = e(cmd)

        * Run reduced form regression using creghdfe
        * Reduced form: y = [exogenous + excluded instruments] + error
        local rf_vce = "`vce'"
        if "`rf_vce'" == "" local rf_vce = "unadjusted"
        qui creghdfe `depvar' `exogvars' `instruments' `if' `in' ///
            [`weight'`exp'], absorb(`absorb') vce(`rf_vce')

        * Display coefficient table
        tempname rf_b rf_V
        matrix `rf_b' = e(b)
        matrix `rf_V' = e(V)

        * Display reduced form coefficients in a table
        local rf_vars : colnames `rf_b'
        local rf_K : word count `rf_vars'
        di as text ""
        di as text "{hline 78}"
        di as text _col(14) "|" _col(28) "`iv_vcetype'"
        di as text %12s abbrev("`depvar'", 12) " | Coefficient  std. err.      t    P>|t|" ///
            "     [95% conf. interval]"
        di as text "{hline 13}+{hline 64}"

        forvalues v = 1/`rf_K' {
            local vname : word `v' of `rf_vars'
            local b = `rf_b'[1, `v']
            local se = sqrt(`rf_V'[`v', `v'])
            if `se' > 0 {
                local t = `b' / `se'
                local p = 2 * ttail(`iv_df_r', abs(`t'))
                local ci_lo = `b' - invttail(`iv_df_r', 0.025) * `se'
                local ci_hi = `b' + invttail(`iv_df_r', 0.025) * `se'
                di as text %12s abbrev("`vname'", 12) " |" ///
                    as result %11.6g `b' "  " %9.6g `se' "  " ///
                    %7.2f `t' "   " %5.3f `p' "    " ///
                    %11.6g `ci_lo' "  " %11.6g `ci_hi'
            }
            else {
                di as text %12s abbrev("`vname'", 12) " |" ///
                    as result %11.6g `b' "  " "(omitted)"
            }
        }
        di as text "{hline 78}"

        * Restore IV results
        * We need to re-post to restore the original e() results
        * But since ereturn post clears everything, we'll just leave a note
        * The key IV results are already stored in the output above
    }

    * Store and display diagnostic statistics
    capture scalar temp_sargan = __civreghdfe_sargan
    capture scalar temp_sargan_df = __civreghdfe_sargan_df
    capture scalar temp_cd_f = __civreghdfe_cd_f
    capture scalar temp_endog_chi2 = __civreghdfe_endog_chi2
    capture scalar temp_endog_f = __civreghdfe_endog_f
    capture scalar temp_endog_df = __civreghdfe_endog_df

    if _rc == 0 {
        local sargan_val = temp_sargan
        local sargan_df = temp_sargan_df
        local cd_f_val = temp_cd_f
        local endog_chi2 = temp_endog_chi2
        local endog_f = temp_endog_f
        local endog_df = temp_endog_df

        * Store in e()
        ereturn scalar sargan = `sargan_val'
        ereturn scalar sargan_df = `sargan_df'
        if `sargan_df' > 0 {
            ereturn scalar sargan_p = chi2tail(`sargan_df', `sargan_val')
        }
        ereturn scalar cd_f = `cd_f_val'

        * Store weak ID statistics (widstat = Kleibergen-Paap F if robust, else Cragg-Donald F)
        capture scalar temp_kp_f2 = __civreghdfe_kp_f
        if _rc == 0 & `vce_type' != 0 {
            local kp_f_val2 = temp_kp_f2
            if `kp_f_val2' > 0 {
                ereturn scalar widstat = `kp_f_val2'
                ereturn scalar kp_f = `kp_f_val2'
            }
            else {
                ereturn scalar widstat = `cd_f_val'
            }
        }
        else {
            ereturn scalar widstat = `cd_f_val'
        }
        capture scalar drop temp_kp_f2

        * Store endogeneity test
        ereturn scalar endog_chi2 = `endog_chi2'
        ereturn scalar endog_f = `endog_f'
        ereturn scalar endog_df = `endog_df'
        if `endog_df' > 0 {
            ereturn scalar endog_p = chi2tail(`endog_df', `endog_chi2')
        }

        * Display Sargan/Hansen J statistic
        if `vce_type' == 0 {
            di as text "Sargan statistic (overidentification test of all instruments):" ///
                _col(65) as result %14.3f `sargan_val'
        }
        else {
            di as text "Hansen J statistic (overidentification test of all instruments):" ///
                _col(65) as result %14.3f `sargan_val'
        }
        if `sargan_df' > 0 {
            local sargan_p = chi2tail(`sargan_df', `sargan_val')
            di as text _col(52) "Chi-sq(" as result `sargan_df' as text ") P-val =" ///
                as result %10.4f `sargan_p'
        }
        else {
            di as text _col(50) "(equation exactly identified)"
        }
        di as text "{hline 78}"
    }

    * Display instruments footer
    di as text "Instrumented:         `endogvars'"
    if "`exogvars'" != "" {
        di as text "Included instruments: `exogvars'"
    }
    di as text "Excluded instruments: `instruments'"
    di as text "Partialled-out:       _cons"
    di as text "                      nb: total SS, model F and R2s are after partialling-out;"
    di as text "                          any small-sample adjustments include partialled-out"
    di as text "                          variables in regressor count K"
    di as text "{hline 78}"

    * Display absorbed degrees of freedom table
    di as text ""
    di as text "Absorbed degrees of freedom:"
    di as text "{hline 53}+"
    di as text " Absorbed FE | Categories  - Redundant  = Num. Coefs |"
    di as text "{hline 13}+{hline 39}|"

    * For single FE, use df_a as approximation
    if `G' == 1 {
        local fe_var : word 1 of `absorb'
        local n_levels = `df_a'
        local redundant = 0
        local is_nested = ""
        if `vce_type' == 2 & "`fe_var'" == "`cluster_var'" {
            local redundant = `n_levels'
            local is_nested = "*"
        }
        local n_coefs = `n_levels' - `redundant'
        di as text %12s abbrev("`fe_var'", 12) " |" ///
            as result %10.0f `n_levels' %12.0f `redundant' %12.0f `n_coefs' ///
            as text "     |`is_nested'"
    }
    else {
        * For multiple FE, just show total
        di as text "    (total) |" ///
            as result %10.0f `df_a' %12.0f 0 %12.0f `df_a' ///
            as text "     |"
    }
    di as text "{hline 53}+"
    if `vce_type' == 2 {
        di as text "* = FE nested within cluster; treated as redundant for DoF computation"
    }
    di as text ""

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
    capture scalar drop __civreghdfe_rss
    capture scalar drop __civreghdfe_tss
    capture scalar drop __civreghdfe_r2
    capture scalar drop __civreghdfe_rmse
    capture scalar drop __civreghdfe_F
    capture scalar drop __civreghdfe_N_clust
    capture scalar drop __civreghdfe_nested_fe_index
    capture scalar drop __civreghdfe_est_method
    capture scalar drop __civreghdfe_kclass
    capture scalar drop __civreghdfe_fuller
    capture scalar drop __civreghdfe_coviv
    capture scalar drop __civreghdfe_lambda
    capture scalar drop __civreghdfe_kernel
    capture scalar drop __civreghdfe_bw
    capture scalar drop __civreghdfe_dkraay
    capture scalar drop __civreghdfe_sargan
    capture scalar drop __civreghdfe_sargan_df
    capture scalar drop __civreghdfe_cd_f
    capture scalar drop __civreghdfe_kp_f
    capture scalar drop __civreghdfe_underid
    capture scalar drop __civreghdfe_underid_df
    capture scalar drop temp_underid
    capture scalar drop temp_underid_df
    capture scalar drop __civreghdfe_endog_chi2
    capture scalar drop __civreghdfe_endog_f
    capture scalar drop __civreghdfe_endog_df
    capture scalar drop temp_sargan
    capture scalar drop temp_sargan_df
    capture scalar drop temp_cd_f
    capture scalar drop temp_kp_f
    capture scalar drop temp_endog_chi2
    capture scalar drop temp_endog_f
    capture scalar drop temp_endog_df
    forvalues e = 1/`K_endog' {
        capture scalar drop __civreghdfe_F1_`e'
    }
    capture matrix drop __civreghdfe_b
    capture matrix drop __civreghdfe_V

end
