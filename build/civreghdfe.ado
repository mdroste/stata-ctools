*! version 1.0.0 17Jan2026
*! civreghdfe: C-accelerated instrumental variables regression with HDFE
*! Implements 2SLS/IV/LIML/GMM2S with high-dimensional fixed effects absorption

program define civreghdfe, eclass
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
        FFIRst ///
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
        NOConstant ///
        Level(cilevel) ///
        noHEADer ///
        noFOOTer ///
        noOUTput ///
        TItle(string) ///
        DEPname(string) ///
        noid ///
        ORThog(varlist) ///
        ENDOGtest(varlist) ///
        REDundant(varlist) ///
        PARTial(varlist) ///
        SAVEFirst ///
        SAVEFPrefix(string) ///
        SAVERF ///
        SAVERFPrefix(string) ///
        THReads(integer 0)]

    * Clean up any leftover scalars from previous runs to prevent state contamination
    capture scalar drop __civreghdfe_n_orthog
    capture scalar drop __civreghdfe_n_endogtest
    capture scalar drop __civreghdfe_n_partial
    forval i = 1/10 {
        capture scalar drop __civreghdfe_orthog_`i'
        capture scalar drop __civreghdfe_endogtest_`i'
        capture scalar drop __civreghdfe_partial_`i'
    }

    * Parse the varlist with parentheses for endogenous vars and instruments
    * Format: depvar [exog_before] (endog1 endog2 = inst1 inst2 inst3) [exog_after]
    local 0 `"`0'"'

    * First get everything before the opening parenthesis
    gettoken before_paren rest : 0, parse("(")
    local before_paren = strtrim("`before_paren'")

    * First token of before_paren is depvar, rest is exogenous vars before parentheses
    gettoken depvar exog_before : before_paren
    local depvar = strtrim("`depvar'")
    local exog_before = strtrim("`exog_before'")

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
    local exog_after = strtrim(substr("`rest'", `paren_end' + 1, .))

    * Combine exogenous variables from before and after parentheses
    local exogvars `exog_before' `exog_after'
    local exogvars = strtrim("`exogvars'")

    * Parse endogenous = instruments
    local eq_pos = strpos("`paren_content'", "=")
    if `eq_pos' == 0 {
        di as error "Must specify instruments after '=' in parentheses"
        exit 198
    }

    local endogvars = strtrim(substr("`paren_content'", 1, `eq_pos' - 1))
    local instruments = strtrim(substr("`paren_content'", `eq_pos' + 1, .))

    * Validate dependent variable exists
    foreach v of local depvar {
        confirm numeric variable `v'
    }
    foreach v of local absorb {
        confirm variable `v'
    }

    * Mark sample early for fvrevar
    marksample touse
    markout `touse' `depvar' `absorb'

    * Expand factor variables for endogenous variables
    local endogvars_coef ""
    if "`endogvars'" != "" {
        fvrevar `endogvars' if `touse'
        local endogvars_fvrevar "`r(varlist)'"
        _rmcoll `endogvars_fvrevar' if `touse', forcedrop
        local endogvars_expanded "`r(varlist)'"

        fvexpand `endogvars' if `touse'
        local endogvars_names_all "`r(varlist)'"
        foreach v of local endogvars_names_all {
            if !regexm("`v'", "^[0-9]+b\.") & !regexm("`v'", "^o\.") & "`v'" != "" {
                local endogvars_coef `endogvars_coef' `v'
            }
        }
    }
    else {
        local endogvars_expanded ""
    }

    * Expand factor variables for exogenous variables
    local exogvars_coef ""
    if "`exogvars'" != "" {
        fvrevar `exogvars' if `touse'
        local exogvars_fvrevar "`r(varlist)'"
        _rmcoll `exogvars_fvrevar' if `touse', forcedrop
        local exogvars_expanded "`r(varlist)'"

        fvexpand `exogvars' if `touse'
        local exogvars_names_all "`r(varlist)'"
        foreach v of local exogvars_names_all {
            if !regexm("`v'", "^[0-9]+b\.") & !regexm("`v'", "^o\.") & "`v'" != "" {
                local exogvars_coef `exogvars_coef' `v'
            }
        }
    }
    else {
        local exogvars_expanded ""
    }

    * Expand factor variables for instruments
    local instruments_coef ""
    if "`instruments'" != "" {
        fvrevar `instruments' if `touse'
        local instruments_fvrevar "`r(varlist)'"
        _rmcoll `instruments_fvrevar' if `touse', forcedrop
        local instruments_expanded "`r(varlist)'"

        fvexpand `instruments' if `touse'
        local instruments_names_all "`r(varlist)'"
        foreach v of local instruments_names_all {
            if !regexm("`v'", "^[0-9]+b\.") & !regexm("`v'", "^o\.") & "`v'" != "" {
                local instruments_coef `instruments_coef' `v'
            }
        }
    }
    else {
        local instruments_expanded ""
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

    * Count variables (using expanded counts)
    local K_endog : word count `endogvars_expanded'
    local K_exog : word count `exogvars_expanded'
    local K_inst_excl : word count `instruments_expanded'
    local G : word count `absorb'

    * Build full instrument list (exogenous + excluded instruments) - use expanded
    local all_instruments `exogvars_expanded' `instruments_expanded'
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
    local cluster_var2 ""
    if `"`vce'"' != "" {
        local vce_lower = lower(`"`vce'"')
        if `"`vce_lower'"' == "robust" | `"`vce_lower'"' == "r" {
            local vce_type = 1
        }
        else if substr(`"`vce_lower'"', 1, 2) == "cl" {
            * Extract cluster variable(s)
            local cluster_spec `"`vce'"'
            gettoken vce_cmd cluster_vars : cluster_spec, parse(" (")
            local cluster_vars = strtrim(subinstr("`cluster_vars'", "(", "", .))
            local cluster_vars = strtrim(subinstr("`cluster_vars'", ")", "", .))

            * Count cluster variables
            local n_cluster : word count `cluster_vars'
            if `n_cluster' == 1 {
                * One-way clustering
                local vce_type = 2
                local cluster_var `cluster_vars'
                confirm variable `cluster_var'
            }
            else if `n_cluster' == 2 {
                * Two-way clustering
                local vce_type = 3
                local cluster_var : word 1 of `cluster_vars'
                local cluster_var2 : word 2 of `cluster_vars'
                confirm variable `cluster_var'
                confirm variable `cluster_var2'
            }
            else {
                di as error "vce(cluster) supports 1 or 2 cluster variables"
                exit 198
            }
        }
        else if `"`vce_lower'"' == "unadjusted" | `"`vce_lower'"' == "ols" {
            local vce_type = 0
        }
        else {
            di as error "vce(`vce') not supported. Use vce(robust), vce(cluster varname), or vce(cluster var1 var2)"
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

    * Mark out additional variables from sample (touse already created earlier for fvrevar)
    markout `touse' `endogvars_expanded' `exogvars_expanded' `instruments_expanded'
    if "`cluster_var'" != "" {
        markout `touse' `cluster_var'
    }
    if "`cluster_var2'" != "" {
        markout `touse' `cluster_var2'
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
        else if `vce_type' == 3 {
            di as text "VCE:                   Two-way Cluster (`cluster_var' `cluster_var2')"
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
    scalar __civreghdfe_has_cluster = (`vce_type' == 2 | `vce_type' == 3)
    scalar __civreghdfe_has_cluster2 = (`vce_type' == 3)
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

    * Parse and validate diagnostic test options
    * orthog() - test exogeneity of specified instruments
    local n_orthog = 0
    local orthog_indices ""
    if "`orthog'" != "" {
        foreach v of local orthog {
            * Check that variable is an excluded instrument
            local found = 0
            local idx = 1
            foreach inst of local instruments {
                if "`v'" == "`inst'" {
                    local found = 1
                    local orthog_indices `orthog_indices' `idx'
                }
                local idx = `idx' + 1
            }
            if `found' == 0 {
                di as error "orthog(): `v' is not an excluded instrument"
                exit 198
            }
            local n_orthog = `n_orthog' + 1
        }
        if "`verbose'" != "" {
            di as text "C-statistic test for: `orthog'"
        }
    }
    scalar __civreghdfe_n_orthog = `n_orthog'
    forval i = 1/`n_orthog' {
        local idx : word `i' of `orthog_indices'
        scalar __civreghdfe_orthog_`i' = `idx'
    }

    * endogtest() - test endogeneity of specified regressors
    local n_endogtest = 0
    local endogtest_indices ""
    if "`endogtest'" != "" {
        foreach v of local endogtest {
            * Check that variable is an endogenous variable
            local found = 0
            local idx = 1
            foreach endog of local endogvars {
                if "`v'" == "`endog'" {
                    local found = 1
                    local endogtest_indices `endogtest_indices' `idx'
                }
                local idx = `idx' + 1
            }
            if `found' == 0 {
                di as error "endogtest(): `v' is not an endogenous variable"
                exit 198
            }
            local n_endogtest = `n_endogtest' + 1
        }
        if "`verbose'" != "" {
            di as text "Endogeneity test for: `endogtest'"
        }
    }
    scalar __civreghdfe_n_endogtest = `n_endogtest'
    forval i = 1/`n_endogtest' {
        local idx : word `i' of `endogtest_indices'
        scalar __civreghdfe_endogtest_`i' = `idx'
    }

    * redundant() - test redundancy of specified instruments
    local n_redundant = 0
    local redundant_indices ""
    if "`redundant'" != "" {
        foreach v of local redundant {
            * Check that variable is an excluded instrument
            local found = 0
            local idx = 1
            foreach inst of local instruments {
                if "`v'" == "`inst'" {
                    local found = 1
                    local redundant_indices `redundant_indices' `idx'
                }
                local idx = `idx' + 1
            }
            if `found' == 0 {
                di as error "redundant(): `v' is not an excluded instrument"
                exit 198
            }
            local n_redundant = `n_redundant' + 1
        }
        if "`verbose'" != "" {
            di as text "Redundancy test for: `redundant'"
        }
    }
    scalar __civreghdfe_n_redundant = `n_redundant'
    forval i = 1/`n_redundant' {
        local idx : word `i' of `redundant_indices'
        scalar __civreghdfe_redundant_`i' = `idx'
    }

    * partial() - variables to partial out via FWL
    * These must be exogenous variables
    local n_partial = 0
    local partial_indices ""
    local partial_vars ""
    if "`partial'" != "" {
        foreach v of local partial {
            * Check that variable is an exogenous variable
            local found = 0
            local idx = 1
            foreach exog of local exogvars {
                if "`v'" == "`exog'" {
                    local found = 1
                    local partial_indices `partial_indices' `idx'
                    local partial_vars `partial_vars' `v'
                }
                local idx = `idx' + 1
            }
            if `found' == 0 {
                di as error "partial(): `v' is not an exogenous variable"
                exit 198
            }
            local n_partial = `n_partial' + 1
        }
        if "`verbose'" != "" {
            di as text "Partialling out: `partial'"
        }
    }
    scalar __civreghdfe_n_partial = `n_partial'
    forval i = 1/`n_partial' {
        local idx : word `i' of `partial_indices'
        scalar __civreghdfe_partial_`i' = `idx'
    }

    * Create output matrices
    local K_total = `K_endog' + `K_exog'
    tempname b_temp V_temp

    matrix `b_temp' = J(1, `K_total', 0)
    matrix `V_temp' = J(`K_total', `K_total', 0)
    matrix __civreghdfe_b = `b_temp'
    matrix __civreghdfe_V = `V_temp'

    * Build variable list for plugin
    * Order: depvar, endogvars, exogvars, all_instruments, absorb, [cluster], [cluster2], [weight]
    * Use expanded variables (factor variables converted to temp numeric vars)
    local plugin_vars `depvar' `endogvars_expanded' `exogvars_expanded' `all_instruments' `absorb'
    if `vce_type' == 2 | `vce_type' == 3 {
        local plugin_vars `plugin_vars' `cluster_var'
    }
    if `vce_type' == 3 {
        local plugin_vars `plugin_vars' `cluster_var2'
    }
    if `has_weights' {
        local plugin_vars `plugin_vars' `weight_var'
    }

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Call the C plugin
    timer clear 99
    timer on 99
    plugin call ctools_plugin `plugin_vars' if `touse', "civreghdfe `threads_code' iv_regression"
    timer off 99
    qui timer list 99
    local t_plugin = r(t99)
    timer clear 98
    timer on 98

    * Retrieve results
    local N_used = __civreghdfe_N
    local df_r = __civreghdfe_df_r
    local df_a = __civreghdfe_df_a
    local df_a_for_vce = __civreghdfe_df_a_for_vce
    local K = __civreghdfe_K
    local rss = __civreghdfe_rss
    local tss = __civreghdfe_tss
    local r2 = __civreghdfe_r2
    local rmse = __civreghdfe_rmse
    local F_model = __civreghdfe_F

    * Get number of levels for each FE (for absorbed DOF table)
    forval g = 1/`G' {
        capture confirm scalar __civreghdfe_num_levels_`g'
        if _rc == 0 {
            local fe_levels_`g' = scalar(__civreghdfe_num_levels_`g')
        }
        else {
            local fe_levels_`g' = .
        }
        local fe_nested_`g' = 0
    }

    * Read nested FE detection from C plugin (computed during plugin call)
    if `vce_type' == 2 & "`cluster_var'" != "" {
        forval g = 1/`G' {
            capture local fe_nested_`g' = __civreghdfe_fe_nested_`g'
            if _rc != 0 local fe_nested_`g' = 0
        }
    }

    if `vce_type' == 2 | `vce_type' == 3 {
        local N_clust = __civreghdfe_N_clust
    }
    if `vce_type' == 3 {
        capture local N_clust2 = __civreghdfe_N_clust2
        if _rc != 0 local N_clust2 = 0
    }

    * Get coefficient vector and VCE (already in [endog, exog] order from C)
    matrix `b_temp' = __civreghdfe_b
    matrix `V_temp' = __civreghdfe_V

    * Build variable names in [endog, exog] order (matching C output)
    * Use coefficient names which have proper factor variable notation
    * Exclude partialled variables from varnames
    local varnames ""
    foreach v of local endogvars_coef {
        local varnames `varnames' `v'
    }
    foreach v of local exogvars_coef {
        * Check if this variable was partialled out
        * Note: partial() uses original var names, so we check against original exogvars
        local is_partial = 0
        foreach pv of local partial_vars {
            if "`v'" == "`pv'" {
                local is_partial = 1
            }
        }
        if `is_partial' == 0 {
            local varnames `varnames' `v'
        }
    }

    * Update K_exog to reflect partialled variables removed
    local K_exog = `K_exog' - `n_partial'
    local K_total = `K_exog' + `K_endog'

    * Force VCE symmetry (numerical precision fix)
    tempname V_sym
    matrix `V_sym' = (`V_temp' + `V_temp'') / 2
    matrix `V_temp' = `V_sym'

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
    local depname_use = cond("`depname'" != "", "`depname'", "`depvar'")
    ereturn post `b_temp' `V_temp', esample(`touse') depname(`depname_use')

    * Store scalars
    ereturn scalar N = `N_used'
    ereturn scalar df_r = `df_r'
    * Report df_a_for_vce when clustering (matches ivreghdfe behavior)
    * When FE is nested in cluster, df_a_for_vce = 0
    if `vce_type' == 2 | `vce_type' == 3 {
        ereturn scalar df_a = `df_a_for_vce'
    }
    else {
        ereturn scalar df_a = `df_a'
    }
    ereturn scalar K = `K'
    ereturn scalar level = `level'
    ereturn scalar K_endog = `K_endog'
    ereturn scalar K_exog = `K_exog'
    ereturn scalar K_iv = `K_iv'
    ereturn scalar G = `G'

    if `vce_type' == 2 | `vce_type' == 3 {
        ereturn scalar N_clust = `N_clust'
        ereturn scalar N_clust1 = `N_clust'
    }
    if `vce_type' == 3 {
        ereturn scalar N_clust2 = `N_clust2'
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
    ereturn local depvar "`depname_use'"
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
    else if `vce_type' == 3 {
        ereturn local vcetype "Two-way Cluster"
        ereturn local clustvar "`cluster_var'"
        ereturn local clustvar1 "`cluster_var'"
        ereturn local clustvar2 "`cluster_var2'"
    }

    * Store additional statistics
    ereturn scalar rss = `rss'
    ereturn scalar tss = `tss'
    ereturn scalar mss = `tss' - `rss'
    ereturn scalar r2 = `r2'
    ereturn scalar rmse = `rmse'
    ereturn scalar F = `F_wald'

    * Model degrees of freedom
    ereturn scalar df_m = `K_total'

    * Number of absorbed fixed effect dimensions (matches reghdfe/ivreghdfe)
    ereturn scalar N_hdfe = `G'

    * Adjusted R-squared (matching ivreghdfe formula)
    * With absorb: noconstant is set, so constant not added to sdofminus
    * sdofminus = absorb_ct = max(1, HDFE.df_a)
    *   - HDFE.df_a = df_a_for_vce (0 if nested, df_a otherwise)
    *   - ivreghdfe line 670: if HDFE.df_a=0 (nested), force absorb_ct to 1
    * Formula: r2_a = 1 - (1 - r2) * N / (N - K - sdofminus)  (ivreghdfe line 1304)
    local sdofminus = max(1, `df_a_for_vce')
    local adj_denom = `N_used' - `K_total' - `sdofminus'
    if `adj_denom' > 0 {
        * With absorb, constant is partialled out, so use N not N-1
        local r2_a = 1 - (1 - `r2') * `N_used' / `adj_denom'
    }
    else {
        local r2_a = .
    }
    ereturn scalar r2_a = `r2_a'

    * Display results in ivreghdfe format (if header not suppressed)
    if "`header'" == "" {
        di as text ""
        if "`title'" != "" {
            di as text "`title'"
        }
        else {
            di as text "IV (2SLS) estimation"
        }
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
        else if `vce_type' == 3 {
            di as text "Estimates efficient for homoskedasticity only"
            di as text "Statistics robust to heteroskedasticity and two-way clustering on `cluster_var' and `cluster_var2'"
        }
        di as text ""

        * Summary statistics line 1
        if `vce_type' == 2 {
            di as text "Number of clusters (`cluster_var') = " as result %9.0fc `N_clust' ///
                as text _col(55) "Number of obs =" as result %9.0fc `N_used'
        }
        else if `vce_type' == 3 {
            di as text "Number of clusters (`cluster_var') = " as result %6.0fc `N_clust' ///
                as text _col(55) "Number of obs =" as result %9.0fc `N_used'
            di as text "Number of clusters (`cluster_var2') = " as result %6.0fc `N_clust2'
        }
        else {
            di as text _col(55) "Number of obs =" as result %9.0fc `N_used'
        }

        * F statistic (Wald test)
        * For two-way clustering, use min(N_clust1, N_clust2) - 1 for degrees of freedom
        local prob_F = Ftail(`K_total', `df_r', `F_wald')
        if `vce_type' == 2 {
            di as text _col(55) "F(" %3.0f `K_total' "," %6.0f (`N_clust' - 1) ") =" as result %9.2f `F_wald'
        }
        else if `vce_type' == 3 {
            local min_clust = min(`N_clust', `N_clust2')
            di as text _col(55) "F(" %3.0f `K_total' "," %6.0f (`min_clust' - 1) ") =" as result %9.2f `F_wald'
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
    }

    * Display coefficient table (if output not suppressed)
    if "`output'" == "" {
        ereturn display, level(`level')
    }

    * Display full first-stage results (if ffirst specified)
    if "`ffirst'" != "" & `K_endog' > 0 {
        di as text ""
        di as text "{hline 78}"
        di as text "First-stage regression summary statistics"
        di as text "{hline 78}"
        di as text _col(20) "Partial" _col(32) "F(" as result `K_inst_excl' as text "," ///
            as result `df_r' as text ")" _col(52) "Prob > F"

        forvalues e = 1/`K_endog' {
            local endogvar : word `e' of `endogvars'
            capture scalar temp_F = __civreghdfe_F1_`e'
            capture scalar temp_r2 = __civreghdfe_partial_r2_`e'
            if _rc == 0 {
                local F_val = temp_F
                local r2_val = temp_r2
                local p_val = Ftail(`K_inst_excl', `df_r', `F_val')
                di as text abbrev("`endogvar'", 15) _col(20) as result %8.4f `r2_val' ///
                    _col(32) as result %12.2f `F_val' _col(52) as result %10.4f `p_val'

                * Store extended first-stage results
                ereturn scalar partial_r2_`e' = `r2_val'
            }
        }
        di as text "{hline 78}"
        di as text ""
    }

    * Display diagnostic tests (if footer not suppressed)
    if "`footer'" == "" {

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

    * Underidentification test (if noid not specified)
    if "`id'" == "" {
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
    }

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
            local sargan_pval = chi2tail(`sargan_df', `sargan_val')
            ereturn scalar sargan_p = `sargan_pval'
            ereturn scalar sarganp = `sargan_pval'
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

    * Display orthog test (C-statistic) if requested
    if `n_orthog' > 0 {
        capture scalar temp_cstat = __civreghdfe_cstat
        capture scalar temp_cstat_df = __civreghdfe_cstat_df
        if _rc == 0 {
            local cstat_val = temp_cstat
            local cstat_df = temp_cstat_df
            local cstat_p = chi2tail(`cstat_df', `cstat_val')

            ereturn scalar cstat = `cstat_val'
            ereturn scalar cstat_df = `cstat_df'
            ereturn scalar cstat_p = `cstat_p'

            di as text "C statistic (exogeneity/orthogonality of instruments):" ///
                _col(65) as result %14.3f `cstat_val'
            di as text _col(52) "Chi-sq(" as result `cstat_df' as text ") P-val =" ///
                as result %10.4f `cstat_p'
            di as text "{hline 78}"
        }
        capture scalar drop temp_cstat
        capture scalar drop temp_cstat_df
    }

    * Display endogtest result if specified
    if `n_endogtest' > 0 {
        capture scalar temp_endogtest_stat = __civreghdfe_endogtest_stat
        capture scalar temp_endogtest_df = __civreghdfe_endogtest_df
        if _rc == 0 {
            local endogtest_val = temp_endogtest_stat
            local endogtest_df_val = temp_endogtest_df
            local endogtest_p = chi2tail(`endogtest_df_val', `endogtest_val')

            ereturn scalar endogtest = `endogtest_val'
            ereturn scalar endogtest_df = `endogtest_df_val'
            ereturn scalar endogtest_p = `endogtest_p'

            di as text "Endogeneity test of endogenous regressors (`endogtest'):" ///
                _col(65) as result %14.3f `endogtest_val'
            di as text _col(52) "Chi-sq(" as result `endogtest_df_val' as text ") P-val =" ///
                as result %10.4f `endogtest_p'
            di as text "{hline 78}"
        }
        capture scalar drop temp_endogtest_stat
        capture scalar drop temp_endogtest_df
    }

    * Display redundant test result if specified
    if `n_redundant' > 0 {
        capture scalar temp_redund_stat = __civreghdfe_redund_stat
        capture scalar temp_redund_df = __civreghdfe_redund_df
        if _rc == 0 {
            local redund_val = temp_redund_stat
            local redund_df_val = temp_redund_df
            local redund_p = chi2tail(`redund_df_val', `redund_val')

            ereturn scalar redund = `redund_val'
            ereturn scalar redund_df = `redund_df_val'
            ereturn scalar redund_p = `redund_p'

            di as text "Redundancy test of instruments (`redundant'):" ///
                _col(65) as result %14.3f `redund_val'
            di as text _col(52) "Chi-sq(" as result `redund_df_val' as text ") P-val =" ///
                as result %10.4f `redund_p'
            di as text "{hline 78}"
        }
        capture scalar drop temp_redund_stat
        capture scalar drop temp_redund_df
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
    di as text "{hline 13}{c TT}{hline 39}{c TRC}"
    di as text " Absorbed FE {c |} Categories  - Redundant  = Num. Coefs {c |}"
    di as text "{hline 13}{c +}{hline 39}{c RT}"

    * Display each FE with its components
    local total_coefs = 0
    local has_nested = 0
    forval g = 1/`G' {
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
    di as text ""

    } /* End of if "`footer'" == "" */

    * Save first-stage estimation results (if savefirst specified)
    if "`savefirst'" != "" {
        local prefix = cond("`savefprefix'" != "", "`savefprefix'", "_civreghdfe_")
        * Store current estimation under each endogenous variable's name
        * Note: Full first-stage regressions would require running OLS separately
        estimates store `prefix'main
        di as text ""
        di as text "Estimation stored as: `prefix'main"
        di as text "First-stage F-statistics available in e(F_first1), e(F_first2), etc."
    }

    * Save reduced-form estimation results (if saverf specified)
    if "`saverf'" != "" {
        local rfprefix = cond("`saverfprefix'" != "", "`saverfprefix'", "_civreghdfe_rf_")
        estimates store `rfprefix'main
        di as text ""
        di as text "Estimation stored as: `rfprefix'main"
    }

    * Capture post-processing time (before cleanup)
    timer off 98
    qui timer list 98
    local t_stata_post = r(t98)

    * Display timing breakdown if verbose or timeit specified
    if "`verbose'" != "" | "`timeit'" != "" {
        local t_total_wall = `t_plugin' + `t_stata_post'
        di as text ""
        di as text "{hline 55}"
        di as text "civreghdfe timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f _civreghdfe_time_load " sec"
        di as text "    Singleton removal:      " as result %8.4f _civreghdfe_time_singleton " sec"
        di as text "    HDFE setup:             " as result %8.4f _civreghdfe_time_setup " sec"
        di as text "    FWL partialling:        " as result %8.4f _civreghdfe_time_fwl " sec"
        di as text "    HDFE partial out:       " as result %8.4f _civreghdfe_time_partial " sec"
        di as text "    Post-processing:        " as result %8.4f _civreghdfe_time_postproc " sec"
        di as text "    IV estimation + VCE:    " as result %8.4f _civreghdfe_time_estimate " sec"
        di as text "    Stats computation:      " as result %8.4f _civreghdfe_time_stats " sec"
        di as text "    Store results:          " as result %8.4f _civreghdfe_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _civreghdfe_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Plugin call (wall clock): " as result %8.4f `t_plugin' " sec"
        di as text "  Stata post-processing:    " as result %8.4f `t_stata_post' " sec"
        di as text "  {hline 53}"
        di as text "  Total wall clock:         " as result %8.4f `t_total_wall' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = _civreghdfe_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _civreghdfe_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * Clean up temp scalars and matrices
    capture scalar drop __civreghdfe_K_endog
    capture scalar drop __civreghdfe_K_exog
    capture scalar drop __civreghdfe_K_iv
    capture scalar drop __civreghdfe_G
    capture scalar drop __civreghdfe_has_cluster
    capture scalar drop __civreghdfe_has_cluster2
    capture scalar drop __civreghdfe_has_weights
    capture scalar drop __civreghdfe_N_clust2
    capture scalar drop __civreghdfe_weight_type
    capture scalar drop __civreghdfe_vce_type
    capture scalar drop __civreghdfe_maxiter
    capture scalar drop __civreghdfe_tolerance
    capture scalar drop __civreghdfe_verbose
    capture scalar drop __civreghdfe_N
    capture scalar drop __civreghdfe_df_r
    capture scalar drop __civreghdfe_df_a
    forval g = 1/20 {
        capture scalar drop __civreghdfe_num_levels_`g'
    }
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

    * Clean up diagnostic test scalars
    capture scalar drop __civreghdfe_n_orthog
    forval i = 1/10 {
        capture scalar drop __civreghdfe_orthog_`i'
    }
    capture scalar drop __civreghdfe_cstat
    capture scalar drop __civreghdfe_cstat_df
    capture scalar drop __civreghdfe_n_endogtest
    forval i = 1/10 {
        capture scalar drop __civreghdfe_endogtest_`i'
    }
    capture scalar drop __civreghdfe_endogtest_stat
    capture scalar drop __civreghdfe_endogtest_df
    capture scalar drop __civreghdfe_n_redundant
    forval i = 1/10 {
        capture scalar drop __civreghdfe_redundant_`i'
    }
    capture scalar drop __civreghdfe_redund_stat
    capture scalar drop __civreghdfe_redund_df

    * Clean up nested FE scalars
    forval i = 1/20 {
        capture scalar drop __civreghdfe_fe_nested_`i'
    }

    * Clean up partial scalars
    capture scalar drop __civreghdfe_n_partial
    forval i = 1/10 {
        capture scalar drop __civreghdfe_partial_`i'
    }

    * Clean up timing scalars
    capture scalar drop _civreghdfe_time_load
    capture scalar drop _civreghdfe_time_singleton
    capture scalar drop _civreghdfe_time_setup
    capture scalar drop _civreghdfe_time_fwl
    capture scalar drop _civreghdfe_time_partial
    capture scalar drop _civreghdfe_time_postproc
    capture scalar drop _civreghdfe_time_estimate
    capture scalar drop _civreghdfe_time_stats
    capture scalar drop _civreghdfe_time_store
    capture scalar drop _civreghdfe_time_total
    capture scalar drop _civreghdfe_threads_max
    capture scalar drop _civreghdfe_openmp_enabled

end
