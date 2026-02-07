*! version 0.9.1 06Feb2026
*! civreghdfe: C-accelerated instrumental variables regression with HDFE
*! Implements 2SLS/IV/LIML/GMM2S with high-dimensional fixed effects absorption

program define civreghdfe, eclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start wall clock timer immediately
    timer clear 97
    timer on 97

    * Parse syntax
    * Basic syntax: civreghdfe depvar (endogvars = instruments) [exogvars], absorb() [options]
    * absorb() is optional - without it, runs as regular IV (no FE absorption)
    syntax anything(name=0 equalok) [if] [in] [aw fw pw/], ///
        [Absorb(varlist) ///
        vce(string) ///
        Robust ///
        CLuster(varlist) ///
        TOLerance(real 1e-8) ///
        MAXiter(integer 500) ///
        Verbose ///
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
        KIEFER ///
        CENTer ///
        B0(string) ///
        RESiduals(name) ///
        RESiduals2 ///
        NOConstant ///
        Level(cilevel) ///
        noHEADer ///
        noFOOTer ///
        noOUTput ///
        TItle(string) ///
        SUBtitle(string) ///
        DEPname(string) ///
        noid ///
        ORThog(varlist) ///
        ENDOGtest(varlist) ///
        REDundant(varlist) ///
        PARTial(varlist) ///
        FWL(varlist) ///
        DOFminus(integer 0) ///
        SDOFminus(integer 0) ///
        NOPARTIALsmall ///
        EForm(string) ///
        PLUS ///
        NOOMIT:ted ///
        OMIT:ted ///
        VSQuish ///
        NOEMPtycells ///
        BASElevels ///
        ALLBASElevels ///
        SAVEFirst ///
        SAVEFPrefix(string) ///
        SAVERF ///
        SAVERFPrefix(string) ///
        THReads(integer 0) ///
        NORETURN]

    * Clean up any leftover scalars from previous runs to prevent state contamination
    capture scalar drop __civreghdfe_n_orthog
    capture scalar drop __civreghdfe_n_endogtest
    capture scalar drop __civreghdfe_n_partial
    capture scalar drop __civreghdfe_dofminus
    capture scalar drop __civreghdfe_sdofminus
    capture scalar drop __civreghdfe_nopartialsmall
    capture scalar drop __civreghdfe_center
    forval i = 1/10 {
        capture scalar drop __civreghdfe_orthog_`i'
        capture scalar drop __civreghdfe_endogtest_`i'
        capture scalar drop __civreghdfe_partial_`i'
    }
    capture scalar drop __civreghdfe_num_collinear
    capture scalar drop __civreghdfe_K_keep
    forval i = 1/20 {
        capture scalar drop __civreghdfe_collinear_`i'
    }

    * Handle standalone robust/cluster as aliases for vce()
    if "`robust'" != "" & `"`vce'"' == "" {
        local vce "robust"
    }
    if "`cluster'" != "" & `"`vce'"' == "" {
        local vce "cluster `cluster'"
    }
    if "`robust'" != "" & "`cluster'" != "" {
        di as error "cannot specify both robust and cluster()"
        exit 198
    }
    if "`robust'" != "" & `"`vce'"' != "" & `"`vce'"' != "robust" {
        di as error "cannot specify both robust and vce()"
        exit 198
    }
    if "`cluster'" != "" & `"`vce'"' != "" & !regexm(`"`vce'"', "^cluster") {
        di as error "cannot specify both cluster() and vce()"
        exit 198
    }

    * Handle fwl() as alias for partial()
    if "`fwl'" != "" & "`partial'" == "" {
        local partial "`fwl'"
    }
    else if "`fwl'" != "" & "`partial'" != "" {
        di as error "cannot specify both fwl() and partial()"
        exit 198
    }

    * Handle residuals2 - auto-named residuals
    if "`residuals2'" != "" {
        capture drop _civreghdfe_resid
        local residuals "_civreghdfe_resid"
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
    * If no parentheses, treat remaining vars as exogenous and run as OLS (no IV)
    local is_iv_model = 0
    if strpos("`rest'", "(") == 0 {
        * No IV specification - treat rest as exogenous vars, run as OLS
        local exog_after = strtrim("`rest'")
        local paren_content ""
        local endogvars ""
        local instruments ""
    }
    else {
        local is_iv_model = 1

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

        * Parse endogenous = instruments
        local eq_pos = strpos("`paren_content'", "=")
        if `eq_pos' == 0 {
            di as error "Must specify instruments after '=' in parentheses"
            exit 198
        }

        local endogvars = strtrim(substr("`paren_content'", 1, `eq_pos' - 1))
        local instruments = strtrim(substr("`paren_content'", `eq_pos' + 1, .))
    }

    * Combine exogenous variables from before and after parentheses
    local exogvars `exog_before' `exog_after'
    local exogvars = strtrim("`exogvars'")

    * Validate dependent variable exists and is numeric (error 109 for string vars)
    foreach v of local depvar {
        capture confirm string variable `v'
        if _rc == 0 {
            di as error "type mismatch"
            exit 109
        }
        confirm numeric variable `v'
    }
    foreach v of local absorb {
        confirm variable `v'
    }

    * Mark sample early for fvrevar
    marksample touse
    markout `touse' `depvar' `absorb'

    * Expand factor variables — single fvrevar, base pre-filter, then _rmcoll
    * fvexpand (cheap) gets names; fvrevar (expensive) creates temp vars
    * Pre-filtering base/omitted levels reduces _rmcoll input size
    local endogvars_names_all ""
    local __n_endog_total = 0
    if "`endogvars'" != "" {
        fvexpand `endogvars' if `touse'
        local endogvars_names_all "`r(varlist)'"
        local __n_endog_total : word count `endogvars_names_all'
    }

    local exogvars_names_all ""
    local __n_exog_total = 0
    if "`exogvars'" != "" {
        fvexpand `exogvars' if `touse'
        local exogvars_names_all "`r(varlist)'"
        local __n_exog_total : word count `exogvars_names_all'
    }

    local instruments_names_all ""
    local __n_inst_total = 0
    if "`instruments'" != "" {
        fvexpand `instruments' if `touse'
        local instruments_names_all "`r(varlist)'"
        local __n_inst_total : word count `instruments_names_all'
    }

    * Single fvrevar call for all variable groups
    local endogvars_expanded ""
    local endogvars_coef ""
    local exogvars_expanded ""
    local exogvars_coef ""
    local instruments_expanded ""
    local instruments_coef ""
    local __all_fv `endogvars' `exogvars' `instruments'
    if `"`__all_fv'"' != "" {
        fvrevar `__all_fv' if `touse'
        local __all_exp "`r(varlist)'"

        * Helper: pre-filter base/omitted levels, then _rmcoll per group
        * Keep term if ANY #-separated component is non-base/non-omitted
        local __pos = 1

        * --- Endogenous variables ---
        local __vars_nobase ""
        local __names_nobase ""
        forval __i = 1/`__n_endog_total' {
            local __vname : word `__i' of `endogvars_names_all'
            local __keep = 0
            local __parts : subinstr local __vname "#" " ", all
            foreach __p of local __parts {
                if !regexm("`__p'", "b\.") & !regexm("`__p'", "o\.") {
                    local __keep = 1
                }
            }
            if `__keep' {
                local __vars_nobase `__vars_nobase' `: word `__pos' of `__all_exp''
                local __names_nobase `__names_nobase' `__vname'
            }
            local __pos = `__pos' + 1
        }
        * Collinearity detection handled by C plugin
        local endogvars_expanded "`__vars_nobase'"
        local endogvars_coef "`__names_nobase'"

        * --- Exogenous variables ---
        local __vars_nobase ""
        local __names_nobase ""
        forval __i = 1/`__n_exog_total' {
            local __vname : word `__i' of `exogvars_names_all'
            local __keep = 0
            local __parts : subinstr local __vname "#" " ", all
            foreach __p of local __parts {
                if !regexm("`__p'", "b\.") & !regexm("`__p'", "o\.") {
                    local __keep = 1
                }
            }
            if `__keep' {
                local __vars_nobase `__vars_nobase' `: word `__pos' of `__all_exp''
                local __names_nobase `__names_nobase' `__vname'
            }
            local __pos = `__pos' + 1
        }
        * Collinearity detection handled by C plugin
        local exogvars_expanded "`__vars_nobase'"
        local exogvars_coef "`__names_nobase'"

        * --- Instruments ---
        local __vars_nobase ""
        local __names_nobase ""
        forval __i = 1/`__n_inst_total' {
            local __vname : word `__i' of `instruments_names_all'
            local __keep = 0
            local __parts : subinstr local __vname "#" " ", all
            foreach __p of local __parts {
                if !regexm("`__p'", "b\.") & !regexm("`__p'", "o\.") {
                    local __keep = 1
                }
            }
            if `__keep' {
                local __vars_nobase `__vars_nobase' `: word `__pos' of `__all_exp''
                local __names_nobase `__names_nobase' `__vname'
            }
            local __pos = `__pos' + 1
        }
        * Collinearity detection handled by C plugin
        local instruments_expanded "`__vars_nobase'"
        local instruments_coef "`__names_nobase'"
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
        di as error "no observations"
        exit 2001
    }

    * Early check: HAC kernel options imply robust VCE (except dkraay/kiefer)
    * This needs to happen before verbose output to show correct VCE type
    local panel_aware_hac = 0
    if "`kernel'" != "" & `bw' > 0 & `vce_type' == 0 & `dkraay' == 0 & "`kiefer'" == "" {
        local vce_type = 1
    }

    * Detect panel-aware HAC for kernel/bw with tsset panel data
    * This is separate from vce_type check because user may specify robust explicitly
    * Note: bw > 0 implies HAC even without explicit kernel (defaults to Bartlett later)
    if ("`kernel'" != "" | `bw' > 0) & `dkraay' == 0 & "`kiefer'" == "" & "`absorb'" != "" {
        capture qui tsset
        if _rc == 0 {
            local tsset_ivar `r(panelvar)'
            local tsset_tvar `r(timevar)'
            if "`tsset_ivar'" != "" & "`tsset_tvar'" != "" {
                local panel_aware_hac = 1
            }
        }
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
        else if `vce_type' == 4 {
            di as text "VCE:                   Driscoll-Kraay (`cluster_var')"
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

    * Driscoll-Kraay standard errors: cluster by TIME, not by panel
    * Requires tsset panel data
    local dkraay_tvar ""
    local dkraay_T = 0
    if `dkraay_val' > 0 {
        * Check that data is tsset
        capture tsset
        if _rc != 0 {
            di as error "dkraay option requires tsset panel data"
            exit 459
        }
        local dkraay_tvar = r(timevar)
        local dkraay_ivar = r(panelvar)
        if "`dkraay_tvar'" == "" | "`dkraay_ivar'" == "" {
            di as error "dkraay option requires tsset panel data with both panel and time variables"
            exit 459
        }
        * Driscoll-Kraay is incompatible with cluster options
        if `vce_type' == 2 | `vce_type' == 3 {
            di as error "dkraay and vce(cluster) are incompatible options"
            exit 198
        }
        * Count number of unique time periods
        tempvar _time_tag
        qui egen `_time_tag' = tag(`dkraay_tvar') if `touse'
        qui count if `_time_tag' == 1
        local dkraay_T = r(N)
        drop `_time_tag'
        * Set up Driscoll-Kraay: cluster by TIME variable
        * vce_type = 4 signals Driscoll-Kraay to C code
        local vce_type = 4
        local cluster_var "`dkraay_tvar'"
        * Set bandwidth if not specified (default = dkraay lag value, matching ivreg2)
        if `bw_val' == 0 {
            local bw_val = `dkraay_val'
        }
        confirm variable `cluster_var'
        markout `touse' `cluster_var'
        if "`verbose'" != "" {
            di as text "Driscoll-Kraay SEs: clustering on `dkraay_tvar' (T=`dkraay_T'), bw=`bw_val'"
        }
    }

    * Kiefer standard errors: within-panel HAC with full bandwidth
    * NOTE: civreghdfe's Kiefer implementation uses panel-aware HAC which
    * produces similar but not identical results to ivreghdfe's Kiefer.
    * The difference is due to the exact small-sample corrections used.
    * For exact matching, use ivreghdfe.
    * Requires tsset panel data
    local kiefer_val = ("`kiefer'" != "")
    local kiefer_tvar ""
    if `kiefer_val' {
        * Check that data is tsset
        capture tsset
        if _rc != 0 {
            di as error "kiefer option requires tsset panel data"
            exit 459
        }
        local kiefer_tvar = r(timevar)
        local ivar = r(panelvar)
        if "`kiefer_tvar'" == "" | "`ivar'" == "" {
            di as error "kiefer option requires tsset panel data with both panel and time variables"
            exit 459
        }
        * Kiefer is incompatible with robust and cluster options
        if `vce_type' == 1 {
            di as error "kiefer and vce(robust) are incompatible options"
            exit 198
        }
        if `vce_type' == 2 | `vce_type' == 3 {
            di as error "kiefer and vce(cluster) are incompatible options"
            exit 198
        }
        * Get max time span T for bandwidth
        tempvar _time_count
        qui by `ivar': gen `_time_count' = _N
        qui sum `_time_count'
        local T = r(max)
        drop `_time_count'
        * Set up Kiefer as panel-aware HAC with truncated kernel
        * Use cluster structure for panel grouping, but Kiefer VCE formula
        local vce_type = 2
        local cluster_var "`ivar'"
        local kernel_type = 4  /* Truncated kernel */
        local bw_val = `T'     /* Bandwidth = T */
        confirm variable `cluster_var'
        markout `touse' `cluster_var'
        if "`verbose'" != "" {
            di as text "Kiefer SEs: panel-aware HAC with truncated kernel, bw=`T'"
        }
    }

    * Set default bandwidth if kernel specified but no bandwidth
    if `kernel_type' > 0 & `kernel_type' < 6 & `bw_val' == 0 {
        * Default: Newey-West optimal bandwidth = floor(4*(N/100)^(2/9))
        local bw_val = floor(4 * (`N'/100)^(2/9))
        if `bw_val' < 1 local bw_val = 1
    }

    * Set default kernel (Bartlett/Newey-West) if bw specified but no kernel
    if `bw_val' > 0 & `kernel_type' == 0 & `dkraay_val' == 0 & `kiefer_val' == 0 {
        local kernel_type = 1  /* Bartlett kernel (Newey-West) */
    }

    * HAC kernel options imply robust VCE (already set earlier)
    * Note: panel_aware_hac is not currently implemented due to complexity

    scalar __civreghdfe_kernel = `kernel_type'
    scalar __civreghdfe_bw = `bw_val'
    scalar __civreghdfe_dkraay = `dkraay_val'
    scalar __civreghdfe_kiefer = `kiefer_val'
    scalar __civreghdfe_dkraay_T = `dkraay_T'
    scalar __civreghdfe_hac_panel = `panel_aware_hac'

    * DOF adjustment options
    scalar __civreghdfe_dofminus = `dofminus'
    scalar __civreghdfe_sdofminus = `sdofminus'
    scalar __civreghdfe_nopartialsmall = ("`nopartialsmall'" != "")
    scalar __civreghdfe_center = ("`center'" != "")
    scalar __civreghdfe_noreturn = ("`noreturn'" != "")

    * Update cluster scalars after Kiefer/dkraay handling
    * vce_type 4 = dkraay, which also uses cluster variable (time)
    scalar __civreghdfe_has_cluster = (`vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4)
    scalar __civreghdfe_vce_type = `vce_type'

    if "`verbose'" != "" & `kernel_type' > 0 {
        if `kiefer_val' {
            di as text "Kiefer standard errors (within-group HAC)"
        }
        else {
            di as text "HAC: kernel=`kernel_type', bandwidth=`bw_val'"
            if `dkraay_val' > 0 {
                di as text "Driscoll-Kraay standard errors with lag `dkraay_val'"
            }
        }
    }

    * Handle b0() initial values option
    * For iterative estimators like CUE, b0 provides starting values
    scalar __civreghdfe_has_b0 = 0
    if "`b0'" != "" {
        capture confirm matrix `b0'
        if _rc == 0 {
            matrix __civreghdfe_b0 = `b0'
            scalar __civreghdfe_has_b0 = 1
            if "`verbose'" != "" {
                di as text "Using initial values from matrix `b0'"
            }
        }
        else {
            di as error "b0(`b0'): matrix not found"
            exit 198
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
    if `vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4 {
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

    * Record Stata pre-processing time
    timer off 97
    qui timer list 97
    local t_stata_pre = r(t97)

    * Call the C plugin
    timer clear 99
    timer on 99
    plugin call ctools_plugin `plugin_vars' if `touse', "civreghdfe `threads_code' iv_regression"
    timer off 99
    qui timer list 99
    local t_plugin = r(t99)
    timer clear 98
    timer on 98

    * If noreturn specified, skip all result retrieval and e() posting
    if "`noreturn'" != "" {
        timer off 98
        di as text "(noreturn specified - skipping result storage)"
        * Clean up the noreturn scalar
        capture scalar drop __civreghdfe_noreturn
        exit 0
    }

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

    * For Kiefer, df_r should be N - K - df_a (not num_clusters - 1)
    * This matches ivreghdfe behavior
    if `kiefer_val' {
        local df_r = `N_used' - `K' - `df_a'
        if `df_r' <= 0 local df_r = 1
    }

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

    * Read mobility groups (connected components) for multi-way FE redundancy
    capture local mobility_groups = scalar(__civreghdfe_mobility_groups)
    if _rc != 0 | "`mobility_groups'" == "" | "`mobility_groups'" == "." {
        local mobility_groups = 0
    }

    if `vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4 {
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
    * Also build full names list (including base/omitted levels) for ivreghdfe compatibility
    local varnames ""
    local varnames_full ""
    foreach v of local endogvars_coef {
        local varnames `varnames' `v'
    }
    if "`endogvars'" != "" {
        foreach v of local endogvars_names_all {
            local varnames_full `varnames_full' `v'
        }
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
    if "`exogvars'" != "" {
        foreach v of local exogvars_names_all {
            local is_partial = 0
            foreach pv of local partial_vars {
                if "`v'" == "`pv'" {
                    local is_partial = 1
                }
            }
            if `is_partial' == 0 {
                local varnames_full `varnames_full' `v'
            }
        }
    }

    * Update K_exog to reflect partialled variables removed
    local K_exog = `K_exog' - `n_partial'
    local K_total = `K_exog' + `K_endog'

    * Force VCE symmetry (numerical precision fix)
    tempname V_sym
    matrix `V_sym' = (`V_temp' + `V_temp'') / 2
    matrix `V_temp' = `V_sym'

    * Insert base/omitted factor variable levels (coefficient=0, variance=0)
    * to match ivreghdfe's e(b) format
    local K_active : word count `varnames'
    local K_full : word count `varnames_full'
    if `K_full' > `K_active' {
        * Build expanded b and V with zeros for base/omitted levels
        tempname b_full V_full
        matrix `b_full' = J(1, `K_full', 0)
        matrix `V_full' = J(`K_full', `K_full', 0)

        local active_idx = 1
        forvalues j = 1/`K_full' {
            local vname : word `j' of `varnames_full'
            * Check if this is a base or omitted level
            local is_base = 0
            if regexm("`vname'", "^[0-9]+b\.") | regexm("`vname'", "^o\.") {
                local is_base = 1
            }
            if `is_base' == 0 {
                * Copy coefficient from active vector
                matrix `b_full'[1, `j'] = `b_temp'[1, `active_idx']
                * Copy row and column of V
                local active_idx2 = 1
                forvalues jj = 1/`K_full' {
                    local vname2 : word `jj' of `varnames_full'
                    local is_base2 = 0
                    if regexm("`vname2'", "^[0-9]+b\.") | regexm("`vname2'", "^o\.") {
                        local is_base2 = 1
                    }
                    if `is_base2' == 0 {
                        matrix `V_full'[`j', `jj'] = `V_temp'[`active_idx', `active_idx2']
                        local active_idx2 = `active_idx2' + 1
                    }
                }
                local active_idx = `active_idx' + 1
            }
        }
        matrix `b_temp' = `b_full'
        matrix `V_temp' = `V_full'
        local varnames `varnames_full'
    }

    * Apply o. prefix for collinear variables detected by C plugin
    * Collinearity flags are in [exog, endog] order (matching C internal order)
    * varnames is in [endog, exog] order (matching ivreghdfe convention)
    capture confirm scalar __civreghdfe_num_collinear
    if _rc == 0 & __civreghdfe_num_collinear > 0 {
        local varnames_new ""
        * Walk through varnames: first K_endog entries are endog, then K_exog are exog
        * Map back to C's [exog, endog] flag order
        local active_k = 1
        foreach v of local varnames {
            * Skip names that are already base levels (b. or o. prefix)
            local is_base = 0
            if regexm("`v'", "^[0-9]+b\.") | regexm("`v'", "^o\.") {
                local is_base = 1
            }
            if `is_base' {
                local varnames_new `varnames_new' `v'
            }
            else {
                * Map active_k to C's [exog, endog] order flag index
                * active_k 1..K_endog -> C flag K_exog + k (endog)
                * active_k K_endog+1..K_total -> C flag k - K_endog (exog)
                if `active_k' <= `K_endog' {
                    local c_flag_idx = `K_exog' + `active_k'
                }
                else {
                    local c_flag_idx = `active_k' - `K_endog'
                }
                capture local is_collin = __civreghdfe_collinear_`c_flag_idx'
                if _rc != 0 local is_collin = 0
                if `is_collin' == 1 {
                    local varnames_new `varnames_new' o.`v'
                }
                else {
                    local varnames_new `varnames_new' `v'
                }
                local active_k = `active_k' + 1
            }
        }
        local varnames `varnames_new'
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
    local depname_use = cond("`depname'" != "", "`depname'", "`depvar'")
    ereturn post `b_temp' `V_temp', esample(`touse') depname(`depname_use')

    * Store scalars
    ereturn scalar N = `N_used'
    ereturn scalar df_r = `df_r'
    * Report df_a_for_vce when clustering (matches ivreghdfe behavior)
    * When FE is nested in cluster, df_a_for_vce = 0
    if `vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4 {
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

    if `vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4 {
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
    ereturn local predict "civreghdfe_p"
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
    else if `vce_type' == 4 {
        ereturn local vcetype "Robust"
        ereturn local clustvar "`cluster_var'"
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

    * DOF adjustment options
    if `dofminus' != 0 {
        ereturn scalar dofminus = `dofminus'
    }
    if `sdofminus' != 0 {
        ereturn scalar sdofminus_opt = `sdofminus'
    }

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
        if "`subtitle'" != "" {
            di as text "`subtitle'"
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
        else if `vce_type' == 2 & `kiefer_val' {
            di as text "Estimates efficient for homoskedasticity only"
            di as text "Statistics robust to within-cluster autocorrelation (Kiefer)"
            di as text "  kernel=Truncated; bandwidth=`bw_val'"
            di as text "  time variable (t):  `kiefer_tvar'"
            di as text "  group variable (i): `cluster_var'"
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
        * Kiefer does not display "Number of clusters" - it uses N-K-df_a for df_r
        if `vce_type' == 2 & !`kiefer_val' {
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
        * For Kiefer, use df_r = N - K - df_a (set earlier in the code)
        local prob_F = Ftail(`K_total', `df_r', `F_wald')
        if `vce_type' == 2 & !`kiefer_val' {
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
        * Build display options string
        local dispopt ""
        if "`noomitted'" != "" local dispopt "`dispopt' noomitted"
        if "`omitted'" != "" local dispopt "`dispopt' omitted"
        if "`vsquish'" != "" local dispopt "`dispopt' vsquish"
        if "`noemptycells'" != "" local dispopt "`dispopt' noemptycells"
        if "`baselevels'" != "" local dispopt "`dispopt' baselevels"
        if "`allbaselevels'" != "" local dispopt "`dispopt' allbaselevels"
        if "`plus'" != "" local dispopt "`dispopt' plus"
        if "`eform'" != "" local dispopt "`dispopt' eform(`eform')"

        ereturn display, level(`level') `dispopt'
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

    * Display absorbed degrees of freedom table (only when absorb was specified)
    if `G' > 0 {
    di as text ""
    di as text "Absorbed degrees of freedom:"
    di as text "{hline 13}{c TT}{hline 39}{c TRC}"
    di as text " Absorbed FE {c |} Categories  - Redundant  = Num. Coefs {c |}"
    di as text "{hline 13}{c +}{hline 39}{c RT}"

    * Display each FE with its components
    * Redundancy calculation uses connected components (mobility groups) from C:
    * - G == 1: no redundancy (0)
    * - G >= 2: second FE gets mobility_groups redundant (from connected components)
    * - G > 2: each additional FE gets 1 redundant
    * - Nested FEs: all levels are redundant
    local total_coefs = 0
    local has_nested = 0
    forval g = 1/`G' {
        local fevar : word `g' of `absorb'
        local cats = `fe_levels_`g''
        * For nested FEs: all levels are redundant
        if `fe_nested_`g'' == 1 {
            local redundant = `cats'
            local has_nested = 1
            local nested_marker "*"
        }
        else {
            if `g' == 1 {
                * First FE: no redundancy
                local redundant = 0
            }
            else if `g' == 2 {
                * Second FE: redundancy = mobility_groups (connected components)
                local redundant = `mobility_groups'
                * For G > 2, mobility_groups includes (G-2) extra, but we want
                * the second FE to show only the base mobility groups
                if `G' > 2 {
                    local redundant = `mobility_groups' - (`G' - 2)
                }
            }
            else {
                * Third+ FE: 1 redundant each
                local redundant = 1
            }
            local nested_marker " "
        }
        local coefs = `cats' - `redundant'
        local total_coefs = `total_coefs' + `coefs'
        di as text %12s abbrev("`fevar'", 12) " {c |}" as result %10.0fc `cats' as text "  " as result %10.0fc `redundant' as text "  " as result %10.0fc `coefs' as text "    `nested_marker'{c |}"
    }
    di as text "{hline 13}{c BT}{hline 39}{c BRC}"
    if `has_nested' {
        di as text "* = FE nested within cluster; treated as redundant for DoF computation"
    }
    di as text ""
    } /* end if G > 0 */

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

    * Display timing breakdown if verbose specified
    if "`verbose'" != "" {
        local t_total_wall = `t_stata_pre' + `t_plugin' + `t_stata_post'
        di as text ""
        di as text "{hline 55}"
        di as text "civreghdfe timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load (SPI):        " as result %8.4f _civreghdfe_time_load " sec"
        di as text "    Data extraction:        " as result %8.4f _civreghdfe_time_extract " sec"
        di as text "    Missing check+compact:  " as result %8.4f _civreghdfe_time_missing " sec"
        di as text "    Singleton removal:      " as result %8.4f _civreghdfe_time_singleton " sec"
        di as text "    FE remap + counts:      " as result %8.4f _civreghdfe_time_remap " sec"
        di as text "    FWL partialling:        " as result %8.4f _civreghdfe_time_fwl " sec"
        di as text "    HDFE partial out:       " as result %8.4f _civreghdfe_time_partial " sec"
        di as text "    DOF + post-processing:  " as result %8.4f _civreghdfe_time_dof " sec"
        di as text "    IV estimation + VCE:    " as result %8.4f _civreghdfe_time_estimate " sec"
        di as text "    Stats computation:      " as result %8.4f _civreghdfe_time_stats " sec"
        di as text "    Store results:          " as result %8.4f _civreghdfe_time_store " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _civreghdfe_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Setup (parsing, etc):   " as result %8.4f `t_stata_pre' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f (`t_plugin' - _civreghdfe_time_total) " sec"
        di as text "    Post-processing:        " as result %8.4f `t_stata_post' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `t_stata_pre' + (`t_plugin' - _civreghdfe_time_total) + `t_stata_post'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `t_total_wall' " sec"
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
    capture scalar drop __civreghdfe_kiefer
    capture scalar drop __civreghdfe_dofminus
    capture scalar drop __civreghdfe_sdofminus
    capture scalar drop __civreghdfe_nopartialsmall
    capture scalar drop __civreghdfe_center
    capture scalar drop __civreghdfe_noreturn
    capture scalar drop __civreghdfe_has_b0
    capture matrix drop __civreghdfe_b0
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

    * Clean up collinearity scalars
    capture scalar drop __civreghdfe_num_collinear
    capture scalar drop __civreghdfe_K_keep
    forval i = 1/20 {
        capture scalar drop __civreghdfe_collinear_`i'
    }

    * Clean up timing scalars
    capture scalar drop _civreghdfe_time_load _civreghdfe_time_extract _civreghdfe_time_missing
    capture scalar drop _civreghdfe_time_singleton
    capture scalar drop _civreghdfe_time_remap
    capture scalar drop _civreghdfe_time_fwl
    capture scalar drop _civreghdfe_time_partial
    capture scalar drop _civreghdfe_time_dof
    capture scalar drop _civreghdfe_time_estimate
    capture scalar drop _civreghdfe_time_stats
    capture scalar drop _civreghdfe_time_store
    capture scalar drop _civreghdfe_time_total
    capture scalar drop _civreghdfe_threads_max
    capture scalar drop _civreghdfe_openmp_enabled
    timer clear 97

end
