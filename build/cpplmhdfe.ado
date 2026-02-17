*! version 1.0.0 11feb2026 github.com/mdroste/stata-ctools

program define cpplmhdfe, eclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start wall clock timer
    local __do_timing = 0
    timer clear 98
    timer on 98

    syntax varlist(min=2 fv ts) [aw fw pw] [if] [in] [, Absorb(string) VCE(string) Verbose ///
        TOLerance(real 1e-8) ITERATE(integer 10000) THReads(integer 0) ///
        EXPosure(varname) OFFset(varname) ///
        IRLSTOLerance(real 1e-12) IRLSMAXiter(integer 1000) ///
        SEPTOLerance(real 1e-8) ///
        DOFadjustments(string)]

    local __do_timing = ("`verbose'" != "")

    * Parse absorb option
    if `"`absorb'"' == "" {
        di as error "cpplmhdfe: absorb() is required"
        exit 198
    }

    * Parse absorb option for suboptions
    local absorb_vars ""
    if strpos(`"`absorb'"', ",") > 0 {
        gettoken absorb_vars absorb_opts : absorb, parse(",")
        local absorb_vars = trim("`absorb_vars'")
    }
    else {
        local absorb_vars "`absorb'"
    }
    local absorb "`absorb_vars'"

    * Parse dofadjustments option
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

    * Mark sample
    marksample touse
    markout `touse' `absorb'

    * Parse weights
    local weight_var ""
    local weight_type = 0
    if "`weight'" != "" {
        local weight_var "`exp'"
        local weight_var = subinstr("`weight_var'", "=", "", .)
        local weight_var = trim("`weight_var'")

        if "`weight'" == "aweight" {
            local weight_type = 1
        }
        else if "`weight'" == "fweight" {
            local weight_type = 2
        }
        else if "`weight'" == "pweight" {
            local weight_type = 3
        }
        markout `touse' `weight_var'
        qui count if `weight_var' <= 0 & `touse'
        if r(N) > 0 {
            di as error "cpplmhdfe: weights must be positive"
            exit 198
        }
    }

    * Handle exposure/offset
    local has_offset = 0
    local offset_var ""
    if "`exposure'" != "" {
        * exposure(var) means offset(log(var))
        tempvar log_exposure
        markout `touse' `exposure'
        qui count if `exposure' <= 0 & `touse'
        if r(N) > 0 {
            di as error "cpplmhdfe: exposure variable must be positive"
            exit 198
        }
        quietly gen double `log_exposure' = ln(`exposure') if `touse'
        local offset_var "`log_exposure'"
        local has_offset = 1
    }
    else if "`offset'" != "" {
        markout `touse' `offset'
        local offset_var "`offset'"
        local has_offset = 1
    }

    * Parse variable list
    gettoken depvar indepvars : varlist
    local depvar_orig "`depvar'"

    * Expand time-series operators
    tsrevar `depvar'
    local depvar "`r(varlist)'"

    * Check y >= 0
    qui count if `depvar' < 0 & `touse'
    if r(N) > 0 {
        di as error "cpplmhdfe: dependent variable must be non-negative"
        exit 198
    }

    * Expand factor variables
    if "`indepvars'" != "" {
        fvrevar `indepvars' if `touse'
        local __indepvars_all "`r(varlist)'"
        fvexpand `indepvars' if `touse'
        local coef_names_all "`r(varlist)'"

        * Filter base/omitted levels
        local __vars_nobase ""
        local __names_nobase ""
        local __n_all : word count `coef_names_all'
        forval __k = 1/`__n_all' {
            local __vname : word `__k' of `coef_names_all'
            local __keep = 0
            local __parts : subinstr local __vname "#" " ", all
            foreach __p of local __parts {
                if !regexm("`__p'", "b\.") & !regexm("`__p'", "o\.") {
                    local __keep = 1
                }
            }
            if `__keep' {
                local __vars_nobase `__vars_nobase' `: word `__k' of `__indepvars_all''
                local __names_nobase `__names_nobase' `__vname'
            }
        }

        local indepvars_expanded "`__vars_nobase'"
        local coef_names "`__names_nobase'"
        local varlist_expanded "`depvar' `indepvars_expanded'"
    }
    else {
        local indepvars_expanded ""
        local varlist_expanded "`depvar'"
        local coef_names ""
    }

    local nvars : word count `varlist_expanded'
    local nfe : word count `absorb'

    * Parse VCE option (default: robust, matching ppmlhdfe)
    local vcetype = 1
    local clustervar ""
    local clustervar_orig ""

    if `"`vce'"' != "" {
        local vce_lower = lower(trim(`"`vce'"'))
        if substr("`vce_lower'", 1, 2) == "cl" {
            local vcetype = 2
            gettoken vce_type clustervar : vce
            local clustervar = trim("`clustervar'")
            if "`clustervar'" == "" {
                di as error "cpplmhdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            local clustervar_orig "`clustervar'"
            capture confirm numeric variable `clustervar'
            if _rc != 0 {
                tempvar clustervar_numeric
                capture which gegen
                if _rc == 0 {
                    quietly gegen `clustervar_numeric' = group(`clustervar')
                }
                else {
                    quietly egen `clustervar_numeric' = group(`clustervar')
                }
                local clustervar "`clustervar_numeric'"
            }
            capture confirm numeric variable `clustervar_orig'
            if _rc == 0 {
                markout `touse' `clustervar_orig'
            }
            else {
                quietly replace `touse' = 0 if `clustervar_orig' == "" & `touse' == 1
            }
        }
        else if "`vce_lower'" == "robust" | "`vce_lower'" == "r" {
            local vcetype = 1
        }
        else {
            di as error "cpplmhdfe: unrecognized vce() option"
            exit 9
        }
    }

    * Count observations
    qui count if `touse'
    local nobs = r(N)
    if `nobs' == 0 {
        di as error "cpplmhdfe: no observations"
        exit 2000
    }

    * Load plugin
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
            di as error "cpplmhdfe: Could not load ctools plugin"
            exit 601
        }
    }

    * Set up scalars for C plugin
    scalar __cpplmhdfe_K = `nvars'
    scalar __cpplmhdfe_G = `nfe'
    scalar __cpplmhdfe_verbose = ("`verbose'" != "")
    scalar __cpplmhdfe_maxiter = `iterate'
    scalar __cpplmhdfe_tolerance = `tolerance'
    scalar __cpplmhdfe_vce_type = `vcetype'
    scalar __cpplmhdfe_has_weights = (`weight_type' > 0)
    scalar __cpplmhdfe_weight_type = `weight_type'
    scalar __cpplmhdfe_has_offset = `has_offset'
    scalar __cpplmhdfe_irls_maxiter = `irlsmaxiter'
    scalar __cpplmhdfe_irls_tol = `irlstolerance'
    scalar __cpplmhdfe_sep_tol = `septolerance'
    scalar __cpplmhdfe_dof_adjust_type = `dof_adjust_type'

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

    * Build varlist for plugin: depvar indepvars fe_vars [cluster_var] [weight_var] [offset_var]
    local plugin_varlist `varlist_expanded' `absorb'
    if `vcetype' == 2 {
        local plugin_varlist `plugin_varlist' `clustervar'
    }
    if `weight_type' > 0 {
        local plugin_varlist `plugin_varlist' `weight_var'
    }
    if `has_offset' {
        local plugin_varlist `plugin_varlist' `offset_var'
    }

    * Create V matrix
    local K_x = `nvars' - 1
    if `K_x' < 1 {
        local K_x = 1
    }
    matrix __cpplmhdfe_V = J(`K_x', `K_x', 0)

    * Call C plugin
    capture noisily plugin call ctools_plugin `plugin_varlist' if `touse', ///
        "cpplmhdfe `threads_code' full_regression"

    local reg_rc = _rc
    if `reg_rc' {
        capture scalar drop __cpplmhdfe_K __cpplmhdfe_G __cpplmhdfe_verbose
        capture scalar drop __cpplmhdfe_maxiter __cpplmhdfe_tolerance
        capture scalar drop __cpplmhdfe_vce_type __cpplmhdfe_has_weights
        capture scalar drop __cpplmhdfe_weight_type __cpplmhdfe_has_offset
        capture scalar drop __cpplmhdfe_irls_maxiter __cpplmhdfe_irls_tol
        capture scalar drop __cpplmhdfe_sep_tol __cpplmhdfe_dof_adjust_type
        capture matrix drop __cpplmhdfe_V
        di as error "Error in PPML regression (rc=`reg_rc')"
        exit `reg_rc'
    }

    * Record plugin call time
    timer off 98
    quietly timer list 98
    local __time_plugin_call = r(t98)
    timer clear 98
    timer on 98

    * Retrieve results
    local N_final = __cpplmhdfe_N
    local num_singletons = __cpplmhdfe_num_singletons
    local num_separated = __cpplmhdfe_num_separated
    local K_keep = __cpplmhdfe_K_keep
    local df_a = __cpplmhdfe_df_a
    local mobility_groups = __cpplmhdfe_mobility_groups
    local deviance = __cpplmhdfe_deviance
    local ll = __cpplmhdfe_ll
    local ll_0 = __cpplmhdfe_ll_0
    local irls_iters = __cpplmhdfe_irls_iterations
    local irls_converged = __cpplmhdfe_irls_converged

    * Get number of levels and nested status per FE
    local df_a_nested = 0
    forval g = 1/`nfe' {
        capture local fe_levels_`g' = __cpplmhdfe_num_levels_`g'
        if _rc != 0 local fe_levels_`g' = .
        capture local fe_nested_`g' = __cpplmhdfe_fe_nested_`g'
        if _rc != 0 local fe_nested_`g' = 0
        if `fe_nested_`g'' == 1 {
            local df_a_nested = `df_a_nested' + `fe_levels_`g''
        }
    }

    local df_a_adjusted = `df_a' - `df_a_nested'
    if `df_a_adjusted' < 0 {
        local df_a_adjusted = 0
    }

    * df_m = number of non-collinear X variables
    local df_m = `K_keep'

    * df_r calculation
    local df_r_ols_base = `N_final' - `K_keep' - `df_a_adjusted'
    if `vcetype' == 2 {
        local num_clusters = __cpplmhdfe_N_clust
        local df_r_cluster = `num_clusters' - 1
        if `df_r_ols_base' < `df_r_cluster' {
            local df_r = `df_r_ols_base'
        }
        else {
            local df_r = `df_r_cluster'
        }
    }
    else {
        local df_r = `df_r_ols_base'
    }
    if `df_r' < 1 {
        local df_r = 1
    }

    * Pseudo R-squared (McFadden)
    local r2_p = 1 - `ll' / `ll_0'

    * Build coefficient vector and names
    tempname b V

    local K_full = `K_x'
    matrix `b' = J(1, `K_full', 0)
    matrix `V' = J(`K_full', `K_full', 0)

    * Fill in coefficients and VCE for non-collinear variables
    if `K_keep' > 0 {
        local kept_idx = 1
        forval k = 1/`K_x' {
            local is_collin = __cpplmhdfe_collinear_`k'
            if `is_collin' == 0 {
                matrix `b'[1, `k'] = __cpplmhdfe_beta_`kept_idx'
                local kept_idx2 = 1
                forval j = 1/`K_x' {
                    local is_collin_j = __cpplmhdfe_collinear_`j'
                    if `is_collin_j' == 0 {
                        matrix `V'[`k', `j'] = __cpplmhdfe_V[`kept_idx', `kept_idx2']
                        local ++kept_idx2
                    }
                }
                local ++kept_idx
            }
        }
    }

    * Build column names
    local colnames ""
    forval k = 1/`K_x' {
        local vname : word `k' of `coef_names'
        if `K_keep' == 0 {
            local colnames `colnames' o.`vname'
        }
        else {
            local is_collin = __cpplmhdfe_collinear_`k'
            if `is_collin' == 1 {
                local colnames `colnames' o.`vname'
            }
            else {
                local colnames `colnames' `vname'
            }
        }
    }

    matrix colnames `b' = `colnames'
    matrix rownames `V' = `colnames'
    matrix colnames `V' = `colnames'

    * Wald F-statistic
    local F = .
    if `K_keep' > 0 {
        tempname b_test V_test Vinv Wald
        matrix `b_test' = `b'[1, 1..`K_keep']
        matrix `V_test' = `V'[1..`K_keep', 1..`K_keep']
        matrix `Vinv' = syminv(`V_test')
        matrix `Wald' = `b_test' * `Vinv' * `b_test''
        local F = `Wald'[1,1] / `df_m'
    }

    * Expand b and V to include base levels for factor variables
    if "`coef_names_all'" != "" {
        local K_all : word count `coef_names_all'
        if `K_all' > `K_x' {
            local K_full_new = `K_all'
            tempname b_exp V_exp
            matrix `b_exp' = J(1, `K_full_new', 0)
            matrix `V_exp' = J(`K_full_new', `K_full_new', 0)

            local __aidx = 1
            forvalues __j = 1/`K_all' {
                local __vn_all : word `__j' of `coef_names_all'
                local __vn_act : word `__aidx' of `coef_names'
                if "`__vn_all'" == "`__vn_act'" {
                    local __bmap_`__j' = `__aidx'
                    local __aidx = `__aidx' + 1
                }
                else {
                    local __bmap_`__j' = 0
                }
            }

            forvalues __j = 1/`K_all' {
                if `__bmap_`__j'' > 0 {
                    local __aj = `__bmap_`__j''
                    matrix `b_exp'[1, `__j'] = `b'[1, `__aj']
                    forvalues __jj = 1/`K_all' {
                        if `__bmap_`__jj'' > 0 {
                            local __ajj = `__bmap_`__jj''
                            matrix `V_exp'[`__j', `__jj'] = `V'[`__aj', `__ajj']
                        }
                    }
                }
            }

            matrix `b' = `b_exp'
            matrix `V' = `V_exp'
            local K_full = `K_full_new'

            local colnames ""
            forvalues __j = 1/`K_all' {
                local __vn : word `__j' of `coef_names_all'
                if `__bmap_`__j'' == 0 {
                    local colnames `colnames' `__vn'
                }
                else {
                    local __aj = `__bmap_`__j''
                    if `K_keep' == 0 {
                        local colnames `colnames' o.`__vn'
                    }
                    else {
                        local is_collin = __cpplmhdfe_collinear_`__aj'
                        if `is_collin' == 1 {
                            local colnames `colnames' o.`__vn'
                        }
                        else {
                            local colnames `colnames' `__vn'
                        }
                    }
                }
            }

            matrix colnames `b' = `colnames'
            matrix rownames `V' = `colnames'
            matrix colnames `V' = `colnames'
        }
    }

    * Record results building time
    timer off 98
    quietly timer list 98
    local __time_build = r(t98)
    timer clear 98
    timer on 98

    * Post results
    ereturn post `b' `V', esample(`touse') depname(`depvar_orig') obs(`N_final')

    * Store e() results
    ereturn scalar N = `N_final'
    ereturn scalar df_m = `df_m'
    ereturn scalar df_r = `df_r'
    ereturn scalar ll = `ll'
    ereturn scalar ll_0 = `ll_0'
    ereturn scalar deviance = `deviance'
    ereturn scalar r2_p = `r2_p'
    ereturn scalar F = `F'
    ereturn scalar rank = `df_m'
    ereturn scalar N_hdfe = `nfe'
    ereturn scalar num_singletons = `num_singletons'
    ereturn scalar num_separated = `num_separated'
    ereturn scalar ic = `irls_iters'
    ereturn scalar converged = `irls_converged'
    ereturn scalar df_a = `df_a_adjusted'
    ereturn scalar df_a_initial = `df_a'
    ereturn scalar df_a_nested = `df_a_nested'

    if `vcetype' == 2 {
        ereturn scalar N_clust = __cpplmhdfe_N_clust
        ereturn local clustvar "`clustervar_orig'"
    }

    ereturn local absorb "`absorb'"
    ereturn local depvar "`depvar_orig'"
    ereturn local indepvars "`indepvars'"
    ereturn local vce = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", "cluster"))
    ereturn local cmd "cpplmhdfe"
    ereturn local cmdline "cpplmhdfe `0'"

    * Record ereturn posting time
    timer off 98
    quietly timer list 98
    local __time_ereturn = r(t98)
    timer clear 98
    timer on 98

    * Display
    if `num_singletons' > 0 {
        di as text "(dropped " as result `num_singletons' as text " singleton observations)"
    }
    if `num_separated' > 0 {
        di as text "(dropped " as result `num_separated' as text " separated observations)"
    }
    if `irls_converged' {
        di as text "(IRLS converged in " as result `irls_iters' as text " iterations)"
    }
    else {
        di as error "(IRLS did NOT converge in " as result `irls_iters' as text " iterations)"
    }

    di as text ""
    di as text "Poisson pseudo-likelihood regression" _col(49) "Number of obs" _col(67) "= " as result %10.0fc `N_final'
    di as text "Absorbing " as result `nfe' as text " HDFE " cond(`nfe'>1, "groups", "group") _col(49) "F(" as result %3.0f `df_m' as text "," as result %8.0f `df_r' as text ")" _col(67) "= " as result %10.2f e(F)
    di as text _col(49) "Prob > F" _col(67) "= " as result %10.4f Ftail(`df_m', `df_r', e(F))
    di as text _col(49) "Pseudo R2" _col(67) "= " as result %10.4f `r2_p'
    di as text _col(49) "Log pseudolikelihood" _col(67) "= " as result %10.4f `ll'

    di as text ""
    ereturn display

    * Absorbed degrees of freedom table
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
        if `fe_nested_`g'' == 1 {
            local redundant = `cats'
            local has_nested = 1
            local nested_marker "*"
        }
        else {
            if `g' == 1 | (`has_nested' == 0 & `g' == 1) {
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
        di as text %12s abbrev("`fevar'", 12) " {c |}" as result %10.0fc `cats' as text "  " as result %10.0fc `redundant' as text "  " as result %10.0fc `coefs' as text "    `nested_marker'{c |}"
    }
    di as text "{hline 13}{c BT}{hline 39}{c BRC}"
    if `has_nested' {
        di as text "* = FE nested within cluster; treated as redundant for DoF computation"
    }

    * Clean up scalars
    capture scalar drop __cpplmhdfe_K __cpplmhdfe_G __cpplmhdfe_verbose
    capture scalar drop __cpplmhdfe_maxiter __cpplmhdfe_tolerance
    capture scalar drop __cpplmhdfe_vce_type __cpplmhdfe_has_weights
    capture scalar drop __cpplmhdfe_weight_type __cpplmhdfe_has_offset
    capture scalar drop __cpplmhdfe_irls_maxiter __cpplmhdfe_irls_tol
    capture scalar drop __cpplmhdfe_sep_tol __cpplmhdfe_dof_adjust_type
    capture scalar drop __cpplmhdfe_N __cpplmhdfe_num_singletons
    capture scalar drop __cpplmhdfe_num_separated __cpplmhdfe_K_keep
    capture scalar drop __cpplmhdfe_df_a __cpplmhdfe_mobility_groups
    capture scalar drop __cpplmhdfe_deviance __cpplmhdfe_ll __cpplmhdfe_ll_0
    capture scalar drop __cpplmhdfe_irls_iterations __cpplmhdfe_irls_converged
    capture scalar drop __cpplmhdfe_N_clust __cpplmhdfe_num_collinear
    forval k = 1/20 {
        capture scalar drop __cpplmhdfe_beta_`k'
        capture scalar drop __cpplmhdfe_collinear_`k'
        capture scalar drop __cpplmhdfe_num_levels_`k'
        capture scalar drop __cpplmhdfe_fe_nested_`k'
    }
    capture matrix drop __cpplmhdfe_V

    * Timing cleanup
    capture scalar drop _cpplmhdfe_time_load _cpplmhdfe_time_remap
    capture scalar drop _cpplmhdfe_time_singleton _cpplmhdfe_time_dof
    capture scalar drop _cpplmhdfe_time_irls _cpplmhdfe_time_vce
    capture scalar drop _cpplmhdfe_time_total
    capture scalar drop _cpplmhdfe_threads_max _cpplmhdfe_openmp_enabled

    * Display timing if verbose
    timer off 98
    quietly timer list 98
    local __time_display = r(t98)
    local __time_total = `__time_setup' + `__time_plugin_call' + `__time_build' + `__time_ereturn' + `__time_display'

    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "cpplmhdfe timing breakdown:"
        di as text "{hline 55}"
        di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
        di as text "  C plugin total:           " as result %8.4f `__time_plugin_call' " sec"
        di as text "  Results building:         " as result %8.4f `__time_build' " sec"
        di as text "  ereturn posting:          " as result %8.4f `__time_ereturn' " sec"
        di as text "  Display output:           " as result %8.4f `__time_display' " sec"
        di as text "{hline 55}"
        di as text "  Wall clock total:         " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"
    }

    timer clear 98

end
