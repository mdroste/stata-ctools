*! version 0.9.0 26Jan2026
*! cqreg: C-accelerated quantile regression
*! Part of the ctools suite
*!
*! Description:
*!   High-performance replacement for qreg using a C plugin
*!   with Frisch-Newton solver and optional HDFE support.
*!
*! Syntax:
*!   cqreg depvar indepvars [if] [in], [options]
*!
*! Options:
*!   quantile(#)        - Quantile to estimate (default: 0.5)
*!   absorb(varlist)    - Fixed effects to absorb (optional)
*!   vce(vcetype)       - Variance estimation: iid (default), robust, cluster
*!   denmethod(method)  - Density estimation: fitted (default), residual
*!   bwmethod(method)   - Bandwidth method: hsheather (default), bofinger, chamberlain
*!   verbose            - Display progress information
*!   tolerance(#)       - Convergence tolerance (default: 1e-12)
*!   maxiter(#)         - Maximum iterations (default: 200)
*!   nopreprocess(#)    - Controls preprocessing algorithm (experimental)
*!                        Default: 0 (preprocessing disabled)
*!                        Set to -1 to enable preprocessing (experimental, may be slow)

program define cqreg, eclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Capture full command line before parsing
    local cmdline "cqreg `0'"

    syntax varlist(min=2 numeric fv) [if] [in], [Quantile(real 0.5) Absorb(varlist) ///
        VCE(string) DENmethod(string) BWmethod(string) Verbose TIMEIT TOLerance(real 1e-12) MAXiter(integer 200) NOPReprocess(integer 0) THReads(integer 0)]

    * Validate quantile
    if `quantile' >= 1 {
        di as error "cqreg: quantile must be between 0 and 1"
        exit 498
    }
    if `quantile' <= 0 {
        di as error "cqreg: quantile must be between 0 and 1"
        exit 198
    }

    * Mark sample
    marksample touse
    if "`absorb'" != "" {
        markout `touse' `absorb'
    }

    * Parse variable list
    gettoken depvar indepvars : varlist

    * Expand factor variables to temporary numeric variables
    * This handles i.varname, c.varname, etc.
    if "`indepvars'" != "" {
        * Step 1: Expand factor variables to temp numeric vars
        fvrevar `indepvars' if `touse'
        local indepvars_fvrevar "`r(varlist)'"

        * Step 2: Get coefficient names from fvexpand (same order as fvrevar)
        * These names include base level notation (e.g., "1b.rep78 2.rep78 3.rep78")
        fvexpand `indepvars' if `touse'
        local coef_names_all "`r(varlist)'"

        * Step 3: Remove collinear variables (including base levels of factors)
        * This matches what Stata's estimation commands do internally
        _rmcoll `indepvars_fvrevar' if `touse', forcedrop
        local indepvars_expanded "`r(varlist)'"

        * Step 4: Build mapping - identify which positions are base levels
        * fvrevar and fvexpand have 1:1 positional correspondence
        * A position is a base level if its temp var was dropped by _rmcoll
        local n_all : word count `coef_names_all'
        local coef_is_base ""
        local coef_names ""
        local plugin_pos = 1
        forval i = 1/`n_all' {
            local fv_var : word `i' of `indepvars_fvrevar'
            local fv_name : word `i' of `coef_names_all'
            * Check if this temp var is in the remaining (non-dropped) list
            local is_kept = 0
            foreach remaining of local indepvars_expanded {
                if "`remaining'" == "`fv_var'" {
                    local is_kept = 1
                }
            }
            if `is_kept' {
                * This variable was kept (not a base level)
                local coef_is_base `coef_is_base' 0
                local coef_names `coef_names' `fv_name'
            }
            else {
                * This variable was dropped (base level or omitted)
                local coef_is_base `coef_is_base' 1
            }
        }

        * Build expanded varlist for plugin (only non-base vars)
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
    local K_x = `nvars' - 1  // number of independent variables (after expansion)
    local nfe = 0
    if "`absorb'" != "" {
        local nfe : word count `absorb'
    }

    * Parse VCE option
    local vcetype = 0   // 0=iid, 1=robust, 2=cluster
    local clustervar ""
    local clustervar_orig ""
    if `"`vce'"' != "" {
        local vce_lower = lower(trim(`"`vce'"'))
        if substr("`vce_lower'", 1, 2) == "cl" {
            local vcetype = 2
            gettoken vce_type clustervar : vce
            local clustervar = trim("`clustervar'")
            if "`clustervar'" == "" {
                di as error "cqreg: cluster variable required with vce(cluster)"
                exit 198
            }
            local clustervar_orig "`clustervar'"
            * Check if cluster variable is string - convert to numeric if so
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
            markout `touse' `clustervar_orig'
        }
        else if "`vce_lower'" == "robust" | "`vce_lower'" == "r" {
            local vcetype = 1
        }
        else if "`vce_lower'" == "iid" {
            local vcetype = 0
        }
        else if substr("`vce_lower'", 1, 4) == "boot" {
            di as error "cqreg: vce(bootstrap) is not currently supported"
            di as error "       Use Stata's qreg for bootstrap standard errors"
            exit 198
        }
        else {
            di as error "cqreg: unrecognized vce() option: `vce'"
            di as error "       Valid options: iid, robust, cluster"
            exit 198
        }
    }

    * Parse bandwidth method
    local bwmethod_num = 0  // 0=hsheather, 1=bofinger, 2=chamberlain
    if "`bwmethod'" != "" {
        local bw_lower = lower(trim("`bwmethod'"))
        if "`bw_lower'" == "bofinger" {
            local bwmethod_num = 1
        }
        else if "`bw_lower'" == "chamberlain" {
            local bwmethod_num = 2
        }
        else if "`bw_lower'" == "hsheather" | "`bw_lower'" == "hall-sheather" {
            local bwmethod_num = 0
        }
        else {
            di as error "cqreg: unrecognized bandwidth method: `bwmethod'"
            di as error "       Valid options: hsheather, bofinger, chamberlain"
            exit 198
        }
    }

    * Parse density estimation method
    local denmethod_num = 1  // 0=residual, 1=fitted (default), 2=kernel
    local denmethod_name "fitted"
    if "`denmethod'" != "" {
        local den_lower = lower(trim("`denmethod'"))
        if "`den_lower'" == "residual" | "`den_lower'" == "res" {
            local denmethod_num = 0
            local denmethod_name "residual"
        }
        else if "`den_lower'" == "fitted" | "`den_lower'" == "fit" {
            local denmethod_num = 1
            local denmethod_name "fitted"
        }
        else if "`den_lower'" == "kernel" {
            local denmethod_num = 2
            local denmethod_name "kernel"
        }
        else {
            di as error "cqreg: unrecognized density method: `denmethod'"
            di as error "       Valid options: fitted (default), residual"
            exit 198
        }
    }

    * Count observations
    qui count if `touse'
    local nobs = r(N)

    if `nobs' == 0 {
        di as error "cqreg: no observations"
        exit 2000
    }

    if `nobs' <= `K_x' + 1 {
        di as error "cqreg: not enough observations"
        exit 2001
    }

    * Note: q_v and sum_rdev are now computed in the C plugin for better performance

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
            di as error "cqreg: Could not load ctools plugin"
            exit 601
        }
    }

    * Pre-compute bandwidth in Stata for exact precision match with qreg.
    * Using Stata's invnormal()/normalden() avoids C approximation errors.
    local __bw_alpha = 0.05
    local __invtau = invnormal(`quantile')
    local __phi_invtau = normalden(`__invtau')
    local __invtau2 = (`__invtau')^2
    if `bwmethod_num' == 0 {
        * Hall-Sheather
        local __z_alpha = invnormal(1 - `__bw_alpha'/2)
        local __bw = (`nobs'^(-1/3)) * (`__z_alpha'^(2/3)) * ///
            (1.5*`__phi_invtau'^2 / (2*`__invtau2'+1))^(1/3)
    }
    else if `bwmethod_num' == 1 {
        * Bofinger
        local __bw = `nobs'^(-1/5) * ((4.5*`__phi_invtau'^4) / ///
            (2*`__invtau2'+1)^2)^(1/5)
    }
    else {
        * Chamberlain
        local __z_alpha = invnormal(1 - `__bw_alpha'/2)
        local __bw = `__z_alpha' * sqrt(`nobs'^(-1)*`quantile'*(1-`quantile'))
    }

    * Set up parameters via Stata scalars (expected by C plugin)
    scalar __cqreg_quantile = `quantile'
    scalar __cqreg_K = `nvars'              // depvar + indepvars
    scalar __cqreg_G = `nfe'                // number of FE groups
    scalar __cqreg_vce_type = `vcetype'
    scalar __cqreg_bw_method = `bwmethod_num'
    scalar __cqreg_density_method = `denmethod_num'
    scalar __cqreg_verbose = ("`verbose'" != "" | "`timeit'" != "")
    scalar __cqreg_tolerance = `tolerance'
    scalar __cqreg_maxiter = `maxiter'
    scalar __cqreg_nopreprocess = `nopreprocess'
    scalar __cqreg_bandwidth = `__bw'

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record start time
    timer clear 99
    timer on 99

    * Build varlist for plugin: depvar indepvars [fe_vars] [cluster_var]
    * Use expanded varlist (factor variables converted to temp numeric vars)
    local plugin_varlist `varlist_expanded'
    if `nfe' > 0 {
        local plugin_varlist `plugin_varlist' `absorb'
    }
    if `vcetype' == 2 {
        local plugin_varlist `plugin_varlist' `clustervar'
    }

    * Create matrix to store VCE (plugin will fill __cqreg_V)
    matrix __cqreg_V = J(`K_x' + 1, `K_x' + 1, 0)

    * Call the C plugin with full_regression subcommand
    capture noisily plugin call ctools_plugin `plugin_varlist' if `touse', ///
        "cqreg `threads_code' full_regression"

    local reg_rc = _rc
    if `reg_rc' {
        * Clean up scalars on error
        capture scalar drop __cqreg_quantile __cqreg_K __cqreg_G
        capture scalar drop __cqreg_vce_type __cqreg_bw_method __cqreg_verbose
        capture scalar drop __cqreg_tolerance __cqreg_maxiter
        capture matrix drop __cqreg_V
        di as error "cqreg: regression failed (rc=`reg_rc')"
        exit `reg_rc'
    }

    timer off 99
    quietly timer list 99
    local elapsed = r(t99)

    * Retrieve results from scalars set by C plugin
    local N_final = __cqreg_N
    local K_keep = __cqreg_K_keep
    local sum_adev = __cqreg_sum_adev
    local sum_rdev = __cqreg_sum_rdev
    local q_v = __cqreg_q_v
    local sparsity = __cqreg_sparsity
    local bandwidth = __cqreg_bandwidth
    local iterations = __cqreg_iterations
    local converged = __cqreg_converged
    local cons = __cqreg_cons

    * Get df_a if HDFE was used
    local df_a = 0
    if `nfe' > 0 {
        capture local df_a = __cqreg_df_a
        if _rc != 0 local df_a = 0
    }

    * Get number of clusters if cluster VCE
    local num_clusters = 0
    if `vcetype' == 2 {
        capture local num_clusters = __cqreg_N_clust
        if _rc != 0 local num_clusters = .
    }

    * Build coefficient vector - including base level columns with zeros
    tempname b V b_full V_full

    * Count total coefficients including base levels
    local n_all : word count `coef_names_all'
    local K_full = `n_all' + 1  // all coefficients plus constant

    * First build the compact b and V from plugin (only non-base vars)
    local K_with_cons = `K_keep' + 1
    matrix `b' = J(1, `K_with_cons', 0)
    forval k = 1/`K_keep' {
        capture scalar define __tmp = __cqreg_beta_`k'
        if _rc == 0 {
            matrix `b'[1, `k'] = __tmp
        }
        capture scalar drop __tmp
    }
    matrix `b'[1, `K_with_cons'] = `cons'
    matrix `V' = __cqreg_V[1..`K_with_cons', 1..`K_with_cons']

    * Now expand to include base level columns
    if `n_all' > 0 & `n_all' > `K_keep' {
        * There are base levels to insert
        matrix `b_full' = J(1, `K_full', 0)
        matrix `V_full' = J(`K_full', `K_full', 0)

        * Map plugin positions to full positions
        local plugin_pos = 1
        forval i = 1/`n_all' {
            local is_base : word `i' of `coef_is_base'
            if `is_base' == 0 {
                * Copy coefficient from plugin
                matrix `b_full'[1, `i'] = `b'[1, `plugin_pos']
                * Copy VCE row and column
                local plugin_pos2 = 1
                forval j = 1/`n_all' {
                    local is_base2 : word `j' of `coef_is_base'
                    if `is_base2' == 0 {
                        matrix `V_full'[`i', `j'] = `V'[`plugin_pos', `plugin_pos2']
                        local plugin_pos2 = `plugin_pos2' + 1
                    }
                }
                * Copy constant's row/col entries
                matrix `V_full'[`i', `K_full'] = `V'[`plugin_pos', `K_with_cons']
                matrix `V_full'[`K_full', `i'] = `V'[`K_with_cons', `plugin_pos']
                local plugin_pos = `plugin_pos' + 1
            }
            * Base levels keep their zeros
        }
        * Copy constant coefficient and VCE diagonal
        matrix `b_full'[1, `K_full'] = `cons'
        matrix `V_full'[`K_full', `K_full'] = `V'[`K_with_cons', `K_with_cons']

        * Use the full matrices
        matrix `b' = `b_full'
        matrix `V' = `V_full'
        local K_with_cons = `K_full'
    }

    * Build column names for matrices - use ALL names including base levels
    local colnames ""
    if `n_all' > 0 {
        foreach v of local coef_names_all {
            local colnames `colnames' `v'
        }
    }
    else if `K_x' > 0 {
        foreach v of local coef_names {
            local colnames `colnames' `v'
        }
    }
    local colnames `colnames' _cons

    * Apply names to matrices
    matrix colnames `b' = `colnames'
    matrix rownames `V' = `colnames'
    matrix colnames `V' = `colnames'

    * Post results
    ereturn post `b' `V', esample(`touse') depname(`depvar') obs(`N_final')

    * Compute pseudo R-squared
    local r2_p = 1 - `sum_adev' / `sum_rdev'

    * Compute density at quantile (f_r = 1/sparsity)
    local f_r = 1 / `sparsity'

    * df_r for quantile regression
    local df_r = `N_final' - `K_keep' - 1
    if `nfe' > 0 {
        local df_r = `df_r' - `df_a'
    }
    if `df_r' < 1 local df_r = 1

    * Store scalar results (matching qreg)
    ereturn scalar rank = `K_with_cons'
    ereturn scalar sparsity = `sparsity'
    ereturn scalar bwidth = `bandwidth'
    ereturn scalar df_m = `K_keep'
    ereturn scalar df_r = `df_r'
    ereturn scalar f_r = `f_r'
    ereturn scalar N = `N_final'
    ereturn scalar sum_w = `N_final'
    ereturn scalar q_v = `q_v'
    ereturn scalar q = `quantile'
    ereturn scalar sum_rdev = `sum_rdev'
    ereturn scalar sum_adev = `sum_adev'
    ereturn scalar convcode = cond(`converged', 0, 1)
    ereturn scalar r2_p = `r2_p'

    if `nfe' > 0 {
        ereturn scalar df_a = `df_a'
    }

    if `vcetype' == 2 {
        ereturn scalar N_clust = `num_clusters'
    }

    * Store macro results (matching qreg)
    ereturn local cmd "cqreg"
    ereturn local cmdline "`cmdline'"
    ereturn local predict "qreg_p"
    ereturn local properties "b V"
    ereturn local marginsnotok "stdp stddp residuals"
    ereturn local vce = cond(`vcetype'==0, "iid", cond(`vcetype'==1, "robust", "cluster"))
    ereturn local denmethod "`denmethod_name'"
    ereturn local bwmethod = cond(`bwmethod_num'==0, "hsheather", cond(`bwmethod_num'==1, "bofinger", "chamberlain"))
    ereturn local depvar "`depvar'"
    if `nfe' > 0 {
        ereturn local absorb "`absorb'"
    }
    if `vcetype' == 2 {
        ereturn local clustvar "`clustervar_orig'"
    }

    * Display results (matching qreg format exactly)
    di as text ""

    * Title line - match qreg exactly
    if `quantile' == 0.5 {
        di as text "Median regression" _col(52) "Number of obs =" _col(68) as result %10.0fc `N_final'
    }
    else {
        di as text "`quantile' Quantile regression" _col(52) "Number of obs =" _col(68) as result %10.0fc `N_final'
    }

    * Raw sum of deviations line - match qreg format
    di as text "  Raw sum of deviations " as result %-8.6g `sum_rdev' as text " (about " as result %-9.7g `q_v' as text ")"

    * Min sum of deviations line
    di as text "  Min sum of deviations " as result %-8.6g `sum_adev' _col(52) as text "Pseudo R2     =" _col(68) as result %10.4f `r2_p'

    if `nfe' > 0 {
        di as text "  Absorbing " as result `nfe' as text " HDFE group(s)" _col(52) "DF absorbed   =" _col(68) as result %10.0fc `df_a'
    }
    di as text ""

    * Display coefficient table
    ereturn display

    * Show timing breakdown if verbose or timeit
    if "`verbose'" != "" | "`timeit'" != "" {
        * Calculate Stata overhead
        local __stata_overhead = `elapsed' - _cqreg_time_total

        di as text ""
        di as text "{hline 55}"
        di as text "cqreg timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f _cqreg_time_load " sec"
        if `nfe' > 0 {
            di as text "    HDFE partial out:       " as result %8.4f _cqreg_time_hdfe " sec"
        }
        di as text "    IPM solver:             " as result %8.4f _cqreg_time_ipm " sec"
        di as text "    VCE computation:        " as result %8.4f _cqreg_time_vce " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _cqreg_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Stata overhead:           " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `elapsed' " sec"
        di as text "{hline 55}"
        di as text ""
        di as text "  Iterations: " as result `iterations' as text ", converged: " as result cond(`converged', "yes", "no")

        * Display thread diagnostics
        capture local __threads_max = _cqreg_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _cqreg_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }

        * Clean up timing scalars
        capture scalar drop _cqreg_time_load _cqreg_time_hdfe
        capture scalar drop _cqreg_time_ipm _cqreg_time_vce _cqreg_time_total
        capture scalar drop _cqreg_threads_max _cqreg_openmp_enabled
    }

    * Clean up scalars
    capture scalar drop __cqreg_quantile __cqreg_K __cqreg_G
    capture scalar drop __cqreg_vce_type __cqreg_bw_method __cqreg_density_method __cqreg_verbose
    capture scalar drop __cqreg_tolerance __cqreg_maxiter __cqreg_nopreprocess
    capture scalar drop __cqreg_N __cqreg_K_keep __cqreg_sum_adev __cqreg_sum_rdev __cqreg_q_v
    capture scalar drop __cqreg_sparsity __cqreg_bandwidth
    capture scalar drop __cqreg_iterations __cqreg_converged __cqreg_cons
    capture scalar drop __cqreg_df_a __cqreg_N_clust
    forval k = 1/`K_keep' {
        capture scalar drop __cqreg_beta_`k'
    }
    capture matrix drop __cqreg_V

end
