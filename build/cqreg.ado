*! version 1.0.0 17Jan2026
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
*!   tolerance(#)       - Convergence tolerance (default: 1e-8)
*!   maxiter(#)         - Maximum iterations (default: 200)
*!   nopreprocess(#)    - Controls preprocessing algorithm (experimental)
*!                        Default: 0 (preprocessing disabled)
*!                        Set to -1 to enable preprocessing (experimental, may be slow)

program define cqreg, eclass
    version 14.0

    * Capture full command line before parsing
    local cmdline "cqreg `0'"

    syntax varlist(min=2 fv) [if] [in], [Quantile(real 0.5) Absorb(varlist) ///
        VCE(string) DENmethod(string) BWmethod(string) Verbose TIMEit TOLerance(real 1e-8) MAXiter(integer 200) NOPReprocess(integer 0)]

    * Validate quantile
    if `quantile' <= 0 | `quantile' >= 1 {
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
    local nvars : word count `varlist'
    local K_x = `nvars' - 1  // number of independent variables
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

    * Set up parameters via Stata scalars (expected by C plugin)
    scalar __cqreg_quantile = `quantile'
    scalar __cqreg_K = `nvars'              // depvar + indepvars
    scalar __cqreg_G = `nfe'                // number of FE groups
    scalar __cqreg_vce_type = `vcetype'
    scalar __cqreg_bw_method = `bwmethod_num'
    scalar __cqreg_density_method = `denmethod_num'
    scalar __cqreg_verbose = ("`verbose'" != "")
    scalar __cqreg_tolerance = `tolerance'
    scalar __cqreg_maxiter = `maxiter'
    scalar __cqreg_nopreprocess = `nopreprocess'

    * Display info if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cqreg: C-Accelerated Quantile Regression"
        di as text "{hline 60}"
        di as text "Quantile:   " as result `quantile'
        di as text "Depvar:     " as result "`depvar'"
        di as text "Indepvars:  " as result "`indepvars'"
        if `nfe' > 0 {
            di as text "Absorb:     " as result "`absorb'"
        }
        di as text "Obs:        " as result `nobs'
        di as text "VCE type:   " as result cond(`vcetype'==0, "iid", cond(`vcetype'==1, "robust", "cluster(`clustervar_orig')"))
        di as text "Density:    " as result "`denmethod_name'"
        di as text "Bandwidth:  " as result cond(`bwmethod_num'==0, "Hall-Sheather", cond(`bwmethod_num'==1, "Bofinger", "Chamberlain"))
        di as text "{hline 60}"
        di ""
    }

    * Record start time
    timer clear 99
    timer on 99

    * Build varlist for plugin: depvar indepvars [fe_vars] [cluster_var]
    local plugin_varlist `varlist'
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
        "cqreg full_regression"

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

    * Build coefficient vector
    tempname b V
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

    * Get VCE matrix
    matrix `V' = __cqreg_V[1..`K_with_cons', 1..`K_with_cons']

    * Build column names for matrices
    local colnames ""
    if `K_x' > 0 {
        foreach v of local indepvars {
            * Handle factor variables
            local vname = subinstr("`v'", ".", "_", .)
            local colnames `colnames' `vname'
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

    * Show timing if verbose or timeit
    if "`verbose'" != "" | "`timeit'" != "" {
        di as text ""
        di as text "Iterations: " as result `iterations' as text ", converged: " as result cond(`converged', "yes", "no")
        di as text "Total time: " as result %9.3f `elapsed' as text " seconds"
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
