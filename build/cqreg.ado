*! version 1.2.1 20261029 github.com/mdroste/stata-ctools

program define cqreg, eclass
    version 14.1

    * Start wall clock timer (must be before any work for accurate rmsg match)
    quietly _ctools_clock __ct_lap99

    * Drop whatever an interrupted run left behind
    _ctools_cleanup cqreg

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Capture full command line before parsing
    local cmdline "cqreg `0'"

    syntax varlist(min=2 numeric fv) [if] [in], [Quantile(real 0.5) Absorb(varlist) ///
        VCE(string) DENmethod(string) BWmethod(string) Verbose TOLerance(real 1e-12) MAXiter(integer 200) NOPReprocess(integer 0) THReads(integer 0)]

    * Least-squares partialling is not valid for the quantile-loss objective.
    if "`absorb'" != "" {
        di as error "cqreg: absorb() is not supported; use explicit factor-variable indicators"
        exit 198
    }
    if `maxiter' < 1 {
        di as error "cqreg: maxiter() must be positive"
        exit 198
    }
    * nopreprocess() is parsed for backward compatibility only
    if `nopreprocess' == -1 {
        di as error "cqreg: the experimental preprocessing solver, nopreprocess(-1), has been removed"
        exit 198
    }

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

    * Parse variable list
    gettoken depvar indepvars : varlist

    * Expand factor variables to temporary numeric variables
    * This handles i.varname, c.varname, etc.
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

    * Parse VCE option
    local vcetype = 0   // 0=iid, 1=robust, 2=cluster
    local clustervar ""
    local clustervar_orig ""
    if `"`vce'"' != "" {
        _ctools_vce vce_kind clustervar : `vce'
        if "`vce_kind'" == "cluster" {
            local vcetype = 2
            if "`clustervar'" == "" {
                di as error "cqreg: cluster variable required with vce(cluster)"
                exit 198
            }
            local clustervar_orig "`clustervar'"
            confirm variable `clustervar', exact
            markout `touse' `clustervar_orig', strok
            * Preserve equality of fractional, large, and string identifiers:
            * those are numbered by _ctools_group. A byte/int/long variable is
            * passed as it is; the plugin numbers its values itself, in the
            * same order.
            local __cl_type : type `clustervar'
            if !inlist("`__cl_type'", "byte", "int", "long") {
                tempvar clustervar_numeric
                quietly _ctools_group `clustervar' if `touse', generate(`clustervar_numeric')
                local clustervar "`clustervar_numeric'"
            }
        }
        else if "`vce_kind'" == "robust" {
            local vcetype = 1
        }
        else if "`vce_kind'" == "iid" {
            local vcetype = 0
        }
        else if "`vce_kind'" == "bootstrap" {
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
    local kernel_name ""
    if "`denmethod'" != "" {
        local den_lower = lower(trim("`denmethod'"))
        local den_kopt = trim(substr("`den_lower'", 7, .))
        if "`den_lower'" == "residual" | "`den_lower'" == "res" {
            local denmethod_num = 0
            local denmethod_name "residual"
        }
        else if "`den_lower'" == "fitted" | "`den_lower'" == "fit" {
            local denmethod_num = 1
            local denmethod_name "fitted"
        }
        else if substr("`den_lower'", 1, 6) == "kernel" & ///
            ("`den_kopt'" == "" | (substr("`den_kopt'", 1, 1) == "(" & ///
            substr("`den_kopt'", -1, 1) == ")")) {
            * kernel or kernel(name), as qreg's vce(, kernel(name))
            if "`den_kopt'" != "" {
                local den_kopt = substr("`den_kopt'", 2, strlen("`den_kopt'") - 2)
            }
            _cqreg_kernel, `den_kopt'
            local kernel_name "`r(kernel)'"
            local denmethod_num = 2
            local denmethod_name "kernel"
        }
        else {
            di as error "cqreg: unrecognized density method: `denmethod'"
            di as error "       Valid options: fitted (default), residual, kernel[(kernel)]"
            exit 198
        }
    }

    * qreg defines the residual method only for the iid VCE
    if `denmethod_num' == 0 & `vcetype' != 0 {
        local __vce_name = cond(`vcetype' == 1, "robust", "cluster")
        di as error "cqreg: denmethod(residual) is invalid with vce(`__vce_name');"
        di as error "       the residual method may only be combined with vce(iid)"
        exit 184
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
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

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
    scalar __cqreg_vce_type = `vcetype'
    scalar __cqreg_bw_method = `bwmethod_num'
    scalar __cqreg_density_method = `denmethod_num'
    scalar __cqreg_verbose = ("`verbose'" != "")
    scalar __cqreg_tolerance = `tolerance'
    scalar __cqreg_maxiter = `maxiter'
    scalar __cqreg_bandwidth = `__bw'

    * Kernel method: qreg's KernelBWidth factor invnormal(q+h) - invnormal(q-h),
    * missing when q-h <= 0 or q+h >= 1 (the plugin then exits with r(498))
    if `denmethod_num' == 2 {
        local __kernels "epanechnikov epan2 biweight cosine gaussian parzen rectangle triangle"
        local __kpos : list posof "`kernel_name'" in __kernels
        scalar __cqreg_kernel = `__kpos' - 1
        scalar __cqreg_kfactor = cond(`quantile' - `__bw' > 0 & ///
            `quantile' + `__bw' < 1, ///
            invnormal(`quantile' + `__bw') - invnormal(`quantile' - `__bw'), .)
    }

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Build varlist for plugin: depvar indepvars [cluster_var]
    * Use expanded varlist (factor variables converted to temp numeric vars)
    local plugin_varlist `varlist_expanded'
    if `vcetype' == 2 {
        local plugin_varlist `plugin_varlist' `clustervar'
    }
    * The sample marker goes last, in place of if `touse': the plugin reads
    * it in parallel, where Stata answers an if one observation at a time
    local plugin_varlist `plugin_varlist' `touse'
    scalar __cqreg_touse = 1

    * Create matrix to store VCE (plugin will fill __cqreg_V)
    matrix __cqreg_V = J(`K_x' + 1, `K_x' + 1, 0)

    * Call the C plugin with full_regression subcommand
    capture noisily plugin call ctools_plugin `plugin_varlist', ///
        "cqreg `threads_code' full_regression"

    local reg_rc = _rc
    if `reg_rc' {
        _ctools_cleanup cqreg
        di as error "cqreg: regression failed (rc=`reg_rc')"
        exit `reg_rc'
    }

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

    * Get number of clusters if cluster VCE
    local num_clusters = 0
    if `vcetype' == 2 {
        capture local num_clusters = __cqreg_N_clust
        if _rc != 0 local num_clusters = .
    }

    * e(b) and e(V) in qreg's layout: every term in order, with zeros for
    * omitted regressors (o. prefix) and for factor-variable base levels
    tempname b V
    _ctools_ebv `b' `V', names(`coef_names_all') base(`coef_is_base') ///
        stub(__cqreg) kkeep(`K_keep') vmatrix(__cqreg_V) cons(__cqreg_cons)
    local K_full = colsof(`b')

    * Post results
    ereturn post `b' `V', esample(`touse') depname(`depvar') obs(`N_final')

    * Compute pseudo R-squared
    local r2_p = 1 - `sum_adev' / `sum_rdev'

    * Compute density at quantile (f_r = 1/sparsity)
    local f_r = 1 / `sparsity'

    * df_r for quantile regression
    local df_r = `N_final' - `K_keep' - 1
    if `df_r' < 1 local df_r = 1

    * Store scalar results (matching qreg)
    _post_vce_rank
    ereturn scalar sparsity = `sparsity'
    ereturn scalar bwidth = `bandwidth'
    if `denmethod_num' == 2 {
        ereturn scalar kbwidth = __cqreg_kbwidth
    }
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
    ereturn scalar iterations = `iterations'
    ereturn scalar r2_p = `r2_p'

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
    if `denmethod_num' == 2 {
        ereturn local kernel "`kernel_name'"
    }
    ereturn local bwmethod = cond(`bwmethod_num'==0, "hsheather", cond(`bwmethod_num'==1, "bofinger", "chamberlain"))
    ereturn local depvar "`depvar'"
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

    di as text ""

    * Display coefficient table
    ereturn display

    * Stop wall clock timer (after all work, before display)
    _ctools_lap elapsed __ct_lap99 `__ct_lap99'

    * Show timing breakdown if verbose
    if "`verbose'" != "" {
        * Calculate Stata overhead
        local __stata_overhead = `elapsed' - _cqreg_time_total

        * Get sub-timers (may not exist in older builds)
        capture local __time_sparsity = _cqreg_time_sparsity
        if _rc != 0 local __time_sparsity = .
        capture local __time_vce_matrix = _cqreg_time_vce_matrix
        if _rc != 0 local __time_vce_matrix = .
        capture local __ipm_iterations = _cqreg_ipm_iterations
        if _rc != 0 local __ipm_iterations = `iterations'

        * IPM sub-timers
        capture local __time_ipm_init = _cqreg_time_ipm_init
        if _rc != 0 local __time_ipm_init = .
        capture local __time_ipm_iterate = _cqreg_time_ipm_iterate
        if _rc != 0 local __time_ipm_iterate = .
        capture local __time_ipm_crossover = _cqreg_time_ipm_crossover
        if _rc != 0 local __time_ipm_crossover = .

        * Sparsity sub-timers
        capture local __time_sparsity_aux = _cqreg_time_sparsity_aux
        if _rc != 0 local __time_sparsity_aux = .

        di as text ""
        di as text "{hline 55}"
        di as text "cqreg timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f _cqreg_time_load " sec"
        di as text "    IPM solver:             " as result %8.4f _cqreg_time_ipm " sec"
        if `__time_ipm_init' != . {
            di as text "      Init (OLS warm start):" as result %8.4f `__time_ipm_init' " sec"
            di as text "      Iterate:              " as result %8.4f `__time_ipm_iterate' " sec"
            if `__ipm_iterations' > 0 {
                local __per_iter = `__time_ipm_iterate' / `__ipm_iterations'
                di as text "        (" as result `__ipm_iterations' as text " iters, " as result %6.4f `__per_iter' as text " sec/iter)"
            }
            di as text "      Crossover:            " as result %8.4f `__time_ipm_crossover' " sec"
        }
        else {
            di as text "      (" as result `__ipm_iterations' as text " iterations)"
        }
        di as text "    VCE total:              " as result %8.4f _cqreg_time_vce " sec"
        if `__time_sparsity' != . {
            di as text "      Sparsity estimation:  " as result %8.4f `__time_sparsity' " sec"
            if `__time_sparsity_aux' != . {
                di as text "        Auxiliary QR solves:" as result %8.4f `__time_sparsity_aux' " sec"
            }
            di as text "      VCE matrix:           " as result %8.4f `__time_vce_matrix' " sec"
        }
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
    }

    _ctools_cleanup cqreg
end

* Kernel of denmethod(kernel(name)), with the names and abbreviations of
* qreg's vce(, kernel(name)); epanechnikov is the default.
program define _cqreg_kernel, rclass
    capture syntax [, EPanechnikov epan2 BIweight COSine GAUssian PARzen RECtangle TRIangle]
    if _rc {
        local 0 : subinstr local 0 "," ""
        local 0 : list retokenize 0
        di as error "cqreg: invalid kernel in denmethod(kernel(`0'))"
        di as error "       Valid kernels: epanechnikov, epan2, biweight, cosine, gaussian,"
        di as error "       parzen, rectangle, triangle"
        exit 198
    }
    local kernel `epanechnikov' `epan2' `biweight' `cosine' `gaussian' `parzen' `rectangle' `triangle'
    local kernel : list retokenize kernel
    if `: word count `kernel'' > 1 {
        di as error "cqreg: denmethod(kernel()) allows only one kernel"
        exit 198
    }
    if "`kernel'" == "" local kernel epanechnikov
    return local kernel `kernel'
end
