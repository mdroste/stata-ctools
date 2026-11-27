*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! original: Stata's qreg (qreg.ado 4.4.2, May 8, 2020)
*! original: see help cqreg for details

* byable(recall) and the sw and mi properties as qreg declares them: by:,
* stepwise, nestreg and mi estimate accept cqreg
program define cqreg, eclass byable(recall) prop(sw mi)
    version 14.1

    if replay() {
        if _by() error 190
        if "`e(cmd)'" != "cqreg" error 301
        syntax [, Level(cilevel) *]
        _get_diopts diopts, `options'
        _cqreg_display, level(`level') `diopts'
        exit
    }

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

    syntax varlist(min=1 numeric fv) [if] [in] [fw iw pw], [Quantile(real 0.5) Absorb(varlist) ///
        VCE(string) DENmethod(string) BWmethod(string) Verbose TOLerance(real 1e-12) MAXiter(integer 200) NOPReprocess(integer 0) THReads(integer 0) ///
        Level(cilevel) NOLOg LOg WLSiter(integer 1) ITERate(string) TRACE *]
    * qreg's display options; log and nolog show or hide the iteration log
    _get_diopts diopts, `options'

    * qreg's wlsiter() and iterate() set the weighted least-squares start
    * and the iteration limit of its simplex solver. They are checked as qreg
    * checks them; the interior-point solver uses neither (its limit is
    * maxiter()).
    if `wlsiter' < 1 {
        di as err "{p}{bf:wlsiter(`wlsiter')} must be a positive integer{p_end}"
        exit 198
    }
    if `"`iterate'"' != "" {
        capture confirm integer number `iterate'
        if _rc | real(`"`iterate'"') <= 0 {
            di as err "{p}{bf:iterate({it:#})} must be a nonnegative integer{p_end}"
            exit 198
        }
    }
    * qreg under version 12 or earlier reports its old standard errors,
    * which cqreg does not reproduce
    if _caller() < 13 {
        di as text "(note: cqreg reports qreg's standard errors of version 13 and later, not those of version " _caller() ")"
    }

    * Least-squares partialling is not valid for the quantile-loss objective.
    if "`absorb'" != "" {
        di as error "cqreg: absorb() is not supported; use explicit factor-variable indicators"
        exit 198
    }
    if `maxiter' < 1 {
        di as error "cqreg: maxiter() must be positive"
        exit 198
    }
    * nopreprocess(#): # > 0 turns off the preprocessing of large samples
    * (the plugin then solves every fit on all observations); 0, the
    * default, leaves it on. nopreprocess(-1) once selected an experimental
    * solver.
    if `nopreprocess' == -1 {
        di as error "cqreg: the experimental preprocessing solver, nopreprocess(-1), has been removed"
        exit 198
    }

    if `nopreprocess' < 0 | missing(`nopreprocess') {
        di as error "cqreg: nopreprocess() must be nonnegative"
        exit 198
    }
    if `tolerance' <= 0 | `tolerance' > 1 | missing(`tolerance') {
        di as error "cqreg: tolerance() must be greater than zero and at most one"
        exit 198
    }

    * Validate quantile; as in qreg, a value of 1 or more is a percentage
    if `quantile' >= 1 local quantile = `quantile' / 100
    if `quantile' <= 0 | `quantile' >= 1 {
        di as error "cqreg: quantile must be between 0 and 1 (or a percentage)"
        exit 198
    }

    * Mark sample
    marksample touse

    * Parse VCE option
    local vcetype = 0   // 0=iid, 1=robust, 2=cluster
    local clustervar ""
    local clustervar_orig ""
    * qreg's vce([type] [, density bandwidth]) syntax. The density method and
    * bandwidth may instead be given in denmethod() and bwmethod().
    local vce_method
    local vce_kernel
    local vce_bw
    if `"`vce'"' != "" {
        _cqreg_vceparse `vce'
        local vce_kind `r(vce)'
        local clustervar `r(clustvar)'
        local vce_method `r(method)'
        local vce_kernel `r(kernel)'
        local vce_bw `r(bwidth)'
        local vce_type_given `vce_kind'
        if "`vce_kind'" == "" local vce_kind iid
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
    }
    * An empty vce() type is the default (iid, robust with pweights), as in
    * qreg: the suboptions alone may be given
    if "`vce_type_given'" == "" local vce

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

    * Bandwidth and density methods, with qreg's names and abbreviations
    * (bofinger, hsheather, chamberlain; fitted, residual or residuals,
    * kernel[(kernel)]), from bwmethod() and denmethod() or from the vce()
    * suboptions. Given in both places, they must agree.
    local bw_name
    if `"`bwmethod'"' != "" {
        _cqreg_bwparse `bwmethod'
        local bw_name `r(bwidth)'
    }
    if "`vce_bw'" != "" {
        if "`bw_name'" != "" & "`bw_name'" != "`vce_bw'" {
            di as error "cqreg: bwmethod(`bwmethod') conflicts with the vce() bandwidth `vce_bw'"
            exit 198
        }
        local bw_name `vce_bw'
    }
    if "`bw_name'" == "" local bw_name hsheather
    local bwmethod_num = cond("`bw_name'" == "hsheather", 0, cond("`bw_name'" == "bofinger", 1, 2))

    local den_name
    local den_kernel
    if `"`denmethod'"' != "" {
        _cqreg_denparse `denmethod'
        local den_name `r(method)'
        local den_kernel `r(kernel)'
    }
    if "`vce_method'" != "" {
        if "`den_name'" != "" & ("`den_name'" != "`vce_method'" | "`den_kernel'" != "`vce_kernel'") {
            di as error "cqreg: denmethod(`denmethod') conflicts with the vce() density method `vce_method'"
            exit 198
        }
        local den_name `vce_method'
        local den_kernel `vce_kernel'
    }
    if "`den_name'" == "" local den_name fitted
    local denmethod_name `den_name'
    local kernel_name `den_kernel'
    local denmethod_num = cond("`den_name'" == "residual", 0, cond("`den_name'" == "fitted", 1, 2))

    * qreg defines the residual method only for the iid VCE
    if `denmethod_num' == 0 & `vcetype' != 0 {
        local __vce_name = cond(`vcetype' == 1, "robust", "cluster")
        di as error "cqreg: denmethod(residual) is invalid with vce(`__vce_name');"
        di as error "       the residual method may only be combined with vce(iid)"
        exit 184
    }

    * Weights as qreg: pweights imply vce(robust) and exclude vce(iid) and
    * the residual method (r(404)); the plugin fits pweights as iweights
    local wtype_num = 0
    if "`weight'" == "pweight" {
        if `"`vce'"' != "" & `vcetype' == 0 {
            di as error "cqreg: option vce(iid) is not possible with pweighted data;"
            di as error "       pweights imply robust"
            exit 404
        }
        if `denmethod_num' == 0 {
            di as error "cqreg: denmethod(residual) is not possible with pweighted data"
            exit 404
        }
        if `"`vce'"' == "" local vcetype = 1
    }
    if "`weight'" == "fweight" local wtype_num = 1
    else if "`weight'" != "" local wtype_num = 2

    * Count observations
    qui count if `touse'
    local nobs = r(N)

    * The weights in double precision, as qreg's GetVCE. marksample has
    * dropped zero and missing weights and rejected negative fweights and
    * pweights; the check loss is not convex with a negative iweight.
    * With fweights, e(N) and the bandwidth count the weights.
    local nobs_w = `nobs'
    local sum_w = `nobs'
    if `wtype_num' {
        tempvar wvar
        quietly gen double `wvar' `exp' if `touse'
        capture assert `wvar' > 0 if `touse'
        if _rc {
            di as error "cqreg: negative weights are not allowed"
            exit 402
        }
        quietly summarize `wvar' if `touse', meanonly
        local sum_w = r(sum)
        if `wtype_num' == 1 local nobs_w = r(sum)
    }

    if `nobs' == 0 {
        di as error "cqreg: no observations"
        exit 2000
    }

    if `nobs' < 2 {
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
    * As qreg, the Hall-Sheather and Chamberlain bandwidths use the
    * confidence level: alpha = 1 - level/100, and N is the sum of the
    * fweights
    local __bw_alpha = 1 - `level'/100
    local __invtau = invnormal(`quantile')
    local __phi_invtau = normalden(`__invtau')
    local __invtau2 = (`__invtau')^2
    if `bwmethod_num' == 0 {
        * Hall-Sheather
        local __z_alpha = invnormal(1 - `__bw_alpha'/2)
        local __bw = (`nobs_w'^(-1/3)) * (`__z_alpha'^(2/3)) * ///
            (1.5*`__phi_invtau'^2 / (2*`__invtau2'+1))^(1/3)
    }
    else if `bwmethod_num' == 1 {
        * Bofinger
        local __bw = `nobs_w'^(-1/5) * ((4.5*`__phi_invtau'^4) / ///
            (2*`__invtau2'+1)^2)^(1/5)
    }
    else {
        * Chamberlain
        local __z_alpha = invnormal(1 - `__bw_alpha'/2)
        local __bw = `__z_alpha' * sqrt(`nobs_w'^(-1)*`quantile'*(1-`quantile'))
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
    scalar __cqreg_nopreprocess = (`nopreprocess' > 0)
    scalar __cqreg_bandwidth = `__bw'
    scalar __cqreg_wtype = `wtype_num'

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

    * Build varlist for plugin: depvar indepvars [cluster_var] [weight]
    * Use expanded varlist (factor variables converted to temp numeric vars)
    local plugin_varlist `varlist_expanded'
    if `vcetype' == 2 {
        local plugin_varlist `plugin_varlist' `clustervar'
    }
    if `wtype_num' {
        local plugin_varlist `plugin_varlist' `wvar'
    }
    * The sample marker goes last, in place of if `touse': the plugin reads
    * it in parallel, where Stata answers an if one observation at a time
    local plugin_varlist `plugin_varlist' `touse'
    scalar __cqreg_touse = 1

    * Create matrix to store VCE (plugin will fill __cqreg_V)
    matrix __cqreg_V = J(`K_x' + 1, `K_x' + 1, 0)
    * and the same covariance in the plugin's scaled response units (rank)
    matrix __cqreg_Vn = J(`K_x' + 1, `K_x' + 1, 0)
    * the observation numbers of the basis (the plugin fills the first
    * __cqreg_basis_n)
    matrix __cqreg_basis = J(1, `K_x' + 1, .)
    scalar __cqreg_basis_n = 0

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
    * e(N) is the sum of the fweights, as in qreg
    local N_final = cond(`wtype_num' == 1, `nobs_w', __cqreg_N)
    local K_keep = __cqreg_K_keep
    local sum_adev = __cqreg_sum_adev
    local sum_rdev = __cqreg_sum_rdev
    local q_v = __cqreg_q_v
    local sparsity = __cqreg_sparsity
    local bandwidth = __cqreg_bandwidth
    local iterations = __cqreg_iterations
    local converged = __cqreg_converged
    local nonunique = __cqreg_nonunique

    * Get number of clusters if cluster VCE
    local num_clusters = 0
    if `vcetype' == 2 {
        capture local num_clusters = __cqreg_N_clust
        if _rc != 0 local num_clusters = .
    }

    * e(b) and e(V) in qreg's layout: every term in order, with zeros for
    * omitted regressors (o. prefix) and for factor-variable base levels
    tempname b V
    if `"`coef_names_all'"' == "" {
        * constant-only model
        matrix `b' = (__cqreg_cons)
        matrix `V' = __cqreg_V
        matrix colnames `b' = _cons
        matrix colnames `V' = _cons
        matrix rownames `V' = _cons
    }
    else {
        _ctools_ebv `b' `V', names(`coef_names_all') base(`coef_is_base') ///
            stub(__cqreg) kkeep(`K_keep') vmatrix(__cqreg_V) cons(__cqreg_cons)
    }
    local K_full = colsof(`b')

    * Assess covariance rank in the plugin's scaled response units (the
    * spread of y near one), so it does not depend on the units of y. At
    * extreme scales Stata's ereturn post can discard even a valid
    * covariance; do not report a successful fit with a zero matrix
    * substituted by that host check.
    tempname rank_V rank_scale rank_inverse posted_scale
    matrix `rank_V' = __cqreg_Vn
    mata: st_numscalar("`rank_scale'", max(abs(st_matrix("`V'"))))
    ereturn post `b' `V', esample(`touse') depname(`depvar') obs(`N_final')
    mata: st_numscalar("`posted_scale'", max(abs(st_matrix("e(V)"))))
    if !missing(`rank_scale') & `rank_scale' > 0 & `posted_scale' == 0 {
        ereturn clear
        _ctools_cleanup cqreg
        di as error "cqreg: Stata cannot post the covariance at this response scale; rescale the dependent variable"
        exit 498
    }

    * Compute pseudo R-squared
    local r2_p = 1 - `sum_adev' / `sum_rdev'

    * Compute density at quantile (f_r = 1/sparsity)
    local f_r = 1 / `sparsity'

    * Store scalar results (matching qreg)
    matrix `rank_inverse' = invsym(`rank_V')
    ereturn scalar rank = rowsof(`rank_V') - diag0cnt(`rank_inverse')

    * df_r as qreg: N - e(rank), so directions the VCE zeroes (no density)
    * do not count
    local df_r = `N_final' - e(rank)
    ereturn scalar sparsity = `sparsity'
    ereturn scalar bwidth = `bandwidth'
    if `denmethod_num' == 2 {
        ereturn scalar kbwidth = __cqreg_kbwidth
    }
    ereturn scalar df_m = `K_keep'
    ereturn scalar df_r = `df_r'
    ereturn scalar f_r = `f_r'
    ereturn scalar N = `N_final'
    ereturn scalar sum_w = `sum_w'
    ereturn scalar q_v = `q_v'
    ereturn scalar q = `quantile'
    ereturn scalar sum_rdev = `sum_rdev'
    ereturn scalar sum_adev = `sum_adev'
    ereturn scalar convcode = cond(`converged', 0, 1)
    ereturn scalar iterations = `iterations'
    ereturn scalar r2_p = `r2_p'
    * 1 when an edge of optimal solutions leaves the reported vertex (other
    * coefficients attain the same minimum), 0 when none does, missing when
    * the solver found no basis
    ereturn scalar nonunique = cond(`nonunique' < 0, ., `nonunique')

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
    if `vcetype' != 0 ereturn local vcetype "Robust"
    ereturn local denmethod "`denmethod_name'"
    if `denmethod_num' == 2 {
        ereturn local kernel "`kernel_name'"
    }
    ereturn local bwmethod = cond(`bwmethod_num'==0, "hsheather", cond(`bwmethod_num'==1, "bofinger", "chamberlain"))
    ereturn local depvar "`depvar'"
    if `vcetype' == 2 {
        ereturn local clustvar "`clustervar_orig'"
    }
    if "`weight'" != "" {
        ereturn local wexp "`exp'"
        ereturn local wtype "`weight'"
    }
    * qreg's hidden results: the observation numbers of the basis of the
    * linear program (the observations the fitted hyperplane passes
    * through), empty when the solver identified none, and the solver's
    * convergence tolerance (qreg's simplex reports its own, 1e-10)
    local basis
    if scalar(__cqreg_basis_n) > 0 & scalar(__cqreg_basis_n) < . {
        forvalues j = 1/`=scalar(__cqreg_basis_n)' {
            local basis `basis' `=__cqreg_basis[1, `j']'
        }
    }
    ereturn hidden local basis `"`basis'"'
    ereturn hidden scalar tolerance = `tolerance'
    * and the legacy globals qreg sets; $S_E_vl is the dependent variable
    * and the coefficients' names (qreg's _rmcoll list)
    local __vl : colnames e(b)
    local __cons _cons
    local __vl : list __vl - __cons
    global S_E_mdf = e(df_m)
    global S_E_tdf = e(df_r)
    global S_E_rd = e(f_r)
    global S_E_nobs = e(N)
    global S_E_q = e(q)
    global S_E_rq = e(q_v)
    global S_E_rsd = e(sum_rdev)
    global S_E_msd = e(sum_adev)
    global S_E_vl "`depvar' `__vl'"
    global S_E_frc = e(convcode)
    global S_E_cmd "cqreg"

    * Iteration log. The interior-point steps have no counterpart in qreg's
    * iterations, so the log is qreg's final line: the iteration count and
    * the minimized sum. nolog (or set iterlog off) hides it, as in qreg.
    if "`trace'" != "" | ("`nolog'" == "" & ("`log'" != "" | "`c(iterlog)'" != "off")) {
        di as text "Iteration `iterations':  Sum of abs. weighted deviations = " ///
            as result %10.0g `sum_adev'
    }
    * trace: qreg prints the coefficients under each iteration's line; here,
    * those of the solution
    if "`trace'" != "" {
        tempname __tb
        matrix `__tb' = e(b)
        di as result _skip(12) %9.0g `__tb'[1, 1] _c
        forvalues j = 2/`=colsof(`__tb')' {
            di as result "  " %9.0g `__tb'[1, `j'] _c
        }
        di
    }

    * Display results (matching qreg format exactly)
    _cqreg_display, level(`level') `diopts'

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

* Header and coefficient table of the posted results (also on replay), as
* qreg's Replay prints them
program define _cqreg_display
    syntax [, Level(cilevel) *]
    di
    if e(q) == 0.5 {
        di in gr "Median regression" _col(53) "Number of obs = " in ye %10.0gc e(N)
    }
    else {
        di in gr e(q) " Quantile regression" _col(53) "Number of obs = " in ye %10.0gc e(N)
    }
    di in gr "  Raw sum of deviations" in ye %9.0g e(sum_rdev) ///
        in gr " (about " in ye e(q_v) in gr ")"
    di in gr "  Min sum of deviations" in ye %9.0g e(sum_adev) _col(53) ///
        in gr "Pseudo R2     = " in ye %10.4f 1 - (e(sum_adev)/e(sum_rdev)) _n

    * Display coefficient table
    _coef_table, level(`level') `options'
    if e(nonunique) == 1 {
        di as text "Note: the solution is not unique; other coefficient vectors attain the"
        di as text "      same minimum sum of deviations."
    }
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

* qreg's vce() syntax: vce([iid | robust | cluster clustvar | excluster
* clustvar] [, fitted | residual | residuals | kernel[(kernel)]
* bofinger | hsheather | chamberlain]), with qreg's abbreviations (vce(r),
* res, fit, ker(), bo, hs, ch). cluster is cqreg's name for qreg's
* excluster. Returns r(vce) (iid, robust, cluster or empty), r(clustvar),
* and r(method), r(kernel) and r(bwidth), each empty unless specified.
program define _cqreg_vceparse, rclass
    syntax [anything(name=vcetype id="vce type")] [, FITted RESidual residuals KERnel(string) *]
    local s_fitted `fitted'
    local s_residual `residual'`residuals'
    local s_kernel `"`kernel'"'
    local s_options `"`options'"'
    local vcetype = strtrim(`"`vcetype'"')
    gettoken word rest : vcetype
    local rest = strtrim(`"`rest'"')
    local lword = lower(`"`word'"')
    local len = strlen(`"`lword'"')
    local vce
    if `"`lword'"' == "" local vce
    else if `"`lword'"' == "iid" local vce iid
    else if `"`lword'"' == substr("robust", 1, `len') local vce robust
    else if `len' >= 2 & `"`lword'"' == substr("cluster", 1, `len') local vce cluster
    else if `len' >= 4 & `"`lword'"' == substr("excluster", 1, `len') local vce cluster
    else if substr(`"`lword'"', 1, 4) == "boot" {
        di as error "cqreg: vce(bootstrap) is not supported; use the bootstrap prefix,"
        di as error "       bootstrap: cqreg ..., or Stata's bsqreg"
        exit 198
    }
    else {
        di as err `"{p}option {bf:vce(`vcetype')} is not allowed{p_end}"'
        exit 198
    }
    if "`vce'" == "cluster" {
        if `"`rest'"' == "" {
            di as error "cqreg: cluster variable required with vce(cluster)"
            exit 198
        }
        local 0 `"`rest'"'
        syntax varname
        return local clustvar `varlist'
    }
    else if `"`rest'"' != "" {
        di as err `"{p}option {bf:vce(`vcetype')} is not allowed{p_end}"'
        exit 198
    }
    return local vce `vce'

    local fitted `s_fitted'
    local residual = cond("`s_residual'" != "", "residual", "")
    local kernel `"`s_kernel'"'
    local options `"`s_options'"'
    local k = ("`fitted'" != "") + ("`residual'" != "") + (`"`kernel'"' != "")
    if `k' > 1 {
        di as err "{p}invalid {bf:vce()} specification; only one of {bf:fitted}, " ///
            "{bf:residual}, or {bf:kernel()} is allowed{p_end}"
        exit 198
    }
    local method `fitted'`residual'
    if !`k' {
        * kernel without its argument
        local 0 `", `options'"'
        syntax [, Kernel *]
        if "`kernel'" != "" local method kernel
        local kernel
    }
    else if `"`kernel'"' != "" local method kernel
    if "`method'" == "kernel" {
        capture _cqreg_kernel, `kernel'
        if _rc {
            local kname = strtrim(`"`kernel'"')
            di as err `"{p}invalid {bf:vce()} specification; suboption {bf:kernel(`kname')} is not allowed{p_end}"'
            exit 198
        }
        return local kernel `r(kernel)'
    }
    return local method `method'

    local 0 `", `options'"'
    capture syntax [, BOfinger HSheather CHamberlain]
    if _rc {
        local opts = strtrim(`"`options'"')
        di as err `"{p}invalid {bf:vce()} specification; suboption {bf:`opts'} is not allowed{p_end}"'
        exit 198
    }
    local bw `bofinger' `hsheather' `chamberlain'
    if `: word count `bw'' > 1 {
        di as err "{p}invalid {bf:vce()} specification; only one of {bf:bofinger}, " ///
            "{bf:hsheather}, or {bf:chamberlain} is allowed{p_end}"
        exit 198
    }
    return local bwidth `bw'
end

* bwmethod(): bofinger, hsheather or chamberlain, with qreg's abbreviations
* (bo, hs, ch), and hall-sheather. Returns r(bwidth).
program define _cqreg_bwparse, rclass
    local spec = lower(strtrim(`"`0'"'))
    if `"`spec'"' == "hall-sheather" local spec hsheather
    local 0 `", `spec'"'
    capture syntax [, BOfinger HSheather CHamberlain]
    local bw `bofinger' `hsheather' `chamberlain'
    if _rc | `: word count `bw'' != 1 {
        di as error `"cqreg: unrecognized bandwidth method: `spec'"'
        di as error "       Valid options: hsheather, bofinger, chamberlain"
        exit 198
    }
    return local bwidth `bw'
end

* denmethod(): fitted, residual (or residuals) or kernel[(kernel)], with
* qreg's abbreviations (fit, res, ker). Returns r(method) and r(kernel).
program define _cqreg_denparse, rclass
    local spec = lower(strtrim(`"`0'"'))
    local 0 `", `spec'"'
    capture syntax [, FITted RESidual residuals KERnel(string) *]
    local bad = _rc
    if !`bad' {
        if "`residuals'" != "" local residual residual
        local method `fitted' `residual'
        if `"`kernel'"' != "" local method `method' kernel
        else if `"`options'"' != "" {
            local 0 `", `options'"'
            capture syntax [, Kernel]
            if _rc local bad = 1
            else if "`kernel'" != "" local method `method' kernel
            local kernel
        }
        if `: word count `method'' != 1 local bad = 1
    }
    if !`bad' & "`method'" == "kernel" {
        capture _cqreg_kernel, `kernel'
        if _rc local bad = 1
        else return local kernel `r(kernel)'
    }
    if `bad' {
        di as error `"cqreg: unrecognized density method: `spec'"'
        di as error "       Valid options: fitted (default), residual, kernel[(kernel)]"
        exit 198
    }
    return local method `method'
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
* ctools: drop the scalars and matrices a command and its plugin leave behind
* _ctools_cleanup <command> [<command> ...] drops every scalar and matrix
* named __<command>_* or _<command>_*.
program define _ctools_cleanup
    version 14.1
    foreach command of local 0 {
        foreach pattern in __`command'_* _`command'_* {
            local found : all scalars "`pattern'"
            if "`found'" != "" capture scalar drop `found'
            local found : all matrices "`pattern'"
            if "`found'" != "" capture matrix drop `found'
        }
    }
end

* Private monotonic timestamp in seconds; no public Stata timer IDs.
program define _ctools_clock
    version 14.1
    args target
    tempname held
    _return hold `held'
    capture noisily {
        * Once registered, this plugin stays loaded for the Stata session.
        * Resolve its path and check its protocol only on first use (or
        * after discard), rather than searching the ado-path at every lap.
        capture plugin call ctools_clock_plugin, "clock"
        if _rc {
            _ctools_load
            capture program ctools_clock_plugin, plugin using("`__ctools_plugin'")
            if _rc != 0 & _rc != 110 exit _rc
            plugin call ctools_clock_plugin, "clock"
        }
    }
    local rc = _rc
    _return restore `held'
    if `rc' exit `rc'
    c_local `target' `ctools_clock_value'
end

* ctools: e(b) and e(V) from a plugin's compact results
* _ctools_ebv <b> <V>, names() vmatrix() [base() stub() kkeep() bmatrix()
*     collinear() cons() kept() colnames()]
*
* names() lists every term of the model in e(b) order (none for a model with
* only the constant) and base() flags the
* base/omitted placeholders among them (_ctools_fvsplit's r(base); default:
* none). The other terms are the plugin's columns, in order. Two layouts of
* the plugin's results are understood:
*
*   stub(<s>) kkeep(#)   the plugin dropped its collinear columns: scalars
*                        <s>_collinear_<k> flag each column, <s>_beta_<i> hold
*                        the kept coefficients and vmatrix() their covariance
*                        in its top-left block. With kkeep(0) every column
*                        counts as dropped.
*   bmatrix(<m>)         the plugin kept a (zero) row and column for every
*                        column: bmatrix() and vmatrix() cover them all and
*                        collinear() lists a 0/1 flag per column.
*
* cons(<scalar>) appends _cons with that coefficient; vmatrix() then carries
* the constant in the row and column after the coefficients.
*
* <b> and <V> are created with zeros for placeholders and dropped columns and
* named: a placeholder keeps its name, a dropped column gets the o. prefix.
* kept() names a matrix for the kept coefficients alone (stub layout, not
* created when there are none) and colnames() a caller local for the names.
program define _ctools_ebv
    version 14.1
    syntax namelist(min=2 max=2), Vmatrix(name) ///
        [NAMES(string) BASE(string) STUB(string) KKeep(integer -1) Bmatrix(name) ///
         COLLinear(string) CONS(name) KEPT(name) COLNames(name)]
    gettoken b V : namelist
    if ("`stub'" == "") == ("`bmatrix'" == "") | ("`stub'" != "" & `kkeep' < 0) {
        di as error "_ctools_ebv: specify stub() with kkeep(), or bmatrix()"
        exit 198
    }
    mata: _ctools_ebv_build("`b'", "`V'")
    matrix colnames `b' = `outnames'
    matrix colnames `V' = `outnames'
    matrix rownames `V' = `outnames'
    if "`colnames'" != "" c_local `colnames' `outnames'
end

capture mata: mata drop _ctools_ebv_build()
capture mata: mata drop _ctools_ebv_scalar()
capture mata: mata drop _ctools_ebv_omit()
mata:
/* The omitted form of a term: o.x, and for an interaction the o operator on
 * every part (co.x#co.x, 1o.g#co.x), which is how Stata names it; an o.
 * prefix on an interaction is not recognized as omitted. A base part keeps
 * its name (3o.race#0b.union, as regress names it): Stata has no "bo"
 * operator, and matrix colnames rejects 0bo.union. */
string scalar _ctools_ebv_omit(string scalar name)
{
    string rowvector parts
    real scalar i, dot

    if (!strpos(name, "#")) return("o." + name)
    parts = tokens(subinstr(name, "#", " "))
    for (i = 1; i <= cols(parts); i++) {
        dot = strpos(parts[i], ".")
        if (dot & (strpos(substr(parts[i], 1, dot - 1), "b") | strpos(substr(parts[i], 1, dot - 1), "o"))) continue
        if (dot) parts[i] = substr(parts[i], 1, dot - 1) + "o" + substr(parts[i], dot, .)
        else parts[i] = "o." + parts[i]
    }
    return(invtokens(parts, "#"))
}

real scalar _ctools_ebv_scalar(string scalar name)
{
    real matrix value
    value = st_numscalar(name)
    if (rows(value) != 1 | cols(value) != 1) {
        errprintf("ctools: plugin result %s not found\n", name)
        exit(498)
    }
    return(value)
}

void _ctools_ebv_build(string scalar bname, string scalar Vname)
{
    string rowvector names
    string scalar stub
    real rowvector base, sent, coll, kept, pos, bsrc, b
    real matrix Vsrc, V
    real scalar T, Ks, Kk, m, n, i, hascons

    names = tokens(st_local("names"))
    T = cols(names)
    base = (st_local("base") == "" ? J(1, T, 0) : strtoreal(tokens(st_local("base"))))
    if (cols(base) != T | missing(base)) {
        errprintf("ctools: %g terms but %g base flags\n", T, cols(base))
        exit(498)
    }
    Ks = T - sum(base :!= 0)
    sent = (Ks ? selectindex(base :== 0) : J(1, 0, .))

    stub = st_local("stub")
    if (stub == "") {
        /* a row and column for every plugin column */
        coll = (st_local("collinear") == "" ? J(1, Ks, 0) : strtoreal(tokens(st_local("collinear"))))
        if (cols(coll) != Ks | missing(coll)) {
            errprintf("ctools: %g plugin columns but %g collinearity flags\n", Ks, cols(coll))
            exit(498)
        }
        kept = (Ks ? 1..Ks : J(1, 0, .))
        bsrc = st_matrix(st_local("bmatrix"))
    }
    else {
        /* collinear columns dropped by the plugin */
        coll = J(1, Ks, 1)
        if (strtoreal(st_local("kkeep")) > 0) {
            for (i = 1; i <= Ks; i++) {
                coll[i] = (_ctools_ebv_scalar(stub + "_collinear_" + strofreal(i)) != 0)
            }
        }
        Kk = Ks - sum(coll)
        kept = (Kk ? selectindex(coll :== 0) : J(1, 0, .))
        bsrc = J(1, Kk, 0)
        for (i = 1; i <= Kk; i++) {
            bsrc[i] = _ctools_ebv_scalar(stub + "_beta_" + strofreal(i))
        }
        if (st_local("kept") != "" & Kk) st_matrix(st_local("kept"), bsrc)
    }

    hascons = (st_local("cons") != "")
    n = T + hascons
    pos = (cols(kept) ? sent[kept] : J(1, 0, .))
    if (hascons) {
        pos = pos, n
        bsrc = bsrc, _ctools_ebv_scalar(st_local("cons"))
    }
    m = cols(pos)
    b = J(1, n, 0)
    V = J(n, n, 0)
    if (m) {
        Vsrc = st_matrix(st_local("vmatrix"))
        if (cols(bsrc) != m | rows(Vsrc) < m | cols(Vsrc) < m) {
            errprintf("ctools: plugin results do not match the %g estimated coefficients\n", m)
            exit(498)
        }
        b[pos] = bsrc
        V[pos, pos] = Vsrc[|1, 1 \ m, m|]
    }
    st_matrix(bname, b)
    st_matrix(Vname, V)

    for (i = 1; i <= Ks; i++) {
        if (coll[i]) names[sent[i]] = _ctools_ebv_omit(names[sent[i]])
    }
    if (hascons) names = names, "_cons"
    st_local("outnames", invtokens(names))
}
end

* ctools: estimated terms and base/omitted placeholders of an fvexpand-ed list
* names() is an fvexpand result. A term is a placeholder, not estimated, only
* when every #-component carries a base (b.) or omitted (o.) marker:
* 1b.a#2.b and 0b.foreign#c.weight are estimable, 1b.a#1b.b and 1b.a#co.x
* are not. vars(), when given, is the fvrevar list for the same terms.
* Returns r(base) (one 0/1 flag per term, 1 = placeholder), r(names) and
* r(vars) (the estimated terms and their variables).
program define _ctools_fvsplit, rclass
    version 14.1
    syntax [, NAMES(string) VARS(string)]
    local has_vars = (`"`vars'"' != "")
    local flags
    local kept_names
    local kept_vars
    foreach name of local names {
        if `has_vars' gettoken var vars : vars
        local base = 1
        local parts : subinstr local name "#" " ", all
        foreach part of local parts {
            if !regexm("`part'", "b\.") & !regexm("`part'", "o\.") local base = 0
        }
        local flags `flags' `base'
        if !`base' {
            local kept_names `kept_names' `name'
            if `has_vars' local kept_vars `kept_vars' `var'
        }
    }
    return local base `flags'
    return local names `kept_names'
    return local vars `kept_vars'
end

* _ctools_group — group identifiers from the ctools plugin, in place of egen group()
* 
* Usage:  _ctools_group varlist [if] [in], generate(newvar) [label]
* _ctools_group varlist [if] [in], count [strict]
* strict returns plugin errors instead of sorting in Stata.
* 
* Numbers the distinct values of varlist 1, 2, ... in sort order, as
* egen group() does: observations excluded by if/in, or with a missing value
* in any variable of varlist (string "" counts as missing), get missing. The
* new variable is a long labeled "group(varlist)". label also defines a value
* label named after the new variable, with egen group(, label)'s text for
* each group. count only counts the groups. Returns r(N_groups), and with
* label r(reps), one observation number of each group in group order.
* 
* The plugin reads, sorts and numbers the keys. Stata numbers them by sorting
* only when the plugin is unavailable or cannot read a key (strL text over
* 2,045 bytes).

program define _ctools_group, rclass
    version 14.1
    syntax varlist [if] [in], [GENerate(name) LABel COUNT STRICT]
    if ("`generate'" == "") == ("`count'" == "") | ("`label'" != "" & "`generate'" == "") {
        di as error "_ctools_group: specify generate() [label] or count"
        exit 198
    }
    if "`generate'" != "" {
        confirm new variable `generate'
        if "`label'" != "" {
            capture label list `generate'
            if !_rc {
                di as error "default label `generate' already defined"
                exit 110
            }
        }
    }
    local nkeys : word count `varlist'

    local rc = 199
    capture _ctools_load
    if !_rc {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if _rc == 0 | _rc == 110 local rc = 0
    }
    capture confirm number 0
    if !`rc' {
        _ctools_strw `varlist'
        local out ""
        if "`generate'" != "" {
            * The plugin writes every observation of the range, so without in
            * the variable is created without initial values.
            if `"`in'"' == "" mata: (void) st_addvar("long", st_local("generate"), 1)
            else quietly generate long `generate' = .
            local out "out"
        }
        local reps = cond("`label'" != "", "reps", "")
        * Stata answers an if one observation at a time. Without one the
        * plugin skips the check; an if that names a numeric variable (the
        * callers' touse) passes that variable last, as a sample marker the
        * plugin reads in parallel. Any other if is marked first, unless in
        * limits the range (marking covers every observation); that marker is
        * the last variable until it is dropped, which is then free.
        local select "skipif"
        local marker ""
        local marked ""
        local ifcond ""
        if `"`if'"' != "" {
            local select "touse"
            local marker = strtrim(substr(`"`if'"', 3, .))
            capture confirm numeric variable `marker', exact
            if _rc | wordcount(`"`marker'"') != 1 {
                local marker ""
                if `"`in'"' == "" {
                    tempvar marker
                    capture mark `marker' `if'
                    local rc = _rc
                    local marked `marker'
                }
                else {
                    local select ""
                    local ifcond `"`if'"'
                }
            }
        }
        local ctools_group_reps ""
        if !`rc' {
            capture plugin call ctools_plugin `varlist' `generate' `marker' `ifcond' `in', ///
                "group nkeys=`nkeys' `out' `reps' `select'"
            local rc = _rc
        }
        if "`marked'" != "" capture drop `marked'
        * Every failure, Break included, drops the unfilled output variable.
        if `rc' & "`generate'" != "" capture drop `generate'
        if `rc' == 1 exit 1
        if `rc' == 134 {
            di as error "too many values for {bf:label define}"
            exit 134
        }
        if !`rc' {
            local n_groups = scalar(_ctools_group_n)
            local group_reps `"`ctools_group_reps'"'
        }
        capture scalar drop _ctools_group_n
    }
    if `rc' & "`strict'" != "" exit `rc'
    if `rc' {
        _ctools_group_sort `varlist' `if' `in', generate(`generate') `count' `label'
        local n_groups = r(N_groups)
        local group_reps `"`r(reps)'"'
    }

    if "`generate'" != "" {
        if udstrlen("group(`varlist')") > 80 {
            note `generate' : group(`varlist')
            label variable `generate' "see notes"
        }
        else {
            label variable `generate' "group(`varlist')"
        }
        if "`label'" != "" {
            if `n_groups' > 0 {
                mata: _ctools_group_label(st_local("varlist"), st_local("group_reps"), ///
                    st_local("generate"))
            }
            capture label list `generate'
            if !_rc label values `generate' `generate'
        }
    }
    return scalar N_groups = `n_groups'
    * label: the observation number of one observation of each group
    if "`label'" != "" return local reps `"`group_reps'"'
end

* Number groups by sorting in Stata (the data and their sort order are restored).
program define _ctools_group_sort, rclass sortpreserve
    syntax varlist [if] [in], [GENerate(name) COUNT LABel]
    marksample touse, strok
    tempvar order g
    quietly generate long `order' = _n
    sort `touse' `varlist' `order'
    quietly by `touse' `varlist': generate long `g' = (_n == 1) if `touse'
    quietly replace `g' = sum(`g') if `touse'
    quietly summarize `g', meanonly
    return scalar N_groups = cond(r(N), r(max), 0)
    if "`label'" != "" {
        mata: _ctools_group_firsts(st_local("g"), st_local("order"))
        return local reps `"`__reps'"'
    }
    if "`generate'" != "" rename `g' `generate'
end

mata:
/* Original observation number of the first observation of each group. */
void _ctools_group_firsts(string scalar g, string scalar order)
{
    real colvector id, obs, previous, first

    id = st_data(., g)
    obs = st_data(., order)
    st_local("__reps", "")
    if (rows(id) == 0) return
    previous = (rows(id) > 1 ? (. \ id[|1 \ rows(id) - 1|]) : .)
    first = selectindex(id :< . :& (id :!= previous))
    if (rows(first)) st_local("__reps", invtokens(strofreal(obs[first], "%12.0g")'))
}

/*
    Value label of egen group(varlist), label for groups 1..G, from one
    observation of each group (reps, in group order). Each variable
    contributes, separated by spaces: a string's value; a numeric value's
    label text when its value label defines one, otherwise the value
    ("%18.0g" for integers, "%9.0g" or "%9,0g" otherwise). Text is cut to 32
    characters and trimmed, and a string's NUL bytes become \0, which is what
    egen group() does when called under version 16 or earlier, as the ctools
    ado files are.
*/
void _ctools_group_label(string scalar varlist, string scalar reps, string scalar lname)
{
    real rowvector v
    real colvector obs, s
    real scalar i, j
    string scalar fmt, vl
    string colvector x, y

    v = st_varindex(tokens(varlist))
    obs = strtoreal(tokens(reps))'
    if (rows(obs) > 65536) {
        errprintf("too many values for {bf:label define}\n")
        exit(134)
    }
    fmt = (st_global("c(dp)") == "comma" ? "%9,0g" : "%9.0g")
    for (j = 1; j <= cols(v); j++) {
        if (st_isnumvar(v[j])) {
            s = st_data(obs, v[j])
            x = J(rows(s), 1, "")
            vl = st_varvaluelabel(v[j])
            if (vl != "") {
                if (st_vlexists(vl)) x = ustrtrim(usubstr(st_vlmap(vl, s), 1, 32))
            }
            for (i = 1; i <= rows(s); i++) {
                if (x[i] == "") x[i] = strofreal(s[i], (s[i] != trunc(s[i]) ? fmt : "%18.0g"))
            }
        }
        else {
            x = st_sdata(obs, v[j])
            x = usubinstr(x, char(0), "\0", .)
            x = ustrtrim(usubstr(x, 1, 32))
        }
        y = (j == 1 ? x : y :+ " " :+ x)
    }
    st_vlmodify(lname, 1::rows(y), y)
}
end

* ctools: seconds since the previous lap of a caller-local wall clock
* _ctools_lap <target> <clock> <clock value> stores the elapsed seconds in the
* caller's local <target> and restarts the clock, a local the caller starts
* with _ctools_clock <clock> and passes back expanded:
*     _ctools_clock __ct_lap
*     ...
*     _ctools_lap __time_setup __ct_lap `__ct_lap'
program define _ctools_lap
    version 14.1
    args target clock start
    quietly _ctools_clock now
    c_local `target' = `now' - `start'
    c_local `clock' `now'
end

* Shared plugin loader and identity contract
program define _ctools_load
    version 14.1
    * Not rclass, and findfile's r() is undone, so the caller's r() survives;
    * identity comes back through c_local.
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
        if `__is_arm' == 0 & !strpos(lower("`__machine'"), "intel") & !strpos(lower("`__machine'"), "x86") {
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

    * The ado/plugin argument protocol, stamped from config/release.json.
    local __ado_protocol "4"
    local __reinstall `"net install ctools, from("https://raw.githubusercontent.com/mdroste/stata-ctools/main/build") replace"'

    tempname __rhold
    _return hold `__rhold'
    capture findfile `__plugin'
    local __found = (_rc == 0)
    local __path `"`r(fn)'"'
    _return restore `__rhold'
    local __wanted "`__plugin'"
    if !`__found' & "`__plugin'" != "ctools.plugin" local __plugin "ctools.plugin"
    capture program ctools_plugin, plugin using("`__plugin'")
    local __load_rc = _rc
    if `__load_rc' != 0 & `__load_rc' != 110 {
        di as error "ctools: could not load the plugin `__plugin'"
        di as error `"  Stata reports c(os) "`__os'" and c(machine_type) "`__machine'"; ctools selected `__wanted'."'
        if `__found' di as error `"  The plugin file is `__path'."'
        else di as error "  `__wanted' was not found on the ado-path; reinstall ctools."
        di as error "  Stata's plugin loader reported:"
        capture noisily program ctools_plugin, plugin using("`__plugin'")
        if "`__os'" == "Unix" & !`__is_mac' {
            di as error "  The distribution Linux plugin needs x86-64 and glibc 2.17 or later; OpenMP is linked statically."
            di as error "  Run ldd -r on the plugin to inspect missing libraries or symbols."
        }
        else if `__is_mac' {
            di as error "  On macOS, check that the plugin is the notarized ctools release and is not quarantined."
        }
        di as error "  See https://github.com/mdroste/stata-ctools/blob/main/docs/INSTALLATION.md#troubleshooting"
        exit `__load_rc'
    }
    local ctools_plugin_hdfe_screen
    local ctools_plugin_hdfe_markout
    local ctools_plugin_hdfe_fvlevels
    capture plugin call ctools_plugin, "version"
    local __id_rc = _rc
    if `__id_rc' | "`ctools_plugin_protocol'" != "`__ado_protocol'" {
        local __has = cond("`ctools_plugin_protocol'" == "", "none", "`ctools_plugin_protocol'")
        di as error "ctools: the plugin `__plugin' does not match the installed ado files."
        local __ident = cond("`ctools_plugin_version'" == "", "", " (plugin version `ctools_plugin_version', revision `ctools_plugin_revision')")
        di as error "  The ado files need plugin protocol `__ado_protocol'; the loaded plugin reports `__has'`__ident'."
        di as error "  Reinstall ctools so its ado files and plugins come from the same release:"
        di as error `"    `__reinstall'"'
        di as error "  Then restart Stata, because a plugin stays loaded until Stata exits."
        exit 9
    }
    c_local __ctools_plugin "`__plugin'"
    c_local __ctools_plugin_version "`ctools_plugin_version'"
    c_local __ctools_plugin_revision "`ctools_plugin_revision'"
    c_local __ctools_plugin_protocol "`ctools_plugin_protocol'"
    c_local __ctools_plugin_hdfe_screen = ("`ctools_plugin_hdfe_screen'" == "1")
    c_local __ctools_plugin_hdfe_markout = ("`ctools_plugin_hdfe_markout'" == "1")
    c_local __ctools_plugin_hdfe_fvlevels = ("`ctools_plugin_hdfe_fvlevels'" == "1")
    capture confirm number 0
end

* _ctools_strw — set string width metadata for ctools flat buffer optimization
* 
* Usage:  _ctools_strw varlist
* 
* Sets local __ctools_strw in CALLER's scope containing comma-separated
* string variable widths (0 = numeric, >0 = str width e.g. 17 for str17).
* The ctools C plugin reads this local automatically via SF_macro_use
* to enable flat buffer string I/O (eliminates per-string strlen, memcpy,
* and arena atomic CAS overhead).
* 
* Example:
* _ctools_strw `allvars'
* plugin call ctools_plugin `allvars', "mycommand args"
* * C plugin automatically detects __ctools_strw local

program define _ctools_strw
    syntax varlist

    local __strw ""
    local __sep ""
    foreach v of local varlist {
        local __type : type `v'
        if substr("`__type'", 1, 3) == "str" {
            local __w = real(substr("`__type'", 4, .))
            if missing(`__w') local __w = 0
        }
        else {
            local __w = 0
        }
        local __strw "`__strw'`__sep'`__w'"
        local __sep ","
    }
    c_local __ctools_strw `__strw'
end
* ctools-helpers:end
