*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! Private helper for cppmlhdfe.
* The suboptions of absorb(absvars, ...), as ppmlhdfe hands them to reghdfe:
* reghdfe's options (savefe, keepsingletons, the partialling tolerance(),
* iterate() or maxiterations(), dofadjustments(), acceleration(),
* transform(), poolsize(), verbose() of the HDFE object, prune, condition,
* rre(), algorithm(), slopemethod()) and, as ppmlhdfe, any of its GLM
* properties (min_ok(), relu_tol()...), by their full names. Sets the
* caller's locals of the same meaning.
program define _cppmlhdfe_abs_options
    syntax [, SAVEfe GENerate KEEPSINgletons TOLerance(string) ITERate(string) MAXITerations(string) ///
        DOFadjustments(string) ACCELeration(string) TRANSform(string) POOLsize(string) ///
        Verbose(numlist min=1 max=1 >=-1 <=5 integer) PRUNE CONDition RRE(varname numeric) ///
        ALGorithm(string) SLOPEmethod(string) *]
    if "`acceleration'"!="" c_local acceleration `acceleration'
    if "`transform'"!="" c_local transform `transform'
    if "`poolsize'"!="" {
        confirm integer number `poolsize'
        if `poolsize'<0 | missing(`poolsize') error 198
        c_local poolsize `poolsize'
    }
    c_local save_all_fe = ("`savefe'`generate'"!="")
    if "`keepsingletons'"!="" c_local keepsingletons keepsingletons
    if "`tolerance'"!="" {
        confirm number `tolerance'
        if `tolerance'<=0 | missing(`tolerance') error 198
        c_local itolerance `tolerance'
    }
    * reghdfe's maxiterations(): its partialling iteration limit, iterate()
    if "`maxiterations'"!="" {
        confirm number `maxiterations'
        local maxiterations = int(`maxiterations')
        if `maxiterations' > 0 local iterate `maxiterations'
    }
    if "`iterate'"!="" {
        confirm integer number `iterate'
        if `iterate'<1 | missing(`iterate') error 198
        c_local iterate `iterate'
    }
    if "`dofadjustments'"!="" c_local dofadjustments `dofadjustments'
    if "`verbose'"!="" c_local hdfe_verbose `verbose'
    foreach option in prune condition {
        if "``option''" != "" c_local `option' `option'
    }
    if "`rre'" != "" c_local rre `rre'
    if `"`options'"' != "" {
        _cppmlhdfe_glm_props, `options'
        * (set here; passed on to cppmlhdfe)
        foreach option in standardize_data min_ok use_exact_solver use_exact_partial use_heuristic_tol ///
            start_inner_tol target_inner_tol remove_collinear_variables use_step_halving ///
            step_halving_memory max_step_halving mu_tol relu_tol relu_zero_tol relu_maxiter ///
            relu_strict relu_accelerate relu_debug relu_report_r2 simplex_tol simplex_maxiter log {
            if `"``option''"' != "" c_local `option' ``option''
        }
    }
end

* ppmlhdfe's GLM properties among the absorb() suboptions, by their exact
* names; any other option is not allowed there, as in ppmlhdfe
program define _cppmlhdfe_glm_props
    syntax [, standardize_data(string) min_ok(string) use_exact_solver(string) ///
        use_exact_partial(string) use_heuristic_tol(string) start_inner_tol(string) ///
        target_inner_tol(string) remove_collinear_variables(string) use_step_halving(string) ///
        step_halving_memory(string) max_step_halving(string) mu_tol(string) relu_tol(string) ///
        relu_zero_tol(string) relu_maxiter(string) relu_strict(string) relu_accelerate(string) ///
        relu_debug(string) relu_report_r2(string) simplex_tol(string) simplex_maxiter(string) ///
        log(string) *]
    if `"`options'"' != "" {
        gettoken name : options, parse(" (")
        di as error "option {bf:`name'} not allowed"
        exit 198
    }
    foreach option in standardize_data min_ok use_exact_solver use_exact_partial use_heuristic_tol ///
        start_inner_tol target_inner_tol remove_collinear_variables use_step_halving ///
        step_halving_memory max_step_halving mu_tol relu_tol relu_zero_tol relu_maxiter ///
        relu_strict relu_accelerate relu_debug relu_report_r2 simplex_tol simplex_maxiter log {
        if `"``option''"' == "" continue
        capture confirm number ``option''
        if _rc {
            di as error "option {bf:`option'} not allowed"
            exit 198
        }
        if "`option'" == "log" {
            if !inlist(`log', 0, 1) error 198
            if !`log' c_local log nolog
        }
        else c_local `option' ``option''
    }
end
