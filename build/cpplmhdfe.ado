*! version 1.1.0 25sep2026 github.com/mdroste/stata-ctools

program define cpplmhdfe, eclass
    version 14.1

    if replay() {
        if "`e(cmd)'" != "cpplmhdfe" error 301
        ereturn display `0'
        exit
    }
    * Accept the legacy bare verbose option as well as ppmlhdfe's verbose(#).
    local 0 = ustrregexra(`"`0'"', "(?<=[ ,])verbose(?=[ ,]|$)", "verbose(1)")

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start wall clock timer
    local __do_timing = 0
    quietly _ctools_clock __ct_lap98

    * Drop whatever an interrupted run left behind
    _ctools_cleanup cpplmhdfe

    syntax varlist(min=1 fv ts) [aw fw pw] [if] [in] [, Absorb(string) VCE(string) Verbose(integer 0) ///
        TOLerance(real 1e-8) ITERATE(integer 10000) THReads(integer 0) ///
        EXPosure(varname) OFFset(varname) ///
        IRLSTOLerance(real -1) IRLSMAXiter(integer 10000) ///
        SEPTOLerance(real 1e-8) ///
        DOFadjustments(string) D(name) D2 NOAbsorb KEEPSINgletons ///
        CLuster(string) MAXITerations(integer 10000) ITOLerance(real -1) ///
        SEParation(string) EForm IRr noLog TIMEit ///
        GUESS(string) Standardize_data(integer 1) Min_ok(integer 1) ///
        Use_exact_solver(integer 0) Use_exact_partial(integer -1) Use_heuristic_tol(integer 1) ///
        Start_inner_tol(real 1e-4) Remove_collinear_variables(integer 1) ///
        Use_step_halving(integer 0) Step_halving_memory(real .9) Max_step_halving(integer 2) ///
        Mu_tol(real 1e-6) Relu_tol(real 1e-4) Relu_zero_tol(real 1e-8) ///
        Relu_maxiter(integer 100) Relu_strict(integer 0) Relu_accelerate(integer 0) ///
        Simplex_tol(real 1e-12) Simplex_maxiter(integer 1000) TAGSEP(name) ZVARname(name) ///
        ACCELeration(string) TRANSform(string) POOLsize(integer 0) BTOL(real 1e-8) CONLIM(real 1e8) *]

    local save_all_fe = 0
    local comma = strpos("`absorb'",",")
    if `comma' {
        local absorb_options = substr("`absorb'",`comma'+1,.)
        local absorb = substr("`absorb'",1,`comma'-1)
        _cpplmhdfe_abs_options, `absorb_options'
    }

    local __do_timing = (`verbose' > 0 | "`timeit'" != "")
    _get_diopts diopts options, `options'
    if "`options'" != "" {
        di as error "unsupported option(s): `options'"
        exit 198
    }
    * Exponentiated coefficients: ereturn display takes the column title
    local eform_opt
    if "`irr'" != "" local eform_opt eform(IRR)
    else if "`eform'" != "" local eform_opt eform(exp(b))
    if "`cluster'" != "" {
        if "`vce'" != "" error 198
        local vce cluster `cluster'
    }
    if `maxiterations' != 10000 local irlsmaxiter `maxiterations'
    if inlist("`separation'", "", "def", "default", "standard", "auto", "on") local separation fe simplex relu
    if inlist("`separation'", "all", "full") local separation fe simplex relu mu
    if inlist("`separation'", "no", "off", "none") local separation
    local separation : subinstr local separation "ir" "relu", word all
    local allowed fe simplex relu mu
    local invalid : list separation - allowed
    if "`invalid'" != "" error 198
    local fe_option fe
    local do_fe : list fe_option in separation
    local relu_option relu
    local simplex_option simplex
    local do_relu : list relu_option in separation
    local do_simplex : list simplex_option in separation
    local mu_option mu
    local do_mu : list mu_option in separation
    * As in ppmlhdfe, tagsep requests a ReLU diagnostic, not estimation.
    if "`tagsep'"!="" {
        local do_fe=0
        local do_simplex=0
        local do_mu=0
        local do_relu=1
    }

    * Exact partialling is the native default, but the reference defaults to
    * inexact partialling and its stopping schedule depends on that choice.
    * Record whether exact partialling was explicitly requested.
    local ref_exact_partial = (`use_exact_partial' == 1)
    if `use_exact_partial' == -1 local use_exact_partial 1
    foreach option in standardize_data use_exact_solver use_exact_partial use_heuristic_tol remove_collinear_variables use_step_halving relu_strict relu_accelerate {
        if !inlist(``option'',0,1) error 198
    }
    if missing(`min_ok',`relu_maxiter',`simplex_maxiter',`max_step_halving') | ///
        `min_ok'<1 | `relu_maxiter'<1 | `simplex_maxiter'<1 | `max_step_halving'<0 error 198
    foreach option in start_inner_tol mu_tol relu_tol relu_zero_tol simplex_tol {
        if missing(``option'') | ``option''<=0 | ``option''>=1 error 198
    }
    if missing(`step_halving_memory') | `step_halving_memory'<=0 | `step_halving_memory'>=1 error 198
    if inlist("`guess'", "", "default") local guess simple
    gettoken guess_method guess_var : guess
    if "`guess_method'"=="var" local guess_method variable
    if !inlist("`guess_method'","simple","ols","variable") error 198
    local guess_code = cond("`guess_method'"=="simple",0,cond("`guess_method'"=="ols",1,2))
    if `guess_code'==2 confirm numeric variable `guess_var'
    else if "`guess_var'"!="" error 198
    if "`tagsep'`zvarname'" != "" _ctools_newvars `tagsep' `zvarname'

    if "`acceleration'"=="" local acceleration cg
    local accel_code = .
    if inlist("`acceleration'","cg","conjugate_gradient") local accel_code=0
    if inlist("`acceleration'","none","no") local accel_code=1
    if inlist("`acceleration'","sd","steepest_descent") local accel_code=2
    if inlist("`acceleration'","a","aitken") local accel_code=3
    if "`acceleration'"=="hybrid" local accel_code=4
    if "`acceleration'"=="lsmr" local accel_code=5
    if missing(`accel_code') error 198
    if "`transform'"=="" local transform symmetric_kaczmarz
    local transform_code=.
    if inlist("`transform'","sym","symmetric_kaczmarz","symmetric") local transform_code=0
    if inlist("`transform'","kac","kaczmarz") local transform_code=1
    if inlist("`transform'","cim","cimmino") local transform_code=2
    if missing(`transform_code',`btol',`conlim',`poolsize') | `poolsize'<0 | `btol'<=0 | `conlim'<=0 error 198

    * Match ppmlhdfe defaults: if IRLS tolerance is not specified, use tolerance()
    if `irlstolerance' < 0 {
        local irlstolerance `tolerance'
    }

    if `iterate' < 1 | `irlsmaxiter' < 1 | `threads'<0 | ///
        missing(`tolerance', `irlstolerance', `septolerance', `iterate', `irlsmaxiter', `threads', `itolerance') | ///
        `tolerance' <= 0 | `irlstolerance' <= 0 | `septolerance' <= 0 {
        di as error "cpplmhdfe: iteration limits and tolerances must be positive"
        exit 198
    }

    * An intercept-only factor implements the no-absorb model in C.
    if "`absorb'" != "" & "`noabsorb'" != "" error 198
    local has_absorb = ("`absorb'" != "")
    if !`has_absorb' {
        tempvar constant_fe
        quietly generate byte `constant_fe' = 1
        local absorb `constant_fe'
        local absorb_orig _cons
    }
    else local absorb_orig "`absorb'"
    if "`d'" != "" & "`d2'" != "" error 198
    if "`d2'" != "" {
        local d _ppmlhdfe_d
        capture drop `d'
    }
    if "`d'" != "" _ctools_newvars `d'

    * Parse independent absorbed-DF adjustments, matching reference semantics.
    local dofmethods = lower(trim("`dofadjustments'"))
    if inlist("`dofmethods'", "", "all") local dofmethods pairwise clusters continuous
    if "`dofmethods'" == "none" local dofmethods
    local allowed_dof firstpair pairwise clusters continuous
    local invalid_dof : list dofmethods - allowed_dof
    if "`invalid_dof'" != "" {
        di as error "invalid dofadjustments(): `invalid_dof'"
        exit 198
    }
    local dof_flags = 0
    local bit = 1
    foreach adjustment in pairwise firstpair clusters continuous {
        local selected : list adjustment in dofmethods
        if `selected' local dof_flags = `dof_flags' + `bit'
        local bit = 2*`bit'
    }
    if mod(`dof_flags',4)==3 error 198
    local dof_adjust_type = cond(`dof_flags'==0,1,cond(mod(`dof_flags',4)==2,2,0))

    * Mark sample
    marksample touse
    * Parse categorical interactions and heterogeneous slopes without Mata
    * estimation or a dependency on reghdfe/ppmlhdfe.
    local absorb_numeric
    local slope_vars
    local expanded_absorb
    local pending `absorb'
    local original_group = 0
    local expanded_index = 0
    while `"`pending'"' != "" {
        gettoken term pending : pending, bind
        local ++original_group
        local target
        local automatic = 0
        local equal = strpos("`term'","=")
        if `equal' {
            local target = substr("`term'",1,`equal'-1)
            local term = substr("`term'",`equal'+1,.)
            confirm name `target'
        }
        else if `save_all_fe' {
            local target __hdfe`original_group'__
            local automatic = 1
        }
        local paren = strpos("`term'", "c.(")
        local slope_number = 0
        if `paren' {
            local prefix = substr("`term'",1,`paren'-1)
            local inside = substr("`term'",`paren'+3,strlen("`term'")-`paren'-3)
            foreach continuous of local inside {
                local ++expanded_index
                local ++slope_number
                local expanded_absorb `expanded_absorb' `prefix'c.`continuous'
                local target_`expanded_index' `target'
                local target_slope_`expanded_index' `slope_number'
                local automatic_`expanded_index' `automatic'
            }
        }
        else {
            local expanded `term'
            if !strpos("`term'", "#") & !strpos("`term'", ".") unab expanded : `term'
            local item_number = 0
            foreach item of local expanded {
                local ++item_number
                if `item_number'>1 {
                    local ++original_group
                    if `automatic' local target __hdfe`original_group'__
                }
                local ++expanded_index
                local expanded_absorb `expanded_absorb' `item'
                local target_`expanded_index' `target'
                local target_slope_`expanded_index' 1
                local automatic_`expanded_index' `automatic'
            }
        }
    }
    local group_count = 0
    local intercept_keys
    local all_fe_keys
    local expanded_index = 0
    local saved_fe_names
    foreach term of local expanded_absorb {
        local ++expanded_index
        local target `target_`expanded_index''
        local target_slope `target_slope_`expanded_index''
        local automatic `automatic_`expanded_index''
        local include_intercept = strpos("`term'", "##") > 0
        local parts : subinstr local term "#" " ", all
        local categories
        local continuous
        foreach part of local parts {
            if substr("`part'",1,2)=="c." {
                local cv = substr("`part'",3,.)
                confirm numeric variable `cv'
                local continuous `continuous' `cv'
            }
            else {
                local cat : subinstr local part "i." "", all
                confirm variable `cat'
                local categories `categories' `cat'
            }
        }
        if "`categories'" == "" error 198
        markout `touse' `categories' `continuous', strok
        tempvar fe_id
        quietly _ctools_group `categories' if `touse', generate(`fe_id')
        local ncontinuous : word count `continuous'
        if `ncontinuous' > 1 {
            di as error "multiple continuous factors in an absorb interaction"
            exit 198
        }
        local key : subinstr local categories " " "#", all
        local all_fe_keys : list all_fe_keys | key
        local seen : list key in intercept_keys
        if (`include_intercept' | !`ncontinuous') & !`seen' {
            local absorb_numeric `absorb_numeric' `fe_id'
            local ++group_count
            local slope_var_`group_count'
            local fe_name_`group_count' `key'
            local fe_target_`group_count' `target'
            local fe_automatic_`group_count' `automatic'
            local intercept_keys `intercept_keys' `key'
        }
        if `ncontinuous' {
            local absorb_numeric `absorb_numeric' `fe_id'
            local ++group_count
            local slope_var_`group_count' `continuous'
            local fe_name_`group_count' `key'#c.`continuous'
            if "`target'"!="" local fe_target_`group_count' `target'Slope`target_slope'
            local fe_automatic_`group_count' `automatic'
        }
    }
    forvalues g=1/`group_count' {
        local name `fe_target_`g''
        if "`name'"!="" {
            if !`fe_automatic_`g'' confirm new variable `name'
            else confirm name `name'
            local saved_fe_names `saved_fe_names' `name'
        }
    }
    local unique_targets : list uniq saved_fe_names
    if `: word count `unique_targets''!=`: word count `saved_fe_names'' error 198
    local all_outputs `saved_fe_names' `d' `tagsep' `zvarname'
    local unique_outputs : list uniq all_outputs
    if `: word count `unique_outputs''!=`: word count `all_outputs'' error 198
    local absorb `absorb_numeric'

    * Parse weights
    local weight_var ""
    local weight_type = 0
    if "`weight'" != "" {
        * Evaluate the weight expression; marks out and validates the weights
        tempvar evaluated_weight
        _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"`exp'"')
        local weight_var "`evaluated_weight'"
        local weight_type = cond("`weight'" == "aweight", 1, cond("`weight'" == "fweight", 2, 3))
    }

    if "`exposure'" != "" & "`offset'" != "" {
        di as error "cpplmhdfe: exposure() and offset() are mutually exclusive"
        exit 198
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

        * Base/omitted placeholders are not sent to the plugin
        _ctools_fvsplit, names(`coef_names_all') vars(`__indepvars_all')
        local indepvars_expanded "`r(vars)'"
        local coef_names "`r(names)'"
        local coef_is_base "`r(base)'"
        local varlist_expanded "`depvar' `indepvars_expanded'"
    }
    else {
        tempvar no_regressors
        quietly generate double `no_regressors' = 0
        local indepvars_expanded `no_regressors'
        local varlist_expanded "`depvar' `no_regressors'"
        local coef_names `no_regressors'
        local coef_names_all `no_regressors'
        local coef_is_base 0
    }

    local nvars : word count `varlist_expanded'
    local nfe : word count `absorb'

    * Parse VCE option (default: robust, matching ppmlhdfe)
    local vcetype = 1
    local clustervar ""
    local clustervar_orig ""
    local cluster_terms
    local num_cluster_vars = 0
    local num_cluster_terms = 0

    if `"`vce'"' != "" {
        _ctools_vce vce_kind clustervar : `vce'
        if "`vce_kind'" == "cluster" {
            local vcetype = 2
            if "`clustervar'" == "" {
                di as error "cpplmhdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            local clustervar_orig "`clustervar'"
            local num_cluster_vars : word count `clustervar'
            if `num_cluster_vars' > 10 {
                di as error "cpplmhdfe: at most ten cluster variables are supported"
                exit 198
            }
            local cluster_numeric
            foreach cl of local clustervar {
                local parts : subinstr local cl "#" " ", all
                local parts : subinstr local parts "i." "", all
                markout `touse' `parts', strok
                * A single numeric variable is passed as it is: the plugin
                * numbers its values in ascending order, the partition and
                * order egen group() gives. Strings and interactions are
                * numbered the same way by _ctools_group.
                local cl_id
                if `: word count `parts'' == 1 {
                    capture confirm numeric variable `parts'
                    if !_rc local cl_id `parts'
                }
                if "`cl_id'" == "" {
                    tempvar cl_id
                    quietly _ctools_group `parts' if `touse', generate(`cl_id')
                }
                local cluster_numeric `cluster_numeric' `cl_id'
            }
            local clustervar `cluster_numeric'
            * Each nonempty subset contributes an inclusion-exclusion term.
            * _ctools_group preserves numeric and string cluster partitions exactly;
            * a subset of one dimension is that dimension's variable.
            local num_cluster_terms = 2^`num_cluster_vars' - 1
            forvalues subset = 1/`num_cluster_terms' {
                local members
                forvalues dim = 1/`num_cluster_vars' {
                    if mod(floor(`subset'/2^(`dim'-1)), 2) {
                        local members `members' `: word `dim' of `clustervar''
                    }
                }
                if `: word count `members'' == 1 {
                    local cluster_id `members'
                }
                else {
                    tempvar cluster_id
                    quietly _ctools_group `members' if `touse', generate(`cluster_id')
                }
                local cluster_terms `cluster_terms' `cluster_id'
            }
        }
        else if inlist("`vce_kind'", "robust", "unadjusted") {
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
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Set up scalars for C plugin
    scalar __cpplmhdfe_api = 0
    foreach option in standardize_data min_ok use_exact_solver use_exact_partial use_heuristic_tol start_inner_tol remove_collinear_variables use_step_halving step_halving_memory max_step_halving mu_tol relu_tol relu_zero_tol relu_maxiter relu_strict relu_accelerate simplex_tol simplex_maxiter {
        local key = cond("`option'"=="remove_collinear_variables","remove_collinear","`option'")
        scalar __cpplmhdfe_`key' = ``option''
    }
    scalar __cpplmhdfe_ref_exact_partial = `ref_exact_partial'
    scalar __cpplmhdfe_acceleration = `accel_code'
    scalar __cpplmhdfe_transform = `transform_code'
    scalar __cpplmhdfe_pool = `poolsize'
    scalar __cpplmhdfe_btol = `btol'
    scalar __cpplmhdfe_conlim = `conlim'
    scalar __cpplmhdfe_guess = `guess_code'
    scalar __cpplmhdfe_diagnostic = ("`tagsep'"!="")
    scalar __cpplmhdfe_check_mu = `do_mu'
    scalar __cpplmhdfe_check_simplex = `do_simplex'
    scalar __cpplmhdfe_guess_idx = 0
    scalar __cpplmhdfe_tag_idx = 0
    scalar __cpplmhdfe_certificate_idx = 0
    scalar __cpplmhdfe_K = `nvars'
    scalar __cpplmhdfe_G = `nfe'
    scalar __cpplmhdfe_verbose = `verbose'
    scalar __cpplmhdfe_maxiter = `iterate'
    scalar __cpplmhdfe_tolerance = `tolerance'
    scalar __cpplmhdfe_vce_type = `vcetype'
    scalar __cpplmhdfe_cluster_terms = `num_cluster_terms'
    scalar __cpplmhdfe_cluster_terms_done = 0
    scalar __cpplmhdfe_has_weights = (`weight_type' > 0)
    scalar __cpplmhdfe_weight_type = `weight_type'
    scalar __cpplmhdfe_has_offset = `has_offset'
    scalar __cpplmhdfe_irls_maxiter = `irlsmaxiter'
    scalar __cpplmhdfe_irls_tol = `irlstolerance'
    scalar __cpplmhdfe_sep_tol = `septolerance'
    scalar __cpplmhdfe_check_fe = `do_fe'
    scalar __cpplmhdfe_check_relu = `do_relu'
    scalar __cpplmhdfe_inner_tol = `itolerance'
    scalar __cpplmhdfe_dof_adjust_type = `dof_adjust_type'
    scalar __cpplmhdfe_dof_flags = `dof_flags'

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record setup time
    _ctools_lap __time_setup __ct_lap98 `__ct_lap98'

    * Build varlist for plugin: depvar indepvars fe_vars [cluster_var] [weight_var] [offset_var]
    local plugin_varlist `varlist_expanded' `absorb'
    if `vcetype' == 2 {
        local plugin_varlist `plugin_varlist' `cluster_terms'
    }
    if `weight_type' > 0 {
        local plugin_varlist `plugin_varlist' `weight_var'
    }
    if `has_offset' {
        local plugin_varlist `plugin_varlist' `offset_var'
    }

    * The plugin returns the final sample after internal row removals.
    tempvar final_sample
    quietly generate byte `final_sample' = 0
    local plugin_varlist `plugin_varlist' `final_sample'
    local sample_idx : word count `plugin_varlist'
    scalar __cpplmhdfe_sample_idx = `sample_idx'

    tempvar saved_d
    quietly generate double `saved_d' = .
    local plugin_varlist `plugin_varlist' `saved_d'
    local d_idx : word count `plugin_varlist'
    scalar __cpplmhdfe_d_idx = `d_idx'

    scalar __cpplmhdfe_keepsingletons = ("`keepsingletons'" != "")
    forvalues g=1/`nfe' {
        scalar __cpplmhdfe_slope_`g' = 0
        if "`slope_var_`g''" != "" {
            local plugin_varlist `plugin_varlist' `slope_var_`g''
            local slope_idx : word count `plugin_varlist'
            scalar __cpplmhdfe_slope_`g' = `slope_idx'
        }
    }

    forvalues g=1/`nfe' {
        scalar __cpplmhdfe_save_fe_`g' = 0
        if "`fe_target_`g''"!="" {
            tempvar fe_output_`g'
            quietly generate double `fe_output_`g'' = .
            local plugin_varlist `plugin_varlist' `fe_output_`g''
            scalar __cpplmhdfe_save_fe_`g' = `: word count `plugin_varlist''
        }
    }
    if `guess_code'==2 {
        local plugin_varlist `plugin_varlist' `guess_var'
        scalar __cpplmhdfe_guess_idx = `: word count `plugin_varlist''
    }
    tempvar separation_tag separation_certificate
    quietly generate byte `separation_tag' = .
    quietly generate double `separation_certificate' = .
    foreach output in tag certificate {
        local variable = cond("`output'"=="tag","`separation_tag'","`separation_certificate'")
        local plugin_varlist `plugin_varlist' `variable'
        scalar __cpplmhdfe_`output'_idx = `: word count `plugin_varlist''
    }

    * Create V matrix
    local K_x = `nvars' - 1
    if `K_x' < 1 {
        local K_x = 1
    }
    matrix __cpplmhdfe_V = J(`K_x', `K_x', 0)
    matrix __cpplmhdfe_V_full = J(`K_x'+1, `K_x'+1, 0)
    matrix __cpplmhdfe_scales = J(1, `K_x'+1, 1)

    * Rebuild the native fit after separation; each pass strictly shrinks the
    * sample and reruns singleton screening and collinearity on retained rows.
    scalar __cpplmhdfe_scale_fixed = 0
    matrix __cpplmhdfe_input_scales = J(1,`nvars',1)
    local total_singletons = 0
    local total_separated = 0
    local refit = 1
    while `refit' {
        scalar __cpplmhdfe_refit = 0
        quietly replace `final_sample' = 0
        capture noisily plugin call ctools_plugin `plugin_varlist' if `touse', ///
            "cpplmhdfe `threads_code' full_regression"
        local reg_rc = _rc
        if !`reg_rc' & __cpplmhdfe_api != 3 {
            di as error "cpplmhdfe: rebuild the ctools plugin alongside the ado files"
            local reg_rc = 601
        }
        if `reg_rc' continue, break
        if "`tagsep'"!="" {
            quietly generate byte `tagsep' = `separation_tag'
            label variable `tagsep' "[cpplmhdfe: ReLU] Obs. is separated"
            format `tagsep' %1.0f
            if "`zvarname'"!="" {
                quietly generate double `zvarname' = `separation_certificate'
                label variable `zvarname' "[cpplmhdfe: ReLU] Certificate of separation"
            }
            _ctools_cleanup cpplmhdfe
            exit
        }
        local total_singletons = `total_singletons' + __cpplmhdfe_num_singletons
        local total_separated = `total_separated' + __cpplmhdfe_num_separated
        local refit = __cpplmhdfe_refit
        if `refit' {
            quietly replace `touse' = `final_sample'
            scalar __cpplmhdfe_scale_fixed = 1
            if `refit'==1 scalar __cpplmhdfe_check_simplex = 0
            if `refit'==2 scalar __cpplmhdfe_check_relu = 0
            if `refit'==3 scalar __cpplmhdfe_check_mu = 0
            scalar __cpplmhdfe_check_fe = 0
        }
    }
    scalar __cpplmhdfe_num_singletons = `total_singletons'
    scalar __cpplmhdfe_num_separated = `total_separated'
    if !`reg_rc' & `num_cluster_terms' > 1 & ///
        __cpplmhdfe_cluster_terms_done != `num_cluster_terms' {
        di as error "cpplmhdfe: rebuild the ctools plugin to enable multiway clustering"
        local reg_rc = 601
    }
    if `reg_rc' {
        _ctools_cleanup cpplmhdfe
        di as error "Error in PPML regression (rc=`reg_rc')"
        exit `reg_rc'
    }

    * Record plugin call time
    _ctools_lap __time_plugin_call __ct_lap98 `__ct_lap98'

    * Match the reference's PSD correction on the full standardized covariance
    * (including its implicit constant), before extracting and rescaling slopes.
    mata: _cpplmhdfe_fix_psd()

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

    if !`has_absorb' local fe_name_1 _cons
    local df_a_nested = __cpplmhdfe_df_nested
    local df_a_adjusted = `df_a'
    tempname dof_table
    matrix `dof_table' = J(`nfe',5,0)
    local fe_names
    forvalues g=1/`nfe' {
        local fe_levels_`g' = __cpplmhdfe_num_levels_`g'
        local fe_nested_`g' = __cpplmhdfe_fe_nested_`g'
        local red = __cpplmhdfe_fe_redundant_`g'
        matrix `dof_table'[`g',1] = `fe_levels_`g''
        matrix `dof_table'[`g',2] = `red'
        matrix `dof_table'[`g',3] = `fe_levels_`g'' - `red'
        matrix `dof_table'[`g',4] = 1 - __cpplmhdfe_fe_exact_`g'
        matrix `dof_table'[`g',5] = `fe_nested_`g''
        local fe_names `fe_names' `fe_name_`g''
    }
    matrix colnames `dof_table' = Categories Redundant Num_Coefs Exact Nested
    matrix rownames `dof_table' = `fe_names'

    * Initialize model degrees of freedom; refine to joint-test rank below.
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

    * Pseudo R-squared (McFadden)
    local r2_p = 1 - `ll' / `ll_0'

    * e(b) and e(V): every term in order, with zeros for omitted regressors
    * (o. prefix) and factor-variable base levels, then the constant. The
    * covariance comes from the PSD-corrected matrix with its implicit constant.
    tempname b V b_kept
    _ctools_ebv `b' `V', names(`coef_names_all') base(`coef_is_base') ///
        stub(__cpplmhdfe) kkeep(`K_keep') vmatrix(__cpplmhdfe_V_full) cons(__cpplmhdfe_cons) kept(`b_kept')

    * Wald statistic on the retained slope coefficients
    local F = .
    if `K_keep' > 0 {
        tempname V_test Vinv Wald
        matrix `V_test' = __cpplmhdfe_V_full[1..`K_keep', 1..`K_keep']
        matrix `Vinv' = syminv(`V_test')
        matrix `Wald' = `b_kept' * `Vinv' * `b_kept''
        if diag0cnt(`Vinv') == 0 local F = `Wald'[1,1] / `df_m'
    }

    * Without regressors only the constant is posted
    if "`indepvars'" == "" {
        local bcols = colsof(`b')
        matrix `b' = `b'[1, `bcols']
        matrix `V' = `V'[`bcols', `bcols']
        matrix colnames `b' = _cons
        matrix colnames `V' = _cons
        matrix rownames `V' = _cons
    }
    if "`d'" != "" {
        quietly generate double `d' = `saved_d'
        label variable `d' "Sum of fixed effects"
    }

    forvalues g=1/`nfe' {
        if "`fe_target_`g''"!="" {
            if `fe_automatic_`g'' capture drop `fe_target_`g''
            quietly generate double `fe_target_`g'' = `fe_output_`g''
            label variable `fe_target_`g'' "[FE] `fe_name_`g''"
        }
    }
    if "`tagsep'"!="" quietly generate byte `tagsep' = `separation_tag'
    if "`zvarname'"!="" quietly generate double `zvarname' = `separation_certificate'

    * Record results building time
    _ctools_lap __time_build __ct_lap98 `__ct_lap98'

    * Post results
    ereturn post `b' `V', esample(`final_sample') depname(`depvar_orig') obs(`N_final')

    * Store e() results
    ereturn scalar N = `N_final'
    ereturn scalar df_m = `df_m'
    ereturn scalar df = `df_r'
    ereturn scalar ll = `ll'
    ereturn scalar ll_0 = `ll_0'
    ereturn scalar deviance = `deviance'
    ereturn scalar r2_p = `r2_p'
    ereturn scalar chi2 = `F' * `df_m'
    ereturn local chi2type Wald
    ereturn local title "HDFE PPML regression"
    ereturn local separation "`separation'"
    ereturn local properties "b V"
    ereturn local marginsok default
    ereturn local dofmethod "`dofmethods'"
    ereturn scalar rank = `K_keep'
    local nfe_original : word count `all_fe_keys'
    ereturn scalar N_hdfe = `nfe_original'
    ereturn scalar num_singletons = `num_singletons'
    ereturn scalar num_separated = `num_separated'
    ereturn scalar ic = `irls_iters'
    ereturn scalar converged = `irls_converged'
    ereturn scalar df_a = `df_a_adjusted'
    ereturn scalar df_a_initial = __cpplmhdfe_df_initial
    ereturn scalar df_a_redundant = __cpplmhdfe_df_redundant
    ereturn scalar N_hdfe_extended = `nfe'
    ereturn scalar N_full = `N_final' + `num_singletons' + `num_separated'
    ereturn scalar drop_singletons = ("`keepsingletons'" == "")
    ereturn matrix dof_table = `dof_table'
    ereturn local extended_absvars "`fe_names'"
    ereturn scalar df_a_nested = `df_a_nested'

    if `vcetype' == 2 {
        ereturn scalar N_clust = __cpplmhdfe_N_clust
        ereturn scalar N_clustervars = `num_cluster_vars'
        forvalues dim = 1/`num_cluster_vars' {
            ereturn scalar N_clust`dim' = __cpplmhdfe_N_clust_`dim'
        }
        ereturn local clustvar "`clustervar_orig'"
    }

    ereturn local absorb "`absorb_orig'"
    ereturn local absvars "`absorb_orig'"
    ereturn local d "`d'"
    ereturn local offset "`exposure'`offset'"
    ereturn local exposure "`exposure'"
    if "`exposure'" != "" ereturn local offset "ln(`exposure')"
    ereturn local depvar "`depvar_orig'"
    ereturn local wtype "`weight'"
    ereturn local wexp "`exp'"
    ereturn local indepvars "`indepvars'"
    ereturn local vce = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", "cluster"))
    ereturn local predict "cpplmhdfe_p"
    ereturn local cmd "cpplmhdfe"
    ereturn local cmdline "cpplmhdfe `0'"

    * Record ereturn posting time
    _ctools_lap __time_ereturn __ct_lap98 `__ct_lap98'

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
    di as text "Absorbing " as result `nfe' as text " HDFE " cond(`nfe'>1, "groups", "group") _col(49) "Wald chi2(" as result %3.0f `df_m' as text ")" _col(67) "= " as result %10.2f e(chi2)
    di as text _col(49) "Prob > chi2" _col(67) "= " as result %10.4f chi2tail(`df_m', e(chi2))
    di as text _col(49) "Pseudo R2" _col(67) "= " as result %10.4f `r2_p'
    di as text _col(49) "Log pseudolikelihood" _col(67) "= " as result %10.4f `ll'

    di as text ""
    ereturn display, `eform_opt' `diopts'

    di as text ""
    di as text "Absorbed degrees of freedom:"
    matrix list e(dof_table), noheader format(%9.0g)

    _ctools_cleanup cpplmhdfe

    * Display timing if verbose
    _ctools_lap __time_display __ct_lap98 `__ct_lap98'
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

end

capture mata: mata drop _cpplmhdfe_fix_psd()
mata:
void _cpplmhdfe_fix_psd()
{
    real matrix V, U
    real rowvector lambda, scales
    real scalar k
    k = st_numscalar("__cpplmhdfe_K_keep")
    V = st_matrix("__cpplmhdfe_V_full")[|1,1 \ k+1,k+1|]
    V = (V + V') / 2
    if (st_numscalar("__cpplmhdfe_cluster_terms") > 1) {
    symeigensystem(V, U=., lambda=.)
    if (min(lambda) < 0) {
        lambda = lambda :* (lambda :>= 0)
        V = quadcross(U', lambda, U')
        printf("{txt}Warning: non-positive-semidefinite clustered covariance; negative eigenvalues set to zero.\n")
    }
    }
    scales = st_matrix("__cpplmhdfe_scales")[1,1..k+1]
    V = V :/ (scales' * scales)
    st_matrix("__cpplmhdfe_V_full", V)
    if (k) st_matrix("__cpplmhdfe_V", V[|1,1 \ k,k|])
}
end
