*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define cppmlhdfe, eclass
    version 14.1

    if replay() {
        * As ppmlhdfe: the standalone options version and reload, else the
        * header, coefficient table and footnote of the last fit
        syntax [, *]
        if inlist(`"`options'"', "version", "reload") {
            _cppmlhdfe_`options'
            exit
        }
        if "`e(cmd)'" != "cppmlhdfe" error 301
        _cppmlhdfe_replay `0'
        exit
    }
    * Accept the legacy bare verbose option as well as ppmlhdfe's verbose(#).
    * Only an option token is rewritten, never a variable named verbose in
    * the varlist, the if condition or an option argument; e(cmdline) keeps
    * the command as typed.
    local cmdline_typed `"`0'"'
    _parse comma __lhs __rhs : 0
    if `"`__rhs'"' != "" {
        gettoken __comma __rhs : __rhs, parse(",")
        local __opts
        while `"`__rhs'"' != "" {
            gettoken __tok __rhs : __rhs, bind quotes
            if `"`__tok'"' == "verbose" local __tok verbose(1)
            local __opts `"`__opts' `__tok'"'
        }
        local 0 `"`__lhs', `__opts'"'
    }

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
    _ctools_cleanup cppmlhdfe

    syntax varlist(min=1 fv ts) [aw fw pw] [if] [in] [, Absorb(string) VCE(string) FEVCE(string) Verbose(integer 0) ///
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
        ACCELeration(string) TRANSform(string) POOLsize(integer 0) BTOL(real 1e-8) CONLIM(real 1e8) ///
        NOWILDtable R2 target_inner_tol(real -1) relu_debug(integer 0) relu_report_r2(integer 0) ///
        PRUNE CONDition RRE(varname numeric) ALGorithm(string) SLOPEmethod(string) *]

    * The absorb() suboptions ppmlhdfe passes to reghdfe, and its GLM
    * properties (min_ok(), relu_tol()...) there
    local save_all_fe = 0
    local hdfe_verbose
    local comma = strpos("`absorb'",",")
    if `comma' {
        local absorb_options = substr("`absorb'",`comma'+1,.)
        local absorb = substr("`absorb'",1,`comma'-1)
        _cppmlhdfe_abs_options, `absorb_options'
    }
    * The partialling options ppmlhdfe hands to reghdfe, where they do not
    * change the estimates: algorithm() and slopemethod() are read and
    * ignored there, rre() reports relative residual errors while
    * partialling, prune removes degree-1 vertices of two fixed effects
    * (a speedup). condition requests reghdfe's finite condition number,
    * which ppmlhdfe posts as -1 (not computed).
    * target_inner_tol() is itolerance(); relu_report_r2() is r2.
    if `target_inner_tol' > 0 local itolerance `target_inner_tol'
    else if `target_inner_tol' != -1 error 198
    foreach option in relu_debug relu_report_r2 {
        if !inlist(``option'', 0, 1) error 198
    }
    if `relu_report_r2' local r2 r2
    * The HDFE object's verbosity (absorb(, verbose())), which shows or hides
    * the note on dropped singletons; else that of the command
    if "`hdfe_verbose'" == "" local hdfe_verbose = cond(`verbose' > 0, `verbose' - 1, `verbose')

    local __do_timing = (`verbose' > 0 | "`timeit'" != "")
    * wildbootstrap[(wildopts)]: the score bootstrap of the coefficients
    _cppmlhdfe_wild_opt, `options'
    _get_diopts diopts options, `__wild_rest'
    if "`options'" != "" {
        di as error "unsupported option(s): `options'"
        exit 198
    }
    * Exponentiated coefficients (eform, or its synonym irr) are shown as
    * exp(b), as in ppmlhdfe
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
        di as text "{bf: (identifying separated observations instead of running regressions)}"
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

    * acceleration() and transform() as reghdfe reads them: any case, cg, sd
    * and off for conjugate_gradient, steepest_descent and none, and any
    * leading part of a name (the last name it begins wins)
    local acceleration = lower(strtrim(`"`acceleration'"'))
    if "`acceleration'"=="" local acceleration cg
    if "`acceleration'"=="cg" local acceleration conjugate_gradient
    if "`acceleration'"=="sd" local acceleration steepest_descent
    if "`acceleration'"=="off" local acceleration none
    foreach name in conjugate_gradient steepest_descent aitken none hybrid lsmr {
        if strpos("`name'", "`acceleration'") == 1 local acceleration `name'
    }
    * the plugin's codes: 0 cg, 1 none, 2 sd, 3 aitken, 4 hybrid, 5 lsmr
    local __accels conjugate_gradient none steepest_descent aitken hybrid lsmr
    local accel_code = `: list posof "`acceleration'" in __accels' - 1
    if `accel_code' < 0 error 198
    local transform = lower(strtrim(`"`transform'"'))
    if "`transform'"=="" local transform symmetric_kaczmarz
    foreach name in cimmino kaczmarz symmetric_kaczmarz {
        if strpos("`name'", "`transform'") == 1 local transform `name'
    }
    local transform_code = cond("`transform'"=="symmetric_kaczmarz", 0, ///
        cond("`transform'"=="kaczmarz", 1, cond("`transform'"=="cimmino", 2, .)))
    if missing(`transform_code',`btol',`conlim',`poolsize') | `poolsize'<0 | `btol'<=0 | `conlim'<=0 error 198

    * Match ppmlhdfe defaults: if IRLS tolerance is not specified, use tolerance()
    if `irlstolerance' < 0 {
        local irlstolerance `tolerance'
    }

    if `iterate' < 1 | `irlsmaxiter' < 1 | `threads'<0 | ///
        missing(`tolerance', `irlstolerance', `septolerance', `iterate', `irlsmaxiter', `threads', `itolerance') | ///
        `tolerance' <= 0 | `irlstolerance' <= 0 | `septolerance' <= 0 {
        di as error "cppmlhdfe: iteration limits and tolerances must be positive"
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
    }
    if "`d'" != "" & "`d2'" == "" _ctools_newvars `d'

    * Parse independent absorbed-DF adjustments, matching reference semantics
    * (abbreviations included); all stands for pairwise clusters continuous.
    local dof_lower = lower(`"`dofadjustments'"')
    _ctools_hdfe_dof, `dof_lower'
    local dof_flags = `s(flags)'
    local dofmethods `s(methods)'
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
        * ppmlhdfe's alternative slope syntax id##(c.x c.z) is id##c.(x z)
        local alt = strpos("`term'", "#(")
        if `alt' {
            if substr("`term'", -1, 1) != ")" error 198
            local inside = substr("`term'", `alt' + 2, strlen("`term'") - `alt' - 2)
            local rebuilt
            foreach item of local inside {
                if substr("`item'", 1, 2) != "c." | strpos(substr("`item'", 3, .), ".") error 198
                local rebuilt `rebuilt' `=substr("`item'", 3, .)'
            }
            local term = substr("`term'", 1, `alt') + "c.(`rebuilt')"
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
    * e(N_hdfe), as reghdfe counts: the absorb() terms (one without absorb())
    local n_absvars = cond(`has_absorb', `original_group', 1)
    if "`condition'" != "" & `n_absvars' != 2 {
        di as error "Computing finite condition number requires two FEs"
        exit 198
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
        local absorb_model_vars `absorb_model_vars' `categories' `continuous'
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
        di as error "cppmlhdfe: exposure() and offset() are mutually exclusive"
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
            di as error "cppmlhdfe: exposure variable must be positive"
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

    * vce(bootstrap): a pairs bootstrap of the coefficients, run by the
    * plugin on the retained fit (docs/development/BOOTSTRAP_PLAN.md, 8). The
    * fit itself uses the robust VCE, posted as e(V_modelbased).
    local __bs = 0
    local bs_unitvar
    local bs_stratavar
    if `"`vce'"' != "" {
        _ctools_hdfe_boot parse, spec(`"`vce'"')
        if "`r(kind)'" == "bootstrap" {
            local __bs = 1
            local bs_reps = r(reps)
            local bs_seed `"`r(seed)'"'
            local bs_cluster `"`r(cluster)'"'
            local bs_strata `"`r(strata)'"'
            local bs_mse = r(mse)
            local bs_level = r(level)
            local __g `r(given)'
            local bs_lgiven : list posof "level" in __g
            local bs_dots = r(dots)
            local bs_draws `r(draws)'
            local bs_bsdraws `r(bsdraws)'
            local bs_save_file `"`r(save_file)'"'
            local bs_save_double `r(save_double)'
            local bs_save_every `r(save_every)'
            local bs_save_replace `r(save_replace)'
            local bs_designopts = r(designopts)
            local vce
        }
    }
    * The command's level() sets vce(bootstrap)'s, as for bootstrap; the
    * table at estimation shows the bootstrap's level
    if `__bs' {
        _ctools_hdfe_boot level, cmd(cppmlhdfe) bslevel(`bs_level') bsgiven(`bs_lgiven') diopts(`"`diopts'"')
        local bs_level = r(level)
        local diopts `"`r(diopts)'"'
    }
    * fevce(bootstrap): bootstrap standard errors of the saved fixed effects
    * and of d (plan 4.3, 8.2). With vce(bootstrap) it shares its replicates;
    * options given to both must agree.
    local __fe = 0
    local fe_any = 0
    local fe_nreps = 0
    local fe_ci = 0
    local fe_dse
    local fe_sestub
    if `"`fevce'"' != "" {
        _ctools_hdfe_boot parse, spec(`"`fevce'"') fe
        if "`r(kind)'" != "bootstrap" {
            di as error "cppmlhdfe: fevce() takes bootstrap[, options]: the fixed effects use the pairs bootstrap"
            exit 198
        }
        forvalues g = 1/`group_count' {
            if "`fe_target_`g''" != "" local fe_any = 1
        }
        if !`fe_any' & "`d'" == "" {
            di as error "cppmlhdfe: fevce(bootstrap) needs saved fixed effects (absorb(..., savefe) or newvar=) or d()"
            exit 198
        }
        if `__bs' & (`"`r(saving)'"' != "" | "`r(bsdraws)'" != "") {
            di as error "cppmlhdfe: with vce(bootstrap), specify saving() and _bsdraws() there"
            exit 198
        }
        if r(ci) & !`fe_any' {
            di as error "cppmlhdfe: fevce(bootstrap, ci) needs saved fixed effects (absorb(..., savefe) or newvar=)"
            exit 198
        }
        local __fe = 1
        local fe_sestub `r(sestub)'
        local fe_nreps = r(nreps)
        local fe_ci = r(ci)
        local fe_dse `r(dse)'
        if `__bs' {
            * level(), mse and dots() of vce(bootstrap) govern the shared
            * replicates: fevce() may repeat them, not change them
            local __given `r(given)'
            if (`: list posof "level" in __given' & r(level) != `bs_level') | ///
                (`: list posof "mse" in __given' & r(mse) != `bs_mse') | ///
                (`: list posof "dots" in __given' & r(dots) != `bs_dots') {
                di as error "cppmlhdfe: with vce(bootstrap), fevce(bootstrap) uses the level(), mse and dots() of vce(); specify them once, in vce()"
                exit 198
            }
        }
        if `__bs' & r(designopts) {
            if r(reps) != `bs_reps' | `"`r(seed)'"' != `"`bs_seed'"' | `"`r(cluster)'"' != `"`bs_cluster'"' | ///
                `"`r(strata)'"' != `"`bs_strata'"' | "`r(draws)'" != "`bs_draws'" {
                di as error "cppmlhdfe: vce(bootstrap) and fevce(bootstrap) share one set of replicates; specify the resampling options once"
                exit 198
            }
        }
        else if !`__bs' {
            local bs_reps = r(reps)
            local bs_seed `"`r(seed)'"'
            local bs_cluster `"`r(cluster)'"'
            local bs_strata `"`r(strata)'"'
            local bs_mse = r(mse)
            local bs_level = r(level)
            local bs_dots = r(dots)
            local bs_draws `r(draws)'
            local bs_bsdraws `r(bsdraws)'
            local bs_save_file `"`r(save_file)'"'
            local bs_save_double `r(save_double)'
            local bs_save_every `r(save_every)'
            local bs_save_replace `r(save_replace)'
        }
    }
    local __bsrun = (`__bs' | `__fe')
    * wildbootstrap(): the score bootstrap (Kline and Santos), boottest's
    * conventions (plan 4.7, 8.4)
    if `__wild' {
        if `__bs' {
            di as error "cppmlhdfe: wildbootstrap() and vce(bootstrap) are two bootstraps of the coefficients; choose one"
            exit 198
        }
        _ctools_hdfe_boot wildparse, spec(`"`__wildspec'"')
        if !r(null) {
            di as error "cppmlhdfe: the score bootstrap imposes the null; nonull is not available"
            exit 198
        }
        local wb_reps = r(reps)
        local wb_seed `"`r(seed)'"'
        local wb_type `r(weighttype)'
        local wb_coefs `"`r(coefs)'"'
        local wb_ptype `r(ptype)'
        local wb_ptype_code = r(ptype_code)
        local wb_enumerate = r(enumerate)
        local wb_ci = r(ci)
        local wb_level = r(level)
        local wb_bootcluster `"`r(bootcluster)'"'
        local wb_dots = r(dots)
    }
    if `__bsrun' | `__wild' {
        local what = cond(`__bs', "vce(bootstrap)", cond(`__fe', "fevce(bootstrap)", "wildbootstrap()"))
        if "`tagsep'`zvarname'" != "" {
            di as error "cppmlhdfe: `what' may not be combined with tagsep() or zvarname()"
            exit 198
        }
        if `: word count `varlist'' < 2 {
            di as error "cppmlhdfe: `what' needs a regressor"
            exit 198
        }
        if `"`bs_cluster'"' != "" {
            _ctools_hdfe_vce, spec(`"cluster `bs_cluster'"') touse(`touse') cmd(cppmlhdfe)
            local bs_cluster_name `r(name_1)'
            tempvar bs_unit
            if r(group_1) quietly _ctools_group `r(vars_1)' if `touse', generate(`bs_unit')
            else quietly generate double `bs_unit' = `r(vars_1)' if `touse'
            local bs_unitvar `bs_unit'
        }
        if "`bs_strata'" != "" {
            markout `touse' `bs_strata', strok
            tempvar bs_str
            quietly _ctools_group `bs_strata' if `touse', generate(`bs_str')
            local bs_stratavar `bs_str'
        }
    }
    * fevce(): output names. Those of automatic outputs (savefe, bare d) are
    * replaced like their effects; named ones must be new.
    local fe_se_names
    local fe_nrep_names
    local fe_lb_names
    local fe_ub_names
    local fe_new_names
    local fe_targets_flags
    if `__fe' {
        forvalues g = 1/`group_count' {
            local fe_targets_flags `fe_targets_flags' `=("`fe_target_`g''" != "")'
            if "`fe_target_`g''" == "" continue
            if "`fe_sestub'" != "" {
                local se_name `fe_sestub'`g'
                local nrep_name `fe_sestub'`g'_nrep
                local lb_name `fe_sestub'`g'_lb
                local ub_name `fe_sestub'`g'_ub
                local fe_new_names `fe_new_names' `se_name' `=cond(`fe_nreps', "`nrep_name'", "")' ///
                    `=cond(`fe_ci', "`lb_name' `ub_name'", "")'
            }
            else if `fe_automatic_`g'' {
                * after savefe's __hdfe#__ and __hdfe#__Slope#
                local se_name : subinstr local fe_target_`g' "__hdfe" "__fese"
                local nrep_name : subinstr local fe_target_`g' "__hdfe" "__fenrep"
                local lb_name : subinstr local fe_target_`g' "__hdfe" "__felb"
                local ub_name : subinstr local fe_target_`g' "__hdfe" "__feub"
            }
            else {
                local se_name `fe_target_`g''_se
                local nrep_name `fe_target_`g''_nrep
                local lb_name `fe_target_`g''_lb
                local ub_name `fe_target_`g''_ub
                local fe_new_names `fe_new_names' `se_name' `=cond(`fe_nreps', "`nrep_name'", "")' ///
                    `=cond(`fe_ci', "`lb_name' `ub_name'", "")'
            }
            local fe_se_names `fe_se_names' `se_name'
            if `fe_nreps' local fe_nrep_names `fe_nrep_names' `nrep_name'
            if `fe_ci' local fe_lb_names `fe_lb_names' `lb_name'
            if `fe_ci' local fe_ub_names `fe_ub_names' `ub_name'
        }
        if "`fe_dse'" != "" local fe_new_names `fe_new_names' `fe_dse'
        else if "`d'" != "" {
            local fe_dse `d'_se
            if "`d2'" == "" local fe_new_names `fe_new_names' `fe_dse'
        }
        if "`fe_new_names'" != "" _ctools_newvars `fe_new_names'
    }

    * Parse VCE option (default: robust, matching ppmlhdfe). Its cluster
    * variables are marked out before the factor variables are expanded, so
    * that their levels come from the estimation sample.
    local vcetype = 1
    local clustervar ""
    local clustervar_orig ""
    local cluster_terms
    local num_cluster_vars = 0
    local num_cluster_terms = 0

    if `"`vce'"' != "" {
        * ppmlhdfe reads vce() with ms_parse_vce: an ols or un prefix and
        * conventional mean unadjusted (which PPML reports as robust), an r
        * prefix robust; iid stays rejected, as there
        gettoken vce_word vce_rest : vce, parse(" ,")
        if substr("`vce_word'", 1, 3) == "ols" | substr("`vce_word'", 1, 2) == "un" | "`vce_word'" == "conventional" {
            local vce_word unadjusted
        }
        else if substr("`vce_word'", 1, 1) == "r" local vce_word robust
        local vce `"`vce_word'`vce_rest'"'
        _ctools_hdfe_vce, spec(`"`vce'"') touse(`touse') cmd(cppmlhdfe)
        local vce_kind `r(kind)'
        if "`vce_kind'" == "cluster" {
            local vcetype = 2
            local num_cluster_vars = r(ndims)
            if `num_cluster_vars' == 0 {
                di as error "cppmlhdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            if `num_cluster_vars' > 10 {
                di as error "cppmlhdfe: at most ten cluster variables are supported"
                exit 198
            }
            forvalues dim = 1/`num_cluster_vars' {
                local cl_name_`dim' `r(name_`dim')'
                local cl_vars_`dim' `r(vars_`dim')'
                local cl_group_`dim' = r(group_`dim')
            }
            * A single numeric variable is passed as it is: the plugin numbers
            * its values in ascending order, the partition and order egen
            * group() gives. Strings and interactions are numbered the same
            * way by _ctools_group.
            local cluster_numeric
            forvalues dim = 1/`num_cluster_vars' {
                local clustervar_orig `clustervar_orig' `cl_name_`dim''
                if `cl_group_`dim'' {
                    tempvar cl_id_`dim'
                    quietly _ctools_group `cl_vars_`dim'' if `touse', generate(`cl_id_`dim'')
                    local cluster_numeric `cluster_numeric' `cl_id_`dim''
                }
                else local cluster_numeric `cluster_numeric' `cl_vars_`dim''
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
            di as error "cppmlhdfe: unrecognized vce() option"
            exit 9
        }
    }

    * wildbootstrap(): one clustering dimension; the bootstrap clusters are
    * the vce() cluster or bootcluster()
    local wb_bootvar
    local wb_bootname
    local wb_vcevar
    if `__wild' {
        if `num_cluster_vars' > 1 {
            di as error "cppmlhdfe: wildbootstrap() supports one clustering dimension"
            exit 198
        }
        if `vcetype' == 2 {
            local wb_vcevar `clustervar'
            local wb_bootname `clustervar_orig'
        }
        if `"`wb_bootcluster'"' != "" {
            _ctools_hdfe_vce, spec(`"cluster `wb_bootcluster'"') touse(`touse') cmd(cppmlhdfe)
            if r(ndims) != 1 {
                di as error "cppmlhdfe: wildbootstrap(bootcluster()) takes one variable"
                exit 198
            }
            local wb_bootname `r(name_1)'
            tempvar wb_boot
            if r(group_1) quietly _ctools_group `r(vars_1)' if `touse', generate(`wb_boot')
            else quietly generate double `wb_boot' = `r(vars_1)' if `touse'
            local wb_bootvar `wb_boot'
        }
        if `vcetype' == 1 & "`wb_bootvar'" == "" & `weight_type' == 2 {
            di as error "cppmlhdfe: wildbootstrap() draws one weight per observation under vce(robust); with fweights a row stands for several observations. Use bootcluster()"
            exit 198
        }
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
        di as error "cppmlhdfe: dependent variable must be non-negative"
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

    * Automatic outputs (d2's _ppmlhdfe_d, savefe's __hdfe#__) replace an
    * existing variable of that name; refuse when the model itself uses it
    local automatic_targets
    if "`d2'" != "" local automatic_targets _ppmlhdfe_d
    forvalues g=1/`nfe' {
        if "`fe_target_`g''" != "" & `fe_automatic_`g'' local automatic_targets `automatic_targets' `fe_target_`g''
    }
    if "`automatic_targets'" != "" {
        local model_vars `absorb_model_vars' `exposure' `offset'
        fvrevar `depvar_orig' `indepvars', list
        local model_vars `model_vars' `r(varlist)'
        local exp_tokens = ustrregexra(`"`exp' `clustervar_orig'"', "[^A-Za-z0-9_]", " ")
        local model_vars `model_vars' `exp_tokens'
        local used_targets : list automatic_targets & model_vars
        if "`used_targets'" != "" {
            local used_target : word 1 of `used_targets'
            di as error "cppmlhdfe: `used_target' is used in the model and would be replaced by an automatic output; rename it or name the output"
            exit 198
        }
    }

    * Count observations
    qui count if `touse'
    local nobs = r(N)
    if `nobs' == 0 {
        di as error "cppmlhdfe: no observations"
        exit 2000
    }

    * Load plugin
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Set up scalars for C plugin
    scalar __cppmlhdfe_api = 0
    foreach option in standardize_data min_ok use_exact_solver use_exact_partial use_heuristic_tol start_inner_tol remove_collinear_variables use_step_halving step_halving_memory max_step_halving mu_tol relu_tol relu_zero_tol relu_maxiter relu_strict relu_accelerate simplex_tol simplex_maxiter {
        local key = cond("`option'"=="remove_collinear_variables","remove_collinear","`option'")
        scalar __cppmlhdfe_`key' = ``option''
    }
    scalar __cppmlhdfe_ref_exact_partial = `ref_exact_partial'
    scalar __cppmlhdfe_acceleration = `accel_code'
    scalar __cppmlhdfe_transform = `transform_code'
    scalar __cppmlhdfe_pool = `poolsize'
    scalar __cppmlhdfe_btol = `btol'
    scalar __cppmlhdfe_conlim = `conlim'
    scalar __cppmlhdfe_guess = `guess_code'
    scalar __cppmlhdfe_diagnostic = ("`tagsep'"!="")
    scalar __cppmlhdfe_check_mu = `do_mu'
    scalar __cppmlhdfe_check_simplex = `do_simplex'
    scalar __cppmlhdfe_guess_idx = 0
    scalar __cppmlhdfe_tag_idx = 0
    scalar __cppmlhdfe_certificate_idx = 0
    scalar __cppmlhdfe_K = `nvars'
    scalar __cppmlhdfe_G = `nfe'
    scalar __cppmlhdfe_verbose = `verbose'
    * ppmlhdfe's iteration log, printed by the plugin unless nolog or
    * verbose(-1)
    scalar __cppmlhdfe_log = ("`log'" != "nolog")
    scalar __cppmlhdfe_maxiter = `iterate'
    scalar __cppmlhdfe_tolerance = `tolerance'
    scalar __cppmlhdfe_vce_type = `vcetype'
    scalar __cppmlhdfe_cluster_terms = `num_cluster_terms'
    scalar __cppmlhdfe_cluster_terms_done = 0
    scalar __cppmlhdfe_has_weights = (`weight_type' > 0)
    scalar __cppmlhdfe_weight_type = `weight_type'
    scalar __cppmlhdfe_has_offset = `has_offset'
    scalar __cppmlhdfe_irls_maxiter = `irlsmaxiter'
    scalar __cppmlhdfe_irls_tol = `irlstolerance'
    scalar __cppmlhdfe_sep_tol = `septolerance'
    scalar __cppmlhdfe_check_fe = `do_fe'
    scalar __cppmlhdfe_check_relu = `do_relu'
    scalar __cppmlhdfe_inner_tol = `itolerance'
    scalar __cppmlhdfe_dof_adjust_type = `dof_adjust_type'
    scalar __cppmlhdfe_dof_flags = `dof_flags'
    * vce(bootstrap): the separation checks each replicate runs (the refit
    * loop below turns off the one that triggered a refit)
    scalar __cppmlhdfe_bs_check_fe = `do_fe'
    scalar __cppmlhdfe_bs_check_simplex = `do_simplex'
    scalar __cppmlhdfe_bs_check_relu = `do_relu'
    scalar __cppmlhdfe_bs_check_mu = `do_mu'
    if `__fe' & `fe_any' {
        matrix __cppmlhdfe_bs_fe_targets = J(1, `nfe', 0)
        forvalues g = 1/`nfe' {
            matrix __cppmlhdfe_bs_fe_targets[1, `g'] = `: word `g' of `fe_targets_flags''
        }
    }

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
    scalar __cppmlhdfe_sample_idx = `sample_idx'

    tempvar saved_d
    quietly generate double `saved_d' = .
    local plugin_varlist `plugin_varlist' `saved_d'
    local d_idx : word count `plugin_varlist'
    scalar __cppmlhdfe_d_idx = `d_idx'

    scalar __cppmlhdfe_keepsingletons = ("`keepsingletons'" != "")
    * The plugin announces dropped singletons unless the HDFE verbosity is -1
    scalar __cppmlhdfe_hdfe_verbose = `hdfe_verbose'
    forvalues g=1/`nfe' {
        scalar __cppmlhdfe_slope_`g' = 0
        if "`slope_var_`g''" != "" {
            local plugin_varlist `plugin_varlist' `slope_var_`g''
            local slope_idx : word count `plugin_varlist'
            scalar __cppmlhdfe_slope_`g' = `slope_idx'
        }
    }

    forvalues g=1/`nfe' {
        scalar __cppmlhdfe_save_fe_`g' = 0
        if "`fe_target_`g''"!="" {
            tempvar fe_output_`g'
            quietly generate double `fe_output_`g'' = .
            local plugin_varlist `plugin_varlist' `fe_output_`g''
            scalar __cppmlhdfe_save_fe_`g' = `: word count `plugin_varlist''
        }
    }
    if `guess_code'==2 {
        local plugin_varlist `plugin_varlist' `guess_var'
        scalar __cppmlhdfe_guess_idx = `: word count `plugin_varlist''
    }
    tempvar separation_tag separation_certificate
    quietly generate byte `separation_tag' = .
    quietly generate double `separation_certificate' = .
    foreach output in tag certificate {
        local variable = cond("`output'"=="tag","`separation_tag'","`separation_certificate'")
        local plugin_varlist `plugin_varlist' `variable'
        scalar __cppmlhdfe_`output'_idx = `: word count `plugin_varlist''
    }

    * vce(bootstrap): the resampling units and strata (numbered in sorted
    * order, as bsample orders them)
    scalar __cppmlhdfe_bs_unit_idx = 0
    scalar __cppmlhdfe_bs_strata_idx = 0
    if "`bs_unitvar'" != "" {
        local plugin_varlist `plugin_varlist' `bs_unitvar'
        scalar __cppmlhdfe_bs_unit_idx = `: word count `plugin_varlist''
    }
    if "`bs_stratavar'" != "" {
        local plugin_varlist `plugin_varlist' `bs_stratavar'
        scalar __cppmlhdfe_bs_strata_idx = `: word count `plugin_varlist''
    }

    * The sample marker goes last, in place of if `touse': the plugin reads
    * it in parallel, where Stata answers an if one observation at a time
    local plugin_varlist `plugin_varlist' `touse'
    scalar __cppmlhdfe_touse = 1

    * Create V matrix
    local K_x = `nvars' - 1
    if `K_x' < 1 {
        local K_x = 1
    }
    matrix __cppmlhdfe_V = J(`K_x', `K_x', 0)
    matrix __cppmlhdfe_V_full = J(`K_x'+1, `K_x'+1, 0)
    matrix __cppmlhdfe_scales = J(1, `K_x'+1, 1)

    * Rebuild the native fit after separation; each pass strictly shrinks the
    * sample and reruns singleton screening and collinearity on retained rows.
    scalar __cppmlhdfe_scale_fixed = 0
    matrix __cppmlhdfe_input_scales = J(1,`nvars',1)
    local total_singletons = 0
    local total_separated = 0
    if `verbose' > -1 & "`keepsingletons'" != "" {
        di as error "warning: keeping singleton groups will keep fixed effects that cause separation"
    }
    local refit = 1
    while `refit' {
        scalar __cppmlhdfe_refit = 0
        quietly replace `final_sample' = 0
        capture noisily plugin call ctools_plugin `plugin_varlist', ///
            "cppmlhdfe `threads_code' full_regression`=cond(`__bsrun' | `__wild', " retain", "")'"
        local reg_rc = _rc
        if !`reg_rc' & __cppmlhdfe_api != 3 {
            di as error "cppmlhdfe: rebuild the ctools plugin alongside the ado files"
            local reg_rc = 601
        }
        if `reg_rc' continue, break
        if "`tagsep'"!="" {
            quietly generate byte `tagsep' = `separation_tag'
            label variable `tagsep' "[cppmlhdfe: ReLU] Obs. is separated"
            format `tagsep' %1.0f
            if "`zvarname'"!="" {
                quietly generate double `zvarname' = `separation_certificate'
                label variable `zvarname' "[cppmlhdfe: ReLU] Certificate of separation"
            }
            _ctools_cleanup cppmlhdfe
            * r2: as ppmlhdfe, regress the certificate on the regressors and
            * the absorbed effects (with creghdfe), whose R-squared of 1
            * verifies it; only when observations were separated
            if "`r2'" != "" & "`zvarname'" != "" {
                quietly count if `tagsep' == 1
                if r(N) {
                    local r2_absorb noabsorb
                    if `has_absorb' {
                        local r2_absvars
                        local r2_rest `"`absorb_orig'"'
                        while `"`r2_rest'"' != "" {
                            gettoken r2_term r2_rest : r2_rest, bind
                            local r2_eq = strpos(`"`r2_term'"', "=")
                            if `r2_eq' local r2_term = substr(`"`r2_term'"', `r2_eq' + 1, .)
                            local r2_absvars `r2_absvars' `r2_term'
                        }
                        local r2_absorb absorb(`r2_absvars')
                    }
                    local r2_cmd creghdfe `zvarname' `indepvars', `r2_absorb'
                    di as input _n "{hline 32} o {hline 32}"
                    di as text "{bf:Verifying certificate of separation:}"
                    di as input `"{bf:. `r2_cmd'}"'
                    `r2_cmd'
                    di as input "{hline 32} o {hline 32}" _n
                }
            }
            exit
        }
        local total_singletons = `total_singletons' + __cppmlhdfe_num_singletons
        local total_separated = `total_separated' + __cppmlhdfe_num_separated
        local refit = __cppmlhdfe_refit
        if `refit' {
            quietly replace `touse' = `final_sample'
            scalar __cppmlhdfe_scale_fixed = 1
            if `refit'==1 scalar __cppmlhdfe_check_simplex = 0
            if `refit'==2 scalar __cppmlhdfe_check_relu = 0
            if `refit'==3 scalar __cppmlhdfe_check_mu = 0
            scalar __cppmlhdfe_check_fe = 0
        }
    }
    scalar __cppmlhdfe_num_singletons = `total_singletons'
    scalar __cppmlhdfe_num_separated = `total_separated'
    if !`reg_rc' & `num_cluster_terms' > 1 & ///
        __cppmlhdfe_cluster_terms_done != `num_cluster_terms' {
        di as error "cppmlhdfe: rebuild the ctools plugin to enable multiway clustering"
        local reg_rc = 601
    }
    if `reg_rc' {
        * a fit the plugin retained for the bootstrap is released too
        if `__bsrun' | `__wild' capture _cppmlhdfe_bs_release
        _ctools_cleanup cppmlhdfe
        di as error "Error in PPML regression (rc=`reg_rc')"
        exit `reg_rc'
    }

    * vce(bootstrap): the replicates on the final fit, in chunks
    if `__wild' {
        capture noisily _cppmlhdfe_wild, reps(`wb_reps') type(`wb_type') ptype(`wb_ptype_code') ///
            enumerate(`wb_enumerate') ci(`wb_ci') level(`wb_level') dots(`wb_dots') ///
            names(`coef_names') coefs(`wb_coefs') seed(`"`wb_seed'"') threads(`"`threads_code'"') ///
            vcevar(`wb_vcevar') bootvar(`wb_bootvar') keep(`__bsrun')
        local wb_rc = _rc
        if `wb_rc' {
            capture _cppmlhdfe_bs_release
            _ctools_cleanup cppmlhdfe
            exit `wb_rc'
        }
        local wb_rngstate `"`r(rngstate)'"'
        local wb_test `r(test)'
        local wb_B = r(reps)
        local wb_enum = r(enumerated)
        local wb_C = r(C)
    }
    local fe_se_tmp
    local fe_nrep_tmp
    local fe_dse_tmp
    * (the retained fit is released if these variables cannot be created)
    capture noisily {
    if `__fe' {
        foreach name of local fe_se_names {
            tempvar t
            quietly generate double `t' = .
            local fe_se_tmp `fe_se_tmp' `t'
        }
        foreach name of local fe_nrep_names {
            tempvar t
            quietly generate long `t' = .
            local fe_nrep_tmp `fe_nrep_tmp' `t'
        }
        foreach name of local fe_lb_names {
            tempvar t
            quietly generate double `t' = .
            local fe_lb_tmp `fe_lb_tmp' `t'
        }
        foreach name of local fe_ub_names {
            tempvar t
            quietly generate double `t' = .
            local fe_ub_tmp `fe_ub_tmp' `t'
        }
        if "`fe_dse'" != "" {
            tempvar fe_dse_tmp
            quietly generate double `fe_dse_tmp' = .
        }
    }
    }
    if _rc {
        local fe_rc = _rc
        capture _cppmlhdfe_bs_release
        _ctools_cleanup cppmlhdfe
        exit `fe_rc'
    }
    if `__bsrun' {
        local bs_savetmp
        if `"`bs_save_file'"' != "" tempfile bs_savetmp
        capture noisily _cppmlhdfe_bootstrap, reps(`bs_reps') seed(`"`bs_seed'"') ///
            draws(`bs_draws') dots(`bs_dots') level(`bs_level') mse(`bs_mse') ///
            threads(`"`threads_code'"') coef(`__bs') fe(`=`__fe' & `fe_any'') ///
            fevars(`fe_se_tmp' `fe_nrep_tmp' `fe_lb_tmp' `fe_ub_tmp' `fe_dse_tmp') nfe(`nfe') ///
            nreps(`=`__fe' & `fe_nreps'') ci(`=`__fe' & `fe_ci'') dse(`=("`fe_dse'" != "")') ///
            bsdraws(`bs_bsdraws') savefile(`"`bs_savetmp'"') every(`bs_save_every') ///
            `=cond("`bs_strata'" != "" & "`bs_unitvar'" == "", "strataobs", "")'
        local bs_rc = _rc
        local bs_rngstate `"`r(rngstate)'"'
        if !`bs_rc' {
            local bs_nunits = scalar(__cppmlhdfe_bs_U)
            local bs_nstrata = scalar(__cppmlhdfe_bs_S)
        }
        if `bs_rc' {
            capture _cppmlhdfe_bs_release
            _ctools_cleanup cppmlhdfe
            exit `bs_rc'
        }
    }

    * Record plugin call time
    _ctools_lap __time_plugin_call __ct_lap98 `__ct_lap98'

    * The end of ppmlhdfe's iteration log, before the covariance warning
    * (as there, verbose(-1) hides it)
    if __cppmlhdfe_irls_converged {
        if `verbose' > -1 & "`log'" != "nolog" di as text "{hline 108}"
        if `verbose' > -1 {
            di as text "Converged in " as result __cppmlhdfe_irls_iterations as text " iterations and " ///
                as result __cppmlhdfe_ic2 as text " HDFE sub-iterations (tol =" %4.0e `irlstolerance' ")"
        }
    }
    else {
        di as error "(IRLS did NOT converge in " as result __cppmlhdfe_irls_iterations as text " iterations)"
    }

    * Match the reference's PSD correction on the full standardized covariance
    * (including its implicit constant), before extracting and rescaling slopes.
    mata: _cppmlhdfe_fix_psd()

    * Retrieve results
    local N_final = __cppmlhdfe_N
    local num_singletons = __cppmlhdfe_num_singletons
    local num_separated = __cppmlhdfe_num_separated
    local K_keep = __cppmlhdfe_K_keep
    local df_a = __cppmlhdfe_df_a
    local mobility_groups = __cppmlhdfe_mobility_groups
    local deviance = __cppmlhdfe_deviance
    local ll = __cppmlhdfe_ll
    local ll_0 = __cppmlhdfe_ll_0
    local irls_iters = __cppmlhdfe_irls_iterations
    local irls_converged = __cppmlhdfe_irls_converged

    if !`has_absorb' local fe_name_1 _cons
    local df_a_nested = __cppmlhdfe_df_nested
    local df_a_adjusted = `df_a'
    * The extended absorbed effects, one per plugin factor (e(dof_table) rows)
    local fe_names
    forvalues g=1/`nfe' {
        local fe_names `fe_names' `fe_name_`g''
    }

    * Initialize model degrees of freedom; refine to joint-test rank below.
    local df_m = `K_keep'

    * df_r calculation
    local df_r_ols_base = `N_final' - `K_keep' - `df_a_adjusted'
    if `vcetype' == 2 {
        local num_clusters = __cppmlhdfe_N_clust
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
        stub(__cppmlhdfe) kkeep(`K_keep') vmatrix(__cppmlhdfe_V_full) cons(__cppmlhdfe_cons) kept(`b_kept')
    if `__bs' {
        * e(V_modelbased) is the fit's own covariance; e(V) and the bootstrap
        * vectors take the same layout (zeros for omitted and base terms)
        tempname V_mb __bx __bd __D
        matrix `V_mb' = `V'
        local __ebv names(`coef_names_all') base(`coef_is_base') stub(__cppmlhdfe) ///
            kkeep(`K_keep') cons(__cppmlhdfe_cons)
        _ctools_ebv `__bx' `V', `__ebv' vmatrix(__cppmlhdfe_bs_V)
        foreach __v in b se bias z0 {
            tempname bs_`__v'
            matrix `__D' = diag(__cppmlhdfe_bs_`__v')
            _ctools_ebv `__bx' `__bd', `__ebv' vmatrix(`__D')
            matrix `bs_`__v'' = vecdiag(`__bd')
        }
        foreach __v in ci_normal ci_percentile ci_bc {
            tempname bs_`__v'
            forvalues __r = 1/2 {
                matrix `__D' = diag(__cppmlhdfe_bs_`__v'[`__r', 1...])
                _ctools_ebv `__bx' `__bd', `__ebv' vmatrix(`__D')
                matrix `bs_`__v'' = nullmat(`bs_`__v'') \ vecdiag(`__bd')
            }
        }
        * As _bs_sum: a term that is zero in every replicate has no
        * bias-corrected interval or z0
        forvalues __j = 1/`=colsof(`b')' {
            if `V'[`__j', `__j'] == 0 & `b'[1, `__j'] == 0 {
                matrix `bs_z0'[1, `__j'] = .
                matrix `bs_ci_bc'[1, `__j'] = .
                matrix `bs_ci_bc'[2, `__j'] = .
            }
        }
    }

    * saving(): the replicate file's column of each e(b) term (0 for terms
    * fixed at zero), by the layout of e(b)
    if `"`bs_save_file'"' != "" {
        tempname bs_map __smv __smx __smV
        local __Kb = `K_keep' + 1
        matrix `__smv' = J(1, `__Kb', 0)
        forvalues __k = 1/`__Kb' {
            matrix `__smv'[1, `__k'] = `__k'
        }
        matrix `__smv' = diag(`__smv')
        _ctools_ebv `__smx' `__smV', names(`coef_names_all') base(`coef_is_base') stub(__cppmlhdfe) ///
            kkeep(`K_keep') cons(__cppmlhdfe_cons) vmatrix(`__smv')
        matrix `bs_map' = vecdiag(`__smV')
    }

    * Wald statistic on the retained slope coefficients
    local F = .
    if `K_keep' > 0 {
        tempname V_test Vinv Wald
        if `__bs' matrix `V_test' = __cppmlhdfe_bs_V[1..`K_keep', 1..`K_keep']
        else matrix `V_test' = __cppmlhdfe_V_full[1..`K_keep', 1..`K_keep']
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
    * Prepare outputs under temporary names. Existing automatic targets
    * remain intact until all estimation, posting and display work succeeds.
    local publish_sources
    local publish_targets
    if "`d'" != "" {
        label variable `saved_d' "Sum of fixed effects"
        local publish_sources `saved_d'
        local publish_targets `d'
    }
    forvalues g=1/`nfe' {
        if "`fe_target_`g''"!="" {
            label variable `fe_output_`g'' "[FE] `fe_name_`g''"
            local publish_sources `publish_sources' `fe_output_`g''
            local publish_targets `publish_targets' `fe_target_`g''
        }
    }
    if `__fe' {
        local j 0
        forvalues g = 1/`nfe' {
            if "`fe_target_`g''" == "" continue
            local ++j
            local t : word `j' of `fe_se_tmp'
            label variable `t' "Bootstrap SE of FE for `fe_name_`g''"
            local publish_sources `publish_sources' `t'
            local publish_targets `publish_targets' `: word `j' of `fe_se_names''
            if `fe_nreps' {
                local t : word `j' of `fe_nrep_tmp'
                label variable `t' "Replicates estimating FE for `fe_name_`g''"
                local publish_sources `publish_sources' `t'
                local publish_targets `publish_targets' `: word `j' of `fe_nrep_names''
            }
            if `fe_ci' {
                local t : word `j' of `fe_lb_tmp'
                label variable `t' "Bootstrap lower bound of FE for `fe_name_`g''"
                local publish_sources `publish_sources' `t'
                local publish_targets `publish_targets' `: word `j' of `fe_lb_names''
                local t : word `j' of `fe_ub_tmp'
                label variable `t' "Bootstrap upper bound of FE for `fe_name_`g''"
                local publish_sources `publish_sources' `t'
                local publish_targets `publish_targets' `: word `j' of `fe_ub_names''
            }
        }
        if "`fe_dse'" != "" {
            label variable `fe_dse_tmp' "Bootstrap SE of the sum of fixed effects"
            local publish_sources `publish_sources' `fe_dse_tmp'
            local publish_targets `publish_targets' `fe_dse'
        }
    }
    if "`tagsep'"!="" {
        local publish_sources `publish_sources' `separation_tag'
        local publish_targets `publish_targets' `tagsep'
    }
    if "`zvarname'"!="" {
        local publish_sources `publish_sources' `separation_certificate'
        local publish_targets `publish_targets' `zvarname'
    }

    * Record results building time
    _ctools_lap __time_build __ct_lap98 `__ct_lap98'

    * Post results. With factor variables, buildfvinfo lets margins check
    * estimability, as after ppmlhdfe.
    * buildfvinfo costs a pass over the data and is only
    * needed when a coefficient other than a base level is omitted; with a
    * full-rank design every margin is estimable.
    local bfv
    _ms_op_info `b'
    if r(fvops) {
        local __bfv_names : colnames `b'
        foreach __bfv_name of local __bfv_names {
            if strpos("`__bfv_name'", "o.") & !strpos("`__bfv_name'", "b.") local bfv buildfvinfo
        }
    }
    * Posting and publication of the automatic outputs are one transaction:
    * a Break cannot separate the new e() from the outputs it refers to
    nobreak {
    ereturn post `b' `V', esample(`final_sample') depname(`depvar_orig') obs(`N_final') `bfv'

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
    ereturn local title = cond(`has_absorb', "HDFE PPML regression", "PPML regression")
    ereturn local marginsnotok "stdp Anscombe Cooksd Deviance Hat Likelihood Pearson Response Score Working ADJusted STAndardized STUdentized MODified"
    ereturn local separation "`separation'"
    ereturn local properties "b V"
    ereturn local marginsok default
    ereturn scalar rank = `K_keep'
    ereturn scalar num_singletons = `num_singletons'
    ereturn scalar num_separated = `num_separated'
    ereturn scalar ic = `irls_iters'
    ereturn scalar converged = `irls_converged'
    ereturn scalar df_a = `df_a_adjusted'
    ereturn scalar N_full = `N_final' + `num_singletons' + `num_separated'
    ereturn scalar drop_singletons = ("`keepsingletons'" == "")
    ereturn scalar df_a_nested = `df_a_nested'

    * The absorbed effects, the absorbed degrees of freedom table, the
    * covariance labels and the clusters, as reghdfe posts them
    local vcename = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", "cluster"))
    local cluster_post
    if `vcetype' == 2 {
        local nclust
        forvalues dim = 1/`num_cluster_vars' {
            local nclust `nclust' `=__cppmlhdfe_N_clust_`dim''
        }
        local cluster_post clustvars(`clustervar_orig') nclust(`nclust')
    }
    _ctools_hdfe_post, stub(__cppmlhdfe) absvars(`absorb_orig') extended(`fe_names') ///
        nhdfe(`n_absvars') vce(`vcename') dofmethod(`dofmethods') `cluster_post'
    * As reghdfe: the heading of a clustered covariance
    if `vcetype' == 2 ereturn local title3 "Statistics robust to heteroskedasticity"
    * condition: ppmlhdfe posts reghdfe's request for the finite condition
    * number, -1, which its PPML fit never replaces by the number
    if "`condition'" != "" ereturn scalar finite_condition = -1

    ereturn local absorb "`absorb_orig'"
    ereturn local d "`d'"
    ereturn local offset "`exposure'`offset'"
    ereturn local exposure "`exposure'"
    if "`exposure'" != "" ereturn local offset "ln(`exposure')"
    ereturn local depvar "`depvar_orig'"
    ereturn local wtype "`weight'"
    ereturn local wexp "`exp'"
    * As ppmlhdfe: the estimated regressors, factor levels in bn form,
    * without base levels and omitted terms, then _cons
    local __iv
    local __bnames : colnames e(b)
    foreach __name of local __bnames {
        if "`__name'" == "_cons" continue
        _ms_parse_parts `__name'
        if r(omit) continue
        _cppmlhdfe_addbn `__name'
        local __iv `__iv' `__bn'
    }
    if "`__iv'" == "" ereturn local indepvars " _cons"
    else ereturn local indepvars `"`__iv'  _cons"'
    ereturn local predict "cppmlhdfe_p"
    ereturn local cmd "cppmlhdfe"
    ereturn local cmdline `"cppmlhdfe `cmdline_typed'"'
    * As ppmlhdfe (from reghdfe's solve of the final working regression): the
    * weighted residual sum of squares and its root mean square on
    * N - df_a - rank - df_a_nested, the partialling subiterations, and the
    * replay helpers (ppmlhdfe names reghdfe's, which refuse other commands)
    local rss_df = `N_final' - `df_a_adjusted' - `K_keep' - `df_a_nested'
    ereturn scalar rss = __cppmlhdfe_rss
    ereturn scalar rmse = cond(`rss_df' == 0, sqrt(__cppmlhdfe_rss), sqrt(__cppmlhdfe_rss / `rss_df'))
    ereturn scalar ic2 = __cppmlhdfe_ic2
    if `has_absorb' ereturn local title2 = "Absorbing `=e(N_hdfe)' HDFE " + cond(e(N_hdfe) > 1, "groups", "group")
    ereturn local footnote "_cppmlhdfe_footnote"
    ereturn local estat_cmd "_cppmlhdfe_estat"
    if `__wild' {
        _cppmlhdfe_wild_post, test(`wb_test') names(`coef_names') type(`wb_type') ///
            ptype(`wb_ptype') reps(`wb_B') enum(`wb_enum') nclust(`wb_C') level(`wb_level') ///
            ci(`wb_ci') bootclust(`wb_bootname') rngstate(`"`wb_rngstate'"')
    }
    if `__fe' {
        * the normalization concerns the absorbed intercepts
        local fe_nint = 0
        forvalues g = 1/`group_count' {
            if "`slope_var_`g''" == "" local ++fe_nint
        }
        _cppmlhdfe_fe_post, sevars(`fe_se_names') nrepvars(`fe_nrep_names') dse(`fe_dse') nint(`fe_nint') ///
            lbvars(`fe_lb_names') ubvars(`fe_ub_names') level(`bs_level') ///
            absvars(`fe_names') table(`fe_any') rngstate(`"`bs_rngstate'"') draws(`bs_draws') ///
            cluster(`bs_cluster_name') strata(`bs_strata')
    }
    if `__bs' {
        _cppmlhdfe_bootstrap_post, vmb(`V_mb') bsb(`bs_b') bsse(`bs_se') bsbias(`bs_bias') ///
            bsz0(`bs_z0') cinormal(`bs_ci_normal') cipercentile(`bs_ci_percentile') ///
            cibc(`bs_ci_bc') level(`bs_level') rngstate(`"`bs_rngstate'"') ///
            cluster(`bs_cluster_name') strata(`bs_strata') draws(`bs_draws') ///
            command(`"cppmlhdfe `cmdline_typed'"') mse(`bs_mse')
    }
    if `"`bs_save_file'"' != "" {
        local __nclust
        local __nstrata
        if "`bs_unitvar'" != "" local __nclust `bs_nunits'
        if "`bs_strata'" != "" local __nstrata `bs_nstrata'
        capture noisily _ctools_hdfe_boot savefinish, tmp(`"`bs_savetmp'"') file(`"`bs_save_file'"') ///
            map(`bs_map') `bs_save_double' `bs_save_replace' command(`"cppmlhdfe `cmdline_typed'"') ///
            seed(`"`bs_rngstate'"') n(`N_final') cluster(`bs_cluster_name') nclust(`__nclust') ///
            strata(`bs_strata') nstrata(`__nstrata') cmdname(cppmlhdfe)
        if _rc {
            local __src = _rc
            ereturn clear
            exit `__src'
        }
    }

    * Publication follows the last posted result and precedes the display
    if "`publish_targets'" != "" {
        capture noisily _cppmlhdfe_publish `publish_sources', names(`publish_targets')
        local publish_rc = _rc
        if `publish_rc' {
            ereturn clear
            exit `publish_rc'
        }
    }
    }

    * Record ereturn posting time
    _ctools_lap __time_ereturn __ct_lap98 `__ct_lap98'

    * Display
    if "$CTOOLS_cppmlhdfe_fault" == "display" error 1
    _cppmlhdfe_replay, `diopts' `eform' `irr' `nowildtable'

    _ctools_cleanup cppmlhdfe

    * Display timing if verbose
    _ctools_lap __time_display __ct_lap98 `__ct_lap98'
    local __time_total = `__time_setup' + `__time_plugin_call' + `__time_build' + `__time_ereturn' + `__time_display'

    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "cppmlhdfe timing breakdown:"
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

* The display of a fit, at estimation and on replay, as ppmlhdfe's Replay:
* the header (unless noheader), the coefficient table (unless notable; eform
* or irr shows exp(b)), the notes of the bootstraps, and the absorbed degrees
* of freedom (not without absorb()). Returns _coef_table's r().
program define _cppmlhdfe_replay, rclass
    version 14.1
    syntax [, EForm IRr noHEADer noTABLE NOWILDtable *]
    _get_diopts diopts, `options'
    if "`irr'" != "" local eform eform
    if "`header'" == "" _cppmlhdfe_header
    if "`table'" == "" {
        _coef_table, `diopts' `eform'
        return add
        if "`e(vce)'" == "bootstrap" & e(N_misreps) > 0 & e(N_misreps) < . {
            di as text "Note: One or more parameters could not be estimated in " ///
                as result e(N_misreps) as text " bootstrap replicates;"
            di as text "      standard-error estimates include only complete replications."
        }
    }
    if "`e(fevce)'" == "bootstrap" & e(fe_norm_warning) == 1 _cppmlhdfe_fe_warning
    if "`e(wb_method)'" != "" & "`nowildtable'" == "" _ctools_wild_table
    _cppmlhdfe_footnote
end

* ppmlhdfe's header (ppmlhdfe_header), and the replications of
* vce(bootstrap)
program define _cppmlhdfe_header
    version 14.1
    di as text ""
    di as text "`e(title)'" _c
    di as text _col(51) "No. of obs" _col(67) "=" _c
    di as text _col(69) as res %10.0gc e(N)
    di as text "`e(title2)'" _c
    di as text _col(51) "Residual df" _col(67) "=" _c
    di as text _col(69) as res %10.0gc e(df)
    di as text "`e(title3)'" _c
    di as text _col(51) "Wald chi2(" as res e(df_m) as text ")" _col(67) "=" _c
    di as text _col(69) as res %10.2f e(chi2)
    di as text "Deviance" _col(22) "=" _c
    di as text _col(24) as res %12.0g e(deviance) _c
    di as text _col(51) "Prob > chi2" _col(67) "=" _c
    di as text _col(69) as res %10.4f chi2tail(e(df_m), e(chi2))
    di as text "Log pseudolikelihood =" _c
    di as text _col(24) as res %12.0g e(ll) _c
    di as text _col(51) "Pseudo R2" _col(67) "=" _c
    di as text _col(69) as res %10.4f e(r2_p)
    if "`e(vce)'" == "bootstrap" {
        di as text _col(51) "Replications" _col(67) "=" _c
        di as text _col(69) as res %10.0fc e(N_reps)
    }
    local N_clustervars = e(N_clustervars)
    if (`N_clustervars' == .) local N_clustervars 0
    if (`N_clustervars' > 0) di as text ""
    forvalues i = 1/`N_clustervars' {
        local cluster "`e(clustvar`i')'"
        local num = "`e(N_clust`i')'"
        di as text "Number of clusters (" as res "`cluster'" as text ")" _col(29) as text "=" _c
        di as text _col(31) as res %10.0fc `num'
    }
end

* cppmlhdfe, version: as ppmlhdfe, version, the command's location and
* version, then the ctools package and plugin it needs
program define _cppmlhdfe_version
    version 14.1
    which cppmlhdfe
    di as text _n "Required packages installed?"
    ctools, version
end

* cppmlhdfe, reload: ppmlhdfe's reload reinstalls ppmlhdfe from its author's
* working copy. ctools is reinstalled with ctools, update; the programs then
* in memory are replaced after program drop _all (or discard).
program define _cppmlhdfe_reload
    version 14.1
    di as text "{bf:cppmlhdfe} is part of ctools: type {stata ctools, update} to reinstall it."
    di as text _n "{bf:Note:} You need to run {stata program drop _all} now."
end

* A factor-variable term in the bn form of ppmlhdfe's e(indepvars)
* (ms_expand_varlist's AddBN): 2.x becomes 2bn.x, 1.a#2.b becomes
* 1bn.a#2bn.b. Returns the caller's local __bn.
program define _cppmlhdfe_addbn
    version 14.1
    local part `0'
    local re "^([0-9]+)b?([.LFSD])"
    if regexm("`part'", `"`re'"') {
        local part = regexr("`part'", "`re'", regexs(1) + "bn" + regexs(2))
    }
    local re "#([0-9]+)b?([.LFSD])"
    local loop = strpos("`part'", "#")
    local old `part'
    while `loop' {
        if regexm("`part'", `"`re'"') {
            local part = regexr("`part'", "`re'", "#" + regexs(1) + "bn" + regexs(2))
        }
        local loop = "`old'" != "`part'"
        local old `part'
    }
    c_local __bn `part'
end

* vce(bootstrap): run the replicates on the fit the plugin retained, in
* chunks of about 50 (dots and Break between chunks; a multiple of the
* replicates the plugin runs at once, __cppmlhdfe_bs_par, so that no thread
* idles at a chunk's end), then summarize them into the
* __cppmlhdfe_bs_* matrices. Returns r(rngstate), the random-number state
* before the first draw (docs/development/BOOTSTRAP_PLAN.md, 5.8).
program define _cppmlhdfe_bootstrap, rclass
    version 14.1
    syntax, REPS(integer) DRAWS(name) DOTS(integer) LEVEL(real) MSE(integer) ///
        [SEED(string) THREADS(string) COEF(integer 1) FE(integer 0) FEVARS(string) ///
        NFE(integer 1) NREPS(integer 0) CI(integer 0) DSE(integer 0) BSDRAWS(string) ///
        SAVEFILE(string) EVERY(string) STRATAOBS]
    local targets = cond(`fe', " fe", "") + cond(`dse', " dse", "") + cond(`ci', " ci", "")
    local Kb = scalar(__cppmlhdfe_K_keep) + 1
    _ctools_hdfe_boot run, cmd(cppmlhdfe) reps(`reps') draws(`draws') dots(`dots') kb(`Kb') ///
        seed(`"`seed'"') threads(`"`threads'"') targets(`"`targets'"') bsdraws(`bsdraws') ///
        savefile(`"`savefile'"') every(`every') `strataobs'
    return local rngstate `"`r(rngstate)'"'
    matrix __cppmlhdfe_bs_V = J(`Kb', `Kb', .)
    foreach m in b se bias z0 {
        matrix __cppmlhdfe_bs_`m' = J(1, `Kb', .)
    }
    foreach m in ci_normal ci_percentile ci_bc {
        matrix __cppmlhdfe_bs_`m' = J(2, `Kb', .)
    }
    matrix __cppmlhdfe_bs_counts = J(1, 5, .)
    if `fe' matrix __cppmlhdfe_bs_fe_table = J(`nfe', 5, .)
    local finish = cond(`coef', " coef", "") + cond(`nreps', " nreps", "") + cond(`ci', " ci", "") + ///
        cond(`mse', " mse", "")
    plugin call ctools_plugin `fevars', ///
        "cppmlhdfe `threads' boot_finish level=`level'`finish'"
end

* wildbootstrap or wildbootstrap(wildopts) among the remaining options:
* caller locals __wild (0/1), __wildspec and __wild_rest
program define _cppmlhdfe_wild_opt
    version 14.1
    * syntax leaves a bare wildbootstrap among the other options (*), so
    * it is looked for there in a second pass
    syntax [, WILDbootstrap(string asis) *]
    local spec = strtrim(`"`wildbootstrap'"')
    local 0 `", `options'"'
    syntax [, WILDbootstrap *]
    c_local __wild = (`"`spec'"' != "" | "`wildbootstrap'" != "")
    c_local __wildspec `"`spec'"'
    c_local __wild_rest `"`options'"'
end

* wildbootstrap(): the score bootstrap on the fit the plugin retained. The
* weights come from Stata's generator as boottest draws them (enumerated in
* the plugin when all 2^C Rademacher vectors fit in reps). The analytic
* standard errors set the scale of the interval search. Results: matrices
* __cppmlhdfe_wb_p, _fail, _ci, _cistat; r(test), r(reps), r(enumerated),
* r(C), r(rngstate) (before the first draw).
program define _cppmlhdfe_wild, rclass
    version 14.1
    syntax, REPS(integer) TYPE(name) PTYPE(integer) ENUMERATE(integer) CI(integer) ///
        LEVEL(real) DOTS(integer) [NAMES(string) COEFS(string) SEED(string) THREADS(string) ///
        VCEVAR(string) BOOTVAR(string) KEEP(integer 0)]
    local K : word count `names'
    local test
    if `"`coefs'"' == "" {
        forvalues k = 1/`K' {
            if scalar(__cppmlhdfe_collinear_`k') != 1 local test `test' `k'
        }
    }
    else {
        foreach c of local coefs {
            local p : list posof `"`c'"' in names
            if !`p' {
                di as error "cppmlhdfe: wildbootstrap(coefs()): `c' is not a regressor"
                exit 198
            }
            if scalar(__cppmlhdfe_collinear_`p') == 1 {
                di as error "cppmlhdfe: wildbootstrap(coefs()): `c' was omitted"
                exit 198
            }
            local test `test' `p'
        }
    }
    local T : word count `test'
    if `T' == 0 {
        di as error "cppmlhdfe: wildbootstrap() has no coefficient to test"
        exit 198
    }
    * the analytic standard error of each tested coefficient (kept position
    * kk among the plugin's columns, in the data's units)
    matrix __cppmlhdfe_wb_test = J(1, `T', .)
    matrix __cppmlhdfe_wb_scale = J(1, `T', .)
    forvalues j = 1/`T' {
        local p : word `j' of `test'
        local kk 0
        forvalues k = 1/`p' {
            if scalar(__cppmlhdfe_collinear_`k') != 1 local ++kk
        }
        matrix __cppmlhdfe_wb_test[1, `j'] = `p'
        matrix __cppmlhdfe_wb_scale[1, `j'] = sqrt(__cppmlhdfe_V_full[`kk', `kk']) / __cppmlhdfe_scales[1, `kk']
    }
    if `"`seed'"' != "" set seed `seed'
    return local rngstate `"`c(rngstate)'"'
    local flags = cond("`vcevar'" != "", " vce", "") + cond("`bootvar'" != "", " boot", "")
    plugin call ctools_plugin `vcevar' `bootvar', "cppmlhdfe `threads' wild_setup`flags'"
    local C = scalar(__cppmlhdfe_wb_C)
    * Any failure or Break while the weights are drawn and passed leaves
    * neither their variables nor the Mata draws behind
    capture noisily {
        _ctools_hdfe_boot wilddraw, clusters(`C') reps(`reps') type(`type') ///
            `=cond(`enumerate', "", "noenumerate")'
        local B = r(reps)
        local enum = r(enumerate)
        * each chunk's weights, packed into variables reused across chunks
        local owned_weights
        if !`enum' {
            forvalues j = 1/`=ceil(min(50, `B') * `C' / _N)' {
                tempvar wild_weight
                quietly generate double `wild_weight' = .
                local owned_weights `owned_weights' `wild_weight'
            }
        }
        local first 0
        while `first' < `B' {
            local n = min(50, `B' - `first')
            if `enum' {
                plugin call ctools_plugin, ///
                    "cppmlhdfe `threads' wild_chunk first=`first' count=`n' total=`B' enum"
            }
            else {
                _ctools_hdfe_boot wildweights, first(`=`first' + 1') count(`n') variables(`owned_weights')
                local wvars `r(varlist)'
                capture noisily plugin call ctools_plugin `wvars', ///
                    "cppmlhdfe `threads' wild_chunk first=`first' count=`n' total=`B'"
                local rc = _rc
                if `rc' {
                    capture drop `owned_weights'
                    _ctools_hdfe_boot wildclear
                    exit `rc'
                }
            }
            local first = `first' + `n'
        }
    }
    local rc = _rc
    if `rc' {
        capture drop `owned_weights'
        capture _ctools_hdfe_boot wildclear
        exit `rc'
    }
    capture drop `owned_weights'
    _ctools_hdfe_boot wildclear
    if `dots' di as text "(score bootstrap: `B' replications, " cond(`ci', "with", "without") " intervals)"
    matrix __cppmlhdfe_wb_p = J(1, `T', .)
    matrix __cppmlhdfe_wb_score = J(1, `T', .)
    matrix __cppmlhdfe_wb_fail = J(1, `T', 0)
    matrix __cppmlhdfe_wb_ci = J(`T', 2, .)
    matrix __cppmlhdfe_wb_cistat = J(`T', 5, 0)
    plugin call ctools_plugin, "cppmlhdfe `threads' wild_finish level=`level' ptype=`ptype'`=cond(`ci', " ci", "")'`=cond(`keep', " keep", "")'"
    return local test `test'
    return scalar reps = `B'
    return scalar enumerated = `enum'
    return scalar C = `C'
end

* wildbootstrap(): its e() results (plan 4.7), with the prefix wb_
program define _cppmlhdfe_wild_post, eclass
    version 14.1
    syntax, TEST(numlist) NAMES(string) TYPE(name) PTYPE(name) REPS(integer) ENUM(integer) ///
        NCLUST(integer) LEVEL(real) CI(integer) [BOOTCLUST(string) RNGSTATE(string)]
    local coefs
    foreach p of local test {
        local coefs `coefs' `: word `p' of `names''
    }
    local T : word count `test'
    tempname t pv fl cim stm
    matrix `t' = J(1, `T', .)
    forvalues j = 1/`T' {
        local c : word `j' of `coefs'
        matrix `t'[1, `j'] = _b[`c'] / _se[`c']
    }
    tempname sc
    matrix `pv' = __cppmlhdfe_wb_p
    matrix `fl' = __cppmlhdfe_wb_fail
    matrix `sc' = __cppmlhdfe_wb_score
    matrix colnames `t' = `coefs'
    matrix colnames `pv' = `coefs'
    matrix colnames `fl' = `coefs'
    matrix colnames `sc' = `coefs'
    ereturn matrix wb_t = `t'
    ereturn matrix wb_p = `pv'
    ereturn matrix wb_fail = `fl'
    * the score statistic the p-values compare the replicates' to (e(wb_t)
    * is the Wald z of e(b) and e(V))
    ereturn matrix wb_score = `sc'
    if `ci' {
        matrix `cim' = __cppmlhdfe_wb_ci
        matrix `stm' = __cppmlhdfe_wb_cistat
        matrix rownames `cim' = `coefs'
        matrix colnames `cim' = ll ul
        matrix rownames `stm' = `coefs'
        matrix colnames `stm' = open_lower open_upper disjoint empty failed
        ereturn matrix wb_ci = `cim'
        ereturn matrix wb_cistatus = `stm'
    }
    ereturn scalar wb_reps = `reps'
    ereturn scalar wb_misreps = scalar(__cppmlhdfe_wb_misreps)
    ereturn scalar wb_enumerated = `enum'
    ereturn scalar wb_N_clust = `nclust'
    ereturn scalar wb_level = `level'
    ereturn local wb_method score
    ereturn local wb_weighttype `type'
    ereturn local wb_ptype `ptype'
    ereturn local wb_coefs `coefs'
    if "`bootclust'" != "" ereturn local wb_bootclust `bootclust'
    ereturn local wb_rngstate `"`rngstate'"'
end

* fevce(bootstrap): e() for the fixed effects' bootstrap (plan 4.4)
program define _cppmlhdfe_fe_post, eclass
    version 14.1
    syntax, ABSVARS(string) TABLE(integer) [SEVARS(string) NREPVARS(string) DSE(string) ///
        LBVARS(string) UBVARS(string) LEVEL(real 95) NINT(integer -1) ///
        RNGSTATE(string) DRAWS(name) CLUSTER(string) STRATA(string)]
    if `nint' < 0 local nint : word count `absvars'
    if `table' {
        tempname t
        matrix `t' = __cppmlhdfe_bs_fe_table
        matrix colnames `t' = levels never always min_reps median_reps
        matrix rownames `t' = `absvars'
        ereturn matrix fe_bs_table = `t'
    }
    ereturn scalar fe_N_reps = scalar(__cppmlhdfe_bs_fe_N_reps)
    ereturn scalar fe_split_reps = scalar(__cppmlhdfe_bs_fe_split_reps)
    * the caveat concerns the per-level SE variables (plan 4.8); d() alone
    * does not depend on the normalization
    ereturn scalar fe_norm_warning = (`nint' >= 2 & "`sevars'" != "")
    ereturn local fevce bootstrap
    if "`sevars'" != "" ereturn local fe_se `sevars'
    if "`nrepvars'" != "" ereturn local fe_nrep `nrepvars'
    if "`lbvars'" != "" {
        ereturn local fe_lb `lbvars'
        ereturn local fe_ub `ubvars'
        ereturn scalar fe_level = `level'
    }
    if "`dse'" != "" ereturn local fe_dse `dse'
    ereturn local fe_rngstate `"`rngstate'"'
    ereturn local fe_draws `draws'
    if "`cluster'" != "" ereturn local fe_cluster_bs `cluster'
    if "`strata'" != "" ereturn local fe_strata_bs `strata'
end

* With two or more absorbed effects, the caveat on per-level fixed-effect
* standard errors (plan 4.8)
program define _cppmlhdfe_fe_warning
    version 14.1
    local G : word count `e(absorb)'
    di as text ""
    di as text "{p 0 9 2}warning: with two or more absorbed effects (`e(absorb)'), individual fixed effects are " ///
        "identified only up to a normalization within each connected set of levels (here the minimum norm). " ///
        "Their bootstrap standard errors" _c
    if "`e(fe_se)'" != "" di as text " in `e(fe_se)'" _c
    di as text " depend on that normalization and can be large where the sets are thinly connected. " _c
    if "`e(fe_dse)'" != "" di as text "The sum of the fixed effects does not depend on it; its " ///
        "standard errors are in `e(fe_dse)'.{p_end}"
    else di as text "The sum of the fixed effects does not depend on it: specify d() with " ///
        "fevce(bootstrap) for its standard errors.{p_end}"
    if e(fe_split_reps) > 0 & e(fe_split_reps) < . {
        di as text "{p 9 9 2}In " e(fe_split_reps) " of " e(fe_N_reps) " replicates a connected set " ///
            "split; levels in a split set were treated as not estimable in those replicates.{p_end}"
    }
end

* vce(bootstrap): the e() results Stata's vce(bootstrap) posts, so that
* estat bootstrap works
program define _cppmlhdfe_bootstrap_post, eclass
    version 14.1
    syntax, VMB(name) BSB(name) BSSE(name) BSBIAS(name) BSZ0(name) CINORMAL(name) ///
        CIPERCENTILE(name) CIBC(name) LEVEL(real) COMMAND(string asis) ///
        [RNGSTATE(string) CLUSTER(string) STRATA(string) DRAWS(name) MSE(integer 0)]
    local names : colnames e(b)
    foreach m in `bsb' `bsse' `bsbias' `bsz0' {
        matrix colnames `m' = `names'
        matrix rownames `m' = y1
    }
    foreach m in `cinormal' `cipercentile' `cibc' {
        matrix colnames `m' = `names'
        matrix rownames `m' = ll ul
    }
    matrix colnames `vmb' = `names'
    matrix rownames `vmb' = `names'
    tempname reps counts
    matrix `reps' = J(1, colsof(e(b)), scalar(__cppmlhdfe_bs_N_reps))
    matrix colnames `reps' = `names'
    matrix rownames `reps' = y1
    matrix `counts' = __cppmlhdfe_bs_counts
    matrix colnames `counts' = ok collinear empty noconv other
    ereturn matrix V_modelbased = `vmb'
    ereturn matrix b_bs = `bsb'
    ereturn matrix reps = `reps'
    ereturn matrix bias = `bsbias'
    ereturn matrix se = `bsse'
    ereturn matrix z0 = `bsz0'
    ereturn matrix ci_normal = `cinormal'
    ereturn matrix ci_percentile = `cipercentile'
    ereturn matrix ci_bc = `cibc'
    ereturn matrix boot_status = `counts'
    ereturn scalar N_reps = scalar(__cppmlhdfe_bs_N_reps)
    ereturn scalar N_misreps = scalar(__cppmlhdfe_bs_N_misreps)
    ereturn scalar level = `level'
    ereturn scalar bs_version = 3
    ereturn scalar k_exp = 0
    ereturn scalar k_eexp = colsof(e(b))
    ereturn scalar k_extra = 0
    ereturn scalar k_eq = 1
    local i 0
    foreach name of local names {
        local ++i
        ereturn local exp`i' "_b[`name']"
    }
    ereturn local command `"`command'"'
    ereturn local cmdname cppmlhdfe
    ereturn local rngstate `"`rngstate'"'
    ereturn local draws `draws'
    if "`cluster'" != "" ereturn local cluster_bs `cluster'
    if "`strata'" != "" ereturn local strata_bs `strata'
    ereturn local prefix bootstrap
    ereturn local vce bootstrap
    * as _bs_sum: with mse, e(V) is about e(b) and labelled so
    ereturn local vcetype = cond(`mse', "Bstrap *", "Bootstrap")
    if `mse' ereturn local mse mse
    * the Wald test of e(chi2) under the bootstrap e(V), as bootstrap posts it
    ereturn scalar p = chi2tail(e(df_m), e(chi2))
end

* Release the fit the plugin retained for vce(bootstrap), after a failure
program define _cppmlhdfe_bs_release
    version 14.1
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    plugin call ctools_plugin, "cppmlhdfe boot_release"
end

* Atomically replace automatic outputs without copying the dataset. Group
* renames validate every destination before making any change; backup names
* retain the old values, metadata and variable positions until commit.
program define _cppmlhdfe_publish
    version 14.1
    syntax varlist, NAMES(string)
    local originals
    local backups
    foreach target of local names {
        capture confirm variable `target', exact
        if !_rc {
            tempvar backup
            local originals `originals' `target'
            local backups `backups' `backup'
        }
    }
    nobreak {
        if "`originals'" != "" rename (`originals') (`backups')
        capture noisily {
            if "$CTOOLS_cppmlhdfe_fault" == "publish" error 459
            rename (`varlist') (`names')
        }
        local rc = _rc
        if `rc' {
            if "`originals'" != "" rename (`backups') (`originals')
            exit `rc'
        }
    }
end

capture mata: mata drop _cppmlhdfe_fix_psd()
mata:
void _cppmlhdfe_fix_psd()
{
    real matrix V, U
    real rowvector lambda, scales
    real scalar k
    k = st_numscalar("__cppmlhdfe_K_keep")
    V = st_matrix("__cppmlhdfe_V_full")[|1,1 \ k+1,k+1|]
    V = (V + V') / 2
    if (st_numscalar("__cppmlhdfe_cluster_terms") > 1) {
    symeigensystem(V, U=., lambda=.)
    if (min(lambda) < 0) {
        lambda = lambda :* (lambda :>= 0)
        V = quadcross(U', lambda, U')
        printf("{txt}Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied.\n")
    }
    }
    scales = st_matrix("__cppmlhdfe_scales")[1,1..k+1]
    V = V :/ (scales' * scales)
    st_matrix("__cppmlhdfe_V_full", V)
    if (k) st_matrix("__cppmlhdfe_V", V[|1,1 \ k,k|])
}
end
