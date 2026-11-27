*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! original: ppmlhdfe 2.3.3 (November 2, 2025) by Sergio Correia, Paulo Guimarães, and Thomas Zylkin
*! original: see help cppmlhdfe for details

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
        quietly mark `constant_fe', noby
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

    * Batch plain numeric effects with the model's missing-value check.
    local __absorb_mark = 0
    capture confirm numeric variable `absorb'
    if !_rc local __absorb_mark = 1
    local __model_varlist `varlist'
    if `__absorb_mark' local varlist `varlist' `absorb'
    * Screen plain inputs after loading, also with a weight that is one
    * numeric variable (marksample checks its rows; the plugin reads it as
    * it is). Weight expressions, offsets and factors use a native streaming
    * pass before weight staging and factor expansion.
    local __wbare = 0
    if "`weight'" != "" {
        local __wexp = strtrim(substr(strtrim(`"`exp'"'), 2, .))
        if `: word count `__wexp'' == 1 {
            capture confirm numeric variable `__wexp', exact
            if !_rc local __wbare = 1
        }
    }
    local __ctools_plugin_hdfe_screen = 0
    local __ctools_plugin_hdfe_markout = 0
    _ctools_load
    local __native_screen = 0
    if "`fevce'" == "" & !`__wild' & ///
        substr(lower(strtrim(`"`vce'"')), 1, 4) != "boot" {
        capture confirm numeric variable `__model_varlist'
        if !_rc & `__absorb_mark' & "`offset'`exposure'" == "" & ("`weight'" == "" | `__wbare') ///
            local __native_screen = `__ctools_plugin_hdfe_screen'
        else if `__ctools_plugin_hdfe_markout' {
            _ctools_hdfe_screen plan `__model_varlist'
            if `__screen_ok' local __native_screen = 2
        }
    }
    * Stata can reuse the old e(sample) storage when it is released before
    * marksample allocates a new marker (posting is about 45 ms faster at
    * 10M rows). Keep previous results when the qualifiers, weights or
    * options can refer to them: e() or _b[], _se[], _coef[] as a token
    * (vce( or tolerance( is not one), and with tagsep(), which posts no
    * fit. A fit that fails after this point leaves e() empty, as a failed
    * creghdfe does.
    if "`tagsep'" == "" & ///
        !ustrregexm(`"`cmdline_typed'"', "(^|[^A-Za-z0-9_])(e[(]|_b\[|_se\[|_coef\[)") {
        ereturn clear
    }
    if `__native_screen' marksample touse, novarlist
    else marksample touse
    local varlist `__model_varlist'
    if `__native_screen' == 2 {
        local __early_inputs `__screen_inputs'
        if `__absorb_mark' local __early_inputs `__early_inputs' `absorb'
        _ctools_hdfe_screen mark `__early_inputs', touse(`touse') plugin(`__ctools_plugin') threads(`threads')
    }
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
        if !`__absorb_mark' markout `touse' `categories' `continuous', strok
        * A numeric variable goes to the plugin as it is (the plugin
        * renumbers every effect on the sample); strings and interactions
        * are numbered first
        local fe_direct = 0
        if `: word count `categories'' == 1 {
            capture confirm numeric variable `categories', exact
            local fe_direct = !_rc
        }
        if `fe_direct' local fe_id `categories'
        else {
            tempvar fe_id
            quietly _ctools_group `categories' if `touse', generate(`fe_id')
        }
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
        _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"`exp'"') marked
        local weight_var "`__ctools_wvar'"
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
        * with the native screen the plugin drops missing cluster values
        _ctools_hdfe_vce, spec(`"`vce'"') touse(`touse') cmd(cppmlhdfe) `=cond(`__native_screen' == 1, "nomark", "")'
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

    * The plugin checks y >= 0 on the loaded sample before singleton or
    * separation removal, returning r(198). Do not scan the data twice.

    * Expand factor variables
    local __fvn_ok 0
    local __fvn_kept
    if "`indepvars'" != "" {
        * Plain i.f terms: the plugin builds the indicators from f and the
        * levels come from one plugin pass (_ctools_fvnative), not from
        * fvrevar and fvexpand (about 0.2 s each per term at 10M rows). The
        * bootstraps and fevce() keep fvrevar's variables.
        if `__ctools_plugin_hdfe_fvlevels' & "`fevce'" == "" & !`__wild' & ///
            substr(lower(strtrim(`"`vce'"')), 1, 4) != "boot" & strpos("`indepvars'", ".") {
            _ctools_fvnative `indepvars', touse(`touse') plugin(`__ctools_plugin') threads(`threads')
        }
        if `__fvn_ok' {
            local __indepvars_all `__fvn_vars'
            local coef_names_all `__fvn_names'
            * the level of each estimated term (-1: not an indicator)
            _ctools_fvsplit, names(`coef_names_all') vars(`__fvn_levels')
            local __fvn_kept "`r(vars)'"
        }
        else {
            fvrevar `indepvars' if `touse'
            local __indepvars_all "`r(varlist)'"
            fvexpand `indepvars' if `touse'
            local coef_names_all "`r(varlist)'"
        }

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

    * Only existence is needed here; stop at the first included row.
    capture assert !`touse', fast
    if !_rc {
        di as error "cppmlhdfe: no observations"
        exit 2000
    }

    * Identity and capabilities were checked before sample preparation.
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
    * i.f indicators the plugin builds: the level of each regressor in the
    * plugin varlist (-1: an ordinary variable)
    if "`__fvn_kept'" != "" {
        local __fvk : word count `__fvn_kept'
        matrix __cppmlhdfe_fvlev = J(1, `__fvk', -1)
        forvalues __c = 1/`__fvk' {
            matrix __cppmlhdfe_fvlev[1, `__c'] = `: word `__c' of `__fvn_kept''
        }
        scalar __cppmlhdfe_fvspec = 1
    }
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

    * The plugin returns the final sample after internal row removals: the
    * sample marker itself (the last plugin variable), in which it clears the
    * rows it removes, so no second marker has to be created and filled
    local final_sample `touse'
    scalar __cppmlhdfe_sample_prefilled = 1

    * Optional outputs get a variable only when requested; index 0 tells the
    * plugin to skip that column's per-observation stores.
    * st_addvar fills these with missing, including all excluded rows.
    tempvar saved_d
    scalar __cppmlhdfe_d_idx = 0
    if "`d'" != "" {
        quietly mata: st_addvar("double", "`saved_d'")
        local plugin_varlist `plugin_varlist' `saved_d'
        scalar __cppmlhdfe_d_idx = `: word count `plugin_varlist''
    }

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
            quietly mata: st_addvar("double", "`fe_output_`g''")
            local plugin_varlist `plugin_varlist' `fe_output_`g''
            scalar __cppmlhdfe_save_fe_`g' = `: word count `plugin_varlist''
        }
    }
    if `guess_code'==2 {
        local plugin_varlist `plugin_varlist' `guess_var'
        scalar __cppmlhdfe_guess_idx = `: word count `plugin_varlist''
    }
    tempvar separation_tag separation_certificate
    scalar __cppmlhdfe_tag_idx = 0
    scalar __cppmlhdfe_certificate_idx = 0
    if "`tagsep'" != "" {
        quietly mata: st_addvar("byte", "`separation_tag'")
        local plugin_varlist `plugin_varlist' `separation_tag'
        scalar __cppmlhdfe_tag_idx = `: word count `plugin_varlist''
    }
    if "`zvarname'" != "" {
        quietly mata: st_addvar("double", "`separation_certificate'")
        local plugin_varlist `plugin_varlist' `separation_certificate'
        scalar __cppmlhdfe_certificate_idx = `: word count `plugin_varlist''
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
    scalar __cppmlhdfe_sample_idx = `: word count `plugin_varlist''
    scalar __cppmlhdfe_touse = 1
    scalar __cppmlhdfe_screen = (`__native_screen' == 1)

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
    local refit_pass = 0
    while `refit' {
        scalar __cppmlhdfe_refit = 0
        capture scalar drop __cppmlhdfe_marked
        local ++refit_pass
        capture noisily plugin call ctools_plugin `plugin_varlist', ///
            "cppmlhdfe `threads_code' full_regression`=cond(`__bsrun' | `__wild', " retain", "")'"
        local reg_rc = _rc
        if !`reg_rc' {
            * a plugin that did not update the marker would leave the
            * dropped rows in e(sample)
            capture confirm scalar __cppmlhdfe_marked
            * and one that did not build the i.f indicators would have
            * fitted the factor variables themselves
            if !_rc & "`__fvn_kept'" != "" capture confirm scalar __cppmlhdfe_fvdone
            if _rc | __cppmlhdfe_api != 3 {
                di as error "cppmlhdfe: rebuild the ctools plugin alongside the ado files"
                local reg_rc = 601
            }
        }
        if `reg_rc' continue, break
        * A separation refit only restricts this complete-case sample.
        * The inputs do not change, so it needs no further missing scan.
        scalar __cppmlhdfe_screen = 0
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
            * the refit reads the marker, from which the plugin has cleared
            * the separated rows
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

    * The plugin's own timers (those of its last pass), read before its
    * scalars are dropped; a timer the plugin did not leave is missing
    if `__do_timing' {
        foreach __t in load screen remap singleton dof setup simplex relu start irls irls_working irls_partial ///
            irls_solve irls_update final vce post total {
            local __pt_`__t' .
            capture confirm scalar _cppmlhdfe_time_`__t'
            if !_rc local __pt_`__t' = scalar(_cppmlhdfe_time_`__t')
        }
        foreach __t in iterations subiterations relusteps {
            local __pt_`__t' .
            capture confirm scalar _cppmlhdfe_irls_`__t'
            if !_rc local __pt_`__t' = scalar(_cppmlhdfe_irls_`__t')
        }
        * a plugin without the breakdown leaves only the totals
        if `__pt_start' >= . local __pt_total .
    }
    _ctools_cleanup cppmlhdfe

    * Display timing if verbose
    _ctools_lap __time_display __ct_lap98 `__ct_lap98'
    local __time_total = `__time_setup' + `__time_plugin_call' + `__time_build' + `__time_ereturn' + `__time_display'

    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "cppmlhdfe timing breakdown:"
        di as text "{hline 55}"
        if `__pt_total' < . {
            di as text "  C plugin internals:"
            di as text "    Data load (SPI):        " as result %8.4f `__pt_load' " sec"
            if `__pt_screen' < . di as text "      of it, missing screen:" as result %8.4f `__pt_screen' " sec"
            di as text "    FE remap + counts:      " as result %8.4f `__pt_remap' " sec"
            di as text "    Singleton removal:      " as result %8.4f `__pt_singleton' " sec"
            di as text "    DOF computation:        " as result %8.4f `__pt_dof' " sec"
            di as text "    Compact + standardize:  " as result %8.4f `__pt_setup' " sec"
            di as text "    Separation (simplex):   " as result %8.4f `__pt_simplex' " sec"
            di as text "    Separation (ReLU):      " as result %8.4f `__pt_relu' " sec" ///
                as text cond(`__pt_relusteps' < . & `__pt_relusteps' > 0, ///
                    "  (`__pt_relusteps' step" + cond(`__pt_relusteps' == 1, "", "s") + ")", "")
            di as text "    Start values, collinear:" as result %8.4f `__pt_start' " sec"
            di as text "    IRLS loop:              " as result %8.4f `__pt_irls' " sec" ///
                as text "  (" as result `__pt_iterations' as text " iterations)"
            di as text "      Working variable:     " as result %8.4f `__pt_irls_working' " sec"
            di as text "      HDFE partial out:     " as result %8.4f `__pt_irls_partial' " sec" ///
                as text "  (" as result `__pt_subiterations' as text " subiterations)"
            di as text "      WLS solve:            " as result %8.4f `__pt_irls_solve' " sec"
            di as text "      Update + deviance:    " as result %8.4f `__pt_irls_update' " sec"
            di as text "    Final projection, solve:" as result %8.4f `__pt_final' " sec"
            di as text "    VCE computation:        " as result %8.4f `__pt_vce' " sec"
            di as text "    Results + stored vars:  " as result %8.4f `__pt_post' " sec"
            di as text "  {hline 53}"
            di as text "    C plugin total:         " as result %8.4f `__pt_total' " sec"
            if `refit_pass' > 1 di as text "    (the last of `refit_pass' passes; the earlier ones are in the call overhead)"
            di as text ""
            di as text "  Stata overhead:"
            di as text "    Setup (parsing, etc):   " as result %8.4f `__time_setup' " sec"
            di as text "    Plugin call overhead:   " as result %8.4f (`__time_plugin_call' - `__pt_total') " sec"
            di as text "    Results building:       " as result %8.4f `__time_build' " sec"
            di as text "    ereturn posting:        " as result %8.4f `__time_ereturn' " sec"
            di as text "    Display output:         " as result %8.4f `__time_display' " sec"
            di as text "  {hline 53}"
            di as text "    Stata overhead total:   " as result %8.4f (`__time_total' - `__pt_total') " sec"
        }
        else {
            di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
            di as text "  C plugin total:           " as result %8.4f `__time_plugin_call' " sec"
            di as text "  Results building:         " as result %8.4f `__time_build' " sec"
            di as text "  ereturn posting:          " as result %8.4f `__time_ereturn' " sec"
            di as text "  Display output:           " as result %8.4f `__time_display' " sec"
        }
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
* working copy. ctools is reinstalled with ctools, update force; the programs
* then in memory are replaced after program drop _all (or discard).
program define _cppmlhdfe_reload
    version 14.1
    di as text "{bf:cppmlhdfe} is part of ctools: type {stata ctools, update force} to reinstall it."
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

* The programs below were separate ado files until ctools 1.3.2. Stata's
* net install takes at most 100 files per package, so helpers that only
* this file calls are kept here as its subprograms.

* Private helper for cppmlhdfe.
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

* ctools: i.f indicators built by the plugin instead of fvrevar and fvexpand
* _ctools_fvnative <terms>, touse(varname) plugin(string) [threads(#)]
*
* terms is an independent-variable list as syntax varlist(fv ts) leaves it.
* When every term is a numeric variable or one of i.v, ib#.v, ibn.v,
* ib(first).v and ib(last).v on a numeric v without fvset base (and set
* fvbase is on), and at least one term is such a factor term, the levels of
* the factor variables on the touse rows come from one plugin pass and the
* caller's locals get
*     __fvn_ok      1 (otherwise 0, and nothing else is set)
*     __fvn_names   the names fvexpand <terms> if touse gives
*     __fvn_vars    one variable per name: the term, or its factor variable
*     __fvn_levels  one value per name: its level, or -1 for a plain term
*     __fvn_counts  one value per term: how many names it gave
* Anything else (interactions, c./time-series operators, i(numlist).v,
* ib(freq).v, noninteger or negative values, levels above 65535) leaves the
* caller on Stata's own fvrevar/fvexpand path, which also gives Stata's
* messages. At 10M rows fvrevar and fvexpand of one i.f cost about 0.2 s
* each; the plugin pass reads each factor variable once.
program define _ctools_fvnative
    version 14.1
    syntax anything(name=terms), TOUSE(varname numeric) PLUGIN(string) [THREADS(integer 0)]
    c_local __fvn_ok 0
    if "`c(fvbase)'" == "off" exit
    local kinds
    local tvars
    local fvars
    foreach t of local terms {
        if ustrregexm("`t'", "^i(b[0-9]+|bn|b\(first\)|b\(last\))?\.([A-Za-z_][A-Za-z0-9_]*)$") {
            local spec = ustrregexs(1)
            local v = ustrregexs(2)
            capture confirm numeric variable `v', exact
            if _rc exit
            if `"`: char `v'[_fv_base]'"' != "" exit
            if "`spec'" == "" | "`spec'" == "b(first)" local kind first
            else if "`spec'" == "b(last)" local kind last
            else if "`spec'" == "bn" local kind bn
            else local kind = substr("`spec'", 2, .)
            local kinds `kinds' `kind'
            local tvars `tvars' `v'
            local fvars `fvars' `v'
        }
        else {
            capture confirm numeric variable `t', exact
            if _rc exit
            local kinds `kinds' plain
            local tvars `tvars' `t'
        }
    }
    local fvars : list uniq fvars
    if "`fvars'" == "" exit
    capture program ctools_plugin, plugin using("`plugin'")
    if _rc != 0 & _rc != 110 exit
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    local nf : word count `fvars'
    forvalues j = 1/`nf' {
        local ctools_fvl_`j'
        local ctools_fvl_ok_`j' 0
    }
    capture plugin call ctools_plugin `fvars' `touse', "fvlevels `threadopt'"
    if _rc exit
    forvalues j = 1/`nf' {
        if "`ctools_fvl_ok_`j''" != "1" exit
    }
    local names
    local vars
    local levels
    local i = 0
    foreach t of local terms {
        local ++i
        local kind : word `i' of `kinds'
        if "`kind'" == "plain" {
            local names `names' `t'
            local vars `vars' `t'
            local levels `levels' -1
            local counts `counts' 1
            continue
        }
        local v : word `i' of `tvars'
        local j : list posof "`v'" in fvars
        local all `ctools_fvl_`j''
        local nl : word count `all'
        local base
        if "`kind'" == "first" local base : word 1 of `all'
        else if "`kind'" == "last" local base : word `nl' of `all'
        else if "`kind'" != "bn" {
            * ib#.v: a base that is not a level is still listed, in order
            local base = real("`kind'")
            if missing(`base') | `base' < 0 | `base' > 65535 | `base' != floor(`base') exit
            local base = strofreal(`base', "%12.0f")
            if !`: list base in all' {
                local merged
                local placed = 0
                foreach l of local all {
                    if !`placed' & `l' > `base' {
                        local merged `merged' `base'
                        local placed = 1
                    }
                    local merged `merged' `l'
                }
                if !`placed' local merged `merged' `base'
                local all `merged'
            }
        }
        local first = 1
        foreach l of local all {
            if "`kind'" == "bn" local name = cond(`first', "`l'bn.`v'", "`l'.`v'")
            else local name = cond("`l'" == "`base'", "`l'b.`v'", "`l'.`v'")
            local first = 0
            local names `names' `name'
            local vars `vars' `v'
            local levels `levels' `l'
        }
        local counts `counts' `: word count `all''
    }
    c_local __fvn_names `names'
    c_local __fvn_vars `vars'
    c_local __fvn_levels `levels'
    c_local __fvn_counts `counts'
    c_local __fvn_ok 1
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

* dofadjustments() of the HDFE estimators.
* _ctools_hdfe_dof [, all none firstpair pairwise clusters continuous]
*
* all (the default), none, or any of firstpair|pairwise, clusters and
* continuous, with reghdfe's abbreviations and exclusions; all stands for
* pairwise clusters continuous. Returns
*   s(flags)     pairwise = 1, firstpair = 2, clusters = 4, continuous = 8
*   s(methods)   the selected adjustments in that order, as e(dofmethod)
*   s(type)      0 = all, 1 = none, 2 = firstpair, 3 = pairwise
*   s(clusters)  1 when the cluster-nesting adjustment applies
program define _ctools_hdfe_dof, sclass
    version 14.1
    syntax [, ALL NONE FIRSTpair PAIRwise CLusters CONTinuous]
    if "`none'`firstpair'`pairwise'`clusters'`continuous'" == "" local all all
    opts_exclusive "`all' `none' `firstpair' `pairwise'" dofadjustments
    opts_exclusive "`all' `none' `clusters'" dofadjustments
    opts_exclusive "`all' `none' `continuous'" dofadjustments
    if "`all'" != "" {
        local pairwise pairwise
        local clusters clusters
        local continuous continuous
    }
    sreturn local flags = ("`pairwise'" != "") + 2 * ("`firstpair'" != "") + ///
        4 * ("`clusters'" != "") + 8 * ("`continuous'" != "")
    sreturn local methods `pairwise' `firstpair' `clusters' `continuous'
    sreturn local type = cond("`all'" != "", 0, cond("`firstpair'" != "", 2, cond("`pairwise'" != "", 3, 1)))
    sreturn local clusters = "`all'`clusters'" != ""
end

* e() results shared by creghdfe and cppmlhdfe.
* _ctools_hdfe_post, stub(name) absvars(string) nhdfe(#) vce(name)
*     [extended(string) clustvars(string) nclust(numlist) dofmethod(string)]
*
* Call after ereturn post, while the plugin's per-effect scalars
* <stub>_num_levels_g, <stub>_fe_redundant_g and <stub>_fe_nested_g (and
* <stub>_fe_exact_g where the plugin knows it) exist, g = 1..the number of
* words in extended(). Posts, as reghdfe does:
*   e(absvars), e(extended_absvars) (default: absvars without i. prefixes),
*   e(N_hdfe) = nhdfe, e(N_hdfe_extended)
*   e(dof_table): Categories, Redundant, Num_Coefs, Inexact (missing where
*   unknown) and Nested per effect; e(df_a_initial) and e(df_a_redundant).
*   Option reghdfenames names its rows and columns exactly as reghdfe does:
*   rows 1.var, 1.var1#1.var2, 1.var#c.slope (_cons without absorb()),
*   columns "Categories" "Redundant" "Num Coefs" "Inexact?" "Nested?"
*   e(dofmethod); e(vce), and e(vcetype) Robust unless vce is unadjusted
*   with clustvars: e(clustvar), e(N_clustervars), e(clustvar#), e(N_clust#)
*   from nclust (one count per variable), e(N_clust) their minimum
program define _ctools_hdfe_post, eclass
    version 14.1
    syntax, STUB(name) ABSvars(string) NHDFE(integer) VCE(name) ///
        [EXTended(string) CLUSTvars(string) NCLUST(numlist integer) DOFmethod(string) REGHDFEnames]
    if `"`extended'"' == "" local extended : subinstr local absvars "i." "", all
    local nfe : word count `extended'
    tempname table
    matrix `table' = J(`nfe', 5, .)
    local initial = 0
    local redundant = 0
    forvalues g = 1/`nfe' {
        matrix `table'[`g', 1] = scalar(`stub'_num_levels_`g')
        matrix `table'[`g', 2] = scalar(`stub'_fe_redundant_`g')
        matrix `table'[`g', 3] = `table'[`g', 1] - `table'[`g', 2]
        capture confirm scalar `stub'_fe_exact_`g'
        if !_rc matrix `table'[`g', 4] = 1 - scalar(`stub'_fe_exact_`g')
        matrix `table'[`g', 5] = scalar(`stub'_fe_nested_`g')
        local initial = `initial' + `table'[`g', 1]
        local redundant = `redundant' + `table'[`g', 2]
    }
    matrix colnames `table' = Categories Redundant Num_Coefs Inexact Nested
    matrix rownames `table' = `extended'

    ereturn local absvars `"`absvars'"'
    ereturn local extended_absvars `"`extended'"'
    ereturn scalar N_hdfe = `nhdfe'
    ereturn scalar N_hdfe_extended = `nfe'
    ereturn scalar df_a_initial = `initial'
    ereturn scalar df_a_redundant = `redundant'
    ereturn matrix dof_table = `table'
    if "`reghdfenames'" != "" {
        local rows
        foreach e of local extended {
            if "`e'" != "_cons" {
                local parts : subinstr local e "#" " ", all
                local e
                foreach p of local parts {
                    if substr("`p'", 1, 2) != "c." local p 1.`p'
                    local e = cond("`e'" == "", "`p'", "`e'#`p'")
                }
            }
            local rows `rows' `e'
        }
        mata: _ctools_hdfe_post_names("`rows'")
    }
    ereturn local dofmethod `"`dofmethod'"'
    ereturn local vce `vce'
    if "`vce'" != "unadjusted" ereturn local vcetype "Robust"

    local ndims : word count `clustvars'
    if `ndims' {
        if `: word count `nclust'' != `ndims' {
            di as error "_ctools_hdfe_post: one cluster count per cluster variable"
            exit 198
        }
        local fewest .
        forvalues d = 1/`ndims' {
            local count : word `d' of `nclust'
            ereturn local clustvar`d' `: word `d' of `clustvars''
            ereturn scalar N_clust`d' = `count'
            local fewest = min(`fewest', `count')
        }
        ereturn local clustvar `"`clustvars'"'
        ereturn scalar N_clustervars = `ndims'
        ereturn scalar N_clust = `fewest'
    }
end

capture mata: mata drop _ctools_hdfe_post_names()
mata:
// e(dof_table) row and column names as reghdfe posts them (with
// st_matrix*stripe, so that 1.var and "Num Coefs" stay as they are)
void _ctools_hdfe_post_names(string scalar rows)
{
    string colvector r, c
    r = tokens(rows)'
    c = ("Categories" \ "Redundant" \ "Num Coefs" \ "Inexact?" \ "Nested?")
    st_matrixrowstripe("e(dof_table)", (J(rows(r), 1, ""), r))
    st_matrixcolstripe("e(dof_table)", (J(5, 1, ""), c))
}
end

* Native sample preparation before factor expansion or weight staging
program define _ctools_hdfe_screen
    version 14.1
    gettoken action 0 : 0
    if "`action'" == "plan" {
        syntax anything(name=model)
        * list strips factor operators without choosing levels or rejecting
        * invalid categories on observations that will be excluded later.
        * Keep time-series operators on the original marksample path. Remove
        * parentheses first so ib(first), ib(last), ib(freq) remain eligible.
        local ops = ustrregexra(`" `model'"', "\([^)]*\)", "")
        local ok = !regexm(`"`ops'"', "[ #][^ .#]*[LlFfDdSs][^ .#]*[.]")
        local inputs
        if `ok' {
            tempname rhold
            _return hold `rhold'
            capture fvrevar `model', list
            local ok = (_rc == 0)
            if `ok' {
                local inputs `r(varlist)'
                capture confirm numeric variable `inputs'
                local ok = (_rc == 0)
            }
            _return restore `rhold'
        }
        c_local __screen_ok = `ok'
        c_local __screen_inputs `inputs'
        exit
    }
    if "`action'" != "mark" error 198
    syntax varlist(numeric), TOUSE(varname numeric) PLUGIN(string) [THREADS(integer 0)]
    local inputs : list uniq varlist
    capture program ctools_plugin, plugin using("`plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    plugin call ctools_plugin `inputs' `touse', "hdfe_markout `threadopt'"
end

* vce() of the HDFE estimators.
* _ctools_hdfe_vce, spec(string) touse(varname) cmd(name)
*
* Classifies spec with _ctools_vce: r(kind) is cluster, robust, unadjusted,
* ols, iid, bootstrap, or empty for anything else. Which kinds and how many
* cluster dimensions a command accepts, and its messages, stay with the
* command. For cluster (vce(cluster a b) or vce(cluster(a b))), each term
* is a dimension: a variable (abbreviations and wildcards expand to one
* dimension per variable), an i.-prefixed variable, or an interaction
* var1#var2 of such variables, as in reghdfe; string variables are allowed.
* Marks out touse for missing cluster values (empty strings included) and
* returns r(ndims) dimensions with, for d = 1..r(ndims):
*   r(name_d)   the dimension as reghdfe names it: turn, turn#foreign
*   r(vars_d)   its variables
*   r(group_d)  1 when they must be numbered with _ctools_group first (a
*               string variable or an interaction), 0 when the numeric
*               variable goes to the plugin as it is
* Option nomark leaves the dimensions with r(group_d) 0 in touse: the
* caller's plugin drops the rows where such a variable is missing (its
* native missing-value screen), so a markout would be a second pass.
program define _ctools_hdfe_vce, rclass
    version 14.1
    syntax, TOUSE(varname) CMD(name) [SPEC(string) NOMARK]
    _ctools_vce kind unused : `spec'
    return local kind `kind'
    local n = 0
    if "`kind'" == "cluster" {
        gettoken first terms : spec, parse(" (")
        local terms : subinstr local terms "(" " ", all
        local terms : subinstr local terms ")" " ", all
        foreach term of local terms {
            if !strpos("`term'", "#") & !strpos("`term'", ".") {
                unab expanded : `term', name(vce(cluster))
                foreach v of local expanded {
                    local ++n
                    capture confirm numeric variable `v'
                    local group = (_rc != 0)
                    return scalar group_`n' = `group'
                    return local name_`n' `v'
                    return local vars_`n' `v'
                    if `group' | "`nomark'" == "" markout `touse' `v', strok
                }
                continue
            }
            local parts : subinstr local term "#" " ", all
            local group = (`: word count `parts'' > 1)
            local name
            local vars
            foreach part of local parts {
                if substr("`part'", 1, 2) == "i." local part = substr("`part'", 3, .)
                if "`part'" == "" | strpos("`part'", ".") {
                    di as error `"`cmd': vce(cluster) term `term' is not supported"'
                    di as error "specify variables, i.var, or interactions such as var1#var2"
                    exit 198
                }
                unab v : `part', max(1) name(vce(cluster))
                capture confirm numeric variable `v'
                if _rc local group = 1
                local vars `vars' `v'
                if "`name'" != "" local name "`name'#"
                local name "`name'`v'"
            }
            local ++n
            return scalar group_`n' = `group'
            return local name_`n' `name'
            return local vars_`n' `vars'
            if `group' | "`nomark'" == "" markout `touse' `vars', strok
        }
    }
    return scalar ndims = `n'
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

* ctools output-name validation
program define _ctools_newvars
    version 14.1
    syntax anything(name=outputs)
    local seen ""
    foreach name of local outputs {
        if `: list name in seen' {
            di as error "ctools: duplicate output variable `name'"
            exit 198
        }
        confirm new variable `name'
        local seen `seen' `name'
    }
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

* ctools: classify a vce() specification
* _ctools_vce <kind> <rest> : <specification> sets two locals of the caller:
* <kind> is cluster (any word starting with cl), robust (robust, r),
* unadjusted, ols, iid, bootstrap (any word starting with boot), or empty for
* anything else, case-insensitively; <rest> is the text after the first word.
* Which kinds a command accepts, and its messages, stay with the command.
program define _ctools_vce
    version 14.1
    gettoken kind_macro 0 : 0
    gettoken rest_macro 0 : 0
    gettoken colon 0 : 0, parse(":")
    local lower = lower(trim(`"`0'"'))
    local kind
    if substr(`"`lower'"', 1, 2) == "cl" local kind cluster
    else if inlist(`"`lower'"', "robust", "r") local kind robust
    else if inlist(`"`lower'"', "unadjusted", "ols", "iid") local kind `lower'
    else if substr(`"`lower'"', 1, 4) == "boot" local kind bootstrap
    gettoken first rest : 0
    local rest = trim(`"`rest'"')
    c_local `kind_macro' `kind'
    c_local `rest_macro' `rest'
end

* ctools weight-expression staging; caller owns generate() tempvar
* Sets the caller's local __ctools_wvar to the variable that holds the
* weights: generate(), filled here and checked. With option marked the
* caller's marksample has seen this weight: the sample then leaves out zero
* and missing weights, and Stata has refused negative weights (r(402)) and
* noninteger frequency weights (r(401)). A weight that is one numeric
* variable is then used as it is, and generate() is not created.
program define _ctools_weight
    version 14.1
    syntax, Generate(name) TOUse(varname) TYPE(string) EXPression(string) [MARKED]
    c_local __ctools_wvar `generate'
    if "`marked'" != "" {
        local wexp = strtrim(substr(strtrim(`"`expression'"'), 2, .))
        if `: word count `wexp'' == 1 {
            capture confirm numeric variable `wexp', exact
            if !_rc {
                c_local __ctools_wvar `wexp'
                exit
            }
        }
    }
    confirm new variable `generate'
    quietly generate double `generate' `expression' if `touse'
    markout `touse' `generate'
    quietly count if `generate' <= 0 & `touse'
    if r(N) {
        di as error "ctools: weights must be positive"
        exit 198
    }
    if "`type'" == "fweight" {
        quietly count if `generate' != floor(`generate') & `touse'
        if r(N) {
            di as error "frequency weights must be integers"
            exit 401
        }
    }
end

* The wild bootstrap table after the coefficient table of an HDFE estimator.
* Design: docs/development/BOOTSTRAP_PLAN.md (section 4.7).
*
* _ctools_wild_table
*     Shows, from the current e(): the tested coefficients with e(b), their
*     analytic t (e(wb_t)), the bootstrap p-values (e(wb_p)) and, when
*     posted, the confidence sets by test inversion (e(wb_ci),
*     e(wb_cistatus)); for the score bootstrap the score statistic
*     (e(wb_score)) in place of e(wb_t); then notes on the p-value type, enumeration, few
*     Rademacher clusters, open, disjoint or empty sets, trial values whose p could
*     not be evaluated, restricted fits that failed
*     (e(wb_fail)) and replicates without a statistic (e(wb_misreps)).
*     e(wb_method): wcr, wcu (creghdfe), wre, wue (civreghdfe) or score
*     (cppmlhdfe).

program define _ctools_wild_table
    version 14.1
    tempname t p ci st fail
    matrix `t' = e(wb_t)
    * the score bootstrap studentizes the score statistic, not the Wald z
    local score_t = 0
    if "`e(wb_method)'" == "score" {
        capture confirm matrix e(wb_score)
        if !_rc {
            matrix `t' = e(wb_score)
            local score_t = 1
        }
    }
    matrix `p' = e(wb_p)
    local hasci = 0
    capture confirm matrix e(wb_ci)
    if !_rc {
        local hasci = 1
        matrix `ci' = e(wb_ci)
        matrix `st' = e(wb_cistatus)
    }
    local hasfail = 0
    capture confirm matrix e(wb_fail)
    if !_rc {
        local hasfail = 1
        matrix `fail' = e(wb_fail)
    }
    if "`e(wb_method)'" == "score" {
        local kind = cond("`e(wb_bootclust)'" != "", "Score cluster bootstrap", "Score bootstrap")
        local how null imposed
    }
    else {
        local kind = cond("`e(wb_bootclust)'" != "", "Wild cluster bootstrap", "Wild bootstrap")
        local how = cond("`e(wb_method)'" == "wcr", "restricted (WCR)", ///
            cond("`e(wb_method)'" == "wre", "restricted efficient (WRE)", ///
            cond("`e(wb_method)'" == "wue", "unrestricted efficient (WUE)", "unrestricted (WCU)")))
    }
    local wt `e(wb_weighttype)'
    if inlist("`wt'", "rademacher", "mammen", "webb") local wt = strproper("`wt'")
    local who = cond("`e(wb_bootclust)'" != "", "(`e(wb_bootclust)')", "(observations)")
    local stat = cond("`e(wb_method)'" == "score", "z", "t")
    local pl = cond("`e(wb_ptype)'" == "lower", "P<`stat'", cond("`e(wb_ptype)'" == "upper", "P>`stat'", "P>|`stat'|"))
    di as text ""
    di as text "`kind': `how', `wt' weights"
    di as text "Replications = " as result strtrim(string(e(wb_reps), "%12.0fc")) as text ///
        _col(40) "Bootstrap clusters = " as result strtrim(string(e(wb_N_clust), "%12.0fc")) ///
        as text " `who'"
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`e(depvar)'", 12) " {c |}" _col(20) "Coef." _col(34) "`stat'" _col(39) "`pl'" ///
        _col(50) "[`=strofreal(e(wb_level))'% Conf. Interval]"
    di as text "{hline 13}{c +}{hline 64}"
    local names : colnames `p'
    local j 0
    local open
    local disjoint
    local empty
    local trialfail
    local failed
    foreach c of local names {
        local ++j
        local lo .
        local hi .
        if `hasci' {
            local lo = `ci'[`j', 1]
            local hi = `ci'[`j', 2]
            if `st'[`j', 1] | `st'[`j', 2] local open `open' `c'
            if `st'[`j', 3] local disjoint `disjoint' `c'
            if `st'[`j', 4] local empty `empty' `c'
            if colsof(`st') >= 5 {
                if `st'[`j', 5] local trialfail `trialfail' `c'
            }
        }
        if `hasfail' {
            if `fail'[1, `j'] local failed `failed' `c'
        }
        di as text %12s abbrev("`c'", 12) " {c |}" as result _col(15) %10.0g _b[`c'] ///
            _col(27) %8.2f `t'[1, `j'] _col(38) %6.4f `p'[1, `j'] _col(49) %10.0g `lo' ///
            _col(62) %10.0g `hi'
    }
    di as text "{hline 13}{c BT}{hline 64}"
    local pdesc = cond("`e(wb_ptype)'" == "symmetric", "symmetric", cond("`e(wb_ptype)'" == "equaltail", "equal-tailed", "`e(wb_ptype)'-tailed"))
    di as text "P-values are `pdesc' bootstrap-`stat'" _c
    if `hasci' di as text "; intervals invert the bootstrap test." _c
    di as text ""
    if `score_t' di as text "z is the score statistic (e(wb_score)); e(wb_t) holds the Wald z."
    if e(wb_enumerated) == 1 {
        di as text "All " as result e(wb_reps) as text " Rademacher sign vectors were enumerated (no random draws)."
    }
    if "`e(wb_weighttype)'" == "rademacher" & e(wb_N_clust) <= 11 {
        di as text "{p 0 6 2}Note: with " e(wb_N_clust) " bootstrap clusters, Rademacher weights allow only " ///
            2^e(wb_N_clust) " distinct draws; consider weighttype(webb).{p_end}"
    }
    if "`e(wb_fe_path)'" == "partialled" {
        di as text "{p 0 6 2}Note: absorbed effects not nested in the bootstrap clusters enter the bootstrap as partialled out (fedgp(partialled)).{p_end}"
    }
    if "`open'" != "" di as text "{p 0 6 2}Note: the confidence set of`open' is unbounded (missing bound) within 1000 standard errors.{p_end}"
    if "`disjoint'" != "" di as text "{p 0 6 2}Note: the confidence set of`disjoint' is not an interval; the bounds enclose it.{p_end}"
    if "`empty'" != "" di as text "{p 0 6 2}Note: no trial value of`empty' was accepted, or its interval could not be searched: no bounds.{p_end}"
    if "`trialfail'" != "" di as text "{p 0 6 2}Note: the p-value could not be evaluated at some trial values for`trialfail'; they were treated as not rejected, so the interval may be wider than the confidence set.{p_end}"
    if "`failed'" != "" di as text "{p 0 6 2}Note: the fit under the null did not converge for`failed'; no bootstrap result.{p_end}"
    if e(wb_misreps) > 0 & e(wb_misreps) < . {
        di as text "Note: " e(wb_misreps) " replicates had no defined statistic and were left out."
    }
end
* ctools-helpers:end
