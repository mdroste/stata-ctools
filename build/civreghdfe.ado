*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define civreghdfe, eclass byable(recall)
    version 14.1

    * Replay: civreghdfe without a variable list redisplays the results,
    * with ivreg2's replay options; by: is refused, as in ivreg2
    if replay() {
        if _by() {
            di as error "invalid syntax - cannot use by with replay"
            exit 601
        }
        if `"`e(cmd)'"' != "civreghdfe" error 301
        syntax [, FIRST FFIRst RF SFIRST dropfirst droprf dropsfirst Level(cilevel) ///
            noHEader noFOoter EForm(string) PLUS NOOMITted vsquish NOEMPtycells ///
            BASElevels ALLBASElevels NOWILDtable]

        * a stored first-stage, reduced-form or system equation (restored
        * with estimates restore) is shown as ivreg2 shows it
        if inlist("`e(model)'", "first", "rf", "sfirst") {
            local dispopt `noomitted' `vsquish' `noemptycells' `baselevels' `allbaselevels' `plus'
            if `"`eform'"' != "" local dispopt `dispopt' eform(`eform')
            _civreghdfe_dispsf, level(`level') dispopt(`dispopt')
            exit
        }
        * the collinearity and duplicates warnings, as ivreg2 repeats them
        if "`e(dups)'" != "" {
            di as text "Warning - duplicate variables detected"
            di as text "Duplicates:" _col(16) "`e(dups)'"
        }
        if "`e(collin)'" != "" {
            di as text "Warning - collinearities detected"
            di as text "Vars dropped:" _col(16) "`e(collin)'"
        }
        local dispopt `noomitted' `vsquish' `noemptycells' `baselevels' `allbaselevels' `plus'
        if `"`eform'"' != "" local dispopt `dispopt' eform(`eform')
        local __stored
        if "`dropsfirst'" == "" local __stored `e(sfirsteq)'
        if "`droprf'" == "" local __stored `__stored' `e(rfeq)'
        if "`dropfirst'" == "" local __stored `__stored' `e(firsteqs)'
        _civreghdfe_disp, `first' `ffirst' `rf' `sfirst' level(`level') `header' `footer' ///
            dispopt(`dispopt') stored(`__stored') `nowildtable'
        * drop*: the stored first-stage, reduced-form and system equations
        * are dropped and their names cleared, as ivreg2 does
        foreach __d in "dropfirst firsteqs" "droprf rfeq" "dropsfirst sfirsteq" {
            gettoken __o __m : __d
            local __m = strtrim("`__m'")
            if "``__o''" != "" & "`e(`__m')'" != "" {
                foreach __eq in `e(`__m')' {
                    capture estimates drop `__eq'
                }
                ereturn local `__m'
            }
        }
        exit
    }

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start wall clock timer immediately
    quietly _ctools_clock __ct_lap97

    * e(cmdline) is the command as typed: syntax below reuses `0' for the
    * variable specification alone
    local cmdline `"civreghdfe `0'"'
    * The options as typed, for the second pass after ivreg2's collinearity
    * rules (below)
    _parse comma __civ_lhs __civ_opts : 0
    * wildbootstrap or wildbootstrap(wildopts) is taken out of the options
    * before syntax (it is given with or without options)
    _civreghdfe_wild_opt `__civ_opts'
    if `__wild' local 0 `"`__civ_lhs' , `__wild_rest'"'
    * reghdfe's verbose(#), which ivreghdfe passes on: a positive level is
    * civreghdfe's verbose
    _civreghdfe_verbose_opt , `__wild_rest'
    if `__verb_given' local 0 `"`__civ_lhs' , `__verb_rest'"'

    * Parse syntax
    * Basic syntax: civreghdfe depvar (endogvars = instruments) [exogvars], absorb() [options]
    * absorb() is optional - without it, runs as regular IV (no FE absorption)
    syntax anything(name=0 equalok) [if] [in] [aw fw pw iw/], ///
        [Absorb(string) ///
        vce(string) ///
        FEVCE(string) ///
        Robust ///
        CLuster(varlist) ///
        TOLerance(real 1e-8) ///
        MAXiter(integer 16000) ///
        ITERATE(integer 0) ///
        Verbose ///
        FIRST ///
        FFIRst ///
        RF ///
        SMall ///
        LIML ///
        FULLER(real 0) ///
        Kclass(real 0) ///
        GMM2s ///
        CUE ///
        COVIV ///
        BW(string) ///
        KERNEL(string) ///
        DKRAAY(integer 0) ///
        KIEFER ///
        CENTer ///
        B0(string) ///
        RESiduals(name) ///
        RESiduals2 ///
        NOConstant ///
        Level(string) ///
        noHEader ///
        noFOoter ///
        noOUTput ///
        TItle(string) ///
        SUBtitle(string) ///
        DEPname(string) ///
        noid ///
        ORThog(varlist fv ts) ///
        ENDOGtest(varlist fv ts) ///
        REDundant(varlist fv ts) ///
        PARTial(string) ///
        FWL(string) ///
        DOFminus(integer 0) ///
        SDOFminus(integer 0) ///
        NOPARTIALsmall ///
        EForm(string) ///
        PLUS ///
        NOOMITted ///
        OMITted ///
        VSQuish ///
        NOEMPtycells ///
        BASElevels ///
        ALLBASElevels ///
        SAVEFirst ///
        SAVEFPrefix(string) ///
        SAVERF ///
        SAVERFPrefix(string) ///
        SAVESFirst ///
        SAVESFPrefix(name) ///
        SMATRIX(string) ///
        WMATRIX(string) ///
        THReads(integer 0) ///
        COEFLegend ///
        CFORMAT(string) ///
        NOCOLLIN ///
        DOFadjustments(string) KEEPSINgletons SFIRST ///
        CIVRERUN CIVCOLLIN(string) CIVECOLLIN(string) CIVDUPS(string) ///
        NORETURN NOWILDtable CIVAUX ///
        PSD0 PSDA SW Tvar(varname) Ivar(varname) BVCLEAN FVALL FVSEP USEQR ///
        TEChnique(string) TRAnsform(string) ACCELeration(string) PREConditioner(string) ///
        PRUNE NOSAMPle POOLsize(integer 10) GROUPVar(name) TIMEit noWARN]

    * Freeze incoming-r() qualifiers before parsing helpers replace r().
    if `"`if'"' != "" {
        tempvar __incoming_sample
        mark `__incoming_sample' `if' `in'
        local if "if `__incoming_sample'"
    }

    * level(): checked as cilevel checks it, and remembered as given (with
    * vce(bootstrap) it also sets the bootstrap's level)
    local __cmdlevel `"`level'"'
    _ctools_hdfe_boot level, cmd(civreghdfe) bslevel(`c(level)') bsgiven(0) cmdlevel(`"`level'"')
    local level = r(level)

    * Drop whatever an earlier run left behind, so no stale state is read
    _ctools_cleanup civreghdfe

    * ivreg2's psd0 and psda: the moment covariance S is made positive
    * semidefinite (negative eigenvalues set to zero, or to their absolute
    * values)
    if "`psd0'" != "" & "`psda'" != "" {
        di as error "cannot use psd0 and psda options together"
        exit 198
    }
    local psd_code = cond("`psd0'" != "", 1, cond("`psda'" != "", 2, 0))
    * sw: ivreg2's Stock-Watson (2008) standard errors, not implemented
    if "`sw'" != "" {
        di as error "civreghdfe: sw (Stock-Watson heteroskedasticity-robust standard errors) is not supported"
        exit 198
    }
    * tvar() and ivar() must name the tsset time and panel variables, which
    * is all ivreg2 accepts
    if "`tvar'" != "" & "`tvar'" != "`_dta[_TStvar]'" {
        di as error "invalid tvar() option - data already -tsset-"
        exit 5
    }
    if "`ivar'" != "" & "`ivar'" != "`_dta[_TSpanel]'" {
        di as error "invalid ivar() option - data already -xtset-"
        exit 5
    }
    * reghdfe's solver options, which ivreghdfe passes on, are checked as
    * reghdfe checks them; civreghdfe's solver does not use them
    _civreghdfe_solver_opts, technique(`technique') transform(`transform') ///
        acceleration(`acceleration') preconditioner(`preconditioner')
    * timeit: the timing breakdown (verbose's)
    local __timeit = ("`timeit'" != "")
    * iweights, as ivreg2 takes them: without absorb() (reghdfe takes none)
    * and with the homoskedastic VCE only
    if "`weight'" == "iweight" {
        if `"`absorb'"' != "" {
            di as error "iweights not allowed"
            exit 101
        }
        if "`robust'`cluster'`gmm2s'`kernel'`bw'`kiefer'" != "" | `dkraay' != 0 | ///
            (`"`vce'"' != "" & !inlist(lower(strtrim(`"`vce'"')), "unadjusted", "ols")) {
            di as error "iweights not allowed with robust or gmm"
            exit 101
        }
    }

    * Plugin statistics are copied into temporary scalars before they are read
    tempname temp_F temp_r2 temp_underid temp_underid_df temp_cd_f temp_kp_f temp_kp_f2
    tempname temp_sargan temp_sargan_df temp_endog_chi2 temp_endog_f temp_endog_df
    tempname temp_cstat temp_cstat_df temp_endogtest_stat temp_endogtest_df
    tempname temp_redund_stat temp_redund_df

    * Handle standalone robust/cluster as aliases for vce()
    if "`cluster'" != "" & `"`vce'"' == "" {
        local vce "cluster `cluster'"
    }
    if "`robust'" != "" & `"`vce'"' == "" {
        local vce "robust"
    }
    * cluster() implies robust, as in ivreg2, so the two may be combined
    if "`robust'" != "" & `"`vce'"' != "" & `"`vce'"' != "robust" & "`cluster'" == "" {
        di as error "cannot specify both robust and vce()"
        exit 198
    }
    if "`cluster'" != "" & `"`vce'"' != "" & !regexm(`"`vce'"', "^cluster") {
        di as error "cannot specify both cluster() and vce()"
        exit 198
    }

    * iterate() is reghdfe's (and so ivreghdfe's) name for maxiter()
    if `iterate' > 0 local maxiter = `iterate'

    * b0(): ivreg2 evaluates the CUE objective at the given coefficients;
    * it takes no other estimator, and implies noid
    if `"`b0'"' != "" {
        if "`gmm2s'`cue'`liml'" != "" | `fuller' != 0 {
            di as error "incompatible options: -b0- and `gmm2s' `cue' `liml'" ///
                cond(`fuller' != 0, " fuller(`fuller')", "")
            exit 198
        }
        if `kclass' != 0 {
            di as error "incompatible options: -b0- and kclass(`kclass')"
            exit 198
        }
        capture confirm matrix `b0'
        if _rc {
            di as error "-b0- option error: supplied b0 columns do not match regressor list"
            exit 198
        }
        local id noid
    }

    * Estimators that cannot be combined, as ivreg2 parses them: liml is not
    * a GMM estimator, and gmm2s and cue are different GMM estimators; the
    * k-class parameter cannot be negative, nor the Fuller parameter with
    * liml. As in ivreg2, fuller() and kclass() take precedence over gmm2s,
    * kclass() over liml (k is the kclass() value) and a negative fuller()
    * is LIML. fuller() with kclass() (ivreg2 fails with r(3200)), and fuller()
    * or kclass() with cue (ivreg2 then runs CUE and labels it LIML or
    * k-class), are refused.
    if "`liml'" != "" & "`gmm2s'`cue'" != "" {
        di as error "GMM estimation not available with LIML or k-class estimators"
        exit 198
    }
    if (`fuller' != 0 | `kclass' != 0) & "`cue'" != "" {
        di as error "GMM estimation not available with LIML or k-class estimators"
        exit 198
    }
    if `fuller' < 0 & "`liml'" != "" {
        di as error "invalid Fuller option"
        exit 198
    }
    if `kclass' < 0 {
        di as error "invalid k-class option"
        exit 198
    }
    if `fuller' != 0 & `kclass' != 0 {
        di as error "cannot use fuller and kclass options together"
        exit 198
    }
    if "`gmm2s'" != "" & "`cue'" != "" {
        di as error "incompatible options: 2-step efficient gmm and cue gmm"
        exit 198
    }

    * partial() and its legacy name fwl(), as ivreg2 parses them: _cons
    * partials out the constant (also alone), _all every included exogenous
    * regressor; the rest is a varlist
    local partial `partial' `fwl'
    local partial : subinstr local partial "_cons" "", all count(local __pcons_given) word
    local partial : list retokenize partial
    if `__pcons_given' > 1 {
        di as error "Error: _cons listed more than once in partial()."
        exit 198
    }
    if `__pcons_given' & "`noconstant'" != "" {
        di as error "Error: _cons listed in partial() but equation specifies -noconstant-."
        exit 198
    }
    if `__pcons_given' & `"`absorb'"' != "" {
        di as error "partial(): _cons cannot be partialled out with absorb(), which absorbs the constant"
        exit 198
    }
    local __partial_all = ("`partial'" == "_all")
    if "`partial'" != "" & !`__partial_all' {
        _civreghdfe_vlcheck `partial'
        local partial `r(varlist)'
    }

    * dofminus() and sdofminus() as given, for the auxiliary regressions
    * (sdofminus is later augmented with the absorbed effects)
    local __user_dofminus = `dofminus'
    local __user_sdofminus = `sdofminus'

    * As in ivreg2, a prefix for the stored equations implies storing them
    if "`savefprefix'" != "" local savefirst savefirst
    if "`saverfprefix'" != "" local saverf saverf
    if "`savesfprefix'" != "" local savesfirst savesfirst
    * nooutput suppresses all output, as in ivreg2; results are still stored
    if "`output'" != "" {
        local header noheader
        local footer nofooter
    }

    * wmatrix() is ivreg2's weighting matrix for IV and GMM; b0() takes none
    if "`b0'" != "" & `"`wmatrix'"' != "" {
        di as error "incompatible options: -b0- and wmatrix"
        exit 198
    }
    * dofadjustments() and keepsingletons, which ivreghdfe passes to reghdfe
    _ctools_hdfe_dof, `=lower(`"`dofadjustments'"')'
    local dof_type = `s(type)'
    local dof_clusters = `s(clusters)'
    local dofmethods `s(methods)'
    local dofmethods : list retokenize dofmethods

    * residuals2: the residuals as _reghdfe_resid, replaced if it exists,
    * as ivreghdfe saves them
    if "`residuals2'" != "" {
        local residuals "_reghdfe_resid"
    }

    if "`residuals'" != "" & "`residuals2'" == "" _ctools_newvars `residuals'
    if "`groupvar'" != "" & `"`absorb'"' != "" confirm new variable `groupvar'

    * Parse the varlist with parentheses for endogenous vars and instruments
    * Format: depvar [exog_before] (endog1 endog2 = inst1 inst2 inst3) [exog_after]
    local 0 `"`0'"'

    * FV groups and lag ranges also contain parentheses. Only a balanced
    * outer group with an equals sign at its own depth is an IV equation.
    gettoken depvar model : 0
    local depvar = strtrim("`depvar'")
    local depth = 0
    local group_eq = 0
    local eq_count = 0
    local iv_start = 0
    local iv_end = 0
    local iv_eq = 0
    local model_length = strlen("`model'")
    if `model_length' > 0 {
        forvalues pos = 1/`model_length' {
            local char = substr("`model'", `pos', 1)
            if "`char'" == "(" {
                if `depth' == 0 {
                    local group_start = `pos'
                    local group_eq = 0
                }
                local ++depth
            }
            else if "`char'" == "=" {
                if `depth' != 1 | `group_eq' != 0 error 198
                local group_eq = `pos'
            }
            else if "`char'" == ")" {
                local --depth
                if `depth' < 0 error 198
                if `depth' == 0 & `group_eq' == 0 {
                    * Without "=", the group must belong to an operator, as in
                    * c.(x1 x2), i.a#(b c) or L(1/2).x; otherwise an IV
                    * equation is missing its equals sign.
                    local before = cond(`group_start' > 1, substr("`model'", `group_start'-1, 1), "")
                    local after = substr("`model'", `pos'+1, 1)
                    if !inlist("`before'", ".", "#") & "`after'" != "." {
                        di as error "Must specify instruments after '=' in parentheses"
                        exit 198
                    }
                }
                if `depth' == 0 & `group_eq' != 0 {
                    local ++eq_count
                    if `eq_count' > 1 error 198
                    local iv_start = `group_start'
                    local iv_end = `pos'
                    local iv_eq = `group_eq'
                }
            }
        }
    }
    if `depth' != 0 {
        di as error "Unmatched parenthesis in variable specification"
        exit 198
    }
    local is_iv_model = (`eq_count' == 1)
    local endogvars
    local instruments
    local exogvars `model'
    if `is_iv_model' {
        local endogvars = strtrim(substr("`model'", `iv_start'+1, `iv_eq'-`iv_start'-1))
        local instruments = strtrim(substr("`model'", `iv_eq'+1, `iv_end'-`iv_eq'-1))
        * (=z) adds excluded instruments without endogenous regressors, as
        * ivreghdfe allows; an equation without instruments is an error.
        if "`instruments'" == "" error 198
        local exog_before = substr("`model'", 1, `iv_start'-1)
        local exog_after = substr("`model'", `iv_end'+1, .)
        local exogvars `exog_before' `exog_after'
    }
    local exogvars = strtrim("`exogvars'")
    * partial(_all): every included exogenous regressor
    if `__partial_all' local partial `exogvars'

    * Validate dependent variable exists and is numeric (error 109 for string
    * vars). A time-series operator (L.y, D.y) is applied by tsrevar: the
    * plugin gets the generated variable, depvar keeps the name as typed.
    capture confirm string variable `depvar'
    if _rc == 0 {
        di as error "type mismatch"
        exit 109
    }
    tsrevar `depvar'
    local depvar_var `r(varlist)'
    confirm numeric variable `depvar_var'
    * Mark sample early for fvrevar
    marksample touse
    markout `touse' `depvar_var'

    * One numeric variable per absorbed effect for the plugin: string
    * variables and interactions are numbered first (the plugin renumbers
    * every effect on the estimation sample). absorb_names holds the effects
    * as reghdfe names them, for e(absorb), the output and the auxiliary
    * first-stage and reduced-form regressions.
    local absorb_names
    * Saved fixed effects (as reghdfe): absorb(..., savefe) names them
    * __hdfe#__, absorb(newvar=absvar) names one; fe_target_# holds the name
    * of effect #, fe_terms the effects without the names
    local fe_any = 0
    local fe_terms
    local any_slopes = 0
    local any_intercept = 0
    local ext_names
    local n_absorb = 0
    if `"`absorb'"' != "" {
        local savefe = 0
        if strpos(`"`absorb'"', ",") {
            gettoken absorb absorb_opts : absorb, parse(",")
            gettoken __comma absorb_opts : absorb_opts, parse(",")
            _civreghdfe_absorb_opts, `absorb_opts'
        }
        * Terms are bound by parentheses (id##c.(x z)); a target may be
        * written with spaces around the equal sign, as reghdfe allows
        local absorb = ustrregexra(`"`absorb'"', "\s*=\s*", "=")
        local __t 0
        local __terms
        while `"`absorb'"' != "" {
            gettoken __term absorb : absorb, bind
            local ++__t
            local fe_target_`__t'
            local __eq = strpos("`__term'", "=")
            if `__eq' {
                local fe_target_`__t' = substr("`__term'", 1, `__eq' - 1)
                local __term = substr("`__term'", `__eq' + 1, .)
                confirm name `fe_target_`__t''
            }
            else if `savefe' local fe_target_`__t' __hdfe`__t'__
            if "`fe_target_`__t''" != "" local fe_any = 1
            local __terms `__terms' `__term'
        }
        local absorb `__terms'
        local fe_terms `__terms'
        * Heterogeneous slopes (id#c.x, id##c.x, id##c.(x z)) are absorbed as
        * reghdfe absorbs them, which ivreghdfe calls: one factor per
        * intercept and per slope (e(extended_absvars))
        _ctools_hdfe_absorb `absorb', touse(`touse') cmd(civreghdfe) slopes
        local n_absorb = r(n)
        if `fe_any' & `n_absorb' != `__t' {
            di as error "civreghdfe: each absorb() term saved with savefe or newvar= must be one effect"
            exit 198
        }
        forvalues g = 1/`n_absorb' {
            local absorb_name_`g' `"`r(name_`g')'"'
            local absorb_vars_`g' `r(vars_`g')'
            local absorb_group_`g' = r(group_`g')
            local absorb_int_`g' = r(intercept_`g')
            local absorb_cvars_`g' `r(cvars_`g')'
        }
        if `fe_any' {
            * A saved effect is its target and, for each slope, the target
            * with Slope# appended (__hdfe#__Slope# with savefe), as reghdfe
            local __fe_names
            forvalues __g = 1/`__t' {
                if "`fe_target_`__g''" == "" continue
                if `absorb_int_`__g'' local __fe_names `__fe_names' `fe_target_`__g''
                local __h = 0
                foreach __cv of local absorb_cvars_`__g' {
                    local ++__h
                    local __fe_names `__fe_names' `fe_target_`__g''Slope`__h'
                }
            }
            _ctools_newvars `__fe_names'
            * creghdfe recovers the effects into __hdfe#__ first
            forvalues __g = 1/`__t' {
                local __own __hdfe`__g'__
                local __h = 0
                foreach __cv of local absorb_cvars_`__g' {
                    local ++__h
                    local __own `__own' __hdfe`__g'__Slope`__h'
                }
                foreach __v of local __own {
                    capture confirm variable `__v', exact
                    if !_rc & "`fe_target_`__g''" != "__hdfe`__g'__" {
                        di as error "civreghdfe: variable `__v' already defined; saving fixed effects needs the name"
                        exit 110
                    }
                }
            }
        }
        * One plugin factor per intercept and per slope, in order: absorb
        * holds the factors' level variables, ext_slope_# the slope
        * variable of factor # (empty for an intercept), ext_label_# its
        * name in e(extended_absvars)
        local absorb
        local ext_names
        local G = 0
        local any_slopes = 0
        local any_intercept = 0
        forvalues g = 1/`n_absorb' {
            local absorb_names `"`absorb_names' `absorb_name_`g''"'
            if `absorb_group_`g'' {
                tempvar absorb_id_`g'
                quietly _ctools_group `absorb_vars_`g'' if `touse', generate(`absorb_id_`g'')
                local fe_var `absorb_id_`g''
            }
            else local fe_var `absorb_vars_`g''
            local ivlabel : subinstr local absorb_vars_`g' " " "#", all
            if `absorb_int_`g'' {
                local ++G
                local any_intercept = 1
                local absorb `absorb' `fe_var'
                local ext_slope_`G'
                local ext_absvar_`G' `g'
                local ext_names `ext_names' `ivlabel'
            }
            foreach cv of local absorb_cvars_`g' {
                local ++G
                local any_slopes = 1
                local absorb `absorb' `fe_var'
                local ext_slope_`G' `cv'
                local ext_absvar_`G' `g'
                local ext_names `ext_names' `ivlabel'#c.`cv'
            }
        }
        local absorb_names : list retokenize absorb_names
        if `any_slopes' & `G' > 10 {
            di as error "civreghdfe: at most ten absorbed intercepts and slopes are supported with heterogeneous slopes"
            exit 198
        }
    }

    * vce(bootstrap) and fevce(bootstrap): pairs bootstraps that the plugin
    * runs on the retained fit (docs/development/BOOTSTRAP_PLAN.md, 9.3).
    * With vce(bootstrap) the fit itself (its diagnostics, the GMM moment
    * covariance, the first-stage and reduced-form regressions) uses a base
    * VCE: cluster on the bootstrap's cluster(), otherwise robust. With
    * fevce() alone it is the user's vce().
    local __bs = 0
    local __fe = 0
    local bs_cluster
    local bs_strata
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
            local vce = cond(`"`bs_cluster'"' != "", `"cluster `bs_cluster'"', "robust")
        }
    }
    * The command's level() sets vce(bootstrap)'s, as for bootstrap; the
    * tables at estimation show the bootstrap's level
    if `__bs' {
        _ctools_hdfe_boot level, cmd(civreghdfe) bslevel(`bs_level') bsgiven(`bs_lgiven') ///
            cmdlevel(`"`__cmdlevel'"')
        local bs_level = r(level)
        local level = r(level)
    }
    local fe_nreps = 0
    local fe_ci = 0
    local fe_dse
    local fe_sestub
    if `"`fevce'"' != "" {
        _ctools_hdfe_boot parse, spec(`"`fevce'"') fe
        if "`r(kind)'" != "bootstrap" {
            di as error "civreghdfe: fevce() takes bootstrap[, options]: the fixed effects use the pairs bootstrap"
            exit 198
        }
        if !`fe_any' {
            di as error "civreghdfe: fevce(bootstrap) needs saved fixed effects (absorb(..., savefe) or newvar=)"
            exit 198
        }
        if `__bs' & (`"`r(saving)'"' != "" | "`r(bsdraws)'" != "") {
            di as error "civreghdfe: with vce(bootstrap), specify saving() and _bsdraws() there"
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
                di as error "civreghdfe: with vce(bootstrap), fevce(bootstrap) uses the level(), mse and dots() of vce(); specify them once, in vce()"
                exit 198
            }
        }
        if `__bs' & r(designopts) {
            if r(reps) != `bs_reps' | `"`r(seed)'"' != `"`bs_seed'"' | `"`r(cluster)'"' != `"`bs_cluster'"' | ///
                `"`r(strata)'"' != `"`bs_strata'"' | "`r(draws)'" != "`bs_draws'" {
                di as error "civreghdfe: vce(bootstrap) and fevce(bootstrap) share one set of replicates; specify the resampling options once"
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
    * wildbootstrap(): the wild restricted efficient (WRE) bootstrap of
    * Davidson and MacKinnon, with boottest's conventions (plan 4.7, 5.9, 9.4)
    if `__wild' {
        if `__bs' {
            di as error "civreghdfe: wildbootstrap() and vce(bootstrap) are two bootstraps of the coefficients; choose one"
            exit 198
        }
        _ctools_hdfe_boot wildparse, spec(`"`__wildspec'"')
        local wb_reps = r(reps)
        local wb_seed `"`r(seed)'"'
        local wb_type `r(weighttype)'
        local wb_null = r(null)
        local wb_coefs `"`r(coefs)'"'
        local wb_ptype `r(ptype)'
        local wb_ptype_code = r(ptype_code)
        local wb_enumerate = r(enumerate)
        local wb_fedgp `r(fedgp)'
        local wb_ci = !r(noci)
        local wb_level = r(level)
        local wb_bootcluster `"`r(bootcluster)'"'
        local wb_dots = r(dots)
        local wb_projected = r(projected)
    }
    * fevce(): the names of the standard-error variables (__fese#__ for
    * __hdfe#__, v_se for a named v, or sestub#), all new
    local fe_se_names
    local fe_nrep_names
    local fe_lb_names
    local fe_ub_names
    local fe_targets_flags
    if `__fe' {
        forvalues __g = 1/`__t' {
            local fe_targets_flags `fe_targets_flags' `=("`fe_target_`__g''" != "")'
            if "`fe_target_`__g''" == "" continue
            if "`fe_sestub'" != "" {
                local __se `fe_sestub'`__g'
                local __nr `fe_sestub'`__g'_nrep
                local __lb `fe_sestub'`__g'_lb
                local __ub `fe_sestub'`__g'_ub
            }
            else if "`fe_target_`__g''" == "__hdfe`__g'__" {
                local __se __fese`__g'__
                local __nr __fenrep`__g'__
                local __lb __felb`__g'__
                local __ub __feub`__g'__
            }
            else {
                local __se `fe_target_`__g''_se
                local __nr `fe_target_`__g''_nrep
                local __lb `fe_target_`__g''_lb
                local __ub `fe_target_`__g''_ub
            }
            local fe_se_names `fe_se_names' `__se'
            if `fe_nreps' local fe_nrep_names `fe_nrep_names' `__nr'
            if `fe_ci' local fe_lb_names `fe_lb_names' `__lb'
            if `fe_ci' local fe_ub_names `fe_ub_names' `__ub'
        }
        _ctools_newvars `fe_se_names' `fe_nrep_names' `fe_lb_names' `fe_ub_names' `fe_dse'
    }

    * Parse VCE. Its cluster variables are marked out before the factor
    * variables are expanded, so that their levels come from the estimation
    * sample. cluster_var and cluster_var2 hold the dimensions as reghdfe
    * names them (an interaction as var1#var2).
    local vce_type = 0
    local cluster_var ""
    local cluster_var2 ""
    if `"`vce'"' != "" {
        _ctools_hdfe_vce, spec(`"`vce'"') touse(`touse') cmd(civreghdfe)
        local vce_kind `r(kind)'
        if "`vce_kind'" == "robust" {
            local vce_type = 1
        }
        else if "`vce_kind'" == "cluster" {
            * One-way (2) or two-way (3) clustering
            if !inlist(r(ndims), 1, 2) {
                di as error "vce(cluster) supports 1 or 2 cluster variables"
                exit 198
            }
            local vce_type = 1 + r(ndims)
            local cluster_var `r(name_1)'
            local cluster_vars_1 `r(vars_1)'
            if r(ndims) == 2 local cluster_var2 `r(name_2)'
        }
        else if inlist("`vce_kind'", "unadjusted", "ols") {
            local vce_type = 0
        }
        else {
            di as error "vce(`vce') not supported. Use vce(robust), vce(cluster varname), or vce(cluster var1 var2)"
            exit 198
        }
    }
    * The bootstrap's resampling units and strata (numbered in sorted order,
    * as bsample orders them). The GMM moment covariance of a replicate
    * clusters on the drawn copies when its clusters are the resampled ones.
    local bs_unitvar
    local bs_stratavar
    local bs_cluster_name
    local bs_s_units = 0
    if `__bsrun' {
        if `"`bs_cluster'"' != "" {
            _ctools_hdfe_vce, spec(`"cluster `bs_cluster'"') touse(`touse') cmd(civreghdfe)
            local bs_cluster_name `r(name_1)'
            tempvar bs_unit
            if r(group_1) quietly _ctools_group `r(vars_1)' if `touse', generate(`bs_unit')
            else quietly generate double `bs_unit' = `r(vars_1)' if `touse'
            local bs_unitvar `bs_unit'
            local bs_s_units = (`vce_type' == 2 & "`cluster_var'" == "`bs_cluster_name'")
        }
        if "`bs_strata'" != "" {
            markout `touse' `bs_strata', strok
            tempvar bs_str
            quietly _ctools_group `bs_strata' if `touse', generate(`bs_str')
            local bs_stratavar `bs_str'
        }
        if `__fe' & !`__bs' & inlist(`vce_type', 2, 3) & `"`bs_cluster'"' == "" {
            di as text "note: fevce(bootstrap) resamples observations, while vce() clusters the coefficients' standard errors"
        }
    }
    * wildbootstrap(): one clustering dimension; the bootstrap clusters are
    * the vce() cluster or bootcluster(), numbered in sorted order as
    * boottest numbers them
    local wb_vcevar
    local wb_bootvar
    local wb_bootname
    if `__wild' {
        if `vce_type' == 3 {
            di as error "civreghdfe: wildbootstrap() supports one clustering dimension"
            exit 198
        }
        if `vce_type' == 2 {
            tempvar wb_vce
            quietly _ctools_group `cluster_vars_1' if `touse', generate(`wb_vce')
            local wb_vcevar `wb_vce'
            local wb_bootname `cluster_var'
        }
        if `"`wb_bootcluster'"' != "" {
            _ctools_hdfe_vce, spec(`"cluster `wb_bootcluster'"') touse(`touse') cmd(civreghdfe)
            if r(ndims) != 1 {
                di as error "civreghdfe: wildbootstrap(bootcluster()) takes one variable"
                exit 198
            }
            local wb_bootname `r(name_1)'
            tempvar wb_boot
            quietly _ctools_group `r(vars_1)' if `touse', generate(`wb_boot')
            local wb_bootvar `wb_boot'
        }
        if `vce_type' != 2 & "`wb_bootvar'" == "" & "`weight'" == "fweight" {
            di as error "civreghdfe: wildbootstrap() draws one weight per observation without clusters; with fweights a row stands for several observations. Use bootcluster()"
            exit 198
        }
    }

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

        * Base/omitted placeholders are not sent to the plugin, which also
        * detects collinearity. The fvrevar variables follow the three name
        * lists in order.
        * A term listed twice is used once; an endogenous term that is also
        * an included exogenous regressor is exogenous, and an excluded
        * instrument that is also a regressor is not an excluded instrument,
        * as in ivreg2 (whose duplicates warning lists them)
        local __dups
        foreach __grp in endogvars exogvars instruments {
            local __grp_vars ""
            local __grp_names ""
            foreach __name of local `__grp'_names_all {
                gettoken __v __all_exp : __all_exp
                local __dup = `: list __name in __grp_names'
                if "`__grp'" == "endogvars" & `: list __name in exogvars_names_all' local __dup = 1
                if "`__grp'" == "instruments" & (`: list __name in exogvars_names_all' | ///
                    `: list __name in endogvars_names_all') local __dup = 1
                if `__dup' {
                    local __dups `__dups' `__name'
                    continue
                }
                local __grp_vars `__grp_vars' `__v'
                local __grp_names `__grp_names' `__name'
            }
            local `__grp'_names_all `__grp_names'
            _ctools_fvsplit, names(``__grp'_names_all') vars(`__grp_vars')
            local `__grp'_expanded "`r(vars)'"
            local `__grp'_coef "`r(names)'"
        }
    }

    * ivreg2's duplicates warning (the duplicates found by the first pass
    * when this is the second)
    local __dups : list uniq __dups
    if "`civrerun'" != "" local __dups `civdups'
    else if "`__dups'" != "" & "`nocollin'" == "" {
        di as text "Warning - duplicate variables detected"
        di as text "Duplicates:" _col(21) "`__dups'"
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * When no absorb specified and noconstant not requested, add explicit constant
    * (With absorb, the constant is implicitly absorbed by the FE)
    local has_cons = 0
    local G : word count `absorb'
    if `G' == 0 & "`noconstant'" == "" {
        tempvar _cons_var
        gen byte `_cons_var' = 1
        local exogvars_expanded `exogvars_expanded' `_cons_var'
        local exogvars_coef `exogvars_coef' _cons
        local exogvars_names_all `exogvars_names_all' _cons
        local has_cons = 1
    }

    * Count variables (using expanded counts)
    local K_endog : word count `endogvars_expanded'
    local K_exog : word count `exogvars_expanded'
    local K_inst_excl : word count `instruments_expanded'

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

    * Handle weights
    local has_weights = 0
    local weight_type = 0
    local weight_var ""
    if "`weight'" != "" {
        local has_weights = 1
        if "`weight'" == "aweight" local weight_type = 1
        else if "`weight'" == "fweight" local weight_type = 2
        else if "`weight'" == "pweight" local weight_type = 3
        * iweights enter as frequencies: N is their sum, as in ivreg2
        else if "`weight'" == "iweight" local weight_type = 2
        else {
            di as error "Weight type `weight' not supported"
            exit 198
        }
        * Evaluate the weight expression; marks out and validates the weights
        tempvar evaluated_weight
        local weight_expr = trim(`"`exp'"')
        if substr(`"`weight_expr'"', 1, 1) == "=" local weight_expr = substr(`"`weight_expr'"', 2, .)
        if "`weight'" == "iweight" {
            quietly gen double `evaluated_weight' = (`weight_expr') if `touse'
            markout `touse' `evaluated_weight'
            quietly count if `touse' & `evaluated_weight' < 0
            if r(N) {
                di as error "civreghdfe: negative iweights are not supported"
                exit 402
            }
        }
        else _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"= `weight_expr'"')
        local weight_var "`evaluated_weight'"
    }

    * Kernel-based VCEs, as in ivreg2: bw() alone gives the AC estimator
    * (autocorrelation-consistent, conditionally homoskedastic errors), bw()
    * with robust HAC, bw() with cluster(<time variable>) or dkraay() the
    * Driscoll-Kraay estimator (clusters = time periods, kernel over periods),
    * bw() with cluster(<panel variable> <time variable>) two-way clustering
    * with Driscoll-Kraay over time, and kiefer the AC estimator with the
    * truncated kernel over the whole panel (HAC with pweights). The data must
    * be tsset: kernel lags are differences of the tsset time variable within
    * panels, so gaps and unbalanced panels pair observations as L. does.
    local kiefer_val = ("`kiefer'" != "")
    local __dk_norob = 0
    * bw(): a real bandwidth, or auto (Newey and West's 1994 selection,
    * which ivreg2 offers for time series); ivreg2's errors for others
    local bw_auto = 0
    local bwchoice ""
    local bw = strtrim(`"`bw'"')
    if `"`bw'"' == "auto" {
        local bw_auto = 1
        local bw_val = 0
        if "`_dta[_TSpanel]'" != "" {
            di as error "Automatic bandwidth selection not available for panel data"
            exit 198
        }
        if `dkraay' > 0 {
            di as error "cannot use dkraay(.) and bw(.) options together"
            exit 198
        }
        if "`kiefer'" != "" {
            di as error "incompatible options: kiefer and bw"
            exit 198
        }
        * ivreg2 selects the bandwidth from the first-step IV/GMM residuals;
        * its LIML, k-class and b0() estimators never run the selection
        if "`liml'`b0'" != "" | `fuller' != 0 | `kclass' != 0 {
            di as error "bw(auto) is not available with LIML, k-class or b0() estimation"
            exit 198
        }
    }
    else if `"`bw'"' != "" {
        capture confirm number `bw'
        if _rc {
            di as error "bandwidth option bw() required for HAC-robust estimation"
            exit 102
        }
        if `bw' <= 0 {
            di as error "invalid bandwidth in option bw() - must be real > 0"
            exit 198
        }
        local bw_val = `bw'
    }
    else local bw_val = 0
    local kernel_type = 0
    local kernel_name ""
    local ts_ivar "`_dta[_TSpanel]'"
    local ts_tvar "`_dta[_TStvar]'"
    local ts_delta = 1
    if "`_dta[_TSdelta]'" != "" local ts_delta = `_dta[_TSdelta]'
    if `bw_val' < 0 | `dkraay' < 0 {
        di as error "invalid bandwidth in option bw() - must be real > 0"
        exit 198
    }
    * dkraay(#): bandwidth #, clusters on the time variable (checked below)
    if `dkraay' > 0 & `bw_val' == 0 local bw_val = `dkraay'
    if `"`kernel'"' != "" | `bw_val' > 0 | `bw_auto' {
        if `bw_val' == 0 & !`bw_auto' {
            di as error "bandwidth option bw() required for HAC-robust estimation"
            exit 102
        }
        local kernel_lower = lower(strtrim(`"`kernel'"'))
        if inlist("`kernel_lower'", "", "bar", "bartlett", "nwest") {
            local kernel_type = 1
            local kernel_name "Bartlett"
        }
        else if inlist("`kernel_lower'", "par", "parzen") {
            local kernel_type = 2
            local kernel_name "Parzen"
        }
        else if inlist("`kernel_lower'", "qua", "qs", "quadratic spectral", "quadratic-spectral", "quadraticspectral") {
            local kernel_type = 3
            local kernel_name "Quadratic Spectral"
        }
        else if inlist("`kernel_lower'", "tru", "trunc", "truncated") {
            local kernel_type = 4
            local kernel_name "Truncated"
        }
        else if inlist("`kernel_lower'", "thann", "tukey", "tukey-hanning") {
            local kernel_type = 5
            local kernel_name "Tukey-Hanning"
        }
        else if inlist("`kernel_lower'", "thamm", "tukey-hamming") {
            local kernel_type = 6
            local kernel_name "Tukey-Hamming"
        }
        * ivreg2 spells the Daniell kernel "Danielle" (its e(kernel))
        else if inlist("`kernel_lower'", "dan", "danielle") {
            local kernel_type = 7
            local kernel_name "Danielle"
        }
        else if inlist("`kernel_lower'", "ten", "tent") {
            local kernel_type = 8
            local kernel_name "Tent"
        }
        else {
            di as error "invalid kernel"
            exit 198
        }
        * bw(auto) is defined for the Bartlett, Parzen and Quadratic
        * Spectral kernels; with another, ivreg2 warns and uses bw=1
        if `bw_auto' & !inlist(`kernel_type', 1, 2, 3) {
            di as error "Error: kernel `kernel_name' not compatible with bw(auto)"
            local bw_auto = 0
            local bw_val = 1
            local bwchoice "Kernel `kernel_name' not compatible with bw(auto); bw=1 (default)"
        }
        else if `bw_auto' local bwchoice "Automatic bw selection according to Newey-West (1994)"
    }
    if `dkraay' > 0 {
        if "`ts_ivar'" == "" | "`ts_tvar'" == "" {
            di as error "invalid use of dkraay option - must use tsset panel data"
            exit 5
        }
        if `bw_val' != `dkraay' {
            di as error "cannot use dkraay(.) and bw(.) options together"
            exit 198
        }
        if `vce_type' == 3 | (`vce_type' == 2 & "`cluster_var'" != "`ts_tvar'") {
            di as error "invalid use of dkraay option - must cluster on `ts_tvar' (or omit cluster option)"
            exit 198
        }
        * (ivreg2's e(vce) says robust only if robust, cluster() or
        * pweights were given too)
        local __dk_norob = (`vce_type' == 0 & "`weight'" != "pweight")
        local vce_type = 2
        local cluster_var "`ts_tvar'"
    }
    if `kiefer_val' {
        if "`ts_ivar'" == "" | "`ts_tvar'" == "" {
            di as error "invalid use of kiefer option - must use tsset panel data"
            exit 5
        }
        if `vce_type' != 0 {
            di as error "incompatible options: kiefer and robust"
            exit 198
        }
        if `kernel_type' != 0 & `kernel_type' != 4 {
            di as error "incompatible options: kiefer and kernel(`kernel_name')"
            exit 198
        }
        if `bw_val' != 0 {
            di as error "incompatible options: kiefer and bw"
            exit 198
        }
        * bandwidth = the time span of the sample (set by the plugin)
        local kernel_type = 4
        local kernel_name "Truncated"
    }
    if `kernel_type' > 0 {
        if "`ts_tvar'" == "" {
            di as error "time variable not set, use tsset varname ..."
            exit 111
        }
        if "`weight'" == "fweight" {
            di as error "fweights not allowed (data are -tsset-)"
            exit 101
        }
        if `bw_val' == 1 & "`bwchoice'" == "" & inlist(`kernel_type', 1, 2, 5, 6) {
            di as text "Note: kernel=`kernel_name' and bw=1 implies zero lags used.  Standard errors and"
            di as text "      test statistics are not autocorrelation-consistent."
        }
    }

    * pweight implies robust VCE (standard Stata behavior)
    if "`weight'" == "pweight" & `vce_type' == 0 {
        local vce_type = 1
    }

    * A kernel with clustering: one-way clustering must be on the time
    * variable (Driscoll-Kraay), two-way clustering on the panel and time
    * variables (panel clusters plus Driscoll-Kraay over time)
    if `kernel_type' > 0 & `vce_type' == 2 {
        if "`cluster_var'" != "`ts_tvar'" {
            di as error "Error: cluster kernel-robust requires clustering on tsset time variable."
            di as error "       tsset time var=`ts_tvar'; cluster var=`cluster_var'"
            exit 198
        }
        local vce_type = 4
    }
    else if `kernel_type' > 0 & `vce_type' == 3 {
        if "`cluster_var'" == "`ts_tvar'" & "`cluster_var2'" == "`ts_ivar'" {
            local cluster_var "`ts_ivar'"
            local cluster_var2 "`ts_tvar'"
        }
        if "`ts_ivar'" == "" | "`cluster_var'" != "`ts_ivar'" | "`cluster_var2'" != "`ts_tvar'" {
            di as error "Error: cluster kernel-robust requires clustering on tsset panel & time vars."
            di as error "       tsset panel var=`ts_ivar'; tsset time var=`ts_tvar'; cluster vars=`cluster_var',`cluster_var2'"
            exit 198
        }
    }

    * tsset variables passed to the plugin: 2 = panel and time, 1 = time only
    local ts_mode = 0
    if `kernel_type' > 0 local ts_mode = cond("`ts_ivar'" != "", 2, 1)

    * Mark out additional variables from sample (touse already created earlier for fvrevar)
    markout `touse' `endogvars_expanded' `exogvars_expanded' `instruments_expanded'
    if `ts_mode' == 2 markout `touse' `ts_ivar' `ts_tvar'
    else if `ts_mode' == 1 markout `touse' `ts_tvar'
    * Count observations
    qui count if `touse'
    local N = r(N)

    if `N' == 0 {
        di as error "no observations"
        exit 2001
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
        di as text "Absorbed FE:           `absorb_names' (`G')"
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
        if `kernel_type' > 0 {
            local __bwtxt = cond(`kiefer_val', "time span", "`bw_val'")
            di as text "Kernel:                `kernel_name', bandwidth `__bwtxt' (" ///
                cond(`vce_type' == 0, "AC", cond(`vce_type' == 1, "HAC", "Driscoll-Kraay")) ")"
            di as text "tsset:                 panel `ts_ivar', time `ts_tvar' (delta `ts_delta')"
        }
        di as text "{hline 60}"
    }

    * Detect if cluster variable is nested within (same as) an absorb variable
    * When cluster var = FE var, that FE shouldn't reduce DOF for VCE
    * Also check for two-way clustering (vce_type == 3)
    local nested_fe_index = 0
    if (`vce_type' == 2 | `vce_type' == 3) & "`cluster_var'" != "" {
        local fe_idx = 1
        foreach fe_var of local absorb {
            * a slope factor is never nested in the clusters
            if ("`fe_var'" == "`cluster_var'" | "`fe_var'" == "`cluster_var2'") & "`ext_slope_`fe_idx''" == "" {
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

    * liml with kclass(): the k-class estimator, labelled LIML (ivreg2)
    local liml_kclass = ("`liml'" != "" & `kclass' != 0)
    if `liml_kclass' {
        local est_method = 3
        local kclass_val = `kclass'
        if "`verbose'" != "" {
            di as text "Estimation method:     k-class (k=`kclass')"
        }
    }
    else if "`liml'" != "" | `fuller' < 0 {
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

    * The bootstrap's restrictions (plan 4.6)
    if `__bsrun' {
        local __what = cond(`__bs', "vce(bootstrap)", "fevce(bootstrap)")
        if "`b0'" != "" {
            di as error "civreghdfe: `__what' is not available with b0()"
            exit 198
        }
        if `"`smatrix'`wmatrix'"' != "" {
            di as error "civreghdfe: `__what' is not available with smatrix() or wmatrix() yet"
            exit 198
        }
        if `__bs' & `kernel_type' > 0 {
            di as error "civreghdfe: vce(bootstrap) is not available with kernel(), bw(), dkraay() or kiefer; resample panels with cluster(panelvar) instead"
            exit 198
        }
        if inlist(`est_method', 4, 5) & (`kernel_type' > 0 | `vce_type' > 2) {
            di as error "civreghdfe: `__what' with gmm2s or cue needs a robust or one-way cluster vce(), whose moment covariance each replicate re-estimates"
            exit 198
        }
        if `__bs' & `est_method' == 5 {
            di as text "note: CUE is re-optimized in every bootstrap replicate"
        }
        if `psd_code' {
            di as error "civreghdfe: `__what' is not available with psd0 or psda"
            exit 198
        }
    }
    if `__wild' {
        if inlist(`est_method', 4, 5) {
            di as error "civreghdfe: wildbootstrap() is the WRE bootstrap of 2SLS, LIML, Fuller and k-class estimates; gmm2s and cue are not supported"
            exit 198
        }
        if `kernel_type' > 0 | `vce_type' > 2 {
            di as error "civreghdfe: wildbootstrap() supports the robust and one-way cluster vce(), not kernel(), bw(), dkraay() or kiefer"
            exit 198
        }
        if `vce_type' == 0 {
            di as error "civreghdfe: wildbootstrap() studentizes with vce(robust) or vce(cluster); specify one"
            exit 198
        }
        if "`b0'" != "" | `"`smatrix'`wmatrix'"' != "" | `psd_code' {
            di as error "civreghdfe: wildbootstrap() is not available with b0(), smatrix(), wmatrix(), psd0 or psda"
            exit 198
        }
        if `K_endog' < 1 {
            di as error "civreghdfe: wildbootstrap() needs an endogenous regressor"
            exit 198
        }
    }

    * bvclean, implied by partial() as in ivreg2 (not with the bootstraps,
    * whose replicate matrices keep the base levels)
    if "`bvclean'" != "" & (`__bsrun' | `__wild') {
        di as error "civreghdfe: bvclean is not available with vce(bootstrap), fevce() or wildbootstrap()"
        exit 198
    }
    local __bvclean = ("`bvclean'" != "" | (("`partial'" != "" | `__pcons_given') & !`__bsrun' & !`__wild'))
    if `any_slopes' & (`__bs' | `__fe' | `__wild') {
        di as error "civreghdfe: vce(bootstrap), fevce() and wildbootstrap() are not available with heterogeneous slopes in absorb()"
        exit 198
    }

    if `maxiter' < 1 | missing(`tolerance') | `tolerance' <= 0 {
        di as error "civreghdfe: maxiter() and tolerance() must be positive"
        exit 198
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

    * Kernel VCE parameters (parsed with the other VCE options above)
    * Kernel types: 0=none, 1=Bartlett, 2=Parzen, 3=Quadratic Spectral,
    * 4=Truncated, 5=Tukey-Hanning. Kiefer: truncated kernel, the plugin sets
    * the bandwidth to the time span of the estimation sample.
    scalar __civreghdfe_kernel = `kernel_type'
    * -1: bw(auto), the plugin selects the bandwidth
    scalar __civreghdfe_bw = cond(`bw_auto', -1, `bw_val')
    scalar __civreghdfe_kiefer = `kiefer_val'
    * tsset structure: the plugin gets the panel and time variables (ts_mode
    * 2) or the time variable (1) after the weight, and the tsset delta
    scalar __civreghdfe_ts_mode = `ts_mode'
    scalar __civreghdfe_tdelta = `ts_delta'

    * DOF adjustment options
    scalar __civreghdfe_small = ("`small'" != "" | "`absorb'" != "")
    scalar __civreghdfe_has_cons = ("`noconstant'" == "")
    scalar __civreghdfe_dofminus = `dofminus'
    scalar __civreghdfe_sdofminus = `sdofminus'
    scalar __civreghdfe_nopartialsmall = ("`nopartialsmall'" != "")
    scalar __civreghdfe_dof_type = `dof_type'
    scalar __civreghdfe_dof_clusters = `dof_clusters'
    scalar __civreghdfe_dof_continuous = strpos(" `dofmethods' ", " continuous ") > 0
    scalar __civreghdfe_keepsingletons = ("`keepsingletons'" != "")
    scalar __civreghdfe_center = ("`center'" != "")
    * psd0 / psda; ivreg2's e(S), e(W), e(condzz), e(condxx) and e(ccev)
    * from the plugin (matrices created below, over the instruments)
    scalar __civreghdfe_psd = `psd_code'
    scalar __civreghdfe_store_sw = 1

    * vce_type 4 = Driscoll-Kraay, whose clusters are the time periods
    scalar __civreghdfe_has_cluster = (`vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4)
    scalar __civreghdfe_vce_type = `vce_type'

    scalar __civreghdfe_has_b0 = 0

    * Parse and validate diagnostic test options and partial().
    * Each option names model terms: its varlist is expanded (factor variables,
    * time-series operators) and matched term by term against the expanded
    * model lists; base levels are skipped and duplicates ignored, as in ivreg2.
    * The plugin gets one 0/1 flag per plugin column (matrix
    * __civreghdfe_<option>_flags), maps the flags through its collinearity
    * checks, and errors if a tested column was dropped.
    * orthog() - test exogeneity of specified instruments, excluded or
    * included (an included one is treated as endogenous in the restricted
    * model, as ivreg2 does). The flags of the included ones are split off
    * once the partialled-out columns are known (below).
    local __orthog_terms
    local __endogtest_terms
    local __redundant_terms
    local n_orthog = 0
    local n_orthogx = 0
    local __on = 0
    local __kxo = 0
    tempname __oflags
    if "`orthog'" != "" {
        local __pcx _cons
        local __onames : list exogvars_names_all - __pcx
        _civreghdfe_optflags "`orthog'" "`__onames' `instruments_names_all'" `touse'
        if "`r(bad)'" != "" {
            di as error "Error: `r(bad)' listed in orthog() but does not appear as exogenous."
            exit 198
        }
        local __orthog_terms `r(terms)'
        if r(n) > 0 matrix `__oflags' = r(flags)
        local __on = r(n)
        local __kxo : word count `: list exogvars_coef - __pcx'
        local __kio : word count `instruments_coef'
        if `__on' > 0 & `__kio' > 0 {
            matrix __civreghdfe_orthog_flags = `__oflags'[1, `__kxo'+1...]
            forvalues __j = 1/`__kio' {
                if __civreghdfe_orthog_flags[1, `__j'] local ++n_orthog
            }
        }
        if "`verbose'" != "" {
            di as text "C-statistic test for: `orthog'"
        }
    }
    scalar __civreghdfe_n_orthog = `n_orthog'

    * endogtest() - test endogeneity of specified regressors
    local n_endogtest = 0
    if "`endogtest'" != "" {
        _civreghdfe_optflags "`endogtest'" "`endogvars_names_all'" `touse'
        if "`r(bad)'" != "" {
            di as error "endogtest(): `r(bad)' is not an endogenous variable"
            exit 198
        }
        local __endogtest_terms `r(terms)'
        local n_endogtest = r(n)
        if `n_endogtest' > 0 matrix __civreghdfe_endogtest_flags = r(flags)
        if "`verbose'" != "" {
            di as text "Endogeneity test for: `endogtest'"
        }
    }
    scalar __civreghdfe_n_endogtest = `n_endogtest'

    * redundant() - test redundancy of specified instruments
    local n_redundant = 0
    if "`redundant'" != "" {
        _civreghdfe_optflags "`redundant'" "`instruments_names_all'" `touse'
        if "`r(bad)'" != "" {
            di as error "redundant(): `r(bad)' is not an excluded instrument"
            exit 198
        }
        local __redundant_terms `r(terms)'
        local n_redundant = r(n)
        if `n_redundant' > 0 matrix __civreghdfe_redundant_flags = r(flags)
        if "`verbose'" != "" {
            di as text "Redundancy test for: `redundant'"
        }
    }
    scalar __civreghdfe_n_redundant = `n_redundant'

    * partial() - variables to partial out via FWL
    * These must be exogenous variables; partial_terms keeps the matched model
    * terms (base levels included) so they are left out of e(b)
    local n_partial = 0
    local partial_terms ""
    * Without absorb(), partial() also partials out the constant, as ivreg2
    * does (with absorb() the constant is absorbed)
    local partial_cons = 0
    local __pc _cons
    if `has_cons' & ("`partial'" != "" | `__pcons_given') local partial_cons = 1
    if "`partial'" != "" {
        _civreghdfe_optflags "`partial'" "`exogvars_names_all'" `touse'
        if "`r(bad)'" != "" {
            di as error "partial(): `r(bad)' is not an exogenous variable"
            exit 198
        }
        local n_partial = r(n)
        local partial_terms `r(terms)'
        if `n_partial' > 0 matrix __civreghdfe_partial_flags = r(flags)
        if "`verbose'" != "" {
            di as text "Partialling out: `partial'"
        }
    }
    if `partial_cons' {
        * the constant is the last exogenous column
        if `n_partial' == 0 matrix __civreghdfe_partial_flags = J(1, `K_exog', 0)
        matrix __civreghdfe_partial_flags[1, `K_exog'] = 1
        local ++n_partial
        local partial_terms `partial_terms' _cons
    }
    scalar __civreghdfe_n_partial = `n_partial'

    * orthog() on included instruments: one flag per exogenous plugin column
    * that is not partialled out (the constant is not a variable)
    if "`orthog'" != "" & `__on' > 0 & `__kxo' > 0 {
        local __xf
        local __k = 0
        foreach __v of local exogvars_coef {
            local ++__k
            local __pf = 0
            if `n_partial' > 0 local __pf = (el(__civreghdfe_partial_flags, 1, `__k') != 0)
            if "`__v'" == "_cons" continue
            local __of = el(`__oflags', 1, `__k')
            if `__pf' & `__of' {
                di as error "Error: `__v' listed in orthog() but does not appear as exogenous."
                exit 198
            }
            if !`__pf' local __xf `__xf' `__of'
        }
        if `has_cons' & !`partial_cons' local __xf `__xf' 0
        local __xf : list retokenize __xf
        local __xf : subinstr local __xf " " ",", all
        if "`__xf'" != "" matrix __civreghdfe_orthogx_flags = (`__xf')
        foreach __f in `: subinstr local __xf "," " ", all' {
            if `__f' local ++n_orthogx
        }
    }
    scalar __civreghdfe_n_orthogx = `n_orthogx'

    * The instruments left after partial(), in the plugin's order (included,
    * then excluded), for smatrix(), wmatrix() and the system of first-stage
    * equations
    local __znames
    local __k = 0
    foreach __v of local exogvars_coef {
        local ++__k
        if `n_partial' > 0 {
            if el(__civreghdfe_partial_flags, 1, `__k') != 0 continue
        }
        local __znames `__znames' `__v'
    }
    local __zx_ct : word count `__znames'
    local __znames `__znames' `instruments_coef'
    local __zn : word count `__znames'
    matrix __civreghdfe_S = J(max(`__zn', 1), max(`__zn', 1), .)
    matrix __civreghdfe_W = J(max(`__zn', 1), max(`__zn', 1), .)
    matrix __civreghdfe_ccev = J(1, max(`K_endog', 1), .)

    * smatrix()/wmatrix(): the user's matrix, its rows and columns picked by
    * the instruments' names (ivreg2's matsort), in the plugin's order
    foreach __mo in smatrix wmatrix {
        if `"``__mo''"' == "" continue
        local __mn = cond("`__mo'" == "smatrix", "S0", "W0")
        tempname __um __uin
        capture matrix `__uin' = ``__mo''
        if _rc {
            local __rc = _rc
            di as error "invalid matrix ``__mo'' in `__mo' option"
            exit `__rc'
        }
        matrix `__um' = J(`__zn', `__zn', 0)
        local __bad = 0
        forvalues __i = 1/`__zn' {
            local __ri = rownumb(`__uin', "`: word `__i' of `__znames''")
            local __ci = colnumb(`__uin', "`: word `__i' of `__znames''")
            if missing(`__ri') | missing(`__ci') local __bad = 1
        }
        if `__bad' {
            di as error "-`__mo'- option error: supplied matrix columns/rows do not match IV list"
            exit 198
        }
        forvalues __j = 1/`__zn' {
            local __cj = colnumb(`__uin', "`: word `__j' of `__znames''")
            forvalues __i = 1/`__zn' {
                local __ri = rownumb(`__uin', "`: word `__i' of `__znames''")
                matrix `__um'[`__i', `__j'] = `__uin'[`__ri', `__cj']
            }
        }
        if !issymmetric(`__um') {
            * (ivreg2's message for both matrices)
            di as error "-smatrix- option error: supplied matrix is not symmetric"
            exit 198
        }
        matrix __civreghdfe_`__mn' = `__um'
    }

    * first, ffirst, sfirst, savefirst, savesfirst: ivreg2's first-stage
    * statistics and the system of first-stage and reduced-form equations
    * (dimensions before collinearity; the plugin fills what it keeps)
    local __want_ff = ("`first'`ffirst'`savefirst'`sfirst'`savesfirst'" != "" & `K_endog' > 0 & `K_inst_excl' > 0)
    * (also for the reduced form alone: its stored equation's e(S) is a
    * block of the system's)
    local __want_sys = (`__want_ff' | ("`rf'`saverf'" != "" & `K_endog' > 0 & `K_inst_excl' > 0))
    scalar __civreghdfe_ffirst = `__want_sys'
    if `__want_sys' {
        matrix __civreghdfe_firstmat = J(21, `K_endog', .)
        local __sfd = (`K_endog' + 1) * `__zn'
        matrix __civreghdfe_sfb = J(1, `__sfd', .)
        matrix __civreghdfe_sfV = J(`__sfd', `__sfd', .)
        matrix __civreghdfe_sfS = J(`__sfd', `__sfd', .)
    }

    * b0(): one coefficient per regressor, matched by name to e(b)'s
    * [endogenous, exogenous] columns (ivreg2's matsort), the partialled-out
    * regressors excluded
    if "`b0'" != "" {
        local __b0_names `endogvars_coef'
        foreach __v of local exogvars_coef {
            if !`: list __v in partial_terms' local __b0_names `__b0_names' `__v'
        }
        local __b0_k : word count `__b0_names'
        local __b0_bad = (colsof(`b0') != `__b0_k')
        matrix __civreghdfe_b0 = J(1, max(`__b0_k', 1), 0)
        local __j = 0
        foreach __v of local __b0_names {
            local ++__j
            local __c = colnumb(`b0', "`__v'")
            if missing(`__c') local __b0_bad = 1
            else matrix __civreghdfe_b0[1, `__j'] = `b0'[1, `__c']
        }
        if `__b0_bad' {
            di as error "-b0- option error: supplied b0 columns do not match regressor list"
            _ctools_cleanup civreghdfe
            exit 198
        }
        scalar __civreghdfe_has_b0 = 1
    }

    * Create output matrices (partialled-out regressors are not in e(b))
    local K_total = `K_endog' + `K_exog' - `n_partial'
    if `K_total' < 1 {
        di as error "partial(): no regressors left to estimate"
        exit 198
    }
    tempname b_temp V_temp

    matrix `b_temp' = J(1, `K_total', 0)
    matrix `V_temp' = J(`K_total', `K_total', 0)
    matrix __civreghdfe_b = `b_temp'
    matrix __civreghdfe_V = `V_temp'

    * Build variable list for plugin
    * Order: depvar, endogvars, exogvars, excluded instruments, absorb, [cluster], [cluster2],
    *        [weight], [tsset panel], [tsset time], residual
    * The exogenous regressors are also instruments (K_iv counts them) but are
    * sent once; the plugin copies them into the instrument matrix.
    * Use expanded variables (factor variables converted to temp numeric vars)
    * The plugin numbers a single cluster variable itself (strings
    * included); an interaction is numbered here first
    forvalues dim = 1/2 {
        local cluster_name = cond(`dim' == 1, "`cluster_var'", "`cluster_var2'")
        local cluster_plugin_`dim' `cluster_name'
        if strpos("`cluster_name'", "#") {
            local cluster_parts : subinstr local cluster_name "#" " ", all
            tempvar cluster_id_`dim'
            quietly _ctools_group `cluster_parts' if `touse', generate(`cluster_id_`dim'')
            local cluster_plugin_`dim' `cluster_id_`dim''
        }
    }
    local plugin_vars `depvar_var' `endogvars_expanded' `exogvars_expanded' `instruments_expanded' `absorb'
    if `vce_type' == 2 | `vce_type' == 3 | `vce_type' == 4 {
        local plugin_vars `plugin_vars' `cluster_plugin_1'
    }
    if `vce_type' == 3 {
        local plugin_vars `plugin_vars' `cluster_plugin_2'
    }
    if `has_weights' {
        local plugin_vars `plugin_vars' `weight_var'
    }
    if `ts_mode' == 2 {
        local plugin_vars `plugin_vars' `ts_ivar' `ts_tvar'
    }
    else if `ts_mode' == 1 {
        local plugin_vars `plugin_vars' `ts_tvar'
    }
    * Heterogeneous slopes: the slope variable of each slope factor follows
    * the tsset columns; __civreghdfe_slope_# flags the slope factors and
    * __civreghdfe_absvar_# names the absorbed effect each factor belongs to
    scalar __civreghdfe_num_slopes = 0
    forvalues g = 1/`G' {
        scalar __civreghdfe_slope_`g' = ("`ext_slope_`g''" != "")
        scalar __civreghdfe_absvar_`g' = `ext_absvar_`g''
        if "`ext_slope_`g''" != "" {
            local plugin_vars `plugin_vars' `ext_slope_`g''
            scalar __civreghdfe_num_slopes = __civreghdfe_num_slopes + 1
        }
    }

    tempvar stored_resid
    quietly gen double `stored_resid' = .
    local plugin_vars `plugin_vars' `stored_resid'
    local resid_idx : word count `plugin_vars'
    scalar __civreghdfe_resid_idx = `resid_idx'

    * vce(bootstrap), fevce(bootstrap): the resampling units and strata
    scalar __civreghdfe_bs_unit_idx = 0
    scalar __civreghdfe_bs_strata_idx = 0
    scalar __civreghdfe_bs_s_units = `bs_s_units'
    if "`bs_unitvar'" != "" {
        local plugin_vars `plugin_vars' `bs_unitvar'
        scalar __civreghdfe_bs_unit_idx = `: word count `plugin_vars''
    }
    if "`bs_stratavar'" != "" {
        local plugin_vars `plugin_vars' `bs_stratavar'
        scalar __civreghdfe_bs_strata_idx = `: word count `plugin_vars''
    }

    * The sample marker goes last, in place of if `touse': the plugin reads
    * it in parallel, where Stata answers an if one observation at a time
    local plugin_vars `plugin_vars' `touse'
    scalar __civreghdfe_touse = 1

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record Stata pre-processing time
    _ctools_lap t_stata_pre __ct_lap97 `__ct_lap97'

    * Call the C plugin
    quietly _ctools_clock __ct_lap99
    capture noisily plugin call ctools_plugin `plugin_vars', ///
        "civreghdfe `threads_code' iv_regression`=cond(`__bsrun' | `__wild', " retain", "")'"
    if _rc {
        local __rc = _rc
        * A variable named in orthog()/endogtest()/redundant() was dropped
        capture confirm scalar __civreghdfe_dropped_opt
        if _rc == 0 {
            local __o = scalar(__civreghdfe_dropped_opt)
            local __opt : word `__o' of orthog endogtest redundant orthog
            local __zexogl : list __znames - instruments_coef
            local __lst : word `__o' of instruments_coef endogvars_coef instruments_coef __zexogl
            local __v : word `=scalar(__civreghdfe_dropped_col)' of ``__lst''
            di as error "`__opt'(): `__v' was dropped (collinear, or absorbed by the fixed effects) and cannot be tested"
        }
        _ctools_cleanup civreghdfe
        exit `__rc'
    }

    * ivreg2 removes collinear variables before it partials out the fixed
    * effects, in its own order (ivreg2's CheckDupsCollin), and an endogenous
    * regressor collinear with the instruments becomes exogenous. The plugin
    * finds collinearity after partialling, which includes the collinearity
    * that the fixed effects create. When it dropped anything, ivreg2's rules
    * are applied to the variables over the estimation sample; if they drop
    * or reclassify a term, the model is estimated again with ivreg2's lists
    * (dropped terms omitted), as ivreg2 then estimates it. nocollin skips
    * the rules, as in ivreg2.
    local __plugin_dropped = 0
    capture local __plugin_dropped = scalar(__civreghdfe_num_collinear) + scalar(__civreghdfe_z_dropped)
    capture local __plugin_dropped = `__plugin_dropped' + scalar(__civreghdfe_endo_in_z)
    if "`civrerun'`nocollin'" == "" & `__plugin_dropped' > 0 & `__plugin_dropped' < . {
        tempvar __civ_sample
        quietly gen byte `__civ_sample' = `touse' & !missing(`stored_resid')
        local __cx_exog : list exogvars_expanded - _cons_var
        local __cx_exogn : list exogvars_coef - __pc
        _civreghdfe_collin, touse(`__civ_sample') endo(`endogvars_coef') endovars(`endogvars_expanded') ///
            inexog(`__cx_exogn') inexogvars(`__cx_exog') exexog(`instruments_coef') ///
            exexogvars(`instruments_expanded') wvar(`weight_var') wtype(`weight') `noconstant'
        if r(changed) {
            local __civ_collin `r(collin)'
            local __civ_ecollin `r(ecollin)'
            foreach __g in endo inexog exexog {
                local __mata_`__g' `r(`__g')'
            }
            * The new lists, with the base levels where the first pass had
            * them; reclassified endogenous regressors leave the endogenous
            * list and come first among the included exogenous regressors,
            * as ivreg2 orders them
            local __n_ec : word count `__civ_ecollin'
            foreach __pair in "endo endogvars" "inexog exogvars" "exexog instruments" {
                gettoken __g __o : __pair
                local __o = strtrim("`__o'")
                local __new `__mata_`__g''
                local __lead
                if "`__g'" == "inexog" {
                    forvalues __j = 1/`__n_ec' {
                        gettoken __nn __new : __new
                        local __lead `__lead' `__nn'
                    }
                }
                local __full
                _ctools_fvsplit, names(``__o'_names_all')
                local __bflags `r(base)'
                foreach __name of local `__o'_names_all {
                    gettoken __bf __bflags : __bflags
                    if "`__name'" == "_cons" continue
                    if "`__g'" == "endo" & `: list __name in __civ_ecollin' continue
                    if `__bf' local __full `__full' `__name'
                    else {
                        gettoken __nn __new : __new
                        local __full `__full' `__nn'
                    }
                }
                local __rerun_`__g' `__lead' `__full'
            }
            if "`__civ_ecollin'" != "" {
                di as text "Warning - endogenous variable(s) collinear with instruments"
                di as text "Vars now exogenous:" _col(21) "`__civ_ecollin'"
            }
            if "`__civ_collin'" != "" {
                di as text "Warning - collinearities detected"
                di as text "Vars dropped:" _col(21) "`__civ_collin'"
            }
            local __wgt
            if "`weight'" != "" local __wgt [`weight'=`weight_expr']
            local __rerun_spec `depvar' `__rerun_inexog'
            if `is_iv_model' local __rerun_spec `__rerun_spec' (`__rerun_endo' = `__rerun_exexog')
            _ctools_cleanup civreghdfe
            if `"`__civ_opts'"' == "" local __civ_opts ","
            civreghdfe `__rerun_spec' if `__civ_sample' `__wgt' `__civ_opts' civrerun ///
                civcollin(`__civ_collin') civecollin(`__civ_ecollin') civdups(`__dups')
            ereturn local cmdline `"`cmdline'"'
            exit
        }
    }
    local __collin `civcollin'
    local __ecollin `civecollin'

    * partial() variables the plugin dropped as collinear or absorbed by the FEs
    if `n_partial' > 0 {
        local __k = 0
        foreach __v of local exogvars_coef {
            local ++__k
            if el(__civreghdfe_partial_flags, 1, `__k') == 2 {
                di as text "note: `__v' omitted from partial() because of collinearity"
            }
        }
    }
    quietly replace `touse' = 0 if missing(`stored_resid')
    if "`residuals'" != "" {
        if "`residuals2'" != "" capture drop `residuals'
        rename `stored_resid' `residuals'
    }
    _ctools_lap t_plugin __ct_lap99 `__ct_lap99'
    quietly _ctools_clock __ct_lap98

    * If noreturn specified, skip all result retrieval and e() posting
    if "`noreturn'" != "" {
        di as text "(noreturn specified - skipping result storage)"
        if `__bsrun' | `__wild' capture _civreghdfe_bs_release
        _ctools_cleanup civreghdfe
        exit 0
    }

    * wildbootstrap(): on the final fit, before fevce()'s pairs replicates.
    * The tested coefficients are named by the plugin's columns: the
    * endogenous, then the exogenous regressors that were not partialled out.
    local wb_names `endogvars_coef'
    foreach __v of local exogvars_coef {
        if !`: list __v in partial_terms' local wb_names `wb_names' `__v'
    }
    if `__wild' {
        capture noisily _civreghdfe_wild, reps(`wb_reps') type(`wb_type') null(`wb_null') ///
            ptype(`wb_ptype_code') enumerate(`wb_enumerate') ///
            partialled(`=("`wb_fedgp'" == "partialled")') ci(`wb_ci') level(`wb_level') ///
            dots(`wb_dots') vcetype(`vce_type') names(`wb_names') coefs(`wb_coefs') ///
            seed(`"`wb_seed'"') threads(`"`threads_code'"') vcevar(`wb_vcevar') ///
            bootvar(`wb_bootvar') keep(`__fe') projected(`wb_projected')
        local wb_rc = _rc
        if `wb_rc' {
            capture _civreghdfe_bs_release
            _ctools_cleanup civreghdfe
            exit `wb_rc'
        }
        local wb_rngstate `"`r(rngstate)'"'
        local wb_test `r(test)'
        local wb_B = r(reps)
        local wb_enum = r(enumerated)
        local wb_C = r(C)
        local wb_path = r(path)
    }

    * vce(bootstrap), fevce(bootstrap): the replicates on the final fit (after
    * ivreg2's collinearity rerun, which bootstraps in its own pass)
    local fe_se_tmp
    local fe_nrep_tmp
    local fe_lb_tmp
    local fe_ub_tmp
    local fe_dse_tmp
    if `__bsrun' {
        if `__fe' {
            foreach __n of local fe_se_names {
                tempvar __t1
                quietly generate double `__t1' = .
                local fe_se_tmp `fe_se_tmp' `__t1'
            }
            foreach __n of local fe_nrep_names {
                tempvar __t1
                quietly generate long `__t1' = .
                local fe_nrep_tmp `fe_nrep_tmp' `__t1'
            }
            foreach __n of local fe_lb_names {
                tempvar __t1
                quietly generate double `__t1' = .
                local fe_lb_tmp `fe_lb_tmp' `__t1'
            }
            foreach __n of local fe_ub_names {
                tempvar __t1
                quietly generate double `__t1' = .
                local fe_ub_tmp `fe_ub_tmp' `__t1'
            }
            if "`fe_dse'" != "" {
                tempvar fe_dse_tmp
                quietly generate double `fe_dse_tmp' = .
            }
            matrix __civreghdfe_bs_fe_targets = J(1, `G', 0)
            forvalues __g = 1/`G' {
                matrix __civreghdfe_bs_fe_targets[1, `__g'] = `: word `__g' of `fe_targets_flags''
            }
        }
        local bs_savetmp
        if `"`bs_save_file'"' != "" tempfile bs_savetmp
        capture noisily _civreghdfe_bootstrap, reps(`bs_reps') seed(`"`bs_seed'"') ///
            draws(`bs_draws') dots(`bs_dots') level(`bs_level') mse(`bs_mse') ///
            threads(`"`threads_code'"') coef(`__bs') fe(`__fe') kb(`=colsof(__civreghdfe_b)') ///
            fevars(`fe_se_tmp' `fe_nrep_tmp' `fe_lb_tmp' `fe_ub_tmp' `fe_dse_tmp') ///
            nfe(`: word count `fe_se_names'') nreps(`=`__fe' & `fe_nreps'') ci(`=`__fe' & `fe_ci'') ///
            dse(`=("`fe_dse'" != "")') bsdraws(`bs_bsdraws') savefile(`"`bs_savetmp'"') ///
            every(`bs_save_every') `=cond("`bs_strata'" != "" & "`bs_unitvar'" == "", "strataobs", "")'
        local bs_rc = _rc
        local bs_rngstate `"`r(rngstate)'"'
        if !`bs_rc' {
            local bs_nunits = scalar(__civreghdfe_bs_U)
            local bs_nstrata = scalar(__civreghdfe_bs_S)
        }
        if `bs_rc' {
            capture _civreghdfe_bs_release
            _ctools_cleanup civreghdfe
            exit `bs_rc'
        }
    }

    * Retrieve results
    local N_used = __civreghdfe_N
    local df_r = __civreghdfe_df_r
    * iweights: N is the integer part of their sum, as ivreg2 posts it, and
    * so are the residual dof
    if "`weight'" == "iweight" {
        local df_r = `df_r' - (`N_used' - floor(`N_used'))
        local N_used = floor(`N_used')
    }
    local df_a = __civreghdfe_df_a
    local df_a_for_vce = __civreghdfe_df_a_for_vce
    local K = __civreghdfe_K
    local rss = __civreghdfe_rss
    local tss = __civreghdfe_tss
    local r2 = __civreghdfe_r2
    local rmse = __civreghdfe_rmse

    * Kiefer's bandwidth is the time span T = tmax - tmin + 1 of the
    * estimation sample (ivreg2's bw for kiefer)
    if `kiefer_val' {
        capture local bw_val = scalar(__civreghdfe_tspan) + 1
    }
    * bw(auto): the bandwidth the plugin selected
    if `bw_auto' {
        capture local bw_val = scalar(__civreghdfe_bw_auto)
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

    * Read nested FE detection from C plugin (computed during plugin call):
    * an effect nested in either cluster dimension (or in the Driscoll-Kraay
    * time clusters)
    if inlist(`vce_type', 2, 3, 4) {
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
        capture local N_clust1 = __civreghdfe_N_clust1
        if _rc != 0 local N_clust1 = `N_clust'
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
        * Skip partialled-out terms (partial_terms holds the expanded model names)
        if !`: list v in partial_terms' {
            local varnames `varnames' `v'
        }
    }
    if "`exogvars'" != "" | `has_cons' {
        foreach v of local exogvars_names_all {
            if !`: list v in partial_terms' {
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

    * e(b) and e(V) in ivreghdfe's layout: zeros for the base/omitted
    * factor-variable levels, which were not sent to the plugin, and the o.
    * prefix for the columns it found collinear.
    _ctools_fvsplit, names(`varnames_full')
    local base_flags `r(base)'
    local K_active : word count `varnames'
    local K_full : word count `varnames_full'
    if `K_full' > `K_active' {
        local K_mapped = 0
        foreach __f of local base_flags {
            if !`__f' local K_mapped = `K_mapped' + 1
        }
        if `K_mapped' != `K_active' | colsof(`b_temp') != `K_active' {
            di as error "civreghdfe: cannot map the estimated coefficients onto the factor-variable terms"
            exit 498
        }
    }

    * Collinearity flags are in the plugin's [exog, endog] order; the columns
    * are in [endog, exog] order (ivreghdfe's convention)
    capture confirm scalar __civreghdfe_num_collinear
    local __any_collinear = (_rc == 0)
    if `__any_collinear' local __any_collinear = (__civreghdfe_num_collinear > 0)
    local collinear_flags ""
    forvalues __k = 1/`K_active' {
        local __is_collin = 0
        if `__any_collinear' {
            local __c = cond(`__k' <= `K_endog', `K_exog' + `__k', `__k' - `K_endog')
            capture local __is_collin = (__civreghdfe_collinear_`__c' == 1)
            if _rc != 0 local __is_collin = 0
        }
        local collinear_flags `collinear_flags' `__is_collin'
    }

    tempname b_full V_full
    _ctools_ebv `b_full' `V_full', names(`varnames_full') base(`base_flags') ///
        bmatrix(`b_temp') vmatrix(`V_temp') collinear(`collinear_flags') colnames(varnames)
    matrix `b_temp' = `b_full'
    matrix `V_temp' = `V_full'
    matrix rownames `b_temp' = `depvar'
    * ivreg2's "0" varlists come from these names (base levels and omitted
    * terms included)
    local __b0names : colnames `b_temp'
    * bvclean (implied by partial(), as in ivreg2): e(b) and e(V) without
    * the base levels and omitted terms; Stata names the factor levels as
    * ivreg2's matrix colnames does
    if `__bvclean' {
        tempname __om
        _ms_omit_info `b_temp'
        matrix `__om' = r(omit)
        local __keep
        local __keepn
        forvalues __j = 1/`=colsof(`b_temp')' {
            if !`__om'[1, `__j'] {
                local __keep `__keep' `__j'
                local __keepn `__keepn' `: word `__j' of `__b0names''
            }
        }
        if "`__keep'" != "" {
            mata: st_matrix(st_local("b_temp"), st_matrix(st_local("b_temp"))[1, strtoreal(tokens(st_local("__keep")))])
            mata: st_matrix(st_local("V_temp"), st_matrix(st_local("V_temp"))[strtoreal(tokens(st_local("__keep"))), strtoreal(tokens(st_local("__keep")))])
            matrix colnames `b_temp' = `__keepn'
            matrix colnames `V_temp' = `__keepn'
            matrix rownames `V_temp' = `__keepn'
            matrix rownames `b_temp' = `depvar'
        }
    }
    * saving(): the replicate file's column of each e(b) term (0 for terms
    * fixed at zero), by the layout of e(b)
    if `"`bs_save_file'"' != "" {
        tempname bs_map __smv __smx __smV
        local __Kb = colsof(__civreghdfe_b)
        matrix `__smv' = J(1, `__Kb', 0)
        forvalues __k = 1/`__Kb' {
            matrix `__smv'[1, `__k'] = `__k'
        }
        matrix `__smv' = diag(`__smv')
        _ctools_ebv `__smx' `__smV', names(`varnames_full') base(`base_flags') ///
            bmatrix(__civreghdfe_b) collinear(`collinear_flags') vmatrix(`__smv')
        matrix `bs_map' = vecdiag(`__smV')
    }

    * Compute Wald F-statistic BEFORE posting (b' * V^{-1} * b / K)
    tempname Vinv b_col wald_temp wald_stat
    * Not `omitted': that local holds the omitted display option
    tempname slope_b slope_V omit_flags
    _ms_omit_info `b_temp'
    matrix `omit_flags' = r(omit)
    local slope_indices
    local names : colnames `b_temp'
    local j = 0
    foreach name of local names {
        local ++j
        if "`name'" != "_cons" & !`omit_flags'[1,`j'] local slope_indices `slope_indices' `j'
    }
    local model_df = 0
    local F_wald = .
    matrix `wald_stat' = J(1,1,.)
    if "`slope_indices'" != "" {
        mata: st_matrix(st_local("slope_b"), st_matrix(st_local("b_temp"))[., strtoreal(tokens(st_local("slope_indices")))])
        mata: st_matrix(st_local("slope_V"), st_matrix(st_local("V_temp"))[strtoreal(tokens(st_local("slope_indices"))), strtoreal(tokens(st_local("slope_indices")))])
        * ivreg2 divides the omnibus Wald by the active slope count,
        * including when clustering makes the covariance rank deficient.
        local model_df : word count `slope_indices'
        matrix `Vinv' = syminv(`slope_V')
        matrix `b_col' = `slope_b''
        matrix `wald_temp' = `Vinv' * `b_col'
        matrix `wald_stat' = `slope_b' * `wald_temp'
        local F_wald = `wald_stat'[1,1] / `model_df'
        * No testable restriction (a perfect fit leaves V at rounding
        * noise): missing, as ivreg2's test gives
        if diag0cnt(`Vinv') == `model_df' local F_wald = .
    }

    * vce(bootstrap): e(V) is the bootstrap covariance in e(b)'s layout, the
    * base VCE's goes to e(V_modelbased), and the Wald test of the slopes
    * is chi2 (plan 4.4)
    local bs_chi2 = .
    if `__bs' {
        tempname V_mb V_bs __bx __bd __D __Vb
        matrix `V_mb' = `V_temp'
        matrix `__Vb' = (__civreghdfe_bs_V + __civreghdfe_bs_V') / 2
        local __ebv names(`varnames_full') base(`base_flags') bmatrix(__civreghdfe_b) ///
            collinear(`collinear_flags')
        _ctools_ebv `__bx' `V_bs', `__ebv' vmatrix(`__Vb')
        foreach __v in b se bias z0 {
            tempname bs_`__v'
            matrix `__D' = diag(__civreghdfe_bs_`__v')
            _ctools_ebv `__bx' `__bd', `__ebv' vmatrix(`__D')
            matrix `bs_`__v'' = vecdiag(`__bd')
        }
        foreach __v in ci_normal ci_percentile ci_bc {
            tempname bs_`__v'
            forvalues __r = 1/2 {
                matrix `__D' = diag(__civreghdfe_bs_`__v'[`__r', 1...])
                _ctools_ebv `__bx' `__bd', `__ebv' vmatrix(`__D')
                matrix `bs_`__v'' = nullmat(`bs_`__v'') \ vecdiag(`__bd')
            }
        }
        * As _bs_sum: a term that is zero in every replicate has no
        * bias-corrected interval or z0
        forvalues __j = 1/`=colsof(`b_temp')' {
            if `V_bs'[`__j', `__j'] == 0 & `b_temp'[1, `__j'] == 0 {
                matrix `bs_z0'[1, `__j'] = .
                matrix `bs_ci_bc'[1, `__j'] = .
                matrix `bs_ci_bc'[2, `__j'] = .
            }
        }
        if "`slope_indices'" != "" {
            tempname __Vs __Vsi __W
            mata: st_matrix(st_local("__Vs"), st_matrix(st_local("V_bs"))[strtoreal(tokens(st_local("slope_indices"))), strtoreal(tokens(st_local("slope_indices")))])
            matrix `__Vsi' = syminv(`__Vs')
            matrix `__W' = `slope_b' * `__Vsi' * `slope_b''
            local bs_chi2 = `__W'[1,1]
            if diag0cnt(`__Vsi') == `model_df' local bs_chi2 = .
        }
        matrix `V_temp' = `V_bs'
    }

    * Post results
    ereturn clear
    local depname_use = cond("`depname'" != "", "`depname'", "`depvar'")
    * No buildfvinfo, as in ivreghdfe: the constant is partialled out, and
    * without _cons in e(b) Stata's estimability check rejects estimable
    * margins such as dydx() over factor interactions.
    ereturn post `b_temp' `V_temp', esample(`touse') depname(`depname_use')
    * The weights and margins permissions, which postestimation commands reuse
    if `has_weights' {
        ereturn local wtype "`weight'"
        ereturn local wexp `"=`weight_expr'"'
    }
    ereturn local marginsok "XB default"
    if `K_endog' > 0 ereturn local marginsnotok "Residuals SCores"

    * Store scalars
    ereturn scalar N = `N_used'
    if "`small'`absorb'" != "" {
        * (z statistics with the bootstrap covariance)
        if !`__bs' ereturn scalar df_r = `df_r'
        ereturn local small "small"
    }
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
        if `vce_type' == 3 {
            ereturn scalar N_clust1 = `N_clust1'
            ereturn scalar N_clust2 = `N_clust2'
        }
        else {
            ereturn scalar N_clust1 = `N_clust'
        }
    }

    * Store first-stage F statistics. The plugin numbers them over the
    * endogenous regressors it kept; e(F_first#) follows the model's
    * endogenous regressors, missing for one dropped as collinear.
    * Positions are those of the endogenous terms, omitted ones included
    * (base levels are not terms)
    local __fs_terms
    local __fs_omit
    _ctools_fvsplit, names(`endogvars_names_all')
    local __bfl `r(base)'
    foreach __n of local endogvars_names_all {
        gettoken __b __bfl : __bfl
        if !`__b' {
            local __fs_terms `__fs_terms' `__n'
            local __fs_omit `__fs_omit' 0
        }
        else if regexm("`__n'", "o\.") & !regexm("`__n'", "b\.") {
            local __fs_terms `__fs_terms' `__n'
            local __fs_omit `__fs_omit' 1
        }
    }
    local __n_fs : word count `__fs_terms'
    local __fs_kept = 0
    local __fs_col = 0
    forvalues e = 1/`__n_fs' {
        local __fs_index_`e' = 0
        if `: word `e' of `__fs_omit'' {
            ereturn scalar F_first`e' = .
            continue
        }
        local ++__fs_col
        if `: word `__fs_col' of `collinear_flags'' {
            ereturn scalar F_first`e' = .
            continue
        }
        local ++__fs_kept
        local __fs_index_`e' = `__fs_kept'
        capture scalar `temp_F' = __civreghdfe_F1_`__fs_kept'
        if _rc == 0 {
            ereturn scalar F_first`e' = `temp_F'
        }
    }

    * Store macros
    ereturn local predict "civreghdfe_p"
    ereturn local resid "`residuals'"
    ereturn local cmd "civreghdfe"
    ereturn local cmdline `"`cmdline'"'
    ereturn local depvar "`depname_use'"
    ereturn local endogvars "`endogvars'"
    ereturn local exogvars "`exogvars'"
    ereturn local instruments "`instruments'"
    ereturn local absorb "`absorb_names'"
    * ivreg2's duplicates, dropped collinear terms and reclassified
    * endogenous regressors
    ereturn local dups "`__dups'"
    ereturn local collin "`__collin'"
    ereturn local ecollin "`__ecollin'"
    ereturn scalar nocollin = ("`nocollin'" != "")
    * partial(): the partialled-out terms, which e(b) leaves out (predict needs them)
    ereturn scalar partial_ct = `n_partial'
    * (without _cons, as ivreg2 lists them)
    if `n_partial' > 0 ereturn local partial "`: list partial_terms - __pc'"

    * ivreg2's e(model): liml (LIML, Fuller's LIML, liml with kclass()),
    * kclass, cue (CUE, b0()), gmm2s, gmmw (wmatrix(), set below), iv, or
    * ols without endogenous regressors. e(title) is set with the display
    * information below.
    if "`b0'" != "" ereturn local model "cue"
    else if `liml_kclass' | inlist(`est_method', 1, 2) ereturn local model "liml"
    else if `est_method' == 3 ereturn local model "kclass"
    else if `est_method' == 4 ereturn local model "gmm2s"
    else if `est_method' == 5 ereturn local model "cue"
    else ereturn local model = cond(`K_endog' == 0, "ols", "iv")
    if `liml_kclass' {
        ereturn scalar kclass = `kclass_val'
        * (ivreg2 computes no LIML lambda when k is given)
        ereturn scalar lambda = .
    }
    else if inlist(`est_method', 1, 2) {
        capture ereturn scalar lambda = __civreghdfe_lambda
        capture ereturn scalar kclass = __civreghdfe_k
        if `est_method' == 2 ereturn scalar fuller = `fuller_val'
    }
    else if `est_method' == 3 ereturn scalar kclass = `kclass_val'

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
        * Both variables, as ivreghdfe posts them
        ereturn local clustvar "`cluster_var' `cluster_var2'"
        ereturn local clustvar1 "`cluster_var'"
        ereturn local clustvar2 "`cluster_var2'"
    }
    else if `vce_type' == 4 {
        ereturn local vcetype "Robust"
        ereturn local clustvar "`cluster_var'"
    }
    if `kernel_type' > 0 {
        ereturn local kernel "`kernel_name'"
        ereturn scalar bw = `bw_val'
        ereturn local tvar "`ts_tvar'"
        if "`ts_ivar'" != "" ereturn local ivar "`ts_ivar'"
    }
    * ivreg2's e(vce): robust, cluster or two-way cluster, hac or ac with
    * the kernel and bandwidth, psd0 or psda (vce(bootstrap) replaces it)
    local __vrob = (`vce_type' != 0 & !`__dk_norob')
    local __vce = cond(`__vrob', "robust", "")
    if inlist(`vce_type', 2, 4) local __vce "`__vce' cluster"
    if `vce_type' == 3 local __vce "`__vce' two-way cluster"
    if `kernel_type' > 0 {
        local __vce = "`__vce' " + cond(`__vrob', "hac", "ac") + " `kernel_name' bw=`bw_val'"
    }
    if `psd_code' local __vce "`__vce' `=cond(`psd_code' == 1, "psd0", "psda")'"
    local __vce : list clean __vce
    local __vce = lower("`__vce'")
    if "`__vce'" != "" ereturn local vce "`__vce'"

    * Store additional statistics
    ereturn scalar rss = `rss'
    ereturn scalar tss = `tss'
    ereturn scalar mss = `tss' - `rss'
    ereturn scalar r2 = `r2'
    ereturn scalar rmse = `rmse'
    if "`small'`absorb'" == "" {
        if `vce_type' < 2 local F_wald = `F_wald' * `df_r' / (`N_used' - `dofminus')
        else local F_wald = `F_wald' * `df_r'/(`df_r'+1) * ///
            (`N_used'-`K'-`sdofminus'-cond("`nopartialsmall'"=="",`n_partial',0))/(`N_used'-1)
    }
    if `__bs' {
        ereturn scalar chi2 = `bs_chi2'
        ereturn scalar p = chi2tail(`model_df', `bs_chi2')
        ereturn local chi2type Wald
        ereturn hidden scalar _civ_wdf = `model_df'
    }
    else {
        ereturn scalar F = `F_wald'
        ereturn scalar Fdf1 = `model_df'
        ereturn scalar Fdf2 = `df_r'
        * ivreg2's p-value of F
        ereturn scalar Fp = Ftail(`model_df', `df_r', `F_wald')
    }
    * ivreg2's centered R2 and variable counts and lists
    * ivreg2's total sums of squares. After absorb() or partial() the
    * constant is partialled out and both are the plugin's. Otherwise the
    * plugin's is the centered one, or with noconstant the uncentered one;
    * the other differs by sum(w) * mean^2 (aweights and pweights are
    * normalized to sum to N, as in ivreg2).
    local __yy = `tss'
    local __yyc = `tss'
    if `G' == 0 & `n_partial' == 0 {
        tempvar __yv
        quietly gen double `__yv' = `depvar' if e(sample)
        local __sw ""
        if inlist("`weight'", "aweight", "pweight") local __sw "[aw=`weight_var']"
        if "`weight'" == "fweight" local __sw "[fw=`weight_var']"
        if "`weight'" == "iweight" local __sw "[iw=`weight_var']"
        quietly summarize `__yv' `__sw' if e(sample)
        local __shift = cond(inlist("`weight'", "fweight", "iweight"), r(sum_w), r(N)) * r(mean)^2
        if `has_cons' local __yy = `tss' + `__shift'
        else local __yyc = `tss' - `__shift'
    }
    ereturn scalar yy = `__yy'
    ereturn scalar yyc = `__yyc'
    ereturn scalar r2u = 1 - `rss' / `__yy'
    ereturn scalar r2c = 1 - `rss' / `__yyc'
    ereturn scalar endog_ct = `K_endog'
    ereturn scalar exexog_ct = `K_inst_excl'
    * (ivreg2's count of the included exogenous regressors not partialled
    * out; its lists with the base levels of factor variables)
    ereturn scalar inexog_ct = `K_exog' - (`has_cons' & !`partial_cons')
    ereturn local instd `endogvars_names_all'
    ereturn local inexog `: list exogvars_names_all - __pc'
    ereturn local exexog `instruments_names_all'
    * (excluded, then included, as ivreg2 lists them)
    ereturn local insts `instruments_names_all' `: list exogvars_names_all - __pc'
    * ivreg2's "0" varlists (base levels and omitted terms included,
    * partialled-out terms not) and "1" varlists (without them)
    local __ne0 : word count `endogvars_names_all'
    local __instd0
    local __inexog0
    local __j = 0
    foreach __v of local __b0names {
        local ++__j
        if `__j' <= `__ne0' local __instd0 `__instd0' `__v'
        else if "`__v'" != "_cons" local __inexog0 `__inexog0' `__v'
    }
    _civreghdfe_list1 `__instd0'
    local __instd1 `r(list)'
    _civreghdfe_list1 `__inexog0'
    local __inexog1 `r(list)'
    local __exexog1
    local __j = 0
    foreach __v of local instruments_coef {
        local ++__j
        local __zd = 0
        capture local __zd = (scalar(__civreghdfe_zcollin_`__j') == 1)
        if !`__zd' local __exexog1 `__exexog1' `__v'
    }
    _civreghdfe_list1 `__exexog1'
    local __exexog1 `r(list)'
    _civreghdfe_list1 `: list partial_terms - __pc'
    local __partial1 `r(list)'
    ereturn local depvar0 `depvar'
    ereturn local depvar1 `depvar'
    ereturn local instd0 `__instd0'
    ereturn local instd1 `__instd1'
    ereturn local inexog0 `__inexog0'
    ereturn local inexog1 `__inexog1'
    ereturn local exexog0 `instruments_names_all'
    ereturn local exexog1 `__exexog1'
    ereturn local insts0 `instruments_names_all' `__inexog0'
    ereturn local insts1 `__exexog1' `__inexog1'
    ereturn local partial1 `__partial1'

    * ivreg2's model degrees of freedom, rankxx - cons + sdofminus -
    * partialcons: the slopes plus sdofminus() and the partialled-out
    * regressors in the small-sample count (always with absorb()), less a
    * partialled-out constant (the absorbed effects with absorb()). Also
    * ivreg2's sdofminus, partialcons and cons.
    local __partialcons = cond(`G' > 0, max(1, `df_a_for_vce'), `partial_cons')
    local __iv_sdof = `sdofminus' + cond(`G' > 0 | "`nopartialsmall'" == "", `n_partial', 0) + ///
        cond(`G' > 0, `__partialcons', 0)
    ereturn scalar df_m = `model_df' + `__iv_sdof' - `__partialcons'
    ereturn scalar sdofminus = `__iv_sdof'
    ereturn scalar partialcons = `__partialcons'
    ereturn scalar cons = (`G' == 0 & "`noconstant'" == "" & !`partial_cons')

    * ivreg2's ranks: regressors and instruments kept, rank of S, nonzero
    * variances; its log likelihood; reghdfe's singleton counts
    ereturn scalar rankxx = `K'
    capture ereturn scalar rankzz = scalar(__civreghdfe_rankzz)
    capture ereturn scalar rankS = scalar(__civreghdfe_rankS)
    tempname __Vd
    matrix `__Vd' = vecdiag(e(V))
    local __rankV = 0
    forvalues __j = 1/`=colsof(`__Vd')' {
        if `__Vd'[1, `__j'] != 0 local ++__rankV
    }
    ereturn scalar rankV = `__rankV'
    ereturn scalar ll = -0.5 * (`N_used' * ln(2 * _pi) + `N_used' * ln(`rss' / `N_used') + `N_used')
    capture ereturn scalar num_singletons = scalar(__civreghdfe_num_singletons)
    ereturn scalar N_full = `N'

    * Number of absorbed fixed effect dimensions (matches reghdfe/ivreghdfe):
    * the effects, and their intercepts and slopes (e(extended_absvars))
    ereturn scalar N_hdfe = cond(`G' > 0, `n_absorb', 0)
    if `G' > 0 {
        ereturn scalar N_hdfe_extended = `G'
        ereturn local extended_absvars "`ext_names'"
    }

    * DOF adjustment options
    * (ivreg2 always posts e(dofminus))
    ereturn scalar dofminus = `dofminus'
    if `sdofminus' != 0 {
        ereturn scalar sdofminus_opt = `sdofminus'
    }

    * Adjusted R-squared (matching ivreghdfe formula)
    * With absorb: noconstant is set, so constant not added to sdofminus
    * sdofminus = absorb_ct = max(1, HDFE.df_a)
    *   - HDFE.df_a = df_a_for_vce (0 if nested, df_a otherwise)
    *   - ivreghdfe line 670: if HDFE.df_a=0 (nested), force absorb_ct to 1
    * Without absorb: constant is already in K_total, so sdofminus = 0
    * Formula: r2_a = 1 - (1 - r2) * N / (N - K - sdofminus)  (ivreghdfe line 1304)
    *   - absorb also implies partialsmall: sdofminus() and the partial()
    *     regressors add to sdofminus, and dofminus() is subtracted too
    if `G' > 0 {
        local sdofminus = max(1, `df_a_for_vce') + `sdofminus' + `n_partial' + `dofminus'
    }
    else {
        local sdofminus = `sdofminus' + `dofminus' + cond("`nopartialsmall'" == "", `n_partial', 0)
    }
    * Use K (plugin's non-collinear count) not K_total, so vars absorbed by FE
    * don't inflate the denominator (matches df_r = N - K - df_a from plugin)
    * (with iweights the sum of the weights, as ivreg2 uses it here)
    local __Nr2 = cond("`weight'" == "iweight", scalar(__civreghdfe_N), `N_used')
    local adj_denom = `__Nr2' - `K' - `sdofminus'
    if `adj_denom' > 0 {
        * With absorb, constant is partialled out, so use N not N-1
        * (N when the constant is partialled out with partial(), as ivreg2)
        local r2_a = 1 - (1 - `r2') * (`__Nr2' - (`G' == 0 & "`noconstant'" == "" & !`partial_cons')) / `adj_denom'
    }
    else {
        local r2_a = .
    }
    ereturn scalar r2_a = `r2_a'

    * Auxiliary regressions (first stage, reduced form): civreghdfe fits
    * without endogenous regressors over the IV estimation sample, as ivreg2
    * runs them. They share the main fit's weights, fixed effects, VCE
    * (kernels, Driscoll-Kraay and Kiefer included), partial() and
    * noconstant, and report ivreg2's small-sample statistics. The IV
    * results are held so they are restored afterwards.
    tempvar aux_sample
    if "`first'`rf'" != "" quietly gen byte `aux_sample' = e(sample)
    local aux_wt ""
    if "`weight'" != "" local aux_wt "[`weight'=`weight_var']"
    local aux_opts small tolerance(`tolerance') maxiter(`maxiter') civaux
    if "`absorb'" != "" local aux_opts `aux_opts' absorb(`absorb_names')
    if `"`vce'"' != "" local aux_opts `aux_opts' vce(`vce')
    if `"`bw'"' != "" local aux_opts `aux_opts' bw(`bw_val')
    if "`kernel'" != "" local aux_opts `aux_opts' kernel(`kernel')
    if `dkraay' != 0 local aux_opts `aux_opts' dkraay(`dkraay')
    if "`kiefer'" != "" local aux_opts `aux_opts' kiefer
    if "`noconstant'" != "" local aux_opts `aux_opts' noconstant
    if "`partial'" != "" | `__pcons_given' local aux_opts `aux_opts' partial(`partial' `=cond(`__pcons_given', "_cons", "")')
    if "`nopartialsmall'" != "" local aux_opts `aux_opts' nopartialsmall
    if `__user_dofminus' != 0 local aux_opts `aux_opts' dofminus(`__user_dofminus')
    if `__user_sdofminus' != 0 local aux_opts `aux_opts' sdofminus(`__user_sdofminus')
    if "`keepsingletons'" != "" local aux_opts `aux_opts' keepsingletons
    if "`dofadjustments'" != "" local aux_opts `aux_opts' dofadjustments(`dofadjustments')
    if `threads' > 0 local aux_opts `aux_opts' threads(`threads')
    if `psd_code' local aux_opts `aux_opts' `=cond(`psd_code' == 1, "psd0", "psda")'

    * Saved fixed effects: d = y - Xb (the actual regressors, the constant
    * kept) absorbed by creghdfe with savefe on e(sample). Its residual is
    * the IV residual on the partialled data, so its effects are those of
    * d = y - _cons - Xb - e, normalized as creghdfe and reghdfe normalize
    * them, with the constant in the first (docs/development/BOOTSTRAP_PLAN.md,
    * 9.1). With partial(), the partialled regressors W, whose coefficients
    * e(b) leaves out, are regressors there: by FWL their coefficients are
    * those of y - Xb on W with the effects absorbed, and the residual is
    * still the IV residual.
    if `fe_any' {
        local __fe_w : list partial_terms - __pc
        tempvar __fe_y __fe_xb __fe_smp
        tempname __fe_hold __fe_b
        quietly gen byte `__fe_smp' = e(sample)
        matrix `__fe_b' = e(b)
        quietly matrix score double `__fe_xb' = `__fe_b' if `__fe_smp'
        quietly gen double `__fe_y' = `depvar_var' - `__fe_xb' if `__fe_smp'
        if colnumb(e(b), "_cons") < . quietly replace `__fe_y' = `__fe_y' + _b[_cons] if `__fe_smp'
        local __fe_wgt
        if "`weight'" != "" local __fe_wgt [`weight'=`weight_expr']
        _estimates hold `__fe_hold', restore
        capture noisily quietly creghdfe `__fe_y' `__fe_w' `__fe_wgt' if `__fe_smp', ///
            absorb(`fe_terms', savefe) keepsingletons
        local __fe_rc = _rc
        local __fe_cons 0
        if !`__fe_rc' local __fe_cons = _b[_cons]
        _estimates unhold `__fe_hold'
        if `__fe_rc' {
            forvalues __g = 1/`n_absorb' {
                capture drop __hdfe`__g'__
                capture drop __hdfe`__g'__Slope*
            }
            di as error "civreghdfe: recovering the fixed effects failed"
            * no partial results: e() is already posted, with e(V) possibly the
            * bootstrap's but e(vce) not yet
            ereturn clear
            _ctools_cleanup civreghdfe
            exit `__fe_rc'
        }
        * As ivreghdfe: without _cons in e(b) (it is partialled out with the
        * effects), the constant goes into the first effect, so that the
        * effects sum to y - Xb - e
        * (the first effect with an intercept: a slope has no constant)
        if colnumb(e(b), "_cons") == . {
            forvalues __g = 1/`n_absorb' {
                if `absorb_int_`__g'' {
                    quietly replace __hdfe`__g'__ = __hdfe`__g'__ + `__fe_cons'
                    continue, break
                }
            }
        }
        local __fe_saved
        forvalues __g = 1/`n_absorb' {
            local __lab : word `__g' of `absorb_names'
            local __tg `fe_target_`__g''
            if `absorb_int_`__g'' {
                if "`__tg'" == "" drop __hdfe`__g'__
                else {
                    if "`__tg'" != "__hdfe`__g'__" rename __hdfe`__g'__ `__tg'
                    label variable `__tg' "FE for `__lab'"
                    local __fe_saved `__fe_saved' `__tg'
                }
            }
            local __h = 0
            foreach __cv of local absorb_cvars_`__g' {
                local ++__h
                if "`__tg'" == "" drop __hdfe`__g'__Slope`__h'
                else {
                    if "`__tg'" != "__hdfe`__g'__" rename __hdfe`__g'__Slope`__h' `__tg'Slope`__h'
                    label variable `__tg'Slope`__h' "FE for `__lab'"
                    local __fe_saved `__fe_saved' `__tg'Slope`__h'
                }
            }
        }
        ereturn local savefe "`__fe_saved'"
    }
    * groupvar(): reghdfe's mobility group of the first two absorbed
    * effects on e(sample) (none with one effect), from creghdfe
    if "`groupvar'" != "" & `n_absorb' >= 2 {
        tempname __gv_hold
        tempvar __gv_smp
        quietly gen byte `__gv_smp' = e(sample)
        _estimates hold `__gv_hold', restore
        capture noisily quietly creghdfe `depvar_var' if `__gv_smp', absorb(`fe_terms') ///
            groupvar(`groupvar') keepsingletons
        local __gv_rc = _rc
        _estimates unhold `__gv_hold'
        if `__gv_rc' {
            capture drop `groupvar'
            di as error "civreghdfe: groupvar() failed"
            ereturn clear
            _ctools_cleanup civreghdfe
            exit `__gv_rc'
        }
    }
    if `__fe' {
        local __j = 0
        forvalues __g = 1/`: word count `fe_terms'' {
            if "`fe_target_`__g''" == "" continue
            local ++__j
            local __n : word `__j' of `fe_se_names'
            rename `: word `__j' of `fe_se_tmp'' `__n'
            label variable `__n' "Bootstrap SE of FE for `: word `__g' of `absorb_names''"
            if `fe_nreps' {
                local __n : word `__j' of `fe_nrep_names'
                rename `: word `__j' of `fe_nrep_tmp'' `__n'
                label variable `__n' "Bootstrap replicates estimating FE for `: word `__g' of `absorb_names''"
            }
            if `fe_ci' {
                local __n : word `__j' of `fe_lb_names'
                rename `: word `__j' of `fe_lb_tmp'' `__n'
                label variable `__n' "Bootstrap lower bound of FE for `: word `__g' of `absorb_names''"
                local __n : word `__j' of `fe_ub_names'
                rename `: word `__j' of `fe_ub_tmp'' `__n'
                label variable `__n' "Bootstrap upper bound of FE for `: word `__g' of `absorb_names''"
            }
        }
        if "`fe_dse'" != "" {
            rename `fe_dse_tmp' `fe_dse'
            label variable `fe_dse' "Bootstrap SE of the sum of fixed effects"
        }
        local __fe_abs
        forvalues __g = 1/`: word count `fe_terms'' {
            if "`fe_target_`__g''" != "" local __fe_abs `__fe_abs' `: word `__g' of `absorb_names''
        }
        _civreghdfe_fe_post, sevars(`fe_se_names') nrepvars(`fe_nrep_names') dse(`fe_dse') ///
            lbvars(`fe_lb_names') ubvars(`fe_ub_names') level(`bs_level') ///
            absvars(`__fe_abs') nabsorb(`G') rngstate(`"`bs_rngstate'"') draws(`bs_draws') ///
            cluster(`bs_cluster_name') strata(`bs_strata')
    }
    if `__bs' {
        _civreghdfe_bootstrap_post, vmb(`V_mb') bsb(`bs_b') bsse(`bs_se') bsbias(`bs_bias') ///
            bsz0(`bs_z0') cinormal(`bs_ci_normal') cipercentile(`bs_ci_percentile') ///
            cibc(`bs_ci_bc') level(`bs_level') rngstate(`"`bs_rngstate'"') ///
            cluster(`bs_cluster_name') strata(`bs_strata') draws(`bs_draws') ///
            command(`"`cmdline'"') mse(`bs_mse')
    }
    if `__wild' {
        _civreghdfe_wild_post, test(`wb_test') names(`wb_names') ///
            method(`=cond(`wb_null', "wre", "wue")') type(`wb_type') ptype(`wb_ptype') ///
            reps(`wb_B') enum(`wb_enum') nclust(`wb_C') level(`wb_level') path(`wb_path') ///
            ci(`wb_ci') bootclust(`wb_bootname') rngstate(`"`wb_rngstate'"')
    }
    if `"`bs_save_file'"' != "" {
        local __nclust
        local __nstrata
        if "`bs_unitvar'" != "" local __nclust `bs_nunits'
        if "`bs_strata'" != "" local __nstrata `bs_nstrata'
        capture noisily _ctools_hdfe_boot savefinish, tmp(`"`bs_savetmp'"') file(`"`bs_save_file'"') ///
            map(`bs_map') `bs_save_double' `bs_save_replace' command(`"`cmdline'"') ///
            seed(`"`bs_rngstate'"') n(`N_used') cluster(`bs_cluster_name') nclust(`__nclust') ///
            strata(`bs_strata') nstrata(`__nstrata') cmdname(civreghdfe)
        if _rc {
            local __src = _rc
            ereturn clear
            * nor the outputs of the failed run
            foreach __v in `__fe_saved' `fe_se_names' `fe_nrep_names' `fe_lb_names' `fe_ub_names' `fe_dse' {
                capture drop `__v'
            }
            _ctools_cleanup civreghdfe
            exit `__src'
        }
    }

    * Display options of the coefficient tables
    local dispopt ""
    if "`noomitted'" != "" local dispopt "`dispopt' noomitted"
    if "`omitted'" != "" local dispopt "`dispopt' omitted"
    if "`vsquish'" != "" local dispopt "`dispopt' vsquish"
    if "`noemptycells'" != "" local dispopt "`dispopt' noemptycells"
    if "`baselevels'" != "" local dispopt "`dispopt' baselevels"
    if "`allbaselevels'" != "" local dispopt "`dispopt' allbaselevels"
    if "`plus'" != "" local dispopt "`dispopt' plus"
    if "`eform'" != "" local dispopt "`dispopt' eform(`eform')"
    if "`coeflegend'" != "" local dispopt "`dispopt' coeflegend"
    if `"`cformat'"' != "" local dispopt `"`dispopt' cformat(`cformat')"'

    * What the statistics are robust to (ivreg2's wording), for the
    * auxiliary and main tables
    local __hac1 = cond(`vce_type' != 0, "heteroskedasticity", "")
    local __hac3 ""
    if `kernel_type' > 0 & inlist(`vce_type', 0, 1) local __hac3 "autocorrelation"
    if `kiefer_val' local __hac3 "within-cluster autocorrelation (Kiefer)"
    if inlist(`vce_type', 2, 4) local __hac3 "clustering on `cluster_var'"
    if `vce_type' == 3 local __hac3 "clustering on `cluster_var' and `cluster_var2'"
    local __hacsep = cond("`__hac1'" != "" & "`__hac3'" != "", " and ", "")
    local __hac "`__hac1'`__hacsep'`__hac3'"
    * ivreg2's stored first-stage equations do not mention Kiefer
    local __hacaux "`__hac'"
    if `kiefer_val' local __hacaux "autocorrelation"
    local __hacV2 ""
    if `kernel_type' > 0 & inlist(`vce_type', 3, 4) {
        local __hacV2 "and kernel-robust to common correlated disturbances (Driscoll-Kraay)"
    }

    * wmatrix(): ivreg2's GMM with a user-supplied weighting matrix (IV and
    * two-step GMM only; LIML, k-class and CUE do not use it)
    local __gmmw = (`"`wmatrix'"' != "" & `est_method' == 0 & "`b0'" == "" & `K_endog' > 0)
    local __gmm2sw = (`"`wmatrix'"' != "" & `est_method' == 4)
    if `__gmmw' ereturn local model "gmmw"

    * ivreg2's title for the estimator
    local __title `"`title'"'
    if `"`__title'"' == "" {
        if "`b0'" != "" local __title "CUE evaluated at user-supplied parameter vector"
        else if `__gmm2sw' local __title "2-Step GMM estimation with user-supplied first-step weighting matrix"
        else if `__gmmw' local __title "GMM estimation with user-supplied weighting matrix"
        else if `K_endog' == 0 local __title "OLS estimation"
        else if `liml_kclass' | inlist(`est_method', 1, 2) local __title "LIML estimation"
        else if `est_method' == 3 local __title "k-class estimation"
        else if `est_method' == 4 local __title "2-Step GMM estimation"
        else if `est_method' == 5 local __title "CUE estimation"
        else local __title "IV (2SLS) estimation"
    }
    if `"`subtitle'"' != "" ereturn local subtitle `"`subtitle'"'
    * what the estimates are efficient for (ivreg2's hacsubtitleB)
    if "`b0'" != "" local __hacB "Estimates based on supplied parameter vector"
    else if "`__hac'" != "" & inlist(`est_method', 4, 5) local __hacB "Estimates efficient for arbitrary `__hac'"
    else if `"`wmatrix'"' != "" local __hacB "Efficiency of estimates dependent on weighting matrix"
    else local __hacB "Estimates efficient for homoskedasticity only"
    * ivreg2's e(title) and the subtitles of its display
    ereturn local title `"`__title'"'
    ereturn local hacsubtitleB "`__hacB'"
    if "`__hac'" != "" ereturn local hacsubtitleV "Statistics robust to `__hac'"
    else ereturn local hacsubtitleV "Statistics consistent for homoskedasticity only"
    if "`__hacV2'" != "" ereturn local hacsubtitleV2 "`__hacV2'"

    * Stock-Yogo critical values, by ivreg2's rules: IV relative bias and
    * size for 2SLS and GMM, LIML size for LIML and CUE, Fuller bias for
    * Fuller's LIML, none for other k-class estimators
    local __sy_model "iv"
    if "`b0'" != "" | `liml_kclass' | inlist(`est_method', 1, 5) local __sy_model "liml"
    else if `est_method' == 2 local __sy_model "fuller"
    else if `est_method' == 3 local __sy_model "kclass"

    * Diagnostic statistics, stored in e() (displayed by _civreghdfe_disp)
    * Underidentification test: without endogenous regressors there are no
    * identification statistics, as in ivreg2, nor with noid
    capture scalar `temp_underid' = __civreghdfe_underid
    capture scalar `temp_underid_df' = __civreghdfe_underid_df
    if _rc == 0 & `K_endog' > 0 & "`id'" == "" {
        ereturn scalar idstat = `temp_underid'
        ereturn scalar iddf = `temp_underid_df'
        ereturn scalar idp = chi2tail(`temp_underid_df', `temp_underid')
    }

    capture scalar `temp_sargan' = __civreghdfe_sargan
    capture scalar `temp_sargan_df' = __civreghdfe_sargan_df
    capture scalar `temp_cd_f' = __civreghdfe_cd_f
    capture scalar `temp_endog_chi2' = __civreghdfe_endog_chi2
    capture scalar `temp_endog_f' = __civreghdfe_endog_f
    capture scalar `temp_endog_df' = __civreghdfe_endog_df
    if _rc == 0 {
        local sargan_val = `temp_sargan'
        local sargan_df = `temp_sargan_df'
        local cd_f_val = `temp_cd_f'
        local endog_chi2 = `temp_endog_chi2'
        local endog_f = `temp_endog_f'
        local endog_df = `temp_endog_df'

        ereturn scalar sargan = `sargan_val'
        ereturn scalar sargan_df = `sargan_df'
        if `sargan_df' > 0 {
            local sargan_pval = chi2tail(`sargan_df', `sargan_val')
            ereturn scalar sargan_p = `sargan_pval'
            ereturn scalar sarganp = `sargan_pval'
        }
        if `K_endog' > 0 ereturn scalar cd_f = `cd_f_val'
        * ivreg2's names for the overidentification and Cragg-Donald
        * statistics
        ereturn scalar j = `sargan_val'
        ereturn scalar jdf = `sargan_df'
        if `sargan_df' > 0 ereturn scalar jp = chi2tail(`sargan_df', `sargan_val')
        if `K_endog' > 0 & "`id'" == "" ereturn scalar cdf = `cd_f_val'

        * widstat = Kleibergen-Paap F if robust, else Cragg-Donald F
        capture scalar `temp_kp_f2' = __civreghdfe_kp_f
        if `K_endog' == 0 | "`id'" != "" {
        }
        else if _rc == 0 & (`vce_type' != 0 | `kernel_type' > 0) {
            local kp_f_val2 = `temp_kp_f2'
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

        ereturn scalar endog_chi2 = `endog_chi2'
        ereturn scalar endog_f = `endog_f'
        ereturn scalar endog_df = `endog_df'
        if `endog_df' > 0 {
            ereturn scalar endog_p = chi2tail(`endog_df', `endog_chi2')
        }
    }

    * C statistic of orthog() (excluded and included instruments)
    if `n_orthog' + `n_orthogx' > 0 {
        capture scalar `temp_cstat' = __civreghdfe_cstat
        capture scalar `temp_cstat_df' = __civreghdfe_cstat_df
        if _rc == 0 {
            ereturn scalar cstat = `temp_cstat'
            ereturn scalar cstat_df = `temp_cstat_df'
            ereturn scalar cstat_p = chi2tail(`temp_cstat_df', `temp_cstat')
        }
    }
    if `n_endogtest' > 0 {
        capture scalar `temp_endogtest_stat' = __civreghdfe_endogtest_stat
        capture scalar `temp_endogtest_df' = __civreghdfe_endogtest_df
        if _rc == 0 {
            ereturn scalar endogtest = `temp_endogtest_stat'
            ereturn scalar endogtest_df = `temp_endogtest_df'
            ereturn scalar endogtest_p = chi2tail(`temp_endogtest_df', `temp_endogtest_stat')
        }
    }
    if `n_redundant' > 0 {
        capture scalar `temp_redund_stat' = __civreghdfe_redund_stat
        capture scalar `temp_redund_df' = __civreghdfe_redund_df
        if _rc == 0 {
            ereturn scalar redund = `temp_redund_stat'
            ereturn scalar redund_df = `temp_redund_df'
            ereturn scalar redund_p = chi2tail(`temp_redund_df', `temp_redund_stat')
        }
    }

    * ivreg2's names and results for the diagnostic tests: the Sargan dof
    * without robust or cluster VCEs; the C statistic (cstatp, cstatdf and
    * clist when it is positive), the endogeneity test (estat...), the
    * redundancy test (redstat...); the Cragg-Donald eigenvalue e(cd), the
    * squared canonical correlations e(ccev) and Cragg-Donald eigenvalues
    * e(cdev), the KP rk Wald F e(rkf); LIML's Anderson-Rubin
    * overidentification statistic (arubin...)
    local __robust = (`vce_type' != 0)
    if !`__robust' & !missing(e(jdf)) ereturn scalar sargandf = e(jdf)
    if "`orthog'" != "" & !missing(e(cstat)) {
        if e(cstat) != 0 {
            ereturn scalar cstatp = e(cstat_p)
            ereturn scalar cstatdf = e(cstat_df)
            _civreghdfe_list1 `__orthog_terms'
            ereturn local clist `r(list)'
        }
    }
    if "`endogtest'" != "" & !missing(e(endogtest_df)) {
        ereturn scalar estat = e(endogtest)
        if e(endogtest) != 0 & !missing(e(endogtest)) {
            ereturn scalar estatp = e(endogtest_p)
            ereturn scalar estatdf = e(endogtest_df)
            _civreghdfe_list1 `__endogtest_terms'
            ereturn local elist `r(list)'
        }
    }
    if "`redundant'" != "" & "`id'" == "" & !missing(e(redund_df)) {
        ereturn scalar redstat = e(redund)
        ereturn scalar redp = e(redund_p)
        ereturn scalar reddf = e(redund_df)
        _civreghdfe_list1 `__redundant_terms'
        ereturn local redlist `r(list)'
    }
    if `K_endog' > 0 & "`id'" == "" {
        tempname __ccev __cdev
        local __kc = 0
        forvalues __j = 1/`=colsof(__civreghdfe_ccev)' {
            if !missing(__civreghdfe_ccev[1, `__j']) local __kc = `__j'
        }
        if `__kc' > 0 {
            matrix `__ccev' = __civreghdfe_ccev[1, 1..`__kc']
            matrix `__cdev' = `__ccev'
            forvalues __j = 1/`__kc' {
                matrix `__cdev'[1, `__j'] = `__ccev'[1, `__j'] / (1 - `__ccev'[1, `__j'])
            }
            ereturn scalar cd = `__cdev'[1, `__kc']
            ereturn matrix ccev = `__ccev'
            ereturn matrix cdev = `__cdev'
        }
        if (`__robust' | `kernel_type' > 0) & !missing(e(kp_f)) ereturn scalar rkf = e(kp_f)
    }
    if `liml_kclass' | inlist(`est_method', 1, 2) {
        local __arn = e(N) - `__user_dofminus'
        ereturn scalar arubin = `__arn' * ln(e(lambda))
        ereturn scalar arubin_lin = `__arn' * (e(lambda) - 1)
        if e(j) != 0 & !missing(e(j)) {
            ereturn scalar arubinp = chi2tail(e(jdf), e(arubin))
            ereturn scalar arubin_linp = chi2tail(e(jdf), e(arubin_lin))
        }
        ereturn scalar arubindf = e(jdf)
    }
    * ivreg2's partialsmall, condition numbers, S and W (the instruments
    * kept, in ivreg2's order: excluded, then included, then _cons)
    if "`nopartialsmall'" == "" ereturn local partialsmall "small"
    capture ereturn scalar condzz = scalar(__civreghdfe_condzz)
    capture ereturn scalar condxx = scalar(__civreghdfe_condxx)
    _civreghdfe_zorder, znames(`__znames') zxct(`__zx_ct')
    local __zp `r(perm)'
    local __zl `r(names)'
    local __zk : word count `__zp'
    if `__zk' > 0 & `__zk' <= rowsof(__civreghdfe_S) & "`civaux'" == "" {
        tempname __Sm __Wm
        mata: st_matrix(st_local("__Sm"), st_matrix("__civreghdfe_S")[strtoreal(tokens(st_local("__zp"))), strtoreal(tokens(st_local("__zp")))])
        if !missing(`__Sm'[1, 1]) {
            mata: st_replacematrix(st_local("__Sm"), makesymmetric(st_matrix(st_local("__Sm"))))
            matrix colnames `__Sm' = `__zl'
            matrix rownames `__Sm' = `__zl'
            ereturn matrix S = `__Sm'
        }
        mata: st_matrix(st_local("__Wm"), st_matrix("__civreghdfe_W")[strtoreal(tokens(st_local("__zp"))), strtoreal(tokens(st_local("__zp")))])
        if !missing(`__Wm'[1, 1]) {
            mata: st_replacematrix(st_local("__Wm"), makesymmetric(st_matrix(st_local("__Wm"))))
            matrix colnames `__Wm' = `__zl'
            matrix rownames `__Wm' = `__zl'
            ereturn matrix W = `__Wm'
        }
    }

    * Partial R2 of each first stage (ivreg2's pr2), e(partial_r2_#) with
    * ffirst and first
    if ("`ffirst'`first'" != "") & `K_endog' > 0 {
        forvalues e = 1/`__n_fs' {
            if `__fs_index_`e'' == 0 {
                ereturn scalar partial_r2_`e' = .
                continue
            }
            capture scalar `temp_r2' = __civreghdfe_partial_r2_`__fs_index_`e''
            if _rc == 0 ereturn scalar partial_r2_`e' = `temp_r2'
        }
    }

    * ivreg2's first-stage statistics: e(first) (one column per endogenous
    * regressor kept), Anderson-Rubin and Stock-Wright tests. The plugin
    * leaves the p-values to Stata's Ftail() and chi2tail().
    local __ff_ok = 0
    local __sys_ok = 0
    if `__want_sys' {
        capture confirm scalar __civreghdfe_ff_K1
        if _rc == 0 local __sys_ok = 1
        local __ff_ok = `__sys_ok' & `__want_ff'
    }
    * the endogenous regressors kept (in e(b)'s order)
    if `__sys_ok' {
        local __k1 = scalar(__civreghdfe_ff_K1)
        local __ffL = scalar(__civreghdfe_ff_L)
        local __endo1
        local __j = 0
        foreach __v of local endogvars_coef {
            local ++__j
            if !`: word `__j' of `collinear_flags'' local __endo1 `__endo1' `__v'
        }
    }
    if `__ff_ok' {
        tempname __fm
        matrix `__fm' = __civreghdfe_firstmat[1..21, 1..`__k1']
        forvalues __c = 1/`__k1' {
            matrix `__fm'[7, `__c'] = Ftail(`__fm'[5, `__c'], `__fm'[6, `__c'], `__fm'[4, `__c'])
            matrix `__fm'[11, `__c'] = Ftail(`__fm'[9, `__c'], `__fm'[10, `__c'], `__fm'[8, `__c'])
            matrix `__fm'[13, `__c'] = chi2tail(`__fm'[9, `__c'], `__fm'[12, `__c'])
            matrix `__fm'[18, `__c'] = Ftail(`__fm'[16, `__c'], `__fm'[17, `__c'], `__fm'[15, `__c'])
            matrix `__fm'[20, `__c'] = chi2tail(`__fm'[16, `__c'], `__fm'[19, `__c'])
        }
        matrix rownames `__fm' = rmse sheapr2 pr2 F df df_r pvalue SWF SWFdf1 SWFdf2 SWFp ///
            SWchi2 SWchi2p SWr2 APF APFdf1 APFdf2 APFp APchi2 APchi2p APr2
        matrix colnames `__fm' = `__endo1'
        ereturn matrix first = `__fm'
        ereturn scalar arf = scalar(__civreghdfe_arf)
        ereturn scalar ardf = scalar(__civreghdfe_ardf)
        ereturn scalar ardf_r = scalar(__civreghdfe_ardf_r)
        ereturn scalar arfp = Ftail(e(ardf), e(ardf_r), e(arf))
        ereturn scalar archi2 = scalar(__civreghdfe_archi2)
        ereturn scalar archi2p = chi2tail(e(ardf), e(archi2))
        ereturn scalar sstat = scalar(__civreghdfe_sstat)
        ereturn scalar sstatdf = scalar(__civreghdfe_sstatdf)
        ereturn scalar sstatp = chi2tail(e(sstatdf), e(sstat))
    }

    * ivreg2's center flag and automatic bandwidth choice
    ereturn scalar center = ("`center'" != "")
    if "`bwchoice'" != "" & `kernel_type' > 0 ereturn local bwchoice "`bwchoice'"

    * Excluded instruments dropped after the fixed effects were partialled
    * out (an included exogenous regressor is an instrument too)
    local __inst_dropped = 0
    capture local __inst_dropped = scalar(__civreghdfe_z_dropped)
    if _rc local __inst_dropped = 0
    forvalues __k = 1/`K_exog' {
        local __kk = `K_endog' + `__k'
        local __inst_dropped = `__inst_dropped' + (`: word `__kk' of `collinear_flags'' == 1)
    }

    * The footer's variable lists, without the terms ivreg2's collinearity
    * rules dropped (listed apart): partialled-out regressors are listed
    * apart from the included instruments, followed by _cons when the
    * constant is partialled out (always with absorb(), as ivreghdfe shows it)
    foreach __l in endogvars exogvars instruments {
        local __disp_`__l'
        foreach __t of local `__l' {
            if !regexm("`__t'", "o\.") local __disp_`__l' `__disp_`__l'' `__t'
        }
    }
    local __disp_partial
    foreach __t of local partial {
        if "`__t'" != "_cons" local __disp_partial `__disp_partial' `__t'
    }
    local __disp_exogvars : list __disp_exogvars - __disp_partial
    local __disp_exogvars : list __disp_exogvars - partial_terms
    if `G' > 0 | `partial_cons' local __disp_partial `__disp_partial' _cons

    * Absorbed degrees of freedom table: categories, redundant, nested
    if `G' > 0 {
        tempname __doft
        matrix `__doft' = J(`G', 4, 0)
        forval g = 1/`G' {
            capture local __red = scalar(__civreghdfe_fe_redundant_`g')
            if _rc local __red = 0
            matrix `__doft'[`g', 1] = `fe_levels_`g''
            matrix `__doft'[`g', 2] = cond(`fe_nested_`g'' == 1, `fe_levels_`g'', `__red')
            matrix `__doft'[`g', 3] = (`fe_nested_`g'' == 1)
            * reghdfe's "?": the redundant count may be higher (a slope, or an
            * intercept past the ones whose count is exact)
            capture local __exact = scalar(__civreghdfe_fe_exact_`g')
            if _rc local __exact = 1
            matrix `__doft'[`g', 4] = !`__exact' & `fe_nested_`g'' != 1
        }
        * reghdfe's footnote results, which ivreghdfe posts: e(dof_table)
        * (rows 1.id, 1.id#1.t, 1.id#c.x), the absorbed dof sums, the dof
        * adjustments, the singleton rule and the title
        tempname __dt
        matrix `__dt' = J(`G', 5, 0)
        local __rows
        local __nest = 0
        forvalues g = 1/`G' {
            matrix `__dt'[`g', 1] = `__doft'[`g', 1]
            matrix `__dt'[`g', 2] = `__doft'[`g', 2]
            matrix `__dt'[`g', 3] = `__doft'[`g', 1] - `__doft'[`g', 2]
            matrix `__dt'[`g', 4] = `__doft'[`g', 4]
            matrix `__dt'[`g', 5] = `__doft'[`g', 3]
            if `__doft'[`g', 3] == 1 local __nest = `__nest' + `__doft'[`g', 2]
            local __r
            local __en : word `g' of `ext_names'
            foreach __p in `: subinstr local __en "#" " ", all' {
                local __r = "`__r'" + cond("`__r'" != "", "#", "") + cond(substr("`__p'", 1, 2) == "c.", "`__p'", "1.`__p'")
            }
            local __rows `__rows' `__r'
        }
        mata: st_matrixrowstripe(st_local("__dt"), (J(`G', 1, ""), tokens(st_local("__rows"))'))
        mata: st_matrixcolstripe(st_local("__dt"), (J(5, 1, ""), ("Categories" \ "Redundant" \ "Num Coefs" \ "Inexact?" \ "Nested?")))
        local __init = 0
        local __red = 0
        forvalues g = 1/`G' {
            local __init = `__init' + `__dt'[`g', 1]
            local __red = `__red' + `__dt'[`g', 2]
        }
        ereturn matrix dof_table = `__dt'
        ereturn scalar df_a_initial = `__init'
        ereturn scalar df_a_redundant = `__red'
        ereturn scalar df_a_nested = `__nest'
        ereturn local dofmethod "`dofmethods'"
        ereturn scalar drop_singletons = ("`keepsingletons'" == "")
        ereturn local absvars `"`fe_terms'"'
        ereturn local title2 "Absorbing `n_absorb' HDFE `=plural(`n_absorb', "group")'"
        ereturn hidden matrix _civ_doft = `__doft'
    }

    * Display information for _civreghdfe_disp (and replay)
    ereturn hidden local _civ_title `"`__title'"'
    ereturn hidden local _civ_hac "`__hac'"
    ereturn hidden local _civ_hacaux "`__hacaux'"
    ereturn hidden local _civ_hacV2 "`__hacV2'"
    ereturn hidden local _civ_hacB "`__hacB'"
    ereturn hidden local _civ_symodel "`__sy_model'"
    ereturn hidden local _civ_depvar "`depvar'"
    ereturn hidden local _civ_cl1 "`cluster_var'"
    ereturn hidden local _civ_cl2 "`cluster_var2'"
    ereturn hidden scalar _civ_vce = `vce_type'
    ereturn hidden scalar _civ_est = `est_method'
    ereturn hidden scalar _civ_lk = `liml_kclass'
    ereturn hidden scalar _civ_b0 = ("`b0'" != "")
    ereturn hidden scalar _civ_noid = ("`id'" != "")
    ereturn hidden scalar _civ_instdrop = `__inst_dropped'
    ereturn hidden scalar _civ_nopsmall = ("`nopartialsmall'" != "")
    ereturn hidden local _civ_endogtest "`endogtest'"
    ereturn hidden local _civ_redundant "`redundant'"
    ereturn hidden local _civ_f_instd "`__disp_endogvars'"
    ereturn hidden local _civ_f_incl "`__disp_exogvars'"
    ereturn hidden local _civ_f_excl "`__disp_instruments'"
    ereturn hidden local _civ_f_part "`__disp_partial'"

    * Stored equations, as ivreg2 stores them: each first-stage regression
    * as <savefprefix><endogenous variable>, the reduced form as
    * <saverfprefix><depvar> (the auxiliary civreghdfe fits of aux_opts,
    * which match ivreg2's stored equations, e(sample) included), and the
    * system of first-stage and reduced-form equations as
    * <savesfprefix>sfirst_<depvar>. first, rf and sfirst store them for the
    * display too; those not requested are dropped afterwards.
    * ivreg2's default prefix for all three is _ivreg2_
    local save_prefix = cond("`savefprefix'" != "", "`savefprefix'", "_ivreg2_")
    local save_rfprefix = cond("`saverfprefix'" != "", "`saverfprefix'", "_ivreg2_")
    local save_sfprefix = cond("`savesfprefix'" != "", "`savesfprefix'", "_ivreg2_")
    * their e(S): the equation's block of the system's S (ivreg2's), over
    * the instruments kept in ivreg2's order; their statistics subtitle
    local __eqhacV = cond("`__hacaux'" != "", "Statistics robust to `__hacaux'", ///
        "Statistics consistent for homoskedasticity only")
    local __eqS_L = 0
    if `__sys_ok' {
        _civreghdfe_zorder, znames(`__znames') zxct(`__zx_ct')
        local __eqS_names `r(names)'
        local __eqS_L : word count `__eqS_names'
        if `__eqS_L' != `__ffL' local __eqS_L = 0
    }
    local save_eqs ""
    if "`savefirst'" != "" | ("`first'" != "" & "`output'" == "") {
        foreach save_dep of local endogvars_coef {
            local save_eqs `save_eqs' first `save_dep'
        }
    }
    if "`saverf'" != "" | ("`rf'" != "" & "`output'" == "") {
        local save_eqs `save_eqs' rf `depvar'
    }
    local __firsteqs
    local __rfeq
    local __sfirsteq
    if "`save_eqs'" != "" | (`__ff_ok' & "`savesfirst'`sfirst'" != "") {
        tempvar save_sample
        quietly gen byte `save_sample' = e(sample)
        * the system of first-stage and reduced-form equations, from the
        * plugin's estimates
        if `__ff_ok' & "`savesfirst'`sfirst'" != "" {
            local __sfcl = inlist(`vce_type', 2, 3, 4)
            local __sfnc = 0
            if inlist(`vce_type', 2, 4) local __sfnc = e(N_clust)
            if `vce_type' == 3 local __sfnc = min(e(N_clust1), e(N_clust2))
            local __sfsdof = e(sdofminus)
            local __sfname = strtoname("`save_sfprefix'sfirst_`depvar'")
            _civreghdfe_sfirst, k1(`__k1') l(`__ffL') endo1(`__endo1') ///
                znames(`__znames') zxct(`__zx_ct') depvar(`depvar') sample(`save_sample') ///
                cluster(`__sfcl') nclust(`__sfnc') dofminus(`__user_dofminus') ///
                sdofminus(`__sfsdof') name(`__sfname') hacv(`"`__eqhacV'"') hacv2(`"`__hacV2'"')

        }
        tempname save_hold
        _estimates hold `save_hold', restore
        while "`save_eqs'" != "" {
            gettoken save_stage save_eqs : save_eqs
            gettoken save_dep save_eqs : save_eqs
            quietly fvrevar `save_dep'
            local save_y `r(varlist)'
            capture quietly _civreghdfe_aux `save_y' `instruments' `exogvars' `aux_wt' ///
                if `save_sample', `aux_opts'
            if _rc {
                di as error "`save_stage' regression of `save_dep' failed (rc=`=_rc'); not stored"
                continue
            }
            * ivreg2 stores the residual dof N - K - sdofminus, also with
            * clustering or Driscoll-Kraay (where the fit's own dof are the
            * number of clusters or periods less one)
            if !missing(e(N_clust)) | `dkraay' > 0 {
                ereturn scalar df_r = e(N) - colsof(e(b)) - e(partial_ct) - ///
                    cond(`G' > 0, e(partialcons), 0)
                if !missing(e(N_clust)) & "`e(clustvar)'" != "" ereturn local vce "robust cluster"
            }
            ereturn local depvar "`save_dep'"
            ereturn local model "`save_stage'"
            * ivreg2's equations count the endogenous regressors
            ereturn scalar k_eq = `K_endog'
            * as ivreg2's: the subtitles, the system's S block, no predict
            * program of the IV model (predict gives the linear prediction)
            ereturn local hacsubtitleV "`__eqhacV'"
            if "`__hacV2'" != "" ereturn local hacsubtitleV2 "`__hacV2'"
            else ereturn local hacsubtitleV2
            ereturn local predict
            if `__eqS_L' > 0 {
                local __vn = 0
                if "`save_stage'" == "first" local __vn : list posof "`save_dep'" in __endo1
                if "`save_stage'" == "rf" | `__vn' > 0 {
                    tempname __eqS
                    local __c0 = 1 + `__vn' * `__eqS_L'
                    local __c1 = (`__vn' + 1) * `__eqS_L'
                    matrix `__eqS' = __civreghdfe_sfS[`__c0'..`__c1', `__c0'..`__c1']
                    mata: st_replacematrix(st_local("__eqS"), makesymmetric(st_matrix(st_local("__eqS"))))
                    matrix colnames `__eqS' = `__eqS_names'
                    matrix rownames `__eqS' = `__eqS_names'
                    ereturn matrix S = `__eqS'
                }
            }
            if "`save_stage'" == "first" {
                ereturn local estimates_title "First-stage regression: `save_dep'"
                local save_name = strtoname("`save_prefix'`save_dep'")
                local __firsteqs `__firsteqs' `save_name'
            }
            else {
                ereturn local estimates_title "Reduced-form regression: `save_dep'"
                local save_name = strtoname("`save_rfprefix'`save_dep'")
                local __rfeq `save_name'
            }
            estimates store `save_name', title(`e(estimates_title)')
        }
        _estimates unhold `save_hold'
    }
    ereturn local firsteqs "`__firsteqs'"
    ereturn local rfeq "`__rfeq'"
    ereturn local sfirsteq "`__sfirsteq'"

    * Display (nooutput shows nothing; results are still stored)
    if "`output'" == "" {
        local __stored
        if "`savesfirst'" != "" local __stored `__sfirsteq'
        if "`saverf'" != "" local __stored `__stored' `__rfeq'
        if "`savefirst'" != "" local __stored `__stored' `__firsteqs'
        _civreghdfe_disp, `first' `ffirst' `rf' `sfirst' level(`level') `header' `footer' ///
            dispopt(`dispopt') stored(`__stored') `nowildtable'
    }

    * Equations stored only for the display are dropped, as in ivreg2
    if "`savefirst'" == "" {
        foreach __eq of local __firsteqs {
            capture estimates drop `__eq'
        }
        ereturn local firsteqs
    }
    if "`saverf'" == "" {
        if "`__rfeq'" != "" capture estimates drop `__rfeq'
        ereturn local rfeq
    }
    if "`savesfirst'" == "" {
        if "`__sfirsteq'" != "" capture estimates drop `__sfirsteq'
        ereturn local sfirsteq
    }

    * Capture post-processing time (before cleanup)
    _ctools_lap t_stata_post __ct_lap98 `__ct_lap98'

    * Display timing breakdown if verbose or timeit specified
    if "`verbose'" != "" | `__timeit' {
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

    _ctools_cleanup civreghdfe
end

* vce(bootstrap), fevce(bootstrap): run the replicates on the fit the plugin
* retained, in chunks of about 50 (dots and Break between chunks; a multiple
* of the replicates the plugin runs at once, __civreghdfe_bs_par), then
* summarize them into the __civreghdfe_bs_* matrices. Returns r(rngstate),
* the random-number state before the first draw (plan 5.8).
program define _civreghdfe_bootstrap, rclass
    version 14.1
    syntax, REPS(integer) DRAWS(name) DOTS(integer) LEVEL(real) MSE(integer) KB(integer) ///
        [SEED(string) THREADS(string) COEF(integer 1) FE(integer 0) FEVARS(string) ///
        NFE(integer 1) NREPS(integer 0) CI(integer 0) DSE(integer 0) BSDRAWS(string) ///
        SAVEFILE(string) EVERY(string) STRATAOBS]
    local targets = cond(`fe', " fe", "") + cond(`dse', " dse", "") + cond(`ci', " ci", "")
    _ctools_hdfe_boot run, cmd(civreghdfe) reps(`reps') draws(`draws') dots(`dots') kb(`kb') ///
        seed(`"`seed'"') threads(`"`threads'"') targets(`"`targets'"') bsdraws(`bsdraws') ///
        savefile(`"`savefile'"') every(`every') `strataobs'
    return local rngstate `"`r(rngstate)'"'
    matrix __civreghdfe_bs_V = J(`kb', `kb', .)
    foreach m in b se bias z0 {
        matrix __civreghdfe_bs_`m' = J(1, `kb', .)
    }
    foreach m in ci_normal ci_percentile ci_bc {
        matrix __civreghdfe_bs_`m' = J(2, `kb', .)
    }
    matrix __civreghdfe_bs_counts = J(1, 5, .)
    if `fe' matrix __civreghdfe_bs_fe_table = J(`nfe', 5, .)
    local finish = cond(`coef', " coef", "") + cond(`nreps', " nreps", "") + cond(`ci', " ci", "") + ///
        cond(`mse', " mse", "")
    plugin call ctools_plugin `fevars', ///
        "civreghdfe `threads' boot_finish level=`level'`finish'"
end

* vce(bootstrap): the e() results Stata's vce(bootstrap) posts, so that
* estat bootstrap works
program define _civreghdfe_bootstrap_post, eclass
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
    matrix `reps' = J(1, colsof(e(b)), scalar(__civreghdfe_bs_N_reps))
    matrix colnames `reps' = `names'
    matrix rownames `reps' = y1
    matrix `counts' = __civreghdfe_bs_counts
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
    ereturn scalar N_reps = scalar(__civreghdfe_bs_N_reps)
    ereturn scalar N_misreps = scalar(__civreghdfe_bs_N_misreps)
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
    ereturn local cmdname civreghdfe
    ereturn local rngstate `"`rngstate'"'
    ereturn local draws `draws'
    if "`cluster'" != "" ereturn local cluster_bs `cluster'
    if "`strata'" != "" ereturn local strata_bs `strata'
    ereturn local prefix bootstrap
    ereturn local vce bootstrap
    * as _bs_sum: with mse, e(V) is about e(b) and labelled so
    ereturn local vcetype = cond(`mse', "Bstrap *", "Bootstrap")
    if `mse' ereturn local mse mse
end

* fevce(bootstrap): e() for the fixed effects' bootstrap (plan 4.4)
program define _civreghdfe_fe_post, eclass
    version 14.1
    syntax, ABSVARS(string) NABSORB(integer) [SEVARS(string) NREPVARS(string) DSE(string) ///
        LBVARS(string) UBVARS(string) LEVEL(real 95) ///
        RNGSTATE(string) DRAWS(name) CLUSTER(string) STRATA(string)]
    tempname t
    matrix `t' = __civreghdfe_bs_fe_table
    matrix colnames `t' = levels never always min_reps median_reps
    matrix rownames `t' = `absvars'
    ereturn matrix fe_bs_table = `t'
    ereturn scalar fe_N_reps = scalar(__civreghdfe_bs_fe_N_reps)
    ereturn scalar fe_split_reps = scalar(__civreghdfe_bs_fe_split_reps)
    ereturn scalar fe_norm_warning = (`nabsorb' >= 2)
    ereturn local fevce bootstrap
    ereturn local fe_se `sevars'
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
program define _civreghdfe_fe_warning
    version 14.1
    local G : word count `e(absorb)'
    di as text ""
    di as text "{p 0 9 2}warning: with `G' absorbed effects (`e(absorb)'), individual fixed effects are " ///
        "identified only up to a normalization within each connected set of levels. The bootstrap " ///
        "standard errors in `e(fe_se)' depend on that normalization and can be large where the sets " ///
        "are thinly connected. " _c
    if "`e(fe_dse)'" != "" di as text "The sum of the fixed effects does not depend on it; its " ///
        "standard errors are in `e(fe_dse)'.{p_end}"
    else di as text "The sum of the fixed effects does not depend on it: specify " ///
        "fevce(bootstrap, dse({it:newvar})) for its standard errors.{p_end}"
    if e(fe_split_reps) > 0 & e(fe_split_reps) < . {
        di as text "{p 9 9 2}In " e(fe_split_reps) " of " e(fe_N_reps) " replicates a connected set " ///
            "split; levels in a split set were treated as not estimable in those replicates (see " ///
            "e(fe_bs_table)).{p_end}"
    }
end

* Release the fit the plugin retained for the bootstrap, after a failure
program define _civreghdfe_bs_release
    version 14.1
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    plugin call ctools_plugin, "civreghdfe boot_release"
end

* A list after a label, as ivreg2's Disp shows it: from column col(),
* wrapped at 80 columns, base levels and omitted terms left out
program define _civreghdfe_displist
    version 14.1
    syntax [anything] [, col(integer 15)]
    local maxlen = 80 - `col'
    local len = 0
    local first = 1
    foreach vn in `anything' {
        _ms_parse_parts `vn'
        if !r(omit) {
            local vnlen : length local vn
            if `len' + `vnlen' > `maxlen' {
                di
                local first = 1
                local len = `vnlen'
            }
            else local len = `len' + `vnlen' + 1
            if `first' {
                local first = 0
                di as text _col(`col') "`vn'" _c
            }
            else di as text " `vn'" _c
        }
    }
    di
end

* ivreg2's "1" varlist of a list of model terms: without base levels and
* omitted terms (r(list))
program define _civreghdfe_list1, rclass
    version 14.1
    local out
    foreach v of local 0 {
        if "`v'" == "_cons" {
            local out `out' `v'
            continue
        }
        _ms_parse_parts `v'
        if !r(omit) local out `out' `v'
    }
    return local list `out'
end

* The instruments kept by the plugin (__znames in the plugin's order:
* the zxct included ones, then the excluded ones; the collinear ones
* dropped), in ivreg2's order (excluded, then included): r(names) and
* r(perm), their positions among the plugin's kept columns
program define _civreghdfe_zorder, rclass
    version 14.1
    syntax, zxct(integer) [znames(string)]
    local kx
    local ke
    local nx = 0
    local j = 0
    foreach v of local znames {
        local ++j
        local drop = 0
        if `j' <= `zxct' {
            capture local drop = (scalar(__civreghdfe_collinear_`j') == 1)
            if !`drop' {
                local ++nx
                local kx `kx' `v'
            }
        }
        else {
            local jj = `j' - `zxct'
            capture local drop = (scalar(__civreghdfe_zcollin_`jj') == 1)
            if !`drop' local ke `ke' `v'
        }
    }
    local perm
    local ne : word count `ke'
    forvalues i = 1/`ne' {
        local perm `perm' `=`nx' + `i''
    }
    forvalues i = 1/`nx' {
        local perm `perm' `i'
    }
    return local names `ke' `kx'
    return local perm `perm'
end

* A varlist option (partial()) checked and unabbreviated as syntax does
program define _civreghdfe_vlcheck, rclass
    version 14.1
    syntax varlist(fv ts)
    return local varlist `varlist'
end

* reghdfe's technique(), transform(), acceleration() and preconditioner(),
* validated as reghdfe validates them (abbreviations included)
program define _civreghdfe_solver_opts
    version 14.1
    syntax [, TECHnique(string) TRAnsform(string) ACCELeration(string) PREConditioner(string)]
    if "`technique'" != "" & !inlist("`technique'", "map", "cg", "lsmr", "lsqr") {
        di as error "invalid technique: `technique'"
        exit 198
    }
    foreach o in transform acceleration preconditioner {
        if `"``o''"' == "" continue
        local v = lower(`"``o''"')
        if "`o'" == "transform" local valid cimmino kaczmarz symmetric_kaczmarz rand_kaczmarz
        if "`o'" == "acceleration" {
            if "`v'" == "cg" local v conjugate_gradient
            if "`v'" == "sd" local v steepest_descent
            if "`v'" == "off" local v none
            local valid conjugate_gradient steepest_descent aitken none
        }
        if "`o'" == "preconditioner" {
            if "`v'" == "off" local v none
            local valid none diagonal block_diagonal
        }
        foreach x of local valid {
            if strpos("`x'", "`v'") == 1 local v `x'
        }
        if !`: list v in valid' {
            di as error "invalid `o': ``o''"
            exit 198
        }
    }
end

* verbose(#) among the options (given with their comma): caller locals
* __verb_given (1 if verbose(#) was given) and __verb_rest (the options,
* with verbose for a positive level and without it otherwise)
program define _civreghdfe_verbose_opt
    version 14.1
    syntax [, Verbose(string) *]
    c_local __verb_given = (`"`verbose'"' != "")
    if `"`verbose'"' == "" exit
    capture confirm integer number `verbose'
    if _rc {
        di as error "option verbose() incorrectly specified"
        exit 198
    }
    c_local __verb_rest `"`options'`=cond(`verbose' > 0, " verbose", "")'"'
end

* wildbootstrap or wildbootstrap(wildopts) among the options (given with
* their comma): caller locals __wild (0/1), __wildspec and __wild_rest (the
* other options)
program define _civreghdfe_wild_opt
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

* wildbootstrap(): the WRE bootstrap on the fit the plugin retained. The
* weights come from Stata's generator as boottest draws them (enumerated in
* the plugin when all 2^C Rademacher vectors fit in reps). The tested
* coefficients are positions among the plugin's coefficients (names); the
* analytic standard errors scale the interval search. Results: matrices
* __civreghdfe_wb_p, _ci, _cistat; r(test), r(reps), r(enumerated), r(C),
* r(path), r(rngstate) (before the first draw).
program define _civreghdfe_wild, rclass
    version 14.1
    syntax, REPS(integer) TYPE(name) NULL(integer) PTYPE(integer) ENUMERATE(integer) ///
        PARTIALLED(integer) CI(integer) LEVEL(real) DOTS(integer) VCETYPE(integer) ///
        [NAMES(string) COEFS(string) SEED(string) THREADS(string) VCEVAR(string) ///
        BOOTVAR(string) KEEP(integer 0) PROJECTED(integer 0)]
    local K : word count `names'
    local test
    if `"`coefs'"' == "" {
        forvalues k = 1/`K' {
            if "`: word `k' of `names''" != "_cons" & __civreghdfe_V[`k', `k'] > 0 ///
                local test `test' `k'
        }
    }
    else {
        foreach c of local coefs {
            local p : list posof `"`c'"' in names
            if !`p' {
                di as error "civreghdfe: wildbootstrap(coefs()): `c' is not a regressor (or was partialled out)"
                exit 198
            }
            if !(__civreghdfe_V[`p', `p'] > 0) {
                di as error "civreghdfe: wildbootstrap(coefs()): `c' was omitted"
                exit 198
            }
            local test `test' `p'
        }
    }
    local T : word count `test'
    if `T' == 0 {
        di as error "civreghdfe: wildbootstrap() has no coefficient to test"
        exit 198
    }
    matrix __civreghdfe_wb_test = J(1, `T', .)
    matrix __civreghdfe_wb_scale = J(1, `T', .)
    forvalues j = 1/`T' {
        local p : word `j' of `test'
        matrix __civreghdfe_wb_test[1, `j'] = `p'
        matrix __civreghdfe_wb_scale[1, `j'] = sqrt(__civreghdfe_V[`p', `p'])
    }
    if `"`seed'"' != "" set seed `seed'
    return local rngstate `"`c(rngstate)'"'
    local flags = " vcetype=`vcetype'" + cond("`vcevar'" != "", " vce", "") + ///
        cond("`bootvar'" != "", " boot", "") + cond(!`null', " nonull", "") + ///
        cond(`partialled', " partialled", "") + cond(`projected', " projected", "")
    plugin call ctools_plugin `vcevar' `bootvar', "civreghdfe `threads' wild_setup`flags'"
    local C = scalar(__civreghdfe_wb_C)
    local path = scalar(__civreghdfe_wb_path)
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
                    "civreghdfe `threads' wild_chunk first=`first' count=`n' total=`B' enum"
            }
            else {
                _ctools_hdfe_boot wildweights, first(`=`first' + 1') count(`n') variables(`owned_weights')
                local wvars `r(varlist)'
                capture noisily plugin call ctools_plugin `wvars', ///
                    "civreghdfe `threads' wild_chunk first=`first' count=`n' total=`B'"
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
    if `dots' di as text "(wild bootstrap: `B' replications" ///
        cond(`ci', ", intervals by test inversion", "") ")"
    matrix __civreghdfe_wb_p = J(1, `T', .)
    matrix __civreghdfe_wb_ci = J(`T', 2, .)
    matrix __civreghdfe_wb_cistat = J(`T', 5, 0)
    plugin call ctools_plugin, "civreghdfe `threads' wild_finish level=`level' ptype=`ptype'`=cond(`ci', " ci", "")'`=cond(`keep', " keep", "")'"
    return local test `test'
    return scalar reps = `B'
    return scalar enumerated = `enum'
    return scalar C = `C'
    return scalar path = `path'
end

* wildbootstrap(): its e() results (plan 4.7), with the prefix wb_
program define _civreghdfe_wild_post, eclass
    version 14.1
    syntax, TEST(numlist) NAMES(string) METHOD(name) TYPE(name) PTYPE(name) REPS(integer) ///
        ENUM(integer) NCLUST(integer) LEVEL(real) PATH(integer) CI(integer) ///
        [BOOTCLUST(string) RNGSTATE(string)]
    local coefs
    foreach p of local test {
        local coefs `coefs' `: word `p' of `names''
    }
    local T : word count `test'
    tempname t p cim stm
    matrix `t' = J(1, `T', .)
    forvalues j = 1/`T' {
        local c : word `j' of `coefs'
        matrix `t'[1, `j'] = _b[`c'] / _se[`c']
    }
    matrix `p' = __civreghdfe_wb_p
    matrix colnames `t' = `coefs'
    matrix colnames `p' = `coefs'
    ereturn matrix wb_t = `t'
    ereturn matrix wb_p = `p'
    if `ci' {
        matrix `cim' = __civreghdfe_wb_ci
        matrix `stm' = __civreghdfe_wb_cistat
        matrix rownames `cim' = `coefs'
        matrix colnames `cim' = ll ul
        matrix rownames `stm' = `coefs'
        matrix colnames `stm' = open_lower open_upper disjoint empty failed
        ereturn matrix wb_ci = `cim'
        ereturn matrix wb_cistatus = `stm'
    }
    ereturn scalar wb_reps = `reps'
    ereturn scalar wb_misreps = scalar(__civreghdfe_wb_misreps)
    ereturn scalar wb_enumerated = `enum'
    ereturn scalar wb_N_clust = `nclust'
    ereturn scalar wb_level = `level'
    ereturn local wb_method `method'
    ereturn local wb_weighttype `type'
    ereturn local wb_ptype `ptype'
    ereturn local wb_coefs `coefs'
    if "`bootclust'" != "" ereturn local wb_bootclust `bootclust'
    ereturn local wb_fe_path = word("nested projected partialled projected", `path' + 1)
    ereturn local wb_rngstate `"`rngstate'"'
end

* Flag the plugin columns named by an orthog()/endogtest()/redundant()/
* partial() varlist. The list is fvexpand-ed and matched term by term against
* the model's expanded terms `all' (base/omitted markers ignored, as ivreg2's
* dictionary match does). Plugin columns are the terms not made only of
* base/omitted components (the rule that built the plugin varlist); matched
* base levels are skipped and duplicates ignored. Returns r(n) flagged
* columns, r(flags) (1 x #plugin columns, if r(n) > 0), r(terms) the matched
* model terms, and r(bad) the first name that matches no term.

* Term name without base/omitted markers (1b.g -> 1.g, 2bn.g -> 2.g,
* o.x -> x, co.x -> c.x), so an option list matches the model's terms
* whatever base each fvexpand call picked

* First-stage and reduced-form regressions (first, rf, savefirst, saverf):
* civreghdfe itself, run on the caller's command line. Every run starts by
* dropping the __civreghdfe_* and _civreghdfe_* scalars and matrices, and
* the caller still reads its own afterwards, so they are set aside and
* restored around the nested run.
program define _civreghdfe_aux, eclass
    version 14.1
    local n = 0
    foreach pattern in __civreghdfe_* _civreghdfe_* {
        foreach kind in scalar matrix {
            if "`kind'" == "scalar" local names : all scalars "`pattern'"
            else local names : all matrices "`pattern'"
            foreach name of local names {
                local ++n
                tempname hold`n'
                `kind' `hold`n'' = `name'
                local kind`n' `kind'
                local name`n' `name'
            }
        }
    }
    capture noisily civreghdfe `0'
    local rc = _rc
    _ctools_cleanup civreghdfe
    forvalues i = 1/`n' {
        `kind`i'' `name`i'' = `hold`i''
    }
    exit `rc'
end

* Stock-Yogo (2005) critical values for the Cragg-Donald F statistic, as
* ivreg2 reports them for each estimator. The tables are those of ivreg2's
* Mata library (s_cdsy() in livreg2.mlib, installed with ivreg2 and with
* ivreghdfe). Without it, the four IV size values for one endogenous
* regressor and one to three instruments, and for two endogenous regressors
* and three instruments, are shown; any other cell is reported as not
* available. Returns r(cdmissing) = 1 when no value was shown.
program define _civreghdfe_cdsy, rclass
    version 14.1
    syntax, model(string) k2(integer) nendog(integer) [col1(integer 42) col2(integer 73) noHEADing]
    if inlist("`model'", "iv") {
        local tables 1 2 3 4 5 6 7 8
        local labels `"" 5% maximal IV relative bias" "10% maximal IV relative bias" "20% maximal IV relative bias" "30% maximal IV relative bias" "10% maximal IV size" "15% maximal IV size" "20% maximal IV size" "25% maximal IV size""'
    }
    else if "`model'" == "liml" {
        local tables 17 18 19 20
        local labels `""10% maximal LIML size" "15% maximal LIML size" "20% maximal LIML size" "25% maximal LIML size""'
    }
    else if "`model'" == "fuller" {
        local tables 9 10 11 12 13 14 15 16
        local labels `"" 5% maximal Fuller rel. bias" "10% maximal Fuller rel. bias" "20% maximal Fuller rel. bias" "30% maximal Fuller rel. bias" " 5% Fuller maximum bias" "10% Fuller maximum bias" "20% Fuller maximum bias" "30% Fuller maximum bias""'
    }
    tempname table
    capture mata: s_cdsy("`table'", 1)
    local have_tables = (_rc == 0)
    local shown = 0
    local i = 0
    foreach t of local tables {
        local ++i
        local label : word `i' of `labels'
        local cv = .
        * ivreg2's limits: up to 100 instruments; 3 endogenous regressors for
        * the IV bias tables, 2 for the others
        local max_endog = cond(`t' <= 4, 3, 2)
        if `have_tables' & `k2' >= 1 & `k2' <= 100 & `nendog' >= 1 & `nendog' <= `max_endog' {
            capture mata: s_cdsy("`table'", `t')
            if !_rc local cv = `table'[`k2', `nendog']
        }
        else if !`have_tables' & inrange(`t', 5, 8) {
            local cell `t'_`k2'_`nendog'
            local cv = cond("`cell'" == "5_1_1", 16.38, cond("`cell'" == "6_1_1", 8.96, ///
                cond("`cell'" == "7_1_1", 6.66, cond("`cell'" == "8_1_1", 5.53, ///
                cond("`cell'" == "5_2_1", 19.93, cond("`cell'" == "6_2_1", 11.59, ///
                cond("`cell'" == "7_2_1", 8.75, cond("`cell'" == "8_2_1", 7.25, ///
                cond("`cell'" == "5_3_1", 22.30, cond("`cell'" == "6_3_1", 12.83, ///
                cond("`cell'" == "7_3_1", 9.54, cond("`cell'" == "8_3_1", 7.80, ///
                cond("`cell'" == "5_3_2", 13.43, cond("`cell'" == "6_3_2", 8.18, ///
                cond("`cell'" == "7_3_2", 6.40, cond("`cell'" == "8_3_2", 5.45, .))))))))))))))))
        }
        if missing(`cv') continue
        if !`shown' & "`heading'" == "" di as text "Stock-Yogo weak ID test critical values:" _continue
        di as text _col(`col1') "`label'" as result _col(`col2') %6.2f `cv'
        local shown = 1
    }
    if !`shown' & "`heading'" != "" {
        if "`model'" == "fuller" di as text "NB: Critical values based on Fuller parameter=1"
    }
    else if !`shown' {
        di as text "Stock-Yogo weak ID test critical values:" _col(64) "<not available>"
        if !`have_tables' & "`model'" != "kclass" {
            di as text "(the full table is read from ivreg2; ssc install ivreg2)"
        }
    }
    else if "`model'" == "fuller" {
        di as text "NB: Critical values based on Fuller parameter=1"
    }
    return scalar cdmissing = !`shown'
end

* absorb(absvars, savefe): reghdfe's only suboption; anything else is an error
program define _civreghdfe_absorb_opts
    version 14.1
    syntax [, SAVEfe]
    c_local savefe = ("`savefe'" != "")
end

* The system of first-stage and reduced-form equations (ivreg2's sfirst),
* posted from the plugin's estimates and stored as estimates `name'. Its
* VCE is the plugin's system VCE with ivreg2's small-sample factor
* (N-d)/(N-L-d-s), or (N-1)/(N-L-s) G/(G-1) with clusters, and the residual
* dof N-L-d-s. The equations are the dependent variable and the endogenous
* regressors kept; the instruments are listed as ivreg2 lists them
* (excluded, then included). The caller's estimation results are kept.
* Sets c_local __sfirsteq to the stored name ("" if it failed).
program define _civreghdfe_sfirst, eclass
    version 14.1
    syntax, k1(integer) l(integer) endo1(string) znames(string) zxct(integer) ///
        depvar(string) sample(varname) cluster(integer) nclust(real) ///
        dofminus(integer) sdofminus(real) name(string) [hacv(string) hacv2(string)]
    c_local __sfirsteq ""
    * the instruments kept, in the plugin's order (included, then excluded)
    local zx
    local ze
    local k_exog : word count `e(inexog)'
    local j = 0
    foreach v of local znames {
        local ++j
        if `j' <= `zxct' {
            local kk = `j'
            local drop = 0
            * included: the collinearity flags of the exogenous columns
            capture local drop = (scalar(__civreghdfe_collinear_`kk') == 1)
            if !`drop' local zx `zx' `v'
        }
        else {
            local jj = `j' - `zxct'
            local drop = 0
            capture local drop = (scalar(__civreghdfe_zcollin_`jj') == 1)
            if !`drop' local ze `ze' `v'
        }
    }
    local zkept `ze' `zx'
    if `: word count `zkept'' != `l' exit
    local D = (`k1' + 1) * `l'
    tempname b V S
    matrix `b' = __civreghdfe_sfb[1, 1..`D']
    matrix `V' = __civreghdfe_sfV[1..`D', 1..`D']
    matrix `S' = __civreghdfe_sfS[1..`D', 1..`D']
    local names
    foreach eq in `depvar' `endo1' {
        foreach z of local zkept {
            local names `names' `eq':`z'
        }
    }
    matrix colnames `b' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'
    matrix colnames `S' = `names'
    matrix rownames `S' = `names'
    local N = e(N)
    if !`cluster' matrix `V' = `V' * (`N' - `dofminus') / (`N' - `l' - `dofminus' - `sdofminus')
    else matrix `V' = `V' * (`N' - 1) / (`N' - `l' - `sdofminus') * `nclust' / (`nclust' - 1)
    local dof = `N' - `l' - `dofminus' - `sdofminus'
    * exactly symmetric, as ivreg2's makesymmetric() leaves it (ereturn post
    * zeroes a V that is indefinite and not exactly symmetric)
    mata: st_replacematrix("`V'", makesymmetric(st_matrix("`V'")))
    * VCE information of the main model
    foreach m in vcetype clustvar clustvar1 clustvar2 kernel tvar ivar {
        local m_`m' `"`e(`m')'"'
    }
    foreach m in N_clust N_clust1 N_clust2 bw {
        local s_`m' = e(`m')
    }
    tempname hold
    _estimates hold `hold', restore
    tempvar es
    quietly gen byte `es' = `sample'
    * (no dependent variable: ivreg2 posts the system without one)
    ereturn post `b' `V', obs(`N') esample(`es') dof(`dof')
    ereturn scalar k_eq = `k1'
    ereturn scalar df_r = `dof'
    ereturn local cmd "civreghdfe"
    ereturn local model "sfirst"
    ereturn matrix S = `S'
    if `"`hacv'"' != "" ereturn local hacsubtitleV `"`hacv'"'
    if `"`hacv2'"' != "" ereturn local hacsubtitleV2 `"`hacv2'"'
    if "`m_vcetype'" != "" ereturn local vcetype "Robust"
    foreach m in clustvar clustvar1 clustvar2 kernel tvar ivar {
        if `"`m_`m''"' != "" ereturn local `m' `"`m_`m''"'
    }
    foreach m in N_clust N_clust1 N_clust2 bw {
        if !missing(`s_`m'') ereturn scalar `m' = `s_`m''
    }
    capture estimates store `name', title("System of first-stage/reduced form regressions")
    if _rc {
        di
        di as result "Unable to store system of first-stage reduced form regressions."
        di
    }
    else c_local __sfirsteq `name'
    _estimates unhold `hold'
end

* Replay of a stored first-stage, reduced-form or system equation (after
* estimates restore), as ivreg2's DispSFirst shows it from the equation's
* own results
program define _civreghdfe_dispsf, eclass
    version 14.1
    syntax [, Level(cilevel) dispopt(string asis)]
    local model `e(model)'
    di
    if "`model'" == "first" di as text "First-stage regression of `e(depvar)':"
    else if "`model'" == "rf" {
        di as text "Reduced-form regression: `e(depvar)'"
        di as text "{hline `=length("`e(depvar)'") + 25'}"
    }
    else {
        di as text "System of first-stage/reduced-form regressions:"
        di as text "{hline 47}"
    }
    if `"`e(hacsubtitleV)'"' != "" di as text _n `"`e(hacsubtitleV)'"'
    if `"`e(hacsubtitleV2)'"' != "" di as text `"`e(hacsubtitleV2)'"'
    di as text "Number of obs = " _col(31) as result %8.0f e(N)
    if "`e(kernel)'" != "" {
        di as text "  kernel=`e(kernel)'; bandwidth=" `e(bw)'
        if "`e(bwchoice)'" != "" di as text "  `e(bwchoice)'"
        di as text "  time variable (t):  " as result e(tvar)
        if "`e(ivar)'" != "" di as text "  group variable (i): " as result e(ivar)
    }
    if "`e(clustvar)'" != "" {
        if "`e(clustvar2)'" == "" {
            local nc = e(N_clust)
            local cv `e(clustvar)'
        }
        else {
            local nc = e(N_clust1)
            local cv `e(clustvar1)'
        }
        di as text "Number of clusters (`cv') = " _col(33) as result %6.0f `nc'
    }
    if "`e(clustvar2)'" != "" {
        di as text "Number of clusters (`e(clustvar2)') = " _col(33) as result %6.0f e(N_clust2)
    }
    local vt `"`e(vcetype)'"'
    local cv `"`e(clustvar)'"'
    if `"`vt'"' != "" ereturn local vcetype "Robust"
    ereturn local clustvar
    capture noisily ereturn display, level(`level') `dispopt'
    local rc = _rc
    ereturn local vcetype `"`vt'"'
    ereturn local clustvar `"`cv'"'
    exit `rc'
end

* One stored first-stage, reduced-form or system equation, displayed as
* ivreg2's DispSFirst does: title, what the statistics are robust to, the
* number of observations, kernel and cluster lines, and the table. The
* VCE information comes from the main results (in memory); the table from
* the stored equation `eq'.
program define _civreghdfe_dispeq, eclass
    version 14.1
    syntax, kind(string) eq(string) [Level(cilevel) dispopt(string asis)]
    c_local __eqdep ""
    local hac `"`e(_civ_hacaux)'"'
    local hacV2 `"`e(_civ_hacV2)'"'
    local kernel `"`e(kernel)'"'
    local bw `"`e(bw)'"'
    local tvar `"`e(tvar)'"'
    local ivar `"`e(ivar)'"'
    local vce = e(_civ_vce)
    local cl1 `"`e(_civ_cl1)'"'
    local cl2 `"`e(_civ_cl2)'"'
    local nc = e(N_clust)
    local nc1 = e(N_clust1)
    local nc2 = e(N_clust2)
    tempname hold
    _estimates hold `hold', copy
    capture estimates restore `eq'
    if _rc {
        _estimates unhold `hold'
        exit 111
    }
    di
    if "`kind'" == "first" di as text "First-stage regression of `e(depvar)':"
    else if "`kind'" == "rf" {
        di as text "Reduced-form regression: `e(depvar)'"
        di as text "{hline `=length("`e(depvar)'") + 25'}"
    }
    else {
        di as text "System of first-stage/reduced-form regressions:"
        di as text "{hline 47}"
    }
    di
    if "`hac'" != "" di as text "Statistics robust to `hac'"
    else di as text "Statistics consistent for homoskedasticity only"
    if "`hacV2'" != "" di as text "`hacV2'"
    di as text "Number of obs = " _col(31) as result %8.0f e(N)
    if "`kernel'" != "" {
        di as text "  kernel=`kernel'; bandwidth=" `bw'
        di as text "  time variable (t):  " as result "`tvar'"
        if "`ivar'" != "" di as text "  group variable (i): " as result "`ivar'"
    }
    if inlist(`vce', 2, 4) {
        di as text "Number of clusters (`cl1') = " _col(33) as result %6.0f `nc'
    }
    else if `vce' == 3 {
        di as text "Number of clusters (`cl1') = " _col(33) as result %6.0f `nc1'
        di as text "Number of clusters (`cl2') = " _col(33) as result %6.0f `nc2'
    }
    * The table, as ivreg2 shows its stored equations: robust standard
    * errors labelled Robust, no cluster-adjustment note; with absorb() the
    * constant is partialled out
    local vt `"`e(vcetype)'"'
    local cv `"`e(clustvar)'"'
    if `vce' != 0 ereturn local vcetype "Robust"
    ereturn local clustvar
    c_local __eqdep `"`e(depvar)'"'
    if "`kind'" == "first" ereturn display
    else ereturn display, level(`level') `dispopt'
    ereturn local vcetype `"`vt'"'
    ereturn local clustvar `"`cv'"'
    _estimates unhold `hold'
end

* civreghdfe's display, at estimation and on replay, from the results in
* e(): the reduced form (rf), the first stages (first), the system of
* first-stage equations (sfirst) and ivreg2's first-stage summary (ffirst,
* also with first and sfirst), then the header, the coefficient table and
* the diagnostics footer, as ivreg2 orders them.
program define _civreghdfe_disp, eclass
    version 14.1
    syntax [, FIRST FFIRst RF SFIRST Level(cilevel) noHEADer noFOOTer dispopt(string asis) ///
        stored(string) NOWILDtable]
    local K_endog = e(K_endog)
    local vce_type = e(_civ_vce)
    local est_method = e(_civ_est)
    local liml_kclass = e(_civ_lk)
    local is_b0 = e(_civ_b0)
    local hac `"`e(_civ_hac)'"'
    local cluster_var `"`e(_civ_cl1)'"'
    local cluster_var2 `"`e(_civ_cl2)'"'
    local kernel_on = ("`e(kernel)'" != "")
    local G = e(G)

    * the saved equations (savesfirst, saverf, savefirst), as ivreg2 lists them
    if "`stored'" != "" {
        di
        di as text "Stored estimation results"
        di as text "{hline 25}" _c
        capture estimates dir `stored'
        if !_rc estimates dir `stored'
        else {
            di
            di as result "Unable to list stored estimations."
            di
        }
    }

    * rf: the stored reduced form
    if "`rf'" != "" {
        if "`e(rfeq)'" == "" {
            di
            di as result "Unable to display stored reduced form estimation."
            di
        }
        else {
            capture noisily _civreghdfe_dispeq, kind(rf) eq(`e(rfeq)') level(`level') dispopt(`dispopt')
            if _rc {
                di
                di as result "Unable to display stored reduced form estimation."
                di
            }
        }
    }

    * first: each stored first stage and its F tests (ivreg2's DispFirst)
    if "`first'" != "" & `K_endog' > 0 {
        tempname fm
        capture matrix `fm' = e(first)
        if _rc {
            di
            di as result "Unable to display first-stage estimates; macro e(first) is missing"
        }
        else {
            di
            di as text "First-stage regressions"
            di as text "{hline 23}"
            di
            local endo1 : colnames `fm'
            local firsteqs `"`e(firsteqs)'"'
            if `: word count `firsteqs'' < `: word count `endo1'' {
                di as result "Unable to display all first-stage regressions."
                di as result "There may be insufficient room to store results using -estimates store-,"
                di as result "or names of endogenous regressors may be too long to store the results."
                di as result "Try dropping one or more estimation results using -estimates drop-,"
                di as result "using the -savefprefix- option, or using shorter variable names."
                di
            }
            foreach eqname of local firsteqs {
                capture noisily _civreghdfe_dispeq, kind(first) eq(`eqname')
                if _rc {
                    di
                    di as result "Unable to list stored estimation `eqname'."
                    di
                    continue
                }
                * the equation's dependent variable is the column of e(first)
                local vn = colnumb(`fm', `"`__eqdep'"')
                if missing(`vn') continue
                di as text "F test of excluded instruments:"
                di as text "  F(" %3.0f el(`fm', 5, `vn') "," %6.0f el(`fm', 6, `vn') ") = " ///
                    as result %8.2f el(`fm', 4, `vn')
                di as text "  Prob > F      = " as result %8.4f el(`fm', 7, `vn')
                di as text "Sanderson-Windmeijer multivariate F test of excluded instruments:"
                di as text "  F(" %3.0f el(`fm', 9, `vn') "," %6.0f el(`fm', 10, `vn') ") = " ///
                    as result %8.2f el(`fm', 8, `vn')
                di as text "  Prob > F      = " as result %8.4f el(`fm', 11, `vn')
                di
            }
        }
    }

    * sfirst: the stored system of equations
    if "`sfirst'" != "" {
        local done = 0
        if "`e(sfirsteq)'" != "" {
            capture noisily _civreghdfe_dispeq, kind(sfirst) eq(`e(sfirsteq)') level(`level') dispopt(`dispopt')
            local done = !_rc
        }
        if !`done' {
            di
            di as result "Unable to display stored first-stage/reduced form estimations."
            di
        }
    }

    * ivreg2's summary of the first-stage and identification statistics
    if "`first'`ffirst'`sfirst'" != "" & `K_endog' > 0 {
        _civreghdfe_dispff
    }

    * Header
    if "`header'" == "" {
        local __title `"`e(_civ_title)'"'
        di as text ""
        di as text `"`__title'"'
        di as text "{hline `=udstrlen(`"`__title'"')'}"
        if `"`e(subtitle)'"' != "" di as text `"`e(subtitle)'"'
        if `K_endog' > 0 & (`liml_kclass' | inlist(`est_method', 1, 2, 3)) & !missing(e(kclass)) {
            di as text "k               =" %7.5f e(kclass)
        }
        if `K_endog' > 0 & (`liml_kclass' | inlist(`est_method', 1, 2)) & !missing(e(lambda)) {
            di as text "lambda          =" %7.5f e(lambda)
        }
        if `K_endog' > 0 & `est_method' == 2 {
            di as text "Fuller parameter=" %-5.0f e(fuller)
        }
        di as text ""
        di as text `"`e(_civ_hacB)'"'
        if "`hac'" != "" di as text "Statistics robust to `hac'"
        else di as text "Statistics consistent for homoskedasticity only"
        if `"`e(_civ_hacV2)'"' != "" di as text `"`e(_civ_hacV2)'"'
        if `kernel_on' {
            di as text "  kernel=`e(kernel)'; bandwidth=" `e(bw)'
            if "`e(bwchoice)'" != "" di as text "  `e(bwchoice)'"
            di as text "  time variable (t):  `e(tvar)'"
            if "`e(ivar)'" != "" di as text "  group variable (i): `e(ivar)'"
        }
        di as text ""

        * Clusters, observations, F, sums of squares and R2s in ivreg2's
        * formats (with clusters, the F test's denominator dof are the
        * clusters less one; two-way: the smaller dimension)
        if inlist(`vce_type', 2, 3, 4) {
            local __ncl = cond(`vce_type' == 3, e(N_clust1), e(N_clust))
            di as text "Number of clusters (`cluster_var') = " _col(33) as result %6.0f `__ncl' _continue
        }
        di as text _col(55) "Number of obs = " as result %8.0f e(N)
        if `vce_type' == 3 {
            di as text "Number of clusters (`cluster_var2') = " _col(33) as result %6.0f e(N_clust2) _continue
        }
        if "`e(vce)'" == "bootstrap" {
            di as text _c _col(55) "Wald chi2(" %3.0f e(_civ_wdf) ") = "
            if e(chi2) < 99999 di as result %8.2f e(chi2)
            else di as result %8.2e e(chi2)
            di as text _col(55) "Prob > chi2   = " as result %8.4f e(p)
            di as text _col(55) "Replications  = " as result %8.0f e(N_reps)
        }
        else {
            local model_df = e(Fdf1)
            local F_wald = e(F)
            local df_r = e(Fdf2)
            di as text _c _col(55) "F(" %3.0f `model_df' "," %6.0f `df_r' ") = "
            if `F_wald' < 99999 di as result %8.2f `F_wald'
            else di as result %8.2e `F_wald'
            di as text _col(55) "Prob > F      = " as result %8.4f e(Fp)
        }
        di as text "Total (centered) SS     = " as result %12.0g e(yyc) _continue
        di as text _col(55) "Centered R2   = " as result %8.4f e(r2c)
        di as text "Total (uncentered) SS   = " as result %12.0g e(yy) _continue
        di as text _col(55) "Uncentered R2 = " as result %8.4f e(r2u)
        di as text "Residual SS             = " as result %12.0g e(rss) _continue
        di as text _col(55) "Root MSE      = " as result %8.4g e(rmse)
        di as text ""
    }

    * Coefficient table: ivreg2 and ivreghdfe label it with the dependent
    * variable itself, whatever depname() says; as in ivreg2, without the
    * note on cluster-adjusted standard errors and with ivreg2's Robust label
    local depname_use `"`e(depvar)'"'
    local __cv `"`e(clustvar)'"'
    local __vt `"`e(vcetype)'"'
    ereturn local depvar `"`e(_civ_depvar)'"'
    ereturn local clustvar
    if `vce_type' != 0 & "`e(vce)'" != "bootstrap" ereturn local vcetype "Robust"
    ereturn display, level(`level') `dispopt'
    ereturn local depvar `"`depname_use'"'
    ereturn local clustvar `"`__cv'"'
    ereturn local vcetype `"`__vt'"'
    * vce(bootstrap): failed replicates, and what the diagnostics below use
    if "`e(vce)'" == "bootstrap" {
        if e(N_misreps) > 0 & e(N_misreps) < . {
            di as text "Note: One or more parameters could not be estimated in " ///
                as result e(N_misreps) as text " bootstrap replicates;"
            di as text "      standard-error estimates include only complete replications."
        }
        local __base = cond(inlist(`vce_type', 2, 3, 4), "cluster-robust", ///
            cond(`vce_type' == 1, "heteroskedasticity-robust", "conventional"))
        di as text "Note: diagnostics use the `__base' covariance; coefficient standard errors are bootstrap."
    }
    if "`e(fevce)'" == "bootstrap" & e(fe_norm_warning) == 1 _civreghdfe_fe_warning
    if "`e(wb_method)'" != "" & "`nowildtable'" == "" _ctools_wild_table

    if "`footer'" != "" exit

    * ivreghdfe's first footer (not for OLS without a J statistic): the
    * identification tests, the redundancy test, the overidentification
    * tests with orthog()'s C statistic, the endogeneity test and LIML's
    * Anderson-Rubin statistic
    local __skip1 = ("`e(model)'" == "ols" & (e(sargan) == 0 | e(j) == 0))
    if !`__skip1' {

    * Underidentification test, and redundant()
    if !e(_civ_noid) & `K_endog' > 0 & !missing(e(idstat)) {
        if `vce_type' == 0 & !`kernel_on' {
            di as text "Underidentification test (Anderson canon. corr. LM statistic):" ///
                _col(65) as result %14.3f e(idstat)
        }
        else {
            di as text "Underidentification test (Kleibergen-Paap rk LM statistic):" ///
                _col(65) as result %14.3f e(idstat)
        }
        di as text _col(52) "Chi-sq(" as result e(iddf) as text ") P-val =" ///
            as result %10.4f e(idp)
        if "`e(redlist)'" != "" {
            di as text "-redundant- option:"
            di as text _c "IV redundancy test (LM test of redundancy of specified instruments):"
            di as result _col(71) %8.3f e(redstat)
            di as text _col(52) "Chi-sq(" as result e(reddf) as text ") P-val =  " ///
                as result _col(73) %6.4f e(redp)
            di as text "Instruments tested: " _c
            _civreghdfe_displist `e(redlist)', col(23)
        }
        di as text "{hline 78}"
    }

    * Weak identification test (with endogenous regressors; not with noid)
    if `K_endog' > 0 & !e(_civ_noid) {
        di as text "Weak identification test (Cragg-Donald Wald F statistic):" ///
            _col(65) as result %14.3f e(cdf)
        if (`vce_type' != 0 | `kernel_on') & e(kp_f) > 0 & !missing(e(kp_f)) {
            di as text "                         (Kleibergen-Paap rk Wald F statistic):" ///
                _col(65) as result %14.3f e(kp_f)
        }
        _civreghdfe_cdsy, model(`e(_civ_symodel)') k2(`=e(exexog_ct)') nendog(`K_endog')
        if !r(cdmissing) {
            di as text "Source: Stock-Yogo (2005).  Reproduced by permission."
            if `vce_type' != 0 | `kernel_on' {
                di as text "NB: Critical values are for Cragg-Donald F statistic and i.i.d. errors."
            }
        }
        di as text "{hline 78}"
    }

    * Sargan/Hansen J statistic with orthog()'s C statistic, or (when the
    * moment covariance matrix is rank deficient, which leaves the J
    * statistic missing) ivreg2's warning
    if !missing(e(sargan_df)) {
        local __hansen = (`vce_type' != 0 | "`e(model)'" == "gmmw")
        local __jname = cond(`__hansen', "Hansen J statistic", "Sargan statistic")
        if !missing(e(j)) {
            if `K_endog' > 0 di as text _c "`__jname' (overidentification test of all instruments):"
            else di as text _c "`__jname' (Lagrange multiplier test of excluded instruments):"
            di as result _col(71) %8.3f e(j)
            if e(jdf) {
                di as text _col(52) "Chi-sq(" as result e(jdf) as text ") P-val =  " ///
                    as result _col(73) %6.4f e(jp)
            }
            else di as text _col(50) "(equation exactly identified)"
            * orthog(): the J statistic of the equation without the tested
            * conditions, and their C statistic
            if !missing(e(cstat)) {
                if e(cstat) > 0 {
                    di as text "-orthog- option:"
                    if `vce_type' != 0 {
                        di as text _c "Hansen J statistic (eqn. excluding suspect orthog. conditions): "
                    }
                    else {
                        di as text _c "Sargan statistic (eqn. excluding suspect orthogonality conditions):"
                    }
                    di as result _col(71) %8.3f e(j) - e(cstat)
                    di as text _col(52) "Chi-sq(" as result e(jdf) - e(cstatdf) as text ") P-val =  " ///
                        as result _col(73) %6.4f chi2tail(e(jdf) - e(cstatdf), e(j) - e(cstat))
                    di as text _c "C statistic (exogeneity/orthogonality of suspect instruments): "
                    di as result _col(71) %8.3f e(cstat)
                    di as text _col(52) "Chi-sq(" as result e(cstatdf) as text ") P-val =  " ///
                        as result _col(73) %6.4f e(cstatp)
                    di as text "Instruments tested:  " _c
                    _civreghdfe_displist `e(clist)', col(23)
                }
                if e(cstat) == 0 {
                    di as text _n "Collinearity/identification problems in eqn. excl. suspect orthog. conditions:"
                    di as text "  C statistic not calculated for -orthog- option"
                }
            }
        }
        else {
            di as error "Warning: estimated covariance matrix of moment conditions not of full rank."
            di as error "         overidentification statistic not reported, and standard errors and"
            di as error "         model tests should be interpreted with caution."
            di as error "Possible causes:"
            if e(nocollin) == 1 {
                di as error "         collinearities in regressors or instruments (with -nocollin- option)"
            }
            if !missing(e(N_clust)) {
                di as error "         number of clusters insufficient to calculate robust covariance matrix"
            }
            if `kernel_on' {
                di as error "         covariance matrix of moment conditions not positive definite"
                di as error "         covariance matrix uses too many lags"
            }
            di as error "         singleton dummy variable (dummy with one 1 and N-1 0s or vice versa)"
            di as error "partial option may address problem."
        }
    }

    * endogtest(): the endogeneity test of the listed regressors
    if !missing(e(estat)) {
        if e(estat) > 0 {
            di as text "-endog- option:"
            di as text _c "Endogeneity test of endogenous regressors: "
            di as result _col(71) %8.3f e(estat)
            di as text _col(52) "Chi-sq(" as result e(estatdf) as text ") P-val =  " ///
                as result _col(73) %6.4f e(estatp)
            di as text "Regressors tested:  " _c
            _civreghdfe_displist `e(elist)', col(23)
        }
        if e(estat) == 0 {
            di as text _n "Collinearity/identification problems in restricted equation:"
            di as text "  Endogeneity test statistic not calculated for -endog- option"
        }
    }
    if !missing(e(sargan_df)) | !missing(e(estat)) di as text "{hline 78}"

    * LIML's Anderson-Rubin overidentification statistic (not robust)
    if "`e(model)'" == "liml" & `vce_type' == 0 & !`kernel_on' & !missing(e(arubin)) {
        if `K_endog' > 0 di as text _c "Anderson-Rubin statistic (overidentification test of all instruments):"
        else di as text _c "Anderson-Rubin statistic (LR test of excluded instruments):"
        di as result _col(72) %7.3f e(arubin)
        if e(arubindf) {
            di as text _col(52) "Chi-sq(" as result e(arubindf) as text ") P-val =  " ///
                as result _col(73) %6.4f e(arubinp)
        }
        else di as text _col(50) "(equation exactly identified)"
        di as text "{hline 78}"
    }

    }

    * Instruments dropped after the fixed effects were partialled out
    if e(_civ_instdrop) > 0 & !missing(e(_civ_instdrop)) {
        di as text "Collinearities detected among instruments: `=e(_civ_instdrop)' instrument(s) dropped"
    }

    * Instruments footer
    if `K_endog' > 0 di as text "Instrumented:         `e(_civ_f_instd)'"
    if "`e(_civ_f_incl)'" != "" di as text "Included instruments: `e(_civ_f_incl)'"
    if "`e(_civ_f_excl)'" != "" di as text "Excluded instruments: `e(_civ_f_excl)'"
    if "`e(_civ_f_part)'" != "" {
        di as text "Partialled-out:       `e(_civ_f_part)'"
        di as text "                      nb: total SS, model F and R2s are after partialling-out;"
        if !e(_civ_nopsmall) {
            di as text "                          any small-sample adjustments include partialled-out"
            di as text "                          variables in regressor count K"
        }
        else {
            di as text "                          any small-sample adjustments do not include"
            di as text "                          partialled-out variables in regressor count K"
        }
    }
    if "`e(dups)'" != "" di as text "Duplicates:" _col(23) "`e(dups)'"
    if "`e(collin)'" != "" di as text "Dropped collinear:" _col(23) "`e(collin)'"
    if "`e(ecollin)'" != "" di as text "Reclassified as exog:" _col(23) "`e(ecollin)'"
    di as text "{hline 78}"

    * Absorbed degrees of freedom table
    if `G' > 0 & !missing(`G') {
        tempname doft
        capture matrix `doft' = e(_civ_doft)
        if !_rc {
            di as text ""
            di as text "Absorbed degrees of freedom:"
            di as text "{hline 13}{c TT}{hline 39}{c TRC}"
            di as text " Absorbed FE {c |} Categories  - Redundant  = Num. Coefs {c |}"
            di as text "{hline 13}{c +}{hline 39}{c RT}"
            local has_nested = 0
            local has_slope = 0
            local __rows `"`e(extended_absvars)'"'
            if `"`__rows'"' == "" local __rows `"`e(absorb)'"'
            forval g = 1/`G' {
                local fevar : word `g' of `__rows'
                local cats = `doft'[`g', 1]
                local redundant = `doft'[`g', 2]
                local nested_marker " "
                if `doft'[`g', 3] == 1 {
                    local has_nested = 1
                    local nested_marker "*"
                }
                else if colsof(`doft') >= 4 {
                    if `doft'[`g', 4] == 1 {
                        local has_slope = 1
                        local nested_marker "?"
                    }
                }
                local coefs = `cats' - `redundant'
                di as text %12s abbrev("`fevar'", 12) " {c |}" as result %10.0fc `cats' as text "  " ///
                    as result %10.0fc `redundant' as text "  " as result %10.0fc `coefs' ///
                    as text "    `nested_marker'{c |}"
            }
            di as text "{hline 13}{c BT}{hline 39}{c BRC}"
            if `has_slope' {
                di as text "? = number of redundant parameters may be higher"
            }
            if `has_nested' {
                di as text "* = FE nested within cluster; treated as redundant for DoF computation"
            }
            di as text ""
        }
    }
end

* ivreg2's summary of the first-stage statistics (its DispFFirst), from
* e(first), the identification statistics and the weak-instrument-robust
* tests in e()
program define _civreghdfe_dispff
    version 14.1
    tempname fm
    capture matrix `fm' = e(first)
    if _rc {
        di
        di as result "Unable to display summary of first-stage estimates; macro e(first) is missing"
        exit
    }
    local endo : colnames `fm'
    local vce_type = e(_civ_vce)
    local kernel_on = ("`e(kernel)'" != "")
    di
    di
    di as text "Summary results for first-stage regressions"
    di as text "{hline 43}"
    di
    di _c _col(44) as text "(Underid)"
    di _col(65) as text "(Weak id)"
    di _c as text "Variable     |"
    di _c _col(16) as text "F("
    di _c as result _col(17) %3.0f el(`fm', 5, 1) as text "," as result %6.0f el(`fm', 6, 1) as text ")  P-val"
    di _c as text _col(37) "|"
    di _c _col(39) as text "SW Chi-sq("
    di _c as result %3.0f el(`fm', 9, 1) as text ") P-val"
    di _c as text _col(60) "|"
    di _c _col(62) as text "SW F("
    di as result _col(67) %3.0f el(`fm', 9, 1) as text "," as result %6.0f el(`fm', 10, 1) as text ")"
    local i = 0
    foreach vn of local endo {
        local ++i
        if udstrlen("`vn'") > 12 local vn = usubstr("`vn'", 1, 12)
        di _c as result %-12s "`vn'" _col(14) as text "|" _col(18) as result %8.2f el(`fm', 4, `i')
        di _c _col(28) as result %8.4f el(`fm', 7, `i')
        di _c _col(37) as text "|" _col(42) as result %8.2f el(`fm', 12, `i') _col(51) as result %8.4f el(`fm', 13, `i')
        di _col(60) as text "|" _col(65) as result %8.2f el(`fm', 8, `i')
    }
    di
    local rtype ""
    if `vce_type' != 0 {
        if inlist(`vce_type', 2, 3, 4) local rtype "cluster-robust"
        else if `kernel_on' local rtype "heteroskedasticity and autocorrelation-robust"
        else local rtype "heteroskedasticity-robust"
    }
    else if `kernel_on' local rtype "autocorrelation-robust"
    if "`rtype'" != "" {
        di as text "NB: first-stage test statistics `rtype'"
        di
    }
    local model `"`e(_civ_symodel)'"'
    di as text "Stock-Yogo weak ID F test critical values for single endogenous regressor:"
    _civreghdfe_cdsy, model(`model') k2(`=e(exexog_ct)') nendog(1) col1(36) col2(67) noheading
    if r(cdmissing) {
        di as text _col(64) "<not available>"
    }
    else {
        di as text "Source: Stock-Yogo (2005).  Reproduced by permission."
        if "`model'" == "iv" & e(_civ_est) == 0 & "`e(model)'" != "gmmw" & `vce_type' == 0 & !`kernel_on' {
            di as text "NB: Critical values are for Sanderson-Windmeijer F statistic."
        }
        else {
            di as text "NB: Critical values are for i.i.d. errors only."
        }
        di
    }
    if e(iddf) != el(`fm', 9, 1) {
        di as error "Warning: Error in calculating first-stage id statistics above;"
        di as error "         dof of SW statistics is " el(`fm', 9, 1) ", should be L-(K-1)=`=e(iddf)'."
    }
    di as text "Underidentification test"
    di as text "Ho: matrix of reduced form coefficients has rank=K1-1 (underidentified)"
    di as text "Ha: matrix has rank=K1 (identified)"
    if `vce_type' == 0 & !`kernel_on' di as result "Anderson canon. corr. LM statistic" _c
    else di as result "Kleibergen-Paap rk LM statistic" _c
    di as text _col(42) "Chi-sq(" as result e(iddf) as text ")=" %-7.2f as result e(idstat) ///
        _col(61) as text "P-val=" %6.4f as result e(idp)
    di
    di as text "Weak identification test"
    di as text "Ho: equation is weakly identified"
    di as result "Cragg-Donald Wald F statistic" _col(65) %8.2f e(cdf)
    if `vce_type' != 0 | `kernel_on' {
        di as result "Kleibergen-Paap Wald rk F statistic" _col(65) %8.2f e(kp_f)
    }
    di
    di as text "Stock-Yogo weak ID test critical values for K1=`=e(endog_ct)' and L1=`=e(exexog_ct)':"
    _civreghdfe_cdsy, model(`model') k2(`=e(exexog_ct)') nendog(`=e(endog_ct)') col1(36) col2(67) noheading
    if r(cdmissing) {
        di as text _col(64) "<not available>"
    }
    else {
        di as text "Source: Stock-Yogo (2005).  Reproduced by permission."
        if `vce_type' != 0 | `kernel_on' {
            di as text "NB: Critical values are for Cragg-Donald F statistic and i.i.d. errors."
        }
    }
    di
    di as text "Weak-instrument-robust inference"
    di as text "Tests of joint significance of endogenous regressors B1 in main equation"
    di as text "Ho: B1=0 and orthogonality conditions are valid"
    di as result _c "Anderson-Rubin Wald test"
    di as text _col(36) "F(" as result e(ardf) as text "," as result e(ardf_r) as text ")=" ///
        _col(49) as result %7.2f e(arf) _col(61) as text "P-val=" as result %6.4f e(arfp)
    di as result _c "Anderson-Rubin Wald test"
    di as text _col(36) "Chi-sq(" as result e(ardf) as text ")=" ///
        _col(49) as result %7.2f e(archi2) _col(61) as text "P-val=" as result %6.4f e(archi2p)
    di as result _c "Stock-Wright LM S statistic"
    di as text _col(36) "Chi-sq(" as result e(sstatdf) as text ")=" ///
        _col(49) as result %7.2f e(sstat) _col(61) as text "P-val=" as result %6.4f e(sstatp)
    di
    if "`rtype'" != "" {
        di as text "NB: Underidentification, weak identification and weak-identification-robust"
        di as text "    test statistics `rtype'"
        di
    }
    if inlist(`vce_type', 2, 4) {
        di as text "Number of clusters             N_clust  = " as result %10.0f e(N_clust)
    }
    else if `vce_type' == 3 {
        di as text "Number of clusters (1)         N_clust1 = " as result %10.0f e(N_clust1)
        di as text "Number of clusters (2)         N_clust2 = " as result %10.0f e(N_clust2)
    }
    di as text "Number of observations               N  = " as result %10.0f e(N)
    di as text "Number of regressors                 K  = " as result %10.0f e(rankxx)
    di as text "Number of endogenous regressors      K1 = " as result %10.0f e(endog_ct)
    di as text "Number of instruments                L  = " as result %10.0f e(rankzz)
    di as text "Number of excluded instruments       L1 = " as result %10.0f e(ardf)
    if "`e(partial)'" != "" {
        di as text "Number of partialled-out regressors/IVs = " as result %10.0f e(partial_ct)
        di as text "NB: K & L do not included partialled-out variables"
    }
end
