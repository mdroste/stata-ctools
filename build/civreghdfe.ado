*! version 1.2.4 20261103 github.com/mdroste/stata-ctools

program define civreghdfe, eclass byable(recall)
    version 14.1

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

    * Parse syntax
    * Basic syntax: civreghdfe depvar (endogvars = instruments) [exogvars], absorb() [options]
    * absorb() is optional - without it, runs as regular IV (no FE absorption)
    syntax anything(name=0 equalok) [if] [in] [aw fw pw/], ///
        [Absorb(string) ///
        vce(string) ///
        Robust ///
        CLuster(varlist) ///
        TOLerance(real 1e-8) ///
        MAXiter(integer 16000) ///
        ITERATE(integer 0) ///
        Verbose ///
        FIRST ///
        FFIRst ///
        RF ///
        SMALL ///
        LIML ///
        FULLER(real 0) ///
        Kclass(real 0) ///
        GMM2s ///
        CUE ///
        COVIV ///
        BW(real 0) ///
        KERNEL(string) ///
        DKRAAY(integer 0) ///
        KIEFER ///
        CENTer ///
        B0(string) ///
        RESiduals(name) ///
        RESiduals2 ///
        NOConstant ///
        Level(cilevel) ///
        noHEADer ///
        noFOOTer ///
        noOUTput ///
        TItle(string) ///
        SUBtitle(string) ///
        DEPname(string) ///
        noid ///
        ORThog(varlist fv ts) ///
        ENDOGtest(varlist fv ts) ///
        REDundant(varlist fv ts) ///
        PARTial(varlist fv ts) ///
        FWL(varlist fv ts) ///
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
        THReads(integer 0) ///
        COEFLegend ///
        CFORMAT(string) ///
        NOCOLLIN ///
        DOFadjustments(string) KEEPSINgletons SFIRST ///
        CIVRERUN CIVCOLLIN(string) CIVECOLLIN(string) CIVDUPS(string) ///
        NORETURN]

    * Drop whatever an earlier run left behind, so no stale state is read
    _ctools_cleanup civreghdfe

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

    * Handle fwl() as alias for partial()
    if "`fwl'" != "" & "`partial'" == "" {
        local partial "`fwl'"
    }
    else if "`fwl'" != "" & "`partial'" != "" {
        di as error "cannot specify both fwl() and partial()"
        exit 198
    }

    if "`center'" != "" {
        di as error "civreghdfe: center is not supported"
        exit 198
    }
    * sfirst: ivreg2's system of first-stage and reduced-form equations
    * (joint VCE) and its Sanderson-Windmeijer and Anderson-Rubin statistics
    * are not computed
    if "`sfirst'" != "" {
        di as error "civreghdfe: sfirst is not supported"
        exit 198
    }
    * dofadjustments() and keepsingletons, which ivreghdfe passes to reghdfe
    _ctools_hdfe_dof, `=lower(`"`dofadjustments'"')'
    local dof_type = `s(type)'
    local dof_clusters = `s(clusters)'

    * Handle residuals2 - auto-named residuals
    if "`residuals2'" != "" {
        local residuals "_civreghdfe_resid"
    }

    if "`residuals'" != "" & "`residuals2'" == "" _ctools_newvars `residuals'

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
    * creghdfe regressions.
    local absorb_names
    if `"`absorb'"' != "" {
        if strpos(`"`absorb'"', ",") {
            di as error "civreghdfe: absorb() takes no suboptions"
            exit 198
        }
        _ctools_hdfe_absorb `absorb', touse(`touse') cmd(civreghdfe)
        local n_absorb = r(n)
        forvalues g = 1/`n_absorb' {
            local absorb_name_`g' `r(name_`g')'
            local absorb_vars_`g' `r(vars_`g')'
            local absorb_group_`g' = r(group_`g')
        }
        local absorb
        forvalues g = 1/`n_absorb' {
            local absorb_names `absorb_names' `absorb_name_`g''
            if `absorb_group_`g'' {
                tempvar absorb_id_`g'
                quietly _ctools_group `absorb_vars_`g'' if `touse', generate(`absorb_id_`g'')
                local absorb `absorb' `absorb_id_`g''
            }
            else local absorb `absorb' `absorb_vars_`g''
        }
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
    * The first-stage and reduced-form regressions run through creghdfe,
    * which clusters in one dimension only: refuse before estimating rather
    * than fail after the main results are shown
    if `vce_type' == 3 & "`first'`rf'`savefirst'`saverf'" != "" {
        di as error "civreghdfe: first, rf, savefirst and saverf are not available with two-way clustering"
        exit 198
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
        else {
            di as error "Weight type `weight' not supported"
            exit 198
        }
        * Evaluate the weight expression; marks out and validates the weights
        tempvar evaluated_weight
        local weight_expr = trim(`"`exp'"')
        if substr(`"`weight_expr'"', 1, 1) == "=" local weight_expr = substr(`"`weight_expr'"', 2, .)
        _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"= `weight_expr'"')
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
    local bw_val = `bw'
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
    if `"`kernel'"' != "" | `bw_val' > 0 {
        if `bw_val' == 0 {
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
        else {
            di as error `"kernel(`kernel') not supported"'
            di as error "Valid kernels: bartlett (bar), parzen (par), quadratic spectral (qs), truncated (tru), tukey-hanning (thann)"
            exit 198
        }
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
        if `bw_val' == 1 & inlist(`kernel_type', 1, 2, 5) {
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
            if "`fe_var'" == "`cluster_var'" | "`fe_var'" == "`cluster_var2'" {
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
    scalar __civreghdfe_bw = `bw_val'
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
    scalar __civreghdfe_keepsingletons = ("`keepsingletons'" != "")
    scalar __civreghdfe_center = ("`center'" != "")

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
    * orthog() - test exogeneity of specified instruments
    local n_orthog = 0
    if "`orthog'" != "" {
        _civreghdfe_optflags "`orthog'" "`instruments_names_all'" `touse'
        if "`r(bad)'" != "" {
            local __bad `r(bad)'
            if `: list __bad in exogvars_names_all' {
                di as error "orthog(): `__bad' is an included exogenous regressor; civreghdfe tests excluded instruments only"
            }
            else di as error "orthog(): `__bad' is not an excluded instrument"
            exit 198
        }
        local n_orthog = r(n)
        if `n_orthog' > 0 matrix __civreghdfe_orthog_flags = r(flags)
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
    if `has_cons' & "`partial'" != "" local partial_cons = 1
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

    tempvar stored_resid
    quietly gen double `stored_resid' = .
    local plugin_vars `plugin_vars' `stored_resid'
    local resid_idx : word count `plugin_vars'
    scalar __civreghdfe_resid_idx = `resid_idx'

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
    capture noisily plugin call ctools_plugin `plugin_vars', "civreghdfe `threads_code' iv_regression"
    if _rc {
        local __rc = _rc
        * A variable named in orthog()/endogtest()/redundant() was dropped
        capture confirm scalar __civreghdfe_dropped_opt
        if _rc == 0 {
            local __o = scalar(__civreghdfe_dropped_opt)
            local __opt : word `__o' of orthog endogtest redundant
            local __lst : word `__o' of instruments_coef endogvars_coef instruments_coef
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
        _ctools_cleanup civreghdfe
        exit 0
    }

    * Retrieve results
    local N_used = __civreghdfe_N
    local df_r = __civreghdfe_df_r
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
        ereturn scalar df_r = `df_r'
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
    if `n_partial' > 0 ereturn local partial "`partial_terms'"

    * Set title and estimation method based on est_method
    if "`b0'" != "" {
        ereturn local title "CUE evaluated at user-supplied parameter vector"
        ereturn local model "cue"
    }
    else if `liml_kclass' {
        ereturn local title "LIML regression with HDFE"
        ereturn local model "liml"
        ereturn scalar kclass = `kclass_val'
    }
    else if `est_method' == 1 {
        ereturn local title "LIML regression with HDFE"
        ereturn local model "liml"
        capture ereturn scalar lambda = __civreghdfe_lambda
        capture ereturn scalar kclass = __civreghdfe_k
    }
    else if `est_method' == 2 {
        ereturn local title "Fuller LIML regression with HDFE"
        ereturn local model "fuller"
        ereturn scalar fuller = `fuller_val'
        capture ereturn scalar lambda = __civreghdfe_lambda
        capture ereturn scalar kclass = __civreghdfe_k
    }
    else if `est_method' == 3 {
        ereturn local title "k-class regression with HDFE"
        ereturn local model "kclass"
        ereturn scalar kclass = `kclass_val'
    }
    else if `est_method' == 4 {
        ereturn local title "2-step GMM with HDFE"
        ereturn local model "gmm2s"
    }
    else if `est_method' == 5 {
        ereturn local title "CUE regression with HDFE"
        ereturn local model "cue"
    }
    else {
        ereturn local title "IV (2SLS) regression with HDFE"
        ereturn local model "2sls"
    }

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
    ereturn scalar F = `F_wald'
    ereturn scalar Fdf1 = `model_df'
    ereturn scalar Fdf2 = `df_r'
    * ivreg2's p-value of F, centered R2 and variable counts and lists
    ereturn scalar Fp = Ftail(`model_df', `df_r', `F_wald')
    ereturn scalar r2c = `r2'
    ereturn scalar endog_ct = `K_endog'
    ereturn scalar exexog_ct = `K_inst_excl'
    ereturn scalar inexog_ct = `K_exog' - `has_cons' + `n_partial'
    ereturn local instd "`endogvars_coef'"
    ereturn local inexog "`: list exogvars_coef - __pc'"
    ereturn local exexog "`instruments_coef'"
    ereturn local insts "`: list exogvars_coef - __pc' `instruments_coef'"

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

    * Number of absorbed fixed effect dimensions (matches reghdfe/ivreghdfe)
    ereturn scalar N_hdfe = `G'

    * DOF adjustment options
    if `dofminus' != 0 {
        ereturn scalar dofminus = `dofminus'
    }
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
    local adj_denom = `N_used' - `K' - `sdofminus'
    if `adj_denom' > 0 {
        * With absorb, constant is partialled out, so use N not N-1
        * (N when the constant is partialled out with partial(), as ivreg2)
        local r2_a = 1 - (1 - `r2') * (`N_used' - (`G' == 0 & "`noconstant'" == "" & !`partial_cons')) / `adj_denom'
    }
    else {
        local r2_a = .
    }
    ereturn scalar r2_a = `r2_a'

    * Auxiliary regressions (first stage, reduced form) run by creghdfe over
    * the IV estimation sample, with the same weights and VCE, while the IV
    * results are held so they are restored afterwards.
    tempvar aux_sample
    if "`first'`rf'" != "" quietly gen byte `aux_sample' = e(sample)
    local aux_wt ""
    if "`weight'" != "" local aux_wt "[`weight'=`weight_var']"
    * The VCE as resolved above: pweights imply robust, as in the main fit
    local aux_vce = cond(`"`vce'"' == "", cond(`vce_type' == 1, "robust", "unadjusted"), `"`vce'"')
    local aux_absorb ""
    if "`absorb'" != "" local aux_absorb "absorb(`absorb_names')"

    * first: as in ivreg2, display each first-stage regression ahead of the
    * main results, and the first-stage summary statistics (ffirst)
    if "`first'" != "" local ffirst "ffirst"
    if "`first'" != "" & "`output'" == "" & `K_endog' > 0 {
        di as text ""
        di as text "First-stage regressions"
        di as text "{hline 23}"
        tempname first_hold
        _estimates hold `first_hold', restore
        foreach fs_dep of local endogvars_coef {
            quietly fvrevar `fs_dep'
            local fs_y `r(varlist)'
            capture quietly creghdfe `fs_y' `instruments' `exogvars' `aux_wt' ///
                if `aux_sample', `aux_absorb' vce(`aux_vce')
            local fs_rc = _rc
            di as text ""
            di as text "First-stage regression of `fs_dep':"
            if `fs_rc' {
                di as error "  first-stage regression failed (rc=`fs_rc')"
                continue
            }
            di as text ""
            if "`e(clustvar)'" != "" {
                di as text "Statistics robust to heteroskedasticity and clustering on `e(clustvar)'"
            }
            else if "`e(vcetype)'" != "" {
                di as text "Statistics robust to heteroskedasticity"
            }
            else {
                di as text "Statistics consistent for homoskedasticity only"
            }
            di as text "Number of obs = " as result %22.0f e(N)
            capture quietly testparm `instruments'
            local fs_test = !_rc
            if `fs_test' {
                local fs_F = r(F)
                local fs_df1 = r(df)
                local fs_df2 = r(df_r)
                local fs_p = r(p)
            }
            * The constant is partialled out with the fixed effects, so it is
            * not reported (creghdfe stores it as the last coefficient)
            tempname fs_b fs_V
            matrix `fs_b' = e(b)
            matrix `fs_V' = e(V)
            local fs_k = colsof(`fs_b')
            if "`absorb'" != "" & `fs_k' > 1 & colnumb(`fs_b', "_cons") == `fs_k' {
                matrix `fs_b' = `fs_b'[1, 1..`fs_k'-1]
                matrix `fs_V' = `fs_V'[1..`fs_k'-1, 1..`fs_k'-1]
            }
            local fs_N = e(N)
            local fs_dof = e(df_r)
            local fs_vcetype `e(vcetype)'
            ereturn post `fs_b' `fs_V', depname(`fs_dep') obs(`fs_N') dof(`fs_dof')
            if "`fs_vcetype'" != "" ereturn local vcetype "`fs_vcetype'"
            ereturn display, level(`level')
            if `fs_test' {
                di as text "F test of excluded instruments:"
                di as text "  F(" %3.0f `fs_df1' "," %6.0f `fs_df2' ") = " as result %8.2f `fs_F'
                di as text "  Prob > F      = " as result %8.4f `fs_p'
            }
        }
        _estimates unhold `first_hold'
    }

    * Display results in ivreghdfe format (if header not suppressed)
    if "`header'" == "" {
        di as text ""
        if "`title'" != "" {
            di as text "`title'"
        }
        else {
            di as text "IV (2SLS) estimation"
        }
        if "`subtitle'" != "" {
            di as text "`subtitle'"
        }
        di as text "{hline 20}"
        di as text ""

        * Efficiency/consistency notes (ivreg2's wording)
        local __hac1 = cond(`vce_type' != 0, "heteroskedasticity", "")
        local __hac3 ""
        if `kernel_type' > 0 & inlist(`vce_type', 0, 1) local __hac3 "autocorrelation"
        if `kiefer_val' local __hac3 "within-cluster autocorrelation (Kiefer)"
        if inlist(`vce_type', 2, 4) local __hac3 "clustering on `cluster_var'"
        if `vce_type' == 3 local __hac3 "clustering on `cluster_var' and `cluster_var2'"
        local __hacsep = cond("`__hac1'" != "" & "`__hac3'" != "", " and ", "")
        local __hac "`__hac1'`__hacsep'`__hac3'"
        if "`b0'" != "" {
            di as text "Estimates based on supplied parameter vector"
        }
        else if "`__hac'" != "" & inlist(`est_method', 4, 5) {
            di as text "Estimates efficient for arbitrary `__hac'"
        }
        else {
            di as text "Estimates efficient for homoskedasticity only"
        }
        if "`__hac'" != "" {
            di as text "Statistics robust to `__hac'"
        }
        else {
            di as text "Statistics consistent for homoskedasticity only"
        }
        if `kernel_type' > 0 & inlist(`vce_type', 3, 4) {
            di as text "and kernel-robust to common correlated disturbances (Driscoll-Kraay)"
        }
        if `kernel_type' > 0 {
            di as text "  kernel=`kernel_name'; bandwidth=`bw_val'"
            di as text "  time variable (t):  `ts_tvar'"
            if "`ts_ivar'" != "" di as text "  group variable (i): `ts_ivar'"
        }
        di as text ""

        * Summary statistics line 1
        if `vce_type' == 2 | `vce_type' == 4 {
            di as text "Number of clusters (`cluster_var') = " as result %9.0fc `N_clust' ///
                as text _col(55) "Number of obs =" as result %9.0fc `N_used'
        }
        else if `vce_type' == 3 {
            di as text "Number of clusters (`cluster_var') = " as result %6.0fc `N_clust1' ///
                as text _col(55) "Number of obs =" as result %9.0fc `N_used'
            di as text "Number of clusters (`cluster_var2') = " as result %6.0fc `N_clust2'
        }
        else {
            di as text _col(55) "Number of obs =" as result %9.0fc `N_used'
        }

        * F statistic (Wald test). Its numerator df is K, the regressors left
        * after collinearity removal, which is also the divisor of F_wald.
        * For two-way clustering, use min(N_clust1, N_clust2) - 1 for degrees of freedom
        local prob_F = Ftail(`model_df', `df_r', `F_wald')
        if `vce_type' == 2 {
            di as text _col(55) "F(" %3.0f `model_df' "," %6.0f (`N_clust' - 1) ") =" as result %9.2f `F_wald'
        }
        else if `vce_type' == 3 {
            local min_clust = min(`N_clust', `N_clust2')
            di as text _col(55) "F(" %3.0f `model_df' "," %6.0f (`min_clust' - 1) ") =" as result %9.2f `F_wald'
        }
        else {
            di as text _col(55) "F(" %3.0f `model_df' "," %6.0f `df_r' ") =" as result %9.2f `F_wald'
        }
        di as text _col(55) "Prob > F      =" as result %9.4f `prob_F'

        * R-squared, sums of squares and root MSE in ivreg2's formats
        di as text "Total (centered) SS     = " as result %12.0g `tss' ///
            as text _col(55) "Centered R2   =" as result %9.4f `r2'
        di as text "Total (uncentered) SS   = " as result %12.0g `tss' ///
            as text _col(55) "Uncentered R2 =" as result %9.4f `r2'
        di as text "Residual SS             = " as result %12.0g `rss' ///
            as text _col(55) "Root MSE      =" as result %9.4g `rmse'
        di as text ""
    }

    * Display coefficient table (if output not suppressed)
    if "`output'" == "" {
        * Build display options string
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

        * ivreg2 and ivreghdfe label the table with the dependent variable
        * itself, whatever depname() says; e(depvar) keeps the depname().
        ereturn local depvar "`depvar'"
        ereturn display, level(`level') `dispopt'
        ereturn local depvar "`depname_use'"
    }

    * Display full first-stage results (if ffirst specified)
    if "`ffirst'" != "" & `K_endog' > 0 {
        di as text ""
        di as text "{hline 78}"
        di as text "First-stage regression summary statistics"
        di as text "{hline 78}"
        * denominator df of the first-stage F (G-1 with clusters, as ivreg2)
        local fs_df_r = `df_r'
        capture confirm scalar __civreghdfe_F1df_r_1
        if _rc == 0 local fs_df_r = scalar(__civreghdfe_F1df_r_1)
        di as text _col(20) "Partial" _col(32) "F(" as result `K_inst_excl' as text "," ///
            as result `fs_df_r' as text ")" _col(52) "Prob > F"

        * One row per endogenous regressor the plugin kept, labelled with its
        * term (factor variables expanded); the statistics of the k-th kept
        * one are the plugin's k-th
        forvalues e = 1/`__n_fs' {
            local endogvar : word `e' of `__fs_terms'
            if `__fs_index_`e'' == 0 {
                ereturn scalar partial_r2_`e' = .
                continue
            }
            local __fsk = `__fs_index_`e''
            capture confirm scalar __civreghdfe_F1df_r_`__fsk'
            if _rc == 0 local fs_df_r = scalar(__civreghdfe_F1df_r_`__fsk')
            capture scalar `temp_F' = __civreghdfe_F1_`__fsk'
            capture scalar `temp_r2' = __civreghdfe_partial_r2_`__fsk'
            if _rc == 0 {
                local F_val = `temp_F'
                local r2_val = `temp_r2'
                local p_val = Ftail(`K_inst_excl', `fs_df_r', `F_val')
                di as text abbrev("`endogvar'", 15) _col(20) as result %8.4f `r2_val' ///
                    _col(32) as result %12.2f `F_val' _col(52) as result %10.4f `p_val'

                * Store extended first-stage results
                ereturn scalar partial_r2_`e' = `r2_val'
            }
        }
        di as text "{hline 78}"
        di as text ""
    }

    * Diagnostic tests: always stored in e(); nofooter suppresses only their
    * display (and that of the reduced form and the footer tables), so the
    * block runs quietly then
    local __footer_block = cond("`footer'" == "", "if 1", "quietly")
    `__footer_block' {

    * Retrieve underidentification test
    capture scalar `temp_underid' = __civreghdfe_underid
    capture scalar `temp_underid_df' = __civreghdfe_underid_df
    * Without endogenous regressors there are no identification statistics,
    * as in ivreg2
    if _rc == 0 & `K_endog' > 0 {
        local underid_stat = `temp_underid'
        local underid_df = `temp_underid_df'
        ereturn scalar idstat = `underid_stat'
        ereturn scalar iddf = `underid_df'
        ereturn scalar idp = chi2tail(`underid_df', `underid_stat')
    }

    * Underidentification test (if noid not specified)
    if "`id'" == "" & `K_endog' > 0 {
        if `vce_type' == 0 & `kernel_type' == 0 {
            di as text "Underidentification test (Anderson canon. corr. LM statistic):" ///
                _col(65) as result %14.3f `underid_stat'
        }
        else {
            di as text "Underidentification test (Kleibergen-Paap rk LM statistic):" ///
                _col(65) as result %14.3f `underid_stat'
        }
        local underid_p = chi2tail(`underid_df', `underid_stat')
        di as text _col(52) "Chi-sq(" as result `underid_df' as text ") P-val =" ///
            as result %10.4f `underid_p'
        di as text "{hline 78}"
    }

    * Weak identification test
    capture scalar `temp_cd_f' = __civreghdfe_cd_f
    capture scalar `temp_kp_f' = __civreghdfe_kp_f
    local cd_f_val = `temp_cd_f'
    local kp_f_val = `temp_kp_f'

    * Always display Cragg-Donald F (with endogenous regressors)
    if `K_endog' > 0 {
    di as text "Weak identification test (Cragg-Donald Wald F statistic):" ///
        _col(65) as result %14.3f `cd_f_val'

    * For robust/cluster/kernel VCEs, also display Kleibergen-Paap rk Wald F
    if (`vce_type' != 0 | `kernel_type' > 0) & `kp_f_val' > 0 {
        di as text "                         (Kleibergen-Paap rk Wald F statistic):" ///
            _col(65) as result %14.3f `kp_f_val'
    }

    * Stock-Yogo critical values (for 1 endogenous variable, varies by # excluded instruments)
    * Using values from Stock-Yogo (2005)
    if `K_endog' == 1 {
        if `K_inst_excl' == 1 {
            * 1 excluded instrument (just-identified)
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 16.38
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 8.96
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 6.66
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 5.53
        }
        else if `K_inst_excl' == 2 {
            * 2 excluded instruments
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 19.93
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 11.59
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 8.75
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 7.25
        }
        else if `K_inst_excl' == 3 {
            * 3 excluded instruments
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 22.30
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 12.83
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 9.54
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 7.80
        }
        else {
            * More than 3 excluded instruments - show approximate values
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 16.38
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 8.96
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 6.66
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 5.53
        }
    }
    else if `K_endog' == 2 {
        * Two endogenous - use appropriate critical values
        if `K_inst_excl' == 3 {
            di as text "Stock-Yogo weak ID test critical values: 10% maximal IV size" _col(65) as result %14.2f 13.43
            di as text "                                         15% maximal IV size" _col(65) as result %14.2f 8.18
            di as text "                                         20% maximal IV size" _col(65) as result %14.2f 6.40
            di as text "                                         25% maximal IV size" _col(65) as result %14.2f 5.45
        }
        else {
            di as text "(Stock-Yogo critical values: see Stock-Yogo (2005) for `K_endog' endogenous, `K_inst_excl' instruments)"
        }
    }
    else {
        * More than 2 endogenous - show message
        di as text "(Stock-Yogo critical values: see Stock-Yogo (2005) for `K_endog' endogenous variables)"
    }
    di as text "Source: Stock-Yogo (2005).  Reproduced by permission."
    if `vce_type' != 0 | `kernel_type' > 0 {
        di as text "NB: Critical values are for Cragg-Donald F statistic and i.i.d. errors."
    }
    di as text "{hline 78}"
    }

    * Display reduced form if requested; a failed reduced form is reported
    * and the remaining statistics are still stored
    local rf_rc = 0
    if "`rf'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "Reduced-form regression: `depvar' on instruments"
        di as text "{hline 60}"

        * Save current IV results
        tempname iv_b iv_V
        matrix `iv_b' = e(b)
        matrix `iv_V' = e(V)
        local iv_N = e(N)
        local iv_df_r = e(df_r)
        local iv_df_a = e(df_a)
        local iv_K = e(K)
        local iv_K_endog = e(K_endog)
        local iv_K_exog = e(K_exog)
        local iv_K_iv = e(K_iv)
        local iv_G = e(G)
        local iv_title = e(title)
        local iv_model = e(model)
        local iv_vcetype = e(vcetype)
        local iv_cmd = e(cmd)

        * Run reduced form regression using creghdfe
        * Reduced form: y = [exogenous + excluded instruments] + error
        tempname rf_hold
        _estimates hold `rf_hold', restore
        capture qui creghdfe `depvar' `exogvars' `instruments' `aux_wt' ///
            if `aux_sample', `aux_absorb' vce(`aux_vce')
        local rf_rc = _rc
    }
    if "`rf'" != "" & `rf_rc' {
        di as error "  reduced-form regression failed (rc=`rf_rc')"
        _estimates unhold `rf_hold'
    }
    else if "`rf'" != "" {

        * Display coefficient table
        tempname rf_b rf_V
        matrix `rf_b' = e(b)
        matrix `rf_V' = e(V)

        * Display reduced form coefficients in a table
        local rf_vars : colnames `rf_b'
        local rf_K : word count `rf_vars'
        di as text ""
        di as text "{hline 78}"
        di as text _col(14) "|" _col(28) "`iv_vcetype'"
        di as text %12s abbrev("`depvar'", 12) " | Coefficient  std. err.      t    P>|t|" ///
            "     [95% conf. interval]"
        di as text "{hline 13}+{hline 64}"

        forvalues v = 1/`rf_K' {
            local vname : word `v' of `rf_vars'
            local b = `rf_b'[1, `v']
            local se = sqrt(`rf_V'[`v', `v'])
            if `se' > 0 {
                local t = `b' / `se'
                local p = 2 * ttail(`iv_df_r', abs(`t'))
                local ci_lo = `b' - invttail(`iv_df_r', 0.025) * `se'
                local ci_hi = `b' + invttail(`iv_df_r', 0.025) * `se'
                di as text %12s abbrev("`vname'", 12) " |" ///
                    as result %11.6g `b' "  " %9.6g `se' "  " ///
                    %7.2f `t' "   " %5.3f `p' "    " ///
                    %11.6g `ci_lo' "  " %11.6g `ci_hi'
            }
            else {
                di as text %12s abbrev("`vname'", 12) " |" ///
                    as result %11.6g `b' "  " "(omitted)"
            }
        }
        di as text "{hline 78}"

        * Restore the IV results
        _estimates unhold `rf_hold'
    }

    * Store and display diagnostic statistics
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

        * Store in e()
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
        if `K_endog' > 0 ereturn scalar cdf = `cd_f_val'

        * Store weak ID statistics (widstat = Kleibergen-Paap F if robust, else Cragg-Donald F)
        capture scalar `temp_kp_f2' = __civreghdfe_kp_f
        if `K_endog' == 0 {
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

        * Store endogeneity test
        ereturn scalar endog_chi2 = `endog_chi2'
        ereturn scalar endog_f = `endog_f'
        ereturn scalar endog_df = `endog_df'
        if `endog_df' > 0 {
            ereturn scalar endog_p = chi2tail(`endog_df', `endog_chi2')
        }

        * Display Sargan/Hansen J statistic (not reported, as in ivreg2, when
        * the moment covariance matrix is rank deficient)
        if missing(`sargan_val') {
            di as error "Warning: estimated covariance matrix of moment conditions not of full rank."
            di as error "         overidentification statistic not reported, and standard errors and"
            di as error "         model tests should be interpreted with caution."
            di as error "Possible causes:"
            di as error "         number of clusters insufficient to calculate robust covariance matrix"
            di as error "         singleton dummy variable (dummy with one 1 and N-1 0s or vice versa)"
        }
        else if `vce_type' == 0 {
            di as text "Sargan statistic (overidentification test of all instruments):" ///
                _col(65) as result %14.3f `sargan_val'
        }
        else {
            di as text "Hansen J statistic (overidentification test of all instruments):" ///
                _col(65) as result %14.3f `sargan_val'
        }
        if !missing(`sargan_val') & `sargan_df' > 0 {
            local sargan_p = chi2tail(`sargan_df', `sargan_val')
            di as text _col(52) "Chi-sq(" as result `sargan_df' as text ") P-val =" ///
                as result %10.4f `sargan_p'
        }
        else if !missing(`sargan_val') {
            di as text _col(50) "(equation exactly identified)"
        }
        di as text "{hline 78}"
    }

    * Display orthog test (C-statistic) if requested
    if `n_orthog' > 0 {
        capture scalar `temp_cstat' = __civreghdfe_cstat
        capture scalar `temp_cstat_df' = __civreghdfe_cstat_df
        if _rc == 0 {
            local cstat_val = `temp_cstat'
            local cstat_df = `temp_cstat_df'
            local cstat_p = chi2tail(`cstat_df', `cstat_val')

            ereturn scalar cstat = `cstat_val'
            ereturn scalar cstat_df = `cstat_df'
            ereturn scalar cstat_p = `cstat_p'

            di as text "C statistic (exogeneity/orthogonality of instruments):" ///
                _col(65) as result %14.3f `cstat_val'
            di as text _col(52) "Chi-sq(" as result `cstat_df' as text ") P-val =" ///
                as result %10.4f `cstat_p'
            di as text "{hline 78}"
        }
    }

    * Display endogtest result if specified
    if `n_endogtest' > 0 {
        capture scalar `temp_endogtest_stat' = __civreghdfe_endogtest_stat
        capture scalar `temp_endogtest_df' = __civreghdfe_endogtest_df
        if _rc == 0 {
            local endogtest_val = `temp_endogtest_stat'
            local endogtest_df_val = `temp_endogtest_df'
            local endogtest_p = chi2tail(`endogtest_df_val', `endogtest_val')

            ereturn scalar endogtest = `endogtest_val'
            ereturn scalar endogtest_df = `endogtest_df_val'
            ereturn scalar endogtest_p = `endogtest_p'

            di as text "Endogeneity test of endogenous regressors (`endogtest'):" ///
                _col(65) as result %14.3f `endogtest_val'
            di as text _col(52) "Chi-sq(" as result `endogtest_df_val' as text ") P-val =" ///
                as result %10.4f `endogtest_p'
            di as text "{hline 78}"
        }
    }

    * Display redundant test result if specified
    if `n_redundant' > 0 {
        capture scalar `temp_redund_stat' = __civreghdfe_redund_stat
        capture scalar `temp_redund_df' = __civreghdfe_redund_df
        if _rc == 0 {
            local redund_val = `temp_redund_stat'
            local redund_df_val = `temp_redund_df'
            local redund_p = chi2tail(`redund_df_val', `redund_val')

            ereturn scalar redund = `redund_val'
            ereturn scalar redund_df = `redund_df_val'
            ereturn scalar redund_p = `redund_p'

            di as text "Redundancy test of instruments (`redundant'):" ///
                _col(65) as result %14.3f `redund_val'
            di as text _col(52) "Chi-sq(" as result `redund_df_val' as text ") P-val =" ///
                as result %10.4f `redund_p'
            di as text "{hline 78}"
        }
    }

    * Instruments dropped after the fixed effects were partialled out, which
    * ivreghdfe reports in its footer (an included exogenous regressor is an
    * instrument too)
    local __inst_dropped = 0
    capture local __inst_dropped = scalar(__civreghdfe_z_dropped)
    if _rc local __inst_dropped = 0
    forvalues __k = 1/`K_exog' {
        local __kk = `K_endog' + `__k'
        local __inst_dropped = `__inst_dropped' + (`: word `__kk' of `collinear_flags'' == 1)
    }
    if `__inst_dropped' > 0 {
        di as text "Collinearities detected among instruments: `__inst_dropped' instrument(s) dropped"
    }

    * Display instruments footer (without the terms ivreg2's collinearity
    * rules dropped, which are listed below)
    foreach __l in endogvars exogvars instruments {
        local __disp_`__l'
        foreach __t of local `__l' {
            if !regexm("`__t'", "o\.") local __disp_`__l' `__disp_`__l'' `__t'
        }
    }
    di as text "Instrumented:         `__disp_endogvars'"
    if "`__disp_exogvars'" != "" {
        di as text "Included instruments: `__disp_exogvars'"
    }
    di as text "Excluded instruments: `__disp_instruments'"
    di as text "Partialled-out:       _cons"
    di as text "                      nb: total SS, model F and R2s are after partialling-out;"
    di as text "                          any small-sample adjustments include partialled-out"
    di as text "                          variables in regressor count K"
    if "`__dups'" != "" di as text "Duplicates:" _col(23) "`__dups'"
    if "`__collin'" != "" di as text "Dropped collinear:" _col(23) "`__collin'"
    if "`__ecollin'" != "" di as text "Reclassified as exog:" _col(23) "`__ecollin'"
    di as text "{hline 78}"

    * Display absorbed degrees of freedom table (only when absorb was specified)
    if `G' > 0 {
    di as text ""
    di as text "Absorbed degrees of freedom:"
    di as text "{hline 13}{c TT}{hline 39}{c TRC}"
    di as text " Absorbed FE {c |} Categories  - Redundant  = Num. Coefs {c |}"
    di as text "{hline 13}{c +}{hline 39}{c RT}"

    * Display each FE with its components: the redundant coefficients the
    * plugin counted as reghdfe does (nested FEs: all levels are redundant)
    local total_coefs = 0
    local has_nested = 0
    forval g = 1/`G' {
        local fevar : word `g' of `absorb_names'
        local cats = `fe_levels_`g''
        capture local redundant = scalar(__civreghdfe_fe_redundant_`g')
        if _rc local redundant = 0
        if `fe_nested_`g'' == 1 {
            local redundant = `cats'
            local has_nested = 1
            local nested_marker "*"
        }
        else {
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
    di as text ""
    } /* end if G > 0 */

    } /* End of the diagnostic tests (quietly under nofooter) */

    * savefirst/saverf: as ivreg2, store each first-stage regression as
    * <savefprefix><endogenous variable> and the reduced form as
    * <saverfprefix><depvar>: creghdfe over the IV estimation sample with the
    * same weights and VCE, the constant partialled out with the effects.
    * The IV results are held and restored afterwards.
    local save_eqs ""
    if "`savefirst'" != "" {
        local save_prefix = cond("`savefprefix'" != "", "`savefprefix'", "_civreghdfe_")
        foreach save_dep of local endogvars_coef {
            local save_eqs `save_eqs' first `save_dep'
        }
    }
    if "`saverf'" != "" {
        local save_rfprefix = cond("`saverfprefix'" != "", "`saverfprefix'", "_civreghdfe_rf_")
        local save_eqs `save_eqs' rf `depvar'
    }
    if "`save_eqs'" != "" {
        tempvar save_sample
        quietly gen byte `save_sample' = e(sample)
        tempname save_hold save_b save_V
        _estimates hold `save_hold', restore
        while "`save_eqs'" != "" {
            gettoken save_stage save_eqs : save_eqs
            gettoken save_dep save_eqs : save_eqs
            quietly fvrevar `save_dep'
            local save_y `r(varlist)'
            capture quietly creghdfe `save_y' `instruments' `exogvars' `aux_wt' ///
                if `save_sample', `aux_absorb' vce(`aux_vce')
            if _rc {
                di as error "`save_stage' regression of `save_dep' failed (rc=`=_rc'); not stored"
                continue
            }
            matrix `save_b' = e(b)
            matrix `save_V' = e(V)
            * the constant is partialled out with the effects (without
            * absorb() it is a coefficient, as in ivreg2)
            local save_k = colsof(`save_b')
            if "`absorb'" != "" & `save_k' > 1 & colnumb(`save_b', "_cons") == `save_k' {
                matrix `save_b' = `save_b'[1, 1..`save_k'-1]
                matrix `save_V' = `save_V'[1..`save_k'-1, 1..`save_k'-1]
            }
            local save_N = e(N)
            local save_dof = e(df_r)
            local save_rmse = e(rmse)
            local save_vcetype `e(vcetype)'
            local save_clust = e(N_clust)
            local save_clustvar `e(clustvar)'
            * As ivreg2 stores them, a clustered first stage keeps the
            * residual dof N - K - sdofminus (sdofminus: the absorbed
            * effects, at least 1), and the root MSE with it, not the
            * number of clusters less one; its small-sample factor uses the
            * same dof, where reghdfe adds 1 when an effect is nested in the
            * clusters
            if !missing(`save_clust') & "`save_clustvar'" != "" {
                local __sdf = 0
                if "`absorb'" != "" local __sdf = cond(missing(e(df_a)) | e(df_a) < 1, 1, e(df_a))
                local save_dof = `save_N' - colsof(`save_b') - `__sdf'
                local save_rmse = sqrt(e(rss) / `save_dof')
                local __nadj = 0
                capture confirm matrix e(dof_table)
                if !_rc {
                    tempname __dt
                    matrix `__dt' = e(dof_table)
                    forvalues __r = 1/`=rowsof(`__dt')' {
                        if `__dt'[`__r', 5] == 1 local __nadj = 1
                    }
                }
                local __dfa = cond(missing(e(df_a)), 0, e(df_a))
                matrix `save_V' = `save_V' * ((`save_N' - e(df_m) - `__dfa' - `__nadj') / `save_dof')
            }
            * no e(sample), as in ivreg2's stored first-stage/reduced-form results
            ereturn post `save_b' `save_V', depname(`save_dep') obs(`save_N') dof(`save_dof')
            if "`save_vcetype'" != "" ereturn local vcetype "`save_vcetype'"
            if !missing(`save_clust') & "`save_clustvar'" != "" {
                ereturn scalar N_clust = `save_clust'
                ereturn local clustvar "`save_clustvar'"
                ereturn local vce "robust cluster"
            }
            ereturn scalar rmse = `save_rmse'
            ereturn scalar df_m = colsof(e(b))
            ereturn local model "`save_stage'"
            ereturn local predict "civreghdfe_p"
            ereturn local cmd "civreghdfe"
            if "`save_stage'" == "first" {
                ereturn local estimates_title "First-stage regression: `save_dep'"
                local save_name = strtoname("`save_prefix'`save_dep'")
            }
            else {
                ereturn local estimates_title "Reduced-form regression: `save_dep'"
                local save_name = strtoname("`save_rfprefix'`save_dep'")
            }
            estimates store `save_name', title(`e(estimates_title)')
        }
        _estimates unhold `save_hold'
    }

    * Save first-stage estimation results (if savefirst specified)
    if "`savefirst'" != "" {
        local prefix = cond("`savefprefix'" != "", "`savefprefix'", "_civreghdfe_")
        * Store current estimation under each endogenous variable's name
        * Note: Full first-stage regressions would require running OLS separately
        estimates store `prefix'main
        di as text ""
        di as text "Estimation stored as: `prefix'main"
        di as text "First-stage F-statistics available in e(F_first1), e(F_first2), etc."
    }

    * Save reduced-form estimation results (if saverf specified)
    if "`saverf'" != "" {
        local rfprefix = cond("`saverfprefix'" != "", "`saverfprefix'", "_civreghdfe_rf_")
        estimates store `rfprefix'main
        di as text ""
        di as text "Estimation stored as: `rfprefix'main"
    }

    * Capture post-processing time (before cleanup)
    _ctools_lap t_stata_post __ct_lap98 `__ct_lap98'

    * Display timing breakdown if verbose specified
    if "`verbose'" != "" {
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
