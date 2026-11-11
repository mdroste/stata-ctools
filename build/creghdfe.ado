*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define creghdfe, eclass
    version 14.1
    * Replay never enters the estimation transaction or clears existing e().
    if replay() {
        if "`e(cmd)'" != "creghdfe" error 301
        _creghdfe_display `0'
        exit
    }
    * A failed or interrupted regression leaves the data as they were without
    * a snapshot of the dataset: creghdfe only adds variables at the end
    * (residuals, groupvar(), savefe and temporary variables) and writes into
    * nothing else, so _ctools_undo drops them and resets c(changed). Option
    * snapshot preserves the whole dataset instead (on disk above
    * c(max_preservemem)); option fast undoes nothing.
    _ctools_mode `0'
    * Break is recognized only inside _creghdfe_impl.
    nobreak {
        if "`mode'" == "snapshot" preserve
        else if "`mode'" == "" _ctools_undo record CTOOLS_u_creghdfe_
        capture noisily break _creghdfe_impl `0'
        local rc = _rc
        if `rc' {
            if "`mode'" == "snapshot" restore
            else if "`mode'" == "" _ctools_undo rollback CTOOLS_u_creghdfe_
            ereturn clear
            * nothing of the failed run stays: the retained fit, the
            * bootstrap's scalars and matrices, the draws' variables (option
            * fast undoes nothing) and the wild bootstrap's Mata draws
            capture _creghdfe_bs_release
            capture _ctools_cleanup creghdfe
            capture _ctools_hdfe_boot wildclear
            exit `rc'
        }
        if "`mode'" == "snapshot" restore, not
        else if "`mode'" == "" _ctools_undo clear CTOOLS_u_creghdfe_
    }
end

program define _creghdfe_impl, eclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    * Start wall clock timer immediately
    local __do_timing = 0
    quietly _ctools_clock __ct_lap98

    * Drop whatever an interrupted run left behind
    _ctools_cleanup creghdfe

    * A refit on the singleton-free sample (see after the plugin call) passes
    * that sample and the first run's singleton count here
    local __refit_sample : word 1 of $CTOOLS_creghdfe_refit_
    local __refit_singletons : word 2 of $CTOOLS_creghdfe_refit_
    global CTOOLS_creghdfe_refit_
    local __creghdfe_args `"`0'"'

    * reghdfe's options with its abbreviations; its solver and memory
    * options are parsed as reghdfe parses them (see _creghdfe_solver_opts)
    syntax varlist(min=1 fv ts) [aw fw pw] [if] [in] [, Absorb(string) NOABsorb VCE(string) FEVCE(string) ///
        CLuster(string) Verbose ///
        TOLerance(real 1e-8) ITERATE(integer 16000) NOSTANDardize ///
        DOFadjustments(string) GROUPVar(name) THReads(integer 0) QUAD SNAPshot FAST ///
        noCONstant KEEPSINgletons noHEader noTABle noFOOTnote NOWILDtable ///
        NOSAMPle TIMEit noWARN VERSION(integer -1) ///
        TEChnique(string) TRAnsform(string) ACCELeration(string) PREConditioner(string) ///
        PRUNE POOLsize(integer 10) COMPACT PARallel(string asis) FASTREGress ///
        noPARTIALout noREGress KEEPMATA varlist_is_touse ///
        Group_id(string) Individual_id(string) AGgregation(string) *]

    * Freeze incoming-r() qualifiers before parsing helpers replace r().
    if `"`if'"' != "" {
        tempvar __incoming_sample
        mark `__incoming_sample' `if' `in'
        local if "if `__incoming_sample'"
    }

    * version(3|5) runs a historical reghdfe in reghdfe; creghdfe implements
    * reghdfe 6 only
    if `version' != -1 {
        if !inlist(`version', 3, 5) {
            di as error "assert failed: inlist(`version', 3, 5)"
            exit 9
        }
        di as error "creghdfe: version(`version') runs reghdfe `version', a historical version; creghdfe implements reghdfe 6 only"
        exit 198
    }
    * reghdfe's programming options hand its Mata FixedEffects object to
    * other commands; creghdfe, a compiled plugin, has none to hand over
    foreach __o in `partialout' `regress' `keepmata' `varlist_is_touse' {
        di as error "creghdfe: option `__o' is a reghdfe programming option that keeps reghdfe's Mata HDFE object; creghdfe has no such object"
        exit 198
    }

    * residuals(newvar) (also resid2(newvar)), or residuals without a name
    * for _reghdfe_resid, which is replaced, as in reghdfe
    _creghdfe_resid_opt, `options'
    * As in reghdfe: verbose(#) (a positive level gives the timing
    * breakdown), and the display options of estimation commands
    * wildbootstrap[(wildopts)]: the wild bootstrap of the coefficients
    _creghdfe_wild_opt, `__resid_rest'
    _creghdfe_verbose_opt, `__wild_rest'
    local verbose = cond("`verbose'" != "" | `__verbose_level' > 0, "verbose", "")
    _get_diopts diopts __undoc, `__options'
    * reghdfe sets any remaining key(value) on its Mata object; creghdfe
    * takes those its solver has: abort(0|1), and accepts miniter(),
    * min_ok(), accel_start() and accel_freq() (no effect, see the help)
    _creghdfe_undoc_opts `__undoc'
    local __abort = `s(abort)'
    * technique(), transform(), acceleration(), preconditioner() and
    * parallel() checked as reghdfe checks them
    _creghdfe_solver_opts, technique(`technique') transform(`transform') ///
        acceleration(`acceleration') preconditioner(`preconditioner') parallel(`parallel')
    local __technique `s(technique)'
    local __transform_code = `s(transform_code)'
    local __accel_code = `s(accel_code)'
    local __kcg = `s(kcg)'
    local __timeit = ("`timeit'" != "")
    * reghdfe's group() and individual() (groupvar() is spelled out, as in
    * reghdfe): the regression has one row per group (a patent), and with
    * individual() an absorbed individual effect (its inventors) enters each
    * group's row as the mean, or with aggregation(sum) the sum, of the
    * effects of the group's individuals
    _creghdfe_teams_opt, group_id(`group_id') individual_id(`individual_id') aggregation(`aggregation')
    local group_id `s(group_id)'
    local individual_id `s(individual_id)'
    local __teams = ("`group_id'" != "")
    local __indiv = ("`individual_id'" != "")
    local __indiv_mean = ("`s(aggregation)'" == "mean")
    * reghdfe's verbose(-1) suppresses its notes and warnings
    local __quiet = (`__verbose_level' < 0)

    * cluster(varlist) is vce(cluster varlist); vce(conventional) is the
    * unadjusted VCE; noabsorb is the default without absorb()
    if `"`cluster'"' != "" {
        if `"`vce'"' != "" {
            di as error "only one of cluster() and vce() can be specified"
            exit 198
        }
        local vce `"cluster `cluster'"'
    }
    * vce() types as reghdfe reads them (case matters): ols* and un* are
    * unadjusted, as is conventional; r* is robust, cl* cluster, dk*
    * dkraay; vce(bootstrap) is creghdfe's own
    if `"`vce'"' != "" {
        gettoken __vw __vrest : vce, parse(" ,(")
        if lower(substr("`__vw'", 1, 4)) != "boot" {
            if substr("`__vw'", 1, 3) == "ols" | substr("`__vw'", 1, 2) == "un" | "`__vw'" == "conventional" local __vw unadjusted
            else if substr("`__vw'", 1, 1) == "r" local __vw robust
            else if substr("`__vw'", 1, 2) == "cl" local __vw cluster
            else if substr("`__vw'", 1, 2) == "dk" local __vw dkraay
            if !inlist("`__vw'", "unadjusted", "robust", "cluster", "dkraay") {
                di as error "vcetype '`__vw'' not allowed"
                exit 9
            }
            local vce `"`__vw'`__vrest'"'
        }
    }
    if "`keepsingletons'" != "" & !`__quiet' & "`warn'" != "nowarn" {
        di as error `"WARNING: Singleton observations not dropped; statistical significance is biased {browse "http://scorreia.com/reghdfe/nested_within_cluster.pdf":(link)}"'
    }

    if `iterate' < 1 | missing(`tolerance') | `tolerance' <= 0 {
        di as error "creghdfe: iterate() and tolerance() must be positive"
        exit 198
    }
    local __do_timing = ("`verbose'" != "")

    * If absorb is not specified, create a constant FE (absorb the intercept only)
    local __noabsorb = 0
    if `"`absorb'"' == "" {
        local __noabsorb = 1
        tempvar _cons_fe
        quietly gen byte `_cons_fe' = 1
        local absorb "`_cons_fe'"
    }

    * Parse absorb option for suboptions (savefe)
    local savefe = 0
    local absorb_vars "`absorb'"
    local absorb_noisily = 0
    local absorb_replace = 0
    local absorb_badopts
    if strpos(`"`absorb'"', ",") > 0 {
        gettoken absorb_vars absorb_opts : absorb, parse(",")
        _creghdfe_absorb_opts `absorb_opts'
    }

    * Parse dofadjustments option as reghdfe does
    * 0=all (default), 1=none, 2=firstpair, 3=pairwise; all implies clusters
    local dof_lower = lower(`"`dofadjustments'"')
    _ctools_hdfe_dof, `dof_lower'
    local dof_adjust_type = `s(type)'
    local dof_clusters = `s(clusters)'
    local dofmethods `s(methods)'

    * Mark sample - include absorb and cluster variables
    marksample touse
    if "`__refit_sample'" != "" quietly replace `touse' = 0 if !`__refit_sample'

    * One variable per absorbed effect, in order: the FE count and every
    * plugin variable position are computed from this list. A numeric
    * variable goes to the plugin as it is; string variables and
    * interactions are numbered first (the plugin renumbers every effect on
    * the estimation sample). absorb_names holds the effects as reghdfe
    * names them, for e(absorb) and the output.
    capture noisily _ctools_hdfe_absorb `absorb_vars', touse(`touse') cmd(creghdfe) slopes ///
        `=cond(`absorb_replace', "replace", "")'
    if _rc {
        local __rc = _rc
        if `__rc' == 110 di as error "error in option absorb()"
        exit `__rc'
    }
    * reghdfe checks the suboptions of absorb() after its variables
    if `"`absorb_badopts'"' != "" {
        di as error `"invalid suboptions in absorb(): `absorb_badopts'"'
        exit 198
    }
    local nfe = r(n)
    forvalues g = 1/`nfe' {
        local absorb_name_`g' `"`r(name_`g')'"'
        local absorb_vars_`g' `r(vars_`g')'
        local absorb_group_`g' = r(group_`g')
        local absorb_int_`g' = r(intercept_`g')
        local absorb_cvars_`g' `r(cvars_`g')'
        local absorb_target_`g' `r(target_`g')'
    }
    * individual(): its variable is one of the absorbed effects, which the
    * plugin absorbs as the individual effect rather than as a factor;
    * __indiv_g is that effect and __indiv_pos the factors before it (its
    * row in e(dof_table))
    local __indiv_g = 0
    local __indiv_pos = 0
    if `__indiv' {
        forvalues g = 1/`nfe' {
            if "`absorb_vars_`g''" == "`individual_id'" local __indiv_g = `g'
        }
        if !`__indiv_g' | `__noabsorb' {
            di as error "individual id not found on absorb()"
            exit 100
        }
        forvalues g = 1/`nfe' {
            if "`absorb_cvars_`g''" != "" {
                di as error "creghdfe: heterogeneous slopes are not available with individual()"
                exit 198
            }
        }
        if `savefe' | "`absorb_target_`__indiv_g''" != "" {
            di as error "creghdfe: the individual fixed effects cannot be saved (savefe or newvar=), as in reghdfe"
            exit 198
        }
    }
    * Heterogeneous slopes and newvar= targets: the plugin partials out one
    * factor per intercept and per slope (the extended absorbed effects of
    * reghdfe, e(extended_absvars)), a slope centered within its levels when
    * its effect has an intercept. ext_* describe these factors in order.
    local absorb
    local absorb_names
    local next = 0
    local any_slopes = 0
    local any_target = 0
    local any_intercept = 0
    local ext_names
    local ext_outvars
    forvalues g = 1/`nfe' {
        local absorb_names `"`absorb_names' `absorb_name_`g''"'
        if `g' == `__indiv_g' {
            local __indiv_pos = `next'
            local __indiv_name `"`absorb_vars_`g''"'
            continue
        }
        if `absorb_group_`g'' {
            tempvar absorb_id_`g'
            quietly _ctools_group `absorb_vars_`g'' if `touse', generate(`absorb_id_`g'')
            local fe_var `absorb_id_`g''
        }
        else local fe_var `absorb_vars_`g''
        local ivlabel : subinstr local absorb_vars_`g' " " "#", all
        local stub `absorb_target_`g''
        if "`stub'" != "" local any_target = 1
        else if `savefe' local stub __hdfe`g'__
        if `absorb_int_`g'' {
            local ++next
            local any_intercept = 1
            local absorb `absorb' `fe_var'
            local ext_slope_`next'
            local ext_absvar_`next' `g'
            local ext_names `ext_names' `ivlabel'
            if "`absorb_cvars_`g''" == "" local ext_label_`next' `"`absorb_name_`g''"'
            else local ext_label_`next' `ivlabel'
            local ext_out_`next' `stub'
        }
        local h = 0
        foreach cv of local absorb_cvars_`g' {
            local ++h
            local ++next
            local any_slopes = 1
            local absorb `absorb' `fe_var'
            local ext_slope_`next' `cv'
            local ext_absvar_`next' `g'
            local ext_names `ext_names' `ivlabel'#c.`cv'
            local ext_label_`next' `ivlabel'#c.`cv'
            local ext_out_`next' = cond("`stub'" != "", "`stub'Slope`h'", "")
        }
    }
    local absorb_names : list retokenize absorb_names
    * absorb(..., noisily) shows the absorbed effects as parsed, as reghdfe
    if `absorb_noisily' & !`__noabsorb' di as result "{bf:`absorb_names'}"
    * absorb(..., replace): a newvar= target that exists is replaced (its
    * slope variables, newvarSlope#, must still be new), as in reghdfe
    local __replace_targets
    if `absorb_replace' {
        forvalues g = 1/`nfe' {
            if "`absorb_target_`g''" != "" local __replace_targets `__replace_targets' `absorb_target_`g''
        }
    }
    * The individual effect counts as an absorbed intercept: reghdfe reports
    * the constant with it, whatever the aggregation
    if `__indiv' local any_intercept = 1
    * Slopes alone (absorb(id#c.x)) absorb no constant: as reghdfe, the
    * model is then fitted without one (e(report_constant) is 0)
    if !`any_intercept' & `: word count `varlist'' < 2 {
        di as error "creghdfe: a model without an absorbed intercept needs a regressor"
        exit 198
    }
    if `next' > 10 {
        di as error "creghdfe: at most ten absorbed intercepts and slopes are supported"
        exit 198
    }
    * Saved fixed effects: savefe, or a newvar= target on any effect. Every
    * factor has an output variable for the plugin (a temporary one where
    * nothing is saved).
    if `any_target' local savefe = 1
    if `savefe' {
        forvalues g = 1/`next' {
            if "`ext_out_`g''" == "" {
                tempvar ext_tmp_`g'
                local ext_out_`g' `ext_tmp_`g''
            }
            else local ext_saved `ext_saved' `ext_out_`g''
            local ext_outvars `ext_outvars' `ext_out_`g''
        }
    }

    * Parse weights
    local weight_var ""
    local weight_type = 0   // 0=none, 1=aweight, 2=fweight, 3=pweight
    if "`weight'" != "" {
        * Evaluate the weight expression; marks out and validates the weights
        tempvar evaluated_weight
        _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"`exp'"')
        local weight_var "`evaluated_weight'"
        local weight_type = cond("`weight'" == "aweight", 1, cond("`weight'" == "fweight", 2, 3))
    }

    * vce(bootstrap): a pairs bootstrap of the coefficients, run by the
    * plugin on the retained fit (docs/development/BOOTSTRAP_PLAN.md, 7). The
    * fit itself uses the unadjusted VCE, posted as e(V_modelbased). Its
    * cluster and strata variables are marked out here, as the vce() cluster.
    local __bs = 0
    local bs_unitvar
    local bs_stratavar
    local bs_save_file
    local bs_bsdraws
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
            local bs_designopts = r(designopts)
            local bs_bsdraws `r(bsdraws)'
            local bs_save_file `"`r(save_file)'"'
            local bs_save_double `r(save_double)'
            local bs_save_every `r(save_every)'
            local bs_save_replace `r(save_replace)'
            local vce
        }
    }
    * The command's level() sets vce(bootstrap)'s, as for bootstrap; the
    * table at estimation shows the bootstrap's level
    if `__bs' {
        _ctools_hdfe_boot level, cmd(creghdfe) bslevel(`bs_level') bsgiven(`bs_lgiven') diopts(`"`diopts'"')
        local bs_level = r(level)
        local diopts `"`r(diopts)'"'
    }
    * fevce(bootstrap): bootstrap standard errors of the saved fixed effects
    * (plan 4.3, 5.7). With vce(bootstrap) it shares its replicates; options
    * given to both must agree.
    local __fe = 0
    local fe_nreps = 0
    local fe_ci = 0
    local fe_dse
    local fe_sestub
    if `"`fevce'"' != "" {
        _ctools_hdfe_boot parse, spec(`"`fevce'"') fe
        if "`r(kind)'" != "bootstrap" {
            di as error "creghdfe: fevce() takes bootstrap[, options]: the fixed effects use the pairs bootstrap"
            exit 198
        }
        if !`savefe' {
            di as error "creghdfe: fevce(bootstrap) needs saved fixed effects: absorb(..., savefe)"
            exit 198
        }
        if `__bs' & (`"`r(saving)'"' != "" | "`r(bsdraws)'" != "") {
            di as error "creghdfe: with vce(bootstrap), specify saving() and _bsdraws() there"
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
                di as error "creghdfe: with vce(bootstrap), fevce(bootstrap) uses the level(), mse and dots() of vce(); specify them once, in vce()"
                exit 198
            }
        }
        if `__bs' & r(designopts) {
            if r(reps) != `bs_reps' | `"`r(seed)'"' != `"`bs_seed'"' | `"`r(cluster)'"' != `"`bs_cluster'"' | ///
                `"`r(strata)'"' != `"`bs_strata'"' | "`r(draws)'" != "`bs_draws'" {
                di as error "creghdfe: vce(bootstrap) and fevce(bootstrap) share one set of replicates; specify the resampling options once"
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
    * fevce() names its outputs after savefe's
    if `__fe' & `any_target' {
        di as error "creghdfe: fevce(bootstrap) is not supported with newvar= targets in absorb(); use absorb(..., savefe)"
        exit 198
    }
    * wildbootstrap(): p-values and intervals of the coefficients by the wild
    * (cluster) bootstrap, after the analytic fit (plan 4.7, 5.9, 7.5)
    if `__wild' {
        if `__bs' {
            di as error "creghdfe: wildbootstrap() and vce(bootstrap) are two bootstraps of the coefficients; choose one"
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
    if `__bsrun' {
        if `"`bs_cluster'"' != "" {
            _ctools_hdfe_vce, spec(`"cluster `bs_cluster'"') touse(`touse') cmd(creghdfe)
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

    * Parse VCE option. Its cluster variables are marked out before the
    * factor variables are expanded, so that their levels come from the
    * estimation sample.
    local vcetype = 0   // 0=unadjusted, 1=robust, 2=cluster, 3=dkraay
    local clustervar ""
    local clustervar_orig ""
    * Cluster variables passed to the plugin: one per dimension of
    * vce(cluster), or the time variable of vce(dkraay)
    local ncl = 0
    local dkraay_bw = 0

    * pweight forces robust VCE (unless clustering)
    if `weight_type' == 3 {
        local vcetype = 1
    }

    * vce(dkraay [#]): Driscoll-Kraay standard errors over the time variable
    * of tsset or xtset, with bandwidth # (lags + 1); reghdfe's rules
    if `"`vce'"' != "" {
        local vce_word
        local vce_dkrest `"`vce'"'
        gettoken vce_word vce_dkrest : vce_dkrest, parse(" ,")
        if substr(`"`vce_word'"', 1, 2) == "dk" {
            if strpos(`"`vce_dkrest'"', ",") {
                gettoken vce_dkarg vce_dkopts : vce_dkrest, parse(",")
                gettoken comma vce_dkopts : vce_dkopts, parse(",")
                di as error `"VCE options not supported: `=strtrim(`"`vce_dkopts'"')'"'
                exit 9
            }
            local vce_dkrest = strtrim(`"`vce_dkrest'"')
            if `"`vce_dkrest'"' != "" {
                capture confirm integer number `vce_dkrest'
                if _rc {
                    di as error "vce(dkraay #): # must be a positive integer (bandwidth = lags + 1)"
                    exit 9
                }
                if `vce_dkrest' < 1 {
                    di as error "vce(dkraay #): # must be >= 1 (bandwidth = lags + 1)"
                    exit 9
                }
                local dkraay_bw = `vce_dkrest'
            }
            local timevar "`_dta[_TStvar]'"
            if "`timevar'" == "" {
                di as error "vce(dkraay) requires data to be tsset or xtset with a time variable"
                exit 9
            }
            markout `touse' `timevar'
            local vcetype = 3
            local ncl = 1
            local clustervar `timevar'
            local clustervar_orig `timevar'
            local vce
        }
    }

    if `"`vce'"' != "" {
        _ctools_hdfe_vce, spec(`"`vce'"') touse(`touse') cmd(creghdfe)
        local vce_kind `r(kind)'
        if "`vce_kind'" == "cluster" {
            local vcetype = 2
            local ncl = r(ndims)
            if `ncl' == 0 {
                di as error "creghdfe: cluster variable required with vce(cluster)"
                exit 198
            }
            if `ncl' > 10 {
                di as error "creghdfe: at most ten cluster variables are supported"
                exit 198
            }
            * Reported as reghdfe names them; a string variable or an
            * interaction is numbered first, as egen group() would
            forvalues d = 1/`ncl' {
                local cl_name_`d' `r(name_`d')'
                local cl_vars_`d' `r(vars_`d')'
                local cl_group_`d' = r(group_`d')
            }
            forvalues d = 1/`ncl' {
                local clustervar_orig `clustervar_orig' `cl_name_`d''
                if `cl_group_`d'' {
                    tempvar clustervar_numeric_`d'
                    quietly _ctools_group `cl_vars_`d'' if `touse', generate(`clustervar_numeric_`d'')
                    local clustervar `clustervar' `clustervar_numeric_`d''
                }
                else local clustervar `clustervar' `cl_vars_`d''
            }
        }
        else if "`vce_kind'" == "robust" {
            local vcetype = 1
        }
        else if inlist("`vce_kind'", "unadjusted", "ols") {
            * As in reghdfe: pweights need robust or clustered standard errors
            if `weight_type' == 3 {
                di as error "pweights do not work with vce(unadjusted), use a different vce()"
                exit 9
            }
            local vcetype = 0
        }
        else {
            di as error "creghdfe: unrecognized vce() option"
            exit 9
        }
    }

    * wildbootstrap(): it studentizes with the robust or cluster VCE; the
    * bootstrap clusters are the vce() cluster or bootcluster()
    local wb_bootvar
    local wb_bootname
    if `__wild' {
        if `vcetype' == 0 | `vcetype' == 3 | `ncl' > 1 {
            di as error "creghdfe: wildbootstrap() needs vce(robust) or vce(cluster {it:clustvar})"
            exit 198
        }
        if `"`wb_bootcluster'"' != "" {
            _ctools_hdfe_vce, spec(`"cluster `wb_bootcluster'"') touse(`touse') cmd(creghdfe)
            if r(ndims) != 1 {
                di as error "creghdfe: wildbootstrap(bootcluster()) takes one variable"
                exit 198
            }
            local wb_bootname `r(name_1)'
            tempvar wb_boot
            if r(group_1) quietly _ctools_group `r(vars_1)' if `touse', generate(`wb_boot')
            else quietly generate double `wb_boot' = `r(vars_1)' if `touse'
            local wb_bootvar `wb_boot'
        }
        else if `vcetype' == 2 local wb_bootname `clustervar_orig'
        if `vcetype' == 1 & "`wb_bootvar'" == "" & `weight_type' == 2 {
            di as error "creghdfe: wildbootstrap() draws one weight per observation under vce(robust); with fweights a row stands for several observations. Use bootcluster()"
            exit 198
        }
    }

    * Parse variable list
    gettoken depvar indepvars : varlist

    * Save original depvar name for labeling (before ts expansion)
    local depvar_orig "`depvar'"
    * syntax writes a factor-variable depvar 1.foreign as i1.foreign
    if regexm("`depvar_orig'", "^i([0-9]+[bno]*[.].+)$") local depvar_orig = regexs(1)

    * Expand factor-variable and time-series operators on depvar (creates
    * a temporary variable if needed); as in reghdfe, a factor-variable
    * depvar must name one indicator (1.foreign, not i.foreign)
    if strpos("`depvar'", ".") {
        fvexpand `depvar'
        if `: word count `r(varlist)'' > 1 {
            di as error "more than one depvar specified: `r(varlist)'"
            exit 9
        }
    }
    fvrevar `depvar'
    local depvar "`r(varlist)'"

    * Expand factor variables to temporary numeric variables
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

        * Build expanded varlist for plugin
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

    * group(): as reghdfe's ValidateGroups, the regression, cluster and
    * weight variables must be constant within each group, a group and an
    * individual identify an observation, and the regression sample is one
    * row per group (the first in group and individual order). With
    * individual(), __itouse marks every row of those groups and the plugin
    * gets the groups and individuals numbered.
    if `__teams' {
        if `__indiv' {
            if `weight_type' == 2 {
                di as error "fweights are incompatible with individual ids as there cannot be two observations for a given group-individual touple"
                exit 198
            }
            if `__bs' | `__fe' | `__wild' | "`groupvar'" != "" {
                di as error "creghdfe: vce(bootstrap), fevce(), wildbootstrap() and groupvar() are not available with individual()"
                exit 198
            }
        }
        markout `touse' `group_id' `individual_id'
        local __gvars
        if `"`varlist'"' != "" {
            fvrevar `varlist', list
            local __gvars `r(varlist)'
        }
        forvalues d = 1/`ncl' {
            local __gvars `__gvars' `cl_vars_`d''
        }
        if `vcetype' == 3 local __gvars `__gvars' `timevar'
        local __gvars : list uniq __gvars
        tempvar __itouse
        _creghdfe_groups, touse(`touse') itouse(`__itouse') group(`group_id') ///
            individual(`individual_id') vars(`__gvars') weightvar(`weight_var') weightexp(`"`exp'"')
        if `__indiv' {
            tempvar __gnum __inum
            quietly _ctools_group `group_id' if `__itouse', generate(`__gnum')
            quietly _ctools_group `individual_id' if `__itouse', generate(`__inum')
        }
    }

    * Count observations
    qui count if `touse'
    local nobs = r(N)

    if `nobs' == 0 {
        di as error "creghdfe: no observations"
        exit 2000
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Handle resid option - determine residual variable name
    local resid_varname ""
    local compute_resid = 0
    if "`__resid_name'" != "" | `__resid_default' {
        local compute_resid = 1
        if "`__resid_name'" != "" local resid_varname "`__resid_name'"
        else local resid_varname "_reghdfe_resid"
    }

    * Validate the entire output list before creating any destination.
    * As reghdfe, every run drops the variables __hdfe* that savefe saved
    * before, residuals without a name replace _reghdfe_resid, and
    * absorb(, replace) replaces the newvar= targets (their checks below
    * come first, so that a failed check leaves them)
    local output_names `resid_varname' `groupvar'
    if `savefe' {
        local output_names `output_names' `ext_saved'
    }
    local __replaced
    if `__resid_default' local __replaced _reghdfe_resid
    foreach __t of local __replace_targets {
        capture confirm variable `__t', exact
        if !_rc local __replaced `__replaced' `__t'
    }
    local fe_sevars
    local fe_nrepvars
    local fe_lbvars
    local fe_ubvars
    if `__fe' {
        * One output per saved factor (intercept or slope), named after
        * savefe's __hdfe#__ and __hdfe#__Slope#: __fese#__, __fese#__Slope#
        * (or sestub#, sestub#Slope#)
        local fe_nint = 0
        forvalues f = 1/`next' {
            local a = `ext_absvar_`f''
            if "`ext_slope_`f''" == "" local ++fe_nint
            local sfx : subinstr local ext_out_`f' "__hdfe`a'__" ""
            if "`fe_sestub'" != "" {
                local fe_sevars `fe_sevars' `fe_sestub'`a'`sfx'
                if `fe_nreps' local fe_nrepvars `fe_nrepvars' `fe_sestub'`a'`sfx'_nrep
                if `fe_ci' local fe_lbvars `fe_lbvars' `fe_sestub'`a'`sfx'_lb
                if `fe_ci' local fe_ubvars `fe_ubvars' `fe_sestub'`a'`sfx'_ub
            }
            else {
                local fe_sevars `fe_sevars' __fese`a'__`sfx'
                if `fe_nreps' local fe_nrepvars `fe_nrepvars' __fenrep`a'__`sfx'
                if `fe_ci' local fe_lbvars `fe_lbvars' __felb`a'__`sfx'
                if `fe_ci' local fe_ubvars `fe_ubvars' __feub`a'__`sfx'
            }
        }
        local output_names `output_names' `fe_sevars' `fe_nrepvars' `fe_lbvars' `fe_ubvars' `fe_dse'
    }
    if "`output_names'" != "" {
        * duplicates and names that must be new, then what is replaced
        local __checked
        foreach __v of local output_names {
            if !`: list __v in __replaced' & substr("`__v'", 1, 6) != "__hdfe" local __checked `__checked' `__v'
        }
        local __dups : list dups output_names
        if "`__dups'" != "" {
            di as error "ctools: duplicate output variable `: word 1 of `__dups''"
            exit 198
        }
        if "`__checked'" != "" _ctools_newvars `__checked'
    }
    capture drop __hdfe*
    foreach __v of local __replaced {
        capture novarabbrev drop `__v'
    }
    if "`output_names'" != "" _ctools_newvars `output_names'
    if `compute_resid' quietly gen double `resid_varname' = .

    * Set up parameters via Stata scalars (expected by C plugin)
    scalar __creghdfe_K = `nvars'           // depvar + indepvars
    scalar __creghdfe_G = `next'            // number of FE factors (intercepts and slopes)
    scalar __creghdfe_drop_singletons = ("`keepsingletons'" == "")
    scalar __creghdfe_verbose = ("`verbose'" != "")
    scalar __creghdfe_maxiter = `iterate'
    scalar __creghdfe_tolerance = `tolerance'
    scalar __creghdfe_standardize = ("`nostandardize'" == "")
    scalar __creghdfe_vce_type = `vcetype'
    scalar __creghdfe_compute_dof = 1
    scalar __creghdfe_df_a_nested = 0
    scalar __creghdfe_compute_resid = `compute_resid'
    scalar __creghdfe_has_weights = (`weight_type' > 0)
    scalar __creghdfe_weight_type = `weight_type'
    scalar __creghdfe_dof_adjust_type = `dof_adjust_type'
    scalar __creghdfe_dof_clusters = `dof_clusters'
    scalar __creghdfe_savefe = `savefe'
    scalar __creghdfe_compute_groupvar = ("`groupvar'" != "")
    scalar __creghdfe_use_quad = ("`quad'" != "")
    scalar __creghdfe_num_cluster = `ncl'
    scalar __creghdfe_dkraay_bw = `dkraay_bw'
    scalar __creghdfe_dof_continuous = strpos(" `dofmethods' ", " continuous ") > 0
    * transform() and acceleration() of reghdfe's MAP solver, and abort()
    scalar __creghdfe_transform = `__transform_code'
    scalar __creghdfe_acceleration = `__accel_code'
    scalar __creghdfe_abort = `__abort'

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record setup time
    _ctools_lap __time_setup __ct_lap98 `__ct_lap98'

    * Build varlist for plugin: depvar indepvars fe_vars [cluster_var] [weight_var] [resid_var]
    * Use expanded varlist (factor variables converted to temp numeric vars)
    local plugin_varlist `varlist_expanded' `absorb'
    if `vcetype' >= 2 {
        local plugin_varlist `plugin_varlist' `clustervar'
    }
    * Add weight variable to plugin varlist if using weights
    if `weight_type' > 0 {
        local plugin_varlist `plugin_varlist' `weight_var'
    }
    * Add residual variable to the plugin varlist if computing residuals
    if `compute_resid' {
        * Position in varlist = K + G + (1 if cluster) + (1 if weights) + 1
        local resid_var_pos = `nvars' + `next' + `ncl' + (`weight_type' > 0) + 1
        local plugin_varlist `plugin_varlist' `resid_varname'
        scalar __creghdfe_resid_var_idx = `resid_var_pos'
    }
    else {
        scalar __creghdfe_resid_var_idx = 0
    }

    * Handle groupvar - create variable and add to plugin varlist
    local groupvar_var_pos = 0
    if "`groupvar'" != "" {
        quietly gen long `groupvar' = .
        local groupvar_var_pos = `nvars' + `next' + `ncl' + (`weight_type' > 0) + `compute_resid' + 1
        local plugin_varlist `plugin_varlist' `groupvar'
        scalar __creghdfe_groupvar_idx = `groupvar_var_pos'
    }
    else {
        scalar __creghdfe_groupvar_idx = 0
    }

    * Handle savefe - create FE variables and add to plugin varlist
    if `savefe' {
        forval g = 1/`next' {
            quietly gen double `ext_out_`g'' = .
            if `__fe' quietly gen double `: word `g' of `fe_sevars'' = .
            if `__fe' & `fe_nreps' quietly gen long `: word `g' of `fe_nrepvars'' = .
            if `__fe' & `fe_ci' quietly gen double `: word `g' of `fe_lbvars'' = .
            if `__fe' & `fe_ci' quietly gen double `: word `g' of `fe_ubvars'' = .
            local plugin_varlist `plugin_varlist' `ext_out_`g''
        }
        * First savefe var position
        local savefe_var_pos = `nvars' + `next' + `ncl' + (`weight_type' > 0) + `compute_resid' + ("`groupvar'" != "") + 1
        scalar __creghdfe_savefe_idx = `savefe_var_pos'
    }
    else {
        scalar __creghdfe_savefe_idx = 0
    }

    * The plugin returns the final sample after internal row removals.
    tempvar final_sample
    quietly generate byte `final_sample' = 0
    local plugin_varlist `plugin_varlist' `final_sample'
    local sample_idx : word count `plugin_varlist'
    scalar __creghdfe_sample_idx = `sample_idx'

    * vce(bootstrap): the resampling units and strata (numbered in sorted
    * order, as bsample orders them)
    scalar __creghdfe_bs_unit_idx = 0
    scalar __creghdfe_bs_strata_idx = 0
    if "`bs_unitvar'" != "" {
        local plugin_varlist `plugin_varlist' `bs_unitvar'
        scalar __creghdfe_bs_unit_idx = `: word count `plugin_varlist''
    }
    if "`bs_stratavar'" != "" {
        local plugin_varlist `plugin_varlist' `bs_stratavar'
        scalar __creghdfe_bs_strata_idx = `: word count `plugin_varlist''
    }

    * Slope variables of the extended absorbed effects, and the absorbed
    * effect each factor belongs to
    forvalues g = 1/`next' {
        scalar __creghdfe_slope_`g' = 0
        scalar __creghdfe_absvar_`g' = `ext_absvar_`g''
        if "`ext_slope_`g''" != "" {
            local plugin_varlist `plugin_varlist' `ext_slope_`g''
            scalar __creghdfe_slope_`g' = `: word count `plugin_varlist''
        }
    }

    * individual(): the group and individual numbers and the rows of the
    * groups (the plugin reads the individuals of each group from them)
    scalar __creghdfe_indiv = `__indiv'
    if `__indiv' {
        local plugin_varlist `plugin_varlist' `__gnum'
        scalar __creghdfe_indiv_gidx = `: word count `plugin_varlist''
        local plugin_varlist `plugin_varlist' `__inum'
        scalar __creghdfe_indiv_iidx = `: word count `plugin_varlist''
        local plugin_varlist `plugin_varlist' `__itouse'
        scalar __creghdfe_indiv_tidx = `: word count `plugin_varlist''
        scalar __creghdfe_indiv_mean = `__indiv_mean'
    }

    * The sample marker goes last, in place of if `touse': the plugin reads
    * it in parallel, where Stata answers an if one observation at a time
    local plugin_varlist `plugin_varlist' `touse'
    scalar __creghdfe_touse = 1

    * Create matrix to store results (plugin will fill __creghdfe_V)
    local K_x = `nvars' - 1  // number of indepvars (excluding depvar)
    matrix __creghdfe_V = J(`K_x' + 1, `K_x' + 1, 0)
    matrix __creghdfe_V_psd = J(`K_x' + 1, `K_x' + 1, 0)

    * Validation hook: fail as Break would once the outputs exist
    if "$CTOOLS_creghdfe_fault" == "outputs" error 1

    * Call the C plugin with full_regression subcommand
    capture noisily plugin call ctools_plugin `plugin_varlist', ///
        "creghdfe `threads_code' full_regression`=cond(`__bsrun' | `__wild', " retain", "")'"

    local reg_rc = _rc
    if `reg_rc' {
        _ctools_cleanup creghdfe
        di as error "Error in regression (rc=`reg_rc')"
        exit `reg_rc'
    }

    * Record plugin call time
    _ctools_lap __time_plugin_call __ct_lap98 `__ct_lap98'

    * Retrieve results from scalars set by C plugin
    local N_final = __creghdfe_N
    local num_singletons = __creghdfe_num_singletons

    * Factor variables were expanded on the sample before singleton removal.
    * As reghdfe (which expands after its singleton pass), refit on the final
    * sample when that changes the terms, e.g. a level or the base level that
    * occurs only on singleton rows. Singletons are a fixed point, so the
    * refit drops none; it reports the first run's count.
    if `num_singletons' > 0 & "`indepvars'" != "" & "`__refit_sample'" == "" {
        fvexpand `indepvars' if `final_sample'
        if "`r(varlist)'" != "`coef_names_all'" {
            if "`resid_varname'" != "" drop `resid_varname'
            if "`groupvar'" != "" drop `groupvar'
            if `savefe' {
                forvalues g = 1/`next' {
                    drop `ext_out_`g''
                }
            }
            * fevce(bootstrap)'s outputs were created with savefe's; the refit
            * creates them again
            if `__fe' drop `fe_sevars' `fe_nrepvars' `fe_lbvars' `fe_ubvars'
            _ctools_cleanup creghdfe
            * individual(): the refit keeps every row of the groups kept
            local __refit_var `final_sample'
            if `__indiv' {
                tempvar __refit_var
                quietly generate byte `__refit_var' = 0
                mata: _g = st_data(., "`__gnum'"); _s = st_data(., "`final_sample'"); _k = J(max(_g), 1, 0); _k[select(_g, _s :== 1)] = J(sum(_s :== 1), 1, 1); _i = selectindex(_g :< .); _r = J(rows(_g), 1, 0); _r[_i] = _k[_g[_i]]; st_store(., "`__refit_var'", _r)
                mata drop _g _s _k _i _r
            }
            global CTOOLS_creghdfe_refit_ "`__refit_var' `num_singletons'"
            _creghdfe_impl `__creghdfe_args'
            exit
        }
    }

    * wildbootstrap(): before the pairs replicates of fevce(), which then
    * reuse the fit the plugin keeps
    if `__wild' {
        _creghdfe_wild, reps(`wb_reps') type(`wb_type') null(`wb_null') ptype(`wb_ptype_code') ///
            enumerate(`wb_enumerate') partialled(`=("`wb_fedgp'" == "partialled")') ci(`wb_ci') ///
            level(`wb_level') dots(`wb_dots') names(`coef_names') coefs(`wb_coefs') ///
            seed(`"`wb_seed'"') threads(`"`threads_code'"') ///
            vcevar(`=cond(`vcetype' == 2, "`clustervar'", "")') bootvar(`wb_bootvar') ///
            keep(`__fe') projected(`wb_projected')
        local wb_rngstate `"`r(rngstate)'"'
        local wb_test `r(test)'
        local wb_B = r(reps)
        local wb_enum = r(enumerated)
        local wb_C = r(C)
        local wb_path = r(path)
    }

    * vce(bootstrap): the replicates, in chunks, so that the dots advance and
    * Break is honored between them (the plugin keeps the fit meanwhile)
    if `__bsrun' {
        if "`fe_dse'" != "" quietly gen double `fe_dse' = .
        local bs_savetmp
        if `"`bs_save_file'"' != "" tempfile bs_savetmp
        _creghdfe_bootstrap, reps(`bs_reps') seed(`"`bs_seed'"') draws(`bs_draws') ///
            dots(`bs_dots') level(`bs_level') mse(`bs_mse') threads(`"`threads_code'"') ///
            coef(`__bs') fe(`__fe') fevars(`fe_sevars' `fe_nrepvars' `fe_lbvars' `fe_ubvars' `fe_dse') ///
            nfe(`next') nreps(`=0`fe_nreps'') ci(`=0`fe_ci'') dse(`=("`fe_dse'" != "")') ///
            bsdraws(`bs_bsdraws') savefile(`"`bs_savetmp'"') every(`bs_save_every') ///
            `=cond("`bs_strata'" != "" & "`bs_unitvar'" == "", "strataobs", "")'
        local bs_rngstate `"`r(rngstate)'"'
        local bs_nunits = scalar(__creghdfe_bs_U)
        local bs_nstrata = scalar(__creghdfe_bs_S)
    }
    if "`__refit_singletons'" != "" local num_singletons = `__refit_singletons'

    local K_keep = __creghdfe_K_keep
    local df_a = __creghdfe_df_a
    capture local num_iters = __creghdfe_iterations
    if _rc != 0 local num_iters = 1
    * abort(0) after iterate() iterations without convergence: as reghdfe,
    * e(converged) is 0, and its MAP solver then reports no iterations
    capture local __converged = scalar(__creghdfe_converged)
    if _rc local __converged = 1
    if !`__converged' & "`__technique'" == "map" local num_iters = 0

    * Get number of levels and nested status for each FE (computed by C plugin)
    local df_a_nested = 0
    forval g = 1/`next' {
        capture local fe_levels_`g' = __creghdfe_num_levels_`g'
        if _rc != 0 local fe_levels_`g' = .

        * Get nested status from C plugin (replaces slow egen-based check)
        capture local fe_nested_`g' = __creghdfe_fe_nested_`g'
        if _rc != 0 local fe_nested_`g' = 0

        * Accumulate df_a_nested
        if `fe_nested_`g'' == 1 {
            local df_a_nested = `df_a_nested' + `fe_levels_`g''
        }
    }

    * Record post-processing time (nested FE now computed in C plugin)
    _ctools_lap __time_postproc __ct_lap98 `__ct_lap98'

    * Adjust df_a for nested FEs
    local df_a_adjusted = `df_a' - `df_a_nested'
    if `df_a_adjusted' < 0 {
        local df_a_adjusted = 0
    }

    * Calculate R-squared and adjusted R-squared
    * Use scalars for precision-critical computations to avoid local precision loss
    tempname cr_rss cr_tss cr_tss_w cr_mss cr_mss_w
    scalar `cr_rss' = scalar(__creghdfe_rss)
    scalar `cr_tss' = scalar(__creghdfe_tss)
    scalar `cr_tss_w' = scalar(__creghdfe_tss_within)
    local r2 = 1 - scalar(`cr_rss') / scalar(`cr_tss')
    local r2_within = 1 - scalar(`cr_rss') / scalar(`cr_tss_w')
    local df_m = `K_keep'

    * df_r depends on VCE type:
    * - For unadjusted/robust: N - K - df_a
    * - For cluster: min(num_clusters - 1, N - K - df_a) (like reghdfe)
    local df_r_ols_base = `N_final' - `K_keep' - `df_a_adjusted'
    if `vcetype' >= 2 {
        * Cluster VCE: df_r = min(num_clusters - 1, df_r_ols), with the
        * fewest clusters of any dimension (the time periods for dkraay)
        * reghdfe uses the minimum to ensure internal consistency
        local num_clusters = __creghdfe_N_clust
        local nclust_list
        forvalues d = 1/`ncl' {
            local nclust_list `nclust_list' `=scalar(__creghdfe_N_clust_`d')'
        }
        local df_r_cluster = `num_clusters' - 1
        if `df_r_ols_base' < `df_r_cluster' {
            local df_r = `df_r_ols_base'
        }
        else {
            local df_r = `df_r_cluster'
        }
    }
    else {
        * Unadjusted or robust: df_r = N - K - df_a
        local df_r = `df_r_ols_base'
    }
    * A saturated model has no residual degrees of freedom, as in reghdfe
    if `df_r' < 0 {
        local df_r = 0
    }

    * R2_a and RMSE use df_r adjusted for nested FEs (like reghdfe)
    * used_df_r = N - df_a - df_m - df_a_nested
    * (both missing or rss itself without residual degrees of freedom)
    local df_r_ols = `N_final' - `df_a_adjusted' - `K_keep' - `df_a_nested'
    local r2_a = 1 - (scalar(`cr_rss')/`df_r_ols') / (scalar(`cr_tss')/(`N_final'-`any_intercept'))
    local rmse = cond(`df_r_ols' == 0, sqrt(scalar(`cr_rss')), sqrt(scalar(`cr_rss') / `df_r_ols'))
    * reghdfe's r2_a_within, ll and ll_0
    local r2_a_within = 1 - (scalar(`cr_rss')/`df_r_ols') / (scalar(`cr_tss_w')/(`df_r_ols' + `K_keep'))
    local ll = -0.5 * `N_final' * (1 + ln(2*_pi) + ln(scalar(`cr_rss') / `N_final'))
    local ll_0 = -0.5 * `N_final' * (1 + ln(2*_pi) + ln(scalar(`cr_tss_w') / `N_final'))

    * Calculate F-statistic
    * For unadjusted VCE: use traditional MSS/MSE formula with within-group variation
    * For robust/cluster VCE: use Wald test F = b'*V^-1*b / k
    * e(mss) = TSS - RSS (total model sum of squares for reporting)
    * But F uses within-group MSS = tss_within - rss
    scalar `cr_mss' = scalar(`cr_tss') - scalar(`cr_rss')
    scalar `cr_mss_w' = scalar(`cr_tss_w') - scalar(`cr_rss')

    * e(b) and e(V) as reghdfe posts them: every term in order, with zeros for
    * omitted regressors (o. prefix) and for factor-variable base levels
    tempname b V b_kept
    _ctools_ebv `b' `V', names(`coef_names_all') base(`coef_is_base') ///
        stub(__creghdfe) kkeep(`K_keep') vmatrix(__creghdfe_V) cons(__creghdfe_cons) kept(`b_kept')
    * HDFE-27: reghdfe's PSD correction of multi-way clustered and dkraay
    * covariances can leave e(V) indefinite, and Stata's ereturn post then
    * sets it to zero or keeps it depending on rounding. In that case the
    * plugin also returns a positive semidefinite version with the same
    * coefficient block, posted only when reghdfe's would be set to zero.
    local __psd = 0
    capture local __psd = !`__bs' & scalar(__creghdfe_psd_made) == 1
    if `__psd' != 1 local __psd = 0
    if `__psd' {
        tempname V_psd __b_psd
        _ctools_ebv `__b_psd' `V_psd', names(`coef_names_all') base(`coef_is_base') ///
            stub(__creghdfe) kkeep(`K_keep') vmatrix(__creghdfe_V_psd) cons(__creghdfe_cons)
    }
    if `__bs' {
        * e(V_modelbased) is the fit's own covariance; e(V) and the bootstrap
        * vectors take the same layout (zeros for omitted and base terms)
        tempname V_mb
        matrix `V_mb' = `V'
        local __ebv names(`coef_names_all') base(`coef_is_base') stub(__creghdfe) ///
            kkeep(`K_keep') cons(__creghdfe_cons)
        tempname __bx __bd __D
        _ctools_ebv `__bx' `V', `__ebv' vmatrix(__creghdfe_bs_V)
        foreach __v in b se bias z0 {
            tempname bs_`__v'
            matrix `__D' = diag(__creghdfe_bs_`__v')
            _ctools_ebv `__bx' `__bd', `__ebv' vmatrix(`__D')
            matrix `bs_`__v'' = vecdiag(`__bd')
        }
        foreach __v in ci_normal ci_percentile ci_bc {
            tempname bs_`__v'
            forvalues __r = 1/2 {
                matrix `__D' = diag(__creghdfe_bs_`__v'[`__r', 1...])
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
        _ctools_ebv `__smx' `__smV', names(`coef_names_all') base(`coef_is_base') stub(__creghdfe) ///
            kkeep(`K_keep') cons(__creghdfe_cons) vmatrix(`__smv')
        matrix `bs_map' = vecdiag(`__smV')
    }

    * Label the residual variable (already filled by C plugin)
    if "`resid_varname'" != "" {
        label var `resid_varname' "Residuals"
    }

    * Label the groupvar variable (mobility group identifier), as reghdfe
    * labels the groups of the first pair of absorbed effects not nested in
    * the clusters. Where reghdfe forms no such pair (effects nested in the
    * clusters, dofadjustments() without mobility groups) it creates no
    * variable, and neither does creghdfe.
    if "`groupvar'" != "" {
        if scalar(__creghdfe_groupvar_created) == 0 {
            drop `groupvar'
        }
        else if scalar(__creghdfe_groupvar_fe2) > 0 {
            local __fe1 : word `=scalar(__creghdfe_groupvar_fe1)' of `ext_names'
            local __fe2 : word `=scalar(__creghdfe_groupvar_fe2)' of `ext_names'
            label var `groupvar' "Mobility group between `__fe1' and `__fe2'"
        }
        else label var `groupvar' "Mobility group"
    }

    * Label the savefe variables (fixed effect estimates)
    if `savefe' {
        forval g = 1/`next' {
            local fevar `ext_label_`g''
            label var `ext_out_`g'' "FE for `fevar'"
            if `__fe' label var `: word `g' of `fe_sevars'' "Bootstrap SE of FE for `fevar'"
            if `__fe' & `fe_nreps' label var `: word `g' of `fe_nrepvars'' "Replicates estimating FE for `fevar'"
            if `__fe' & `fe_ci' label var `: word `g' of `fe_lbvars'' "Bootstrap lower bound of FE for `fevar'"
            if `__fe' & `fe_ci' label var `: word `g' of `fe_ubvars'' "Bootstrap upper bound of FE for `fevar'"
        }
    }

    * The Wald F-statistic of all kept regressors, computed by the plugin as
    * reghdfe computes it for every VCE (in the units of the standardized
    * regression; missing when V is rank deficient: too few clusters, a
    * perfect fit, a saturated model)
    local F = scalar(__creghdfe_F)

    * noconstant: as reghdfe, _cons is not reported (the constant is still
    * absorbed with the fixed effects); without an absorbed intercept there
    * is no constant
    if ("`constant'" == "noconstant" | !`any_intercept') & colsof(`b') > 1 {
        local __kb = colsof(`b') - 1
        matrix `b' = `b'[1, 1..`__kb']
        matrix `V' = `V'[1..`__kb', 1..`__kb']
        if `__bs' {
            matrix `V_mb' = `V_mb'[1..`__kb', 1..`__kb']
            foreach __v in b se bias z0 ci_normal ci_percentile ci_bc {
                matrix `bs_`__v'' = `bs_`__v''[1..., 1..`__kb']
            }
        }
        if `"`bs_save_file'"' != "" matrix `bs_map' = `bs_map'[1, 1..`__kb']
    }

    * Record results building time
    _ctools_lap __time_build __ct_lap98 `__ct_lap98'

    * Post results. With factor variables, buildfvinfo lets margins check
    * estimability, as after reghdfe.
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
    * HDFE-27: a trial post (no e(sample)) tells whether Stata keeps
    * reghdfe's indefinite e(V); when it would be set to zero, the PSD version
    * is posted instead. Without a constant the coefficient block is PSD.
    if `__psd' {
        if colsof(`V') != colsof(`V_psd') local __psd = 0
    }
    if `__psd' {
        tempname __b_try __V_try
        matrix `__b_try' = `b'
        matrix `__V_try' = `V'
        capture quietly ereturn post `__b_try' `__V_try'
        local __vzero = 0
        if !_rc mata: st_local("__vzero", strofreal(allof(st_matrix("e(V)"), 0) & !allof(st_matrix("`V'"), 0)))
        if `__vzero' matrix `V' = `V_psd'
        ereturn clear
    }
    * As reghdfe: a coefficient with zero b and zero variance is omitted
    _ms_findomitted `b' `V'
    * nosample: no e(sample), as reghdfe
    local __esample esample(`final_sample')
    if "`nosample'" != "" local __esample
    ereturn post `b' `V', `__esample' depname(`depvar_orig') obs(`N_final') `bfv'
    * The weights, which margins and other postestimation commands reuse
    if "`weight'" != "" {
        ereturn local wtype "`weight'"
        ereturn local wexp `"`exp'"'
    }

    * Store additional e() results
    ereturn scalar N = `N_final'
    ereturn scalar df_m = `df_m'
    * Bootstrap inference is z-based, as after vce(bootstrap): no e(df_r)
    if !`__bs' ereturn scalar df_r = `df_r'
    ereturn scalar r2 = `r2'
    ereturn scalar r2_within = `r2_within'
    ereturn scalar r2_a = `r2_a'
    ereturn scalar rmse = `rmse'
    ereturn scalar rss = scalar(`cr_rss')
    ereturn scalar mss = scalar(`cr_mss')
    ereturn scalar tss = scalar(`cr_tss')
    ereturn scalar tss_within = scalar(`cr_tss_w')
    ereturn scalar df_a = `df_a_adjusted'
    ereturn scalar df_a_nested = `df_a_nested'
    if !`__bs' ereturn scalar F = `F'
    ereturn scalar rank = `df_m'
    ereturn scalar num_singletons = `num_singletons'
    ereturn scalar drop_singletons = ("`keepsingletons'" == "")
    ereturn scalar N_full = `N_final' + `num_singletons'
    * Likelihood, fit and solver results reghdfe posts
    ereturn scalar r2_a_within = `r2_a_within'
    ereturn scalar ll = `ll'
    ereturn scalar ll_0 = `ll_0'
    ereturn scalar sumweights = scalar(__creghdfe_sumweights)
    ereturn scalar ic = `num_iters'
    ereturn scalar converged = `__converged'
    ereturn scalar report_constant = ("`constant'" != "noconstant") & `any_intercept'

    * The absorbed effects and their intercepts and slopes, the absorbed
    * degrees of freedom table and the covariance labels, as reghdfe posts
    * them; the model without absorb() absorbs _cons
    local extended
    if `__noabsorb' {
        ereturn local absorb ""
        local absvars _cons
    }
    else {
        ereturn local absorb `"`absorb_names'"'
        local absvars `"`absorb_names'"'
        local extended extended(`ext_names')
    }
    local vcename = cond(`vcetype'==0, "unadjusted", cond(`vcetype'==1, "robust", cond(`vcetype'==2, "cluster", "dkraay")))
    local cluster_post
    if `vcetype' >= 2 {
        local cluster_post clustvars(`clustervar_orig') nclust(`nclust_list')
    }
    * individual(): the individual effect's row of e(dof_table), in its
    * place among the absorbed effects: its individuals, none counted as
    * redundant, and that count marked inexact, as reghdfe reports them
    if `__indiv' {
        forvalues g = `next'(-1)`=`__indiv_pos' + 1' {
            foreach s in num_levels fe_redundant fe_exact fe_nested {
                scalar __creghdfe_`s'_`=`g' + 1' = scalar(__creghdfe_`s'_`g')
            }
        }
        local __p = `__indiv_pos' + 1
        scalar __creghdfe_num_levels_`__p' = scalar(__creghdfe_indiv_levels)
        scalar __creghdfe_fe_redundant_`__p' = 0
        scalar __creghdfe_fe_exact_`__p' = 0
        scalar __creghdfe_fe_nested_`__p' = 0
        local __ext
        forvalues g = 1/`__indiv_pos' {
            local __ext `__ext' `: word `g' of `ext_names''
        }
        local __ext `__ext' `__indiv_name'
        forvalues g = `=`__indiv_pos' + 1'/`next' {
            local __ext `__ext' `: word `g' of `ext_names''
        }
        local extended extended(`__ext')
    }
    _ctools_hdfe_post, stub(__creghdfe) absvars(`absvars') nhdfe(`nfe') vce(`vcename') ///
        `extended' dofmethod(`dofmethods') `cluster_post' reghdfenames
    if `vcetype' >= 2 {
        ereturn local title3 "Statistics robust to heteroskedasticity"
    }
    if `vcetype' == 3 {
        ereturn local title4 "and autocorrelation (Driscoll-Kraay)"
        ereturn scalar dkraay_bw = scalar(__creghdfe_dkraay_bw)
    }
    * reghdfe calls the unadjusted VCE ols
    if `vcetype' == 0 & !`__bs' ereturn local vce "ols"
    if `__fe' {
        if "`fe_dse'" != "" label var `fe_dse' "Bootstrap SE of the summed fixed effects"
        _creghdfe_fe_post, sevars(`fe_sevars') nrepvars(`fe_nrepvars') dse(`fe_dse') ///
            lbvars(`fe_lbvars') ubvars(`fe_ubvars') level(`bs_level') ///
            absvars(`ext_names') nint(`fe_nint') rngstate(`"`bs_rngstate'"') draws(`bs_draws') ///
            cluster(`bs_cluster') strata(`bs_strata')
    }
    if `__bs' {
        _creghdfe_bootstrap_post, vmb(`V_mb') bsb(`bs_b') bsse(`bs_se') bsbias(`bs_bias') ///
            bsz0(`bs_z0') cinormal(`bs_ci_normal') cipercentile(`bs_ci_percentile') ///
            cibc(`bs_ci_bc') level(`bs_level') rngstate(`"`bs_rngstate'"') ///
            cluster(`bs_cluster_name') strata(`bs_strata') draws(`bs_draws') ///
            command(`"creghdfe `__creghdfe_args'"') dfm(`df_m') mse(`bs_mse')
    }
    if `__wild' {
        _creghdfe_wild_post, test(`wb_test') names(`coef_names') method(`=cond(`wb_null', "wcr", "wcu")') ///
            type(`wb_type') ptype(`wb_ptype') reps(`wb_B') enum(`wb_enum') nclust(`wb_C') ///
            level(`wb_level') path(`wb_path') ci(`wb_ci') bootclust(`wb_bootname') ///
            rngstate(`"`wb_rngstate'"')
    }
    if `"`bs_save_file'"' != "" {
        local __nclust
        local __nstrata
        if "`bs_unitvar'" != "" local __nclust `bs_nunits'
        if "`bs_strata'" != "" local __nstrata `bs_nstrata'
        _ctools_hdfe_boot savefinish, tmp(`"`bs_savetmp'"') file(`"`bs_save_file'"') map(`bs_map') ///
            `bs_save_double' `bs_save_replace' command(`"creghdfe `__creghdfe_args'"') ///
            seed(`"`bs_rngstate'"') n(`N_final') cluster(`bs_cluster_name') nclust(`__nclust') ///
            strata(`bs_strata') nstrata(`__nstrata') cmdname(creghdfe)
    }
    ereturn local model "ols"
    ereturn local marginsnotok "Residuals SCore"
    ereturn local title "HDFE Linear regression"
    ereturn local title2 = "Absorbing `nfe' HDFE " + cond(`nfe' > 1, "groups", "group")
    ereturn local depvar "`depvar_orig'"
    * As reghdfe: the regressors partialled out (base levels left out,
    * levels in bn form) and _cons when it is reported
    local __indep
    local __k = 0
    foreach __n of local coef_names_all {
        local ++__k
        if `: word `__k' of `coef_is_base'' continue
        local __parts : subinstr local __n "#" " ", all
        local __t
        foreach __p of local __parts {
            if regexm("`__p'", "^([0-9]+)[bno]*[.](.+)$") local __p = regexs(1) + "bn." + regexs(2)
            local __t = cond("`__t'" == "", "`__p'", "`__t'#`__p'")
        }
        local __indep `__indep' `__t'
    }
    if e(report_constant) local __indep `__indep' _cons
    ereturn local indepvars "`__indep'"
    if "`resid_varname'" != "" {
        ereturn local resid "`resid_varname'"
    }
    ereturn local predict "creghdfe_p"
    * estat (summarize, vce, ic with e(df_m) + 1 degrees of freedom) and
    * the absorbed-degrees-of-freedom footnote of prefix displays, as reghdfe
    ereturn local estat_cmd "_creghdfe_estat"
    ereturn local footnote "_creghdfe_footnote"
    ereturn local cmd "creghdfe"
    ereturn local cmdline "creghdfe `0'"

    * Record ereturn posting time
    _ctools_lap __time_ereturn __ct_lap98 `__ct_lap98'

    * Display singleton and convergence messages (before header, like reghdfe)
    if `num_singletons' > 0 & !`__quiet' {
        di as text "(dropped " as result `num_singletons' as text " singleton observations)"
    }
    * reghdfe warns of Kaczmarz with conjugate gradient (with more than one
    * absorbed effect), and runs it
    if `__kcg' & `nfe' > 1 di as error "(WARNING: convergence is {bf:unlikely} with transform=kaczmarz and accel=CG)"
    * reghdfe's LSMR and LSQR warn when they stop at iterate() under abort(0)
    if !`__converged' & "`__technique'" != "map" & `__verbose_level' >= 0 {
        di as error "convergence not achieved in " strtrim(string(`num_iters', "%8.0fc")) " iterations (stop code=13); try increasing maxiter() or decreasing tol()."
    }
    * reghdfe shows the convergence line only at the default verbosity
    if `__verbose_level' == 0 di as text `"({browse "http://scorreia.com/research/hdfe.pdf":MWFE estimator} converged in "' as result `num_iters' as text " iterations)"
    * reghdfe's notes on the regressors it omits, and its warning when the
    * fixed effects explain the dependent variable (which is then zeroed),
    * after the convergence line, as there
    local __coltol = min(1e-6, `tolerance'/10)
    if scalar(__creghdfe_depvar_explained) == 1 & !`__quiet' {
        di as text "warning: dependent variable {res}`depvar_orig'{txt} is likely perfectly explained by the fixed effects (tol =" %3.1e `__coltol' ")"
    }
    local __k = 0
    foreach __v of local coef_names {
        local ++__k
        if `__quiet' continue
        if scalar(__creghdfe_fe_collinear_`__k') == 1 {
            di as text "note: {res}`__v'{txt} is probably collinear with the fixed effects (all partialled-out values are close to zero; tol =" %3.1e `__coltol' ")"
        }
        else if scalar(__creghdfe_collinear_`__k') == 1 {
            di as text "note: {res}`__v'{txt} omitted because of collinearity"
        }
    }
    * reghdfe's warnings from the covariance step: a dkraay bandwidth beyond
    * the time periods, and the eigenvalue correction of a multi-way
    * clustered covariance
    if `vcetype' == 3 & !`__quiet' {
        if `dkraay_bw' > scalar(__creghdfe_N_clust) {
            di as text "Warning: bandwidth `dkraay_bw' too large for `=scalar(__creghdfe_N_clust)' time periods; reducing to `=scalar(__creghdfe_N_clust)'"
        }
    }
    if `vcetype' == 2 & `ncl' > 1 & !`__quiet' {
        if scalar(__creghdfe_psd_fixed) == 1 di as text "Warning: VCV matrix was non-positive semi-definite; adjustment from Cameron, Gelbach & Miller applied."
    }
    * reghdfe's warning when the F test cannot be computed (a rank-deficient
    * covariance matrix: too few clusters, collinear regressors)
    if !`__bs' & `df_m' > 0 & missing(`F') & !`__quiet' {
        di as text "warning: missing F statistic; dropped variables due to collinearity or too few clusters"
    }

    * Header, coefficient table and absorbed degrees of freedom, as replay
    * shows them
    _creghdfe_display, `header' `table' `footnote' `nowildtable' `diopts'

    * Record display time and compute total
    _ctools_lap __time_display __ct_lap98 `__ct_lap98'
    local __time_total = `__time_setup' + `__time_plugin_call' + `__time_postproc' + `__time_build' + `__time_ereturn' + `__time_display'

    * Display timing if requested
    if `__do_timing' {
        di as text ""
        di as text "{hline 55}"
        di as text "creghdfe timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load (SPI):        " as result %8.4f _creghdfe_time_load " sec"
        di as text "    Data copy + sums:       " as result %8.4f _creghdfe_time_copy " sec"
        di as text "    FE remap + counts:      " as result %8.4f _creghdfe_time_remap " sec"
        di as text "    Singleton removal:      " as result %8.4f _creghdfe_time_singleton " sec"
        di as text "    DOF computation:        " as result %8.4f _creghdfe_time_dof " sec"
        di as text "    HDFE partial out:       " as result %8.4f _creghdfe_time_partial " sec"
        di as text "    OLS estimation:         " as result %8.4f _creghdfe_time_ols " sec"
        di as text "    VCE computation:        " as result %8.4f _creghdfe_time_vce " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f _creghdfe_time_total " sec"
        di as text ""
        di as text "  Stata overhead:"
        di as text "    Setup (parsing, etc):   " as result %8.4f `__time_setup' " sec"
        di as text "    Plugin call overhead:   " as result %8.4f (`__time_plugin_call' - _creghdfe_time_total) " sec"
        di as text "    Post-processing:        " as result %8.4f `__time_postproc' " sec"
        di as text "    Results building:       " as result %8.4f `__time_build' " sec"
        di as text "    ereturn posting:        " as result %8.4f `__time_ereturn' " sec"
        di as text "    Display output:         " as result %8.4f `__time_display' " sec"
        di as text "  {hline 53}"
        local __stata_overhead = `__time_setup' + (`__time_plugin_call' - _creghdfe_time_total) + `__time_postproc' + `__time_build' + `__time_ereturn' + `__time_display'
        di as text "    Stata overhead total:   " as result %8.4f `__stata_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `__time_total' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = _creghdfe_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = _creghdfe_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * timeit: reghdfe's two timer tables, with creghdfe's steps, and their
    * last one in r(timer), as reghdfe leaves it
    if `__timeit' {
        _creghdfe_timeit, setup(`__time_setup') total(`__time_total') ///
            plugin(`__time_plugin_call') post(`=`__time_postproc' + `__time_build' + `__time_ereturn'')
    }
    _ctools_cleanup creghdfe
    * Validation hook: fail after the plugin filled the outputs and e() was posted
    if "$CTOOLS_creghdfe_fault" == "post" error 1
end

* Results as reghdfe shows them, after estimation and on replay: the header
* (reghdfe_header's layout), the coefficient table (display options as for
* any estimation command) and the absorbed degrees of freedom (e(dof_table),
* _creghdfe_footnote), each of which noheader, notable and nofootnote
* suppress
program define _creghdfe_display
    version 14.1
    syntax [, noHEader noTABle noFOOTnote NOWILDtable *]
    _get_diopts diopts, `options'
    if "`header'" == "" _creghdfe_header
    if "`table'" == "" {
        if "`header'" == "" di as text ""
        ereturn display, `diopts'
        if "`e(vce)'" == "bootstrap" & e(N_misreps) > 0 & e(N_misreps) < . {
            di as text "Note: One or more parameters could not be estimated in " ///
                as result e(N_misreps) as text " bootstrap replicates;"
            di as text "      standard-error estimates include only complete replications."
        }
        if "`e(wb_method)'" != "" & "`nowildtable'" == "" _creghdfe_wild_table
    }
    if "`e(fevce)'" == "bootstrap" & e(fe_norm_warning) == 1 _creghdfe_fe_warning
    * The absorbed degrees of freedom (only when absorb() was specified)
    if "`footnote'" == "" _creghdfe_footnote
end

* reghdfe's header: on the right the sample size, the F (or chi2) test and
* the fit statistics; on the left the titles, then (aligned to the bottom)
* one line per cluster variable and the Driscoll-Kraay bandwidth. After
* creghdfe's own vce(bootstrap) the right column also has its replications.
program define _creghdfe_header
    version 14.1
    if !c(noisily) exit
    local C1 _col(1)
    local C2 _col(30)
    local C3 _col(51)
    local C4 _col(67)
    local nr = 0
    local R`++nr' `"`C3' "Number of obs" `C4' "= " as result %10.0fc e(N)"'
    capture confirm matrix e(boot_status)
    if !_rc & "`e(vce)'" == "bootstrap" {
        local R`++nr' `"`C3' "Replications" `C4' "= " as result %10.0fc e(N_reps)"'
    }
    if `"`e(chi2)'"' != "" | "`e(df_r)'" == "" {
        local type `e(chi2type)'
        if `"`type'"' == "" local type Wald
        if !missing(e(chi2)) {
            local R`++nr' `"`C3' "`type' chi2(" as result e(df_m) as text ")" `C4' "= " as result %10.2f e(chi2)"'
            local R`++nr' `"`C3' "Prob > chi2" `C4' "= " as result %10.4f chi2tail(e(df_m), e(chi2))"'
        }
        else {
            local R`++nr' `"`C3' "{help j_robustsingular##|_new:`type' chi2(`e(df_m)')}" as text `C4' "= " as result %10s ".""'
            local R`++nr' `"`C3' "Prob > chi2" `C4' "= " as result %10s ".""'
        }
    }
    else if !missing(e(F)) {
        local R`++nr' `"`C3' "F(" as result %4.0f e(df_m) as text "," as result %7.0f e(df_r) as text ")" `C4' "= " as result %10.2f e(F)"'
        local R`++nr' `"`C3' "Prob > F" `C4' "= " as result %10.4f Ftail(e(df_m), e(df_r), e(F))"'
    }
    else {
        local dfm_l : di %4.0f e(df_m)
        local dfm_l2 : di %7.0f e(df_r)
        local R`++nr' `"`C3' "{help j_robustsingular##|_new:F(`dfm_l',`dfm_l2')}" as text `C4' "= " as result %10s ".""'
        local R`++nr' `"`C3' "Prob > F" `C4' "= " as result %10s ".""'
    }
    if !missing(e(r2)) local R`++nr' `"`C3' "R-squared" `C4' "= " as result %10.4f e(r2)"'
    if !missing(e(r2_a)) local R`++nr' `"`C3' "Adj R-squared" `C4' "= " as result %10.4f e(r2_a)"'
    if !missing(e(r2_within)) local R`++nr' `"`C3' "Within R-sq." `C4' "= " as result %10.4f e(r2_within)"'
    if !missing(e(rmse)) local R`++nr' `"`C3' "Root MSE" `C4' "= " as result %10.4f e(rmse)"'

    * The titles go on the left when they fit, followed by an empty line
    local len = strlen(`"`e(title)'"')
    forvalues i = 2/4 {
        if `"`e(title`i')'"' != "" local len = max(`len', strlen(`"`e(title`i')'"'))
    }
    local nl = 0
    local above = (`len' >= 49)
    if !`above' {
        local L`++nl' `""`e(title)'""'
        forvalues i = 2/4 {
            if `"`e(title`i')'"' != "" local L`++nl' `""`e(title`i')'""'
        }
        local L`++nl'
    }
    local nc = cond(missing(e(N_clustervars)), 0, e(N_clustervars))
    local space = `nr' - `nl' - `nc'
    forvalues i = 1/`space' {
        local L`++nl'
    }
    local clustvars `e(clustvar)'
    forvalues d = 1/`nc' {
        gettoken cv clustvars : clustvars
        local L`++nl' `"`C1' "Number of clusters (" as result "`cv'" as text ") " `C2' as text "= " as result %10.0fc e(N_clust`d')"'
    }
    if "`e(vce)'" == "dkraay" & e(dkraay_bw) < . {
        local L`++nl' `"`C1' "DK Bandwidth" `C2' "= " as result %10.0fc e(dkraay_bw)"'
    }
    local K = max(`nl', `nr')
    di
    if `above' {
        di as text `"`e(title)'"'
        forvalues i = 2/4 {
            if `"`e(title`i')'"' != "" di as text `"`e(title`i')'"'
        }
        if `K' di
    }
    forvalues i = 1/`K' {
        di as text `L`i'' as text `R`i''
    }
end

* wildbootstrap or wildbootstrap(wildopts) among the remaining options:
* caller locals __wild (0/1), __wildspec and __wild_rest (the other options)
program define _creghdfe_wild_opt
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

* verbose(#) among the remaining options, returned in caller locals
program define _creghdfe_verbose_opt
    version 14.1
    syntax [, VERbose(integer 0) *]
    c_local __verbose_level `verbose'
    c_local __options `"`options'"'
end

* absorb(absvars, suboptions): reghdfe's savefe (or its synonym generate),
* replace (newvar= targets that exist are replaced) and noisily (the parsed
* effects are shown); anything else (in caller local absorb_badopts) is an
* error, as in reghdfe
program define _creghdfe_absorb_opts
    version 14.1
    syntax [, SAVEfe Generate NOIsily REPLACE *]
    c_local absorb_badopts `"`options'"'
    c_local savefe = ("`savefe'`generate'" != "")
    c_local absorb_noisily = ("`noisily'" != "")
    c_local absorb_replace = ("`replace'" != "")
end

* residuals(newvar), also resid2(newvar) and reghdfe's abbreviations res()
* and resid(); residuals (or res, resid) without a name, for _reghdfe_resid.
* Caller locals __resid_name (the new variable or empty), __resid_default
* (1 for _reghdfe_resid) and __resid_rest (the other options).
program define _creghdfe_resid_opt
    version 14.1
    * resid2() first: its name is spelled out, so it takes no other spelling
    syntax [, RESID2(name) *]
    local named `resid2'
    local 0 `", `options'"'
    syntax [, RESiduals(name) RESiduals2 *]
    if "`named'" != "" & "`residuals'" != "" {
        di as error "option residuals() invalid: specify residuals() once"
        exit 198
    }
    if "`named'" == "" local named `residuals'
    opts_exclusive "`named' `residuals2'" residuals
    c_local __resid_name `named'
    c_local __resid_default = ("`residuals2'" != "")
    c_local __resid_rest `"`options'"'
end

* reghdfe's undocumented key(value) options, which it assigns to its Mata
* FixedEffects object: s(abort) is abort() (1 by default: stop with r(430)
* when the projection does not converge in iterate() iterations; 0: keep
* it, with e(converged) = 0). miniter(), min_ok(), accel_start() and
* accel_freq() are accepted; they tune parts of reghdfe's solvers that
* creghdfe's does not have (miniter() applies to LSMR and LSQR only, and
* creghdfe's Aitken acceleration starts at iteration 6, every third, as
* reghdfe's defaults). Anything else is an error, as in reghdfe.
program define _creghdfe_undoc_opts, sclass
    version 14.1
    local abort = 1
    local rest `"`0'"'
    if `"`rest'"' != "" {
        local all = strtrim(`"`rest'"')
        while `"`rest'"' != "" {
            gettoken key rest : rest, parse(" (")
            if `"`key'"' == "" continue, break
            local open = substr(`"`rest'"', 1, 1)
            if `"`open'"' != "(" {
                di as error `"options {bf:`all'} not allowed"'
                exit 198
            }
            gettoken value rest : rest, match(paren)
            local known abort miniter min_ok accel_start accel_freq
            if !`: list key in known' {
                di as error `"option {bf:`key'} not allowed"'
                exit 198
            }
            capture confirm number `value'
            if _rc {
                di as error `"option {bf:`key'} not allowed"'
                exit 198
            }
            if "`key'" == "abort" local abort = (`value' != 0)
        }
    }
    sreturn clear
    sreturn local abort `abort'
end

* technique(), transform(), acceleration() and preconditioner() with
* reghdfe's names, abbreviations and errors, and parallel() as reghdfe
* parses it. creghdfe runs reghdfe's MAP solver (symmetric Kaczmarz with
* conjugate gradient by default); transform() and acceleration() select
* its variants, as in reghdfe. technique(lsmr|lsqr) and preconditioner()
* are accepted: creghdfe still solves with MAP, to the same tolerance, as
* it does with individual(). Returns s(technique), s(transform_code) (0
* symmetric Kaczmarz, 1 Kaczmarz, 2 Cimmino) and s(accel_code) (0 none, 1
* conjugate gradient, 2 steepest descent, 3 Aitken).
program define _creghdfe_solver_opts, sclass
    version 14.1
    syntax [, TECHnique(string) TRANSform(string) ACCELeration(string) PRECONDitioner(string) ///
        PARALLEL(string asis)]
    if "`technique'" == "" local technique map
    if !inlist("`technique'", "map", "cg", "lsmr", "lsqr") {
        di as error "invalid technique: `technique'"
        exit 9
    }
    if "`technique'" == "cg" {
        di as error `"invalid algorithm "cg""'
        exit 3498
    }
    if "`transform'" == "" local transform symmetric_kaczmarz
    local transform = lower("`transform'")
    foreach x in cimmino kaczmarz symmetric_kaczmarz rand_kaczmarz {
        if strpos("`x'", "`transform'") == 1 local transform `x'
    }
    if !inlist("`transform'", "cimmino", "kaczmarz", "symmetric_kaczmarz", "rand_kaczmarz") {
        di as error "invalid transform: `transform'"
        exit 9
    }
    if "`acceleration'" == "" local acceleration conjugate_gradient
    local acceleration = lower("`acceleration'")
    if "`acceleration'" == "cg" local acceleration conjugate_gradient
    if "`acceleration'" == "sd" local acceleration steepest_descent
    if "`acceleration'" == "off" local acceleration none
    foreach x in conjugate_gradient steepest_descent aitken none {
        if strpos("`x'", "`acceleration'") == 1 local acceleration `x'
    }
    if !inlist("`acceleration'", "conjugate_gradient", "steepest_descent", "aitken", "none") {
        di as error "invalid acceleration: `acceleration'"
        exit 9
    }
    if "`preconditioner'" == "" local preconditioner block_diagonal
    local preconditioner = lower("`preconditioner'")
    if "`preconditioner'" == "off" local preconditioner none
    foreach x in none diagonal block_diagonal {
        if strpos("`x'", "`preconditioner'") == 1 local preconditioner `x'
    }
    if !inlist("`preconditioner'", "none", "diagonal", "block_diagonal") {
        di as error "invalid preconditioner: `preconditioner'"
        exit 9
    }
    * reghdfe's MAP has no random Kaczmarz transform
    if "`technique'" == "map" & "`transform'" == "rand_kaczmarz" {
        di as error "NULL transform function"
        exit 3498
    }
    * parallel(#[, options]): reghdfe's worker processes; creghdfe runs in
    * threads (threads()) and ignores it
    if `"`parallel'"' != "" {
        local 0 `"`parallel'"'
        syntax anything(name=maxprocesses id="number of worker processes"), ///
            [MAXprocesses2(integer 0) ID(integer 0) TMP_path(string) FORCE *]
        confirm integer number `maxprocesses'
        if !inrange(`maxprocesses', 0, 1000) {
            di as error "assertion is false"
            exit 9
        }
        if `maxprocesses' == 1 & "`force'" == "" {
            di as text "(ignoring parallel(1) as it's slower than parallel(0)"
        }
        else if `maxprocesses' > 1 {
            di as text "(creghdfe ignores parallel(): it runs in threads; see threads())"
        }
    }
    local tcode = cond("`transform'" == "kaczmarz", 1, cond("`transform'" == "cimmino", 2, 0))
    local acode = cond("`acceleration'" == "none", 0, cond("`acceleration'" == "steepest_descent", 2, ///
        cond("`acceleration'" == "aitken", 3, 1)))
    * LSMR and LSQR have no transform or acceleration
    if "`technique'" != "map" {
        local tcode 0
        local acode 1
    }
    sreturn clear
    sreturn local kcg = ("`transform'" == "kaczmarz" & "`acceleration'" == "conjugate_gradient")
    sreturn local technique `technique'
    sreturn local transform_code `tcode'
    sreturn local accel_code `acode'
end

* timeit: reghdfe's timer tables (top level and partialling out) with
* creghdfe's steps, in seconds and as a percentage of the total; r(timer)
* is the second, as reghdfe leaves it
program define _creghdfe_timeit, rclass
    version 14.1
    syntax, SETUP(real) TOTAL(real) PLUGIN(real) POST(real)
    tempname top part
    local load = scalar(_creghdfe_time_load) + scalar(_creghdfe_time_copy)
    local hdfe = scalar(_creghdfe_time_remap) + scalar(_creghdfe_time_singleton)
    local dof = scalar(_creghdfe_time_dof)
    local partial = scalar(_creghdfe_time_partial)
    local ols = scalar(_creghdfe_time_ols) + scalar(_creghdfe_time_vce)
    local other = max(0, `plugin' - scalar(_creghdfe_time_total))
    matrix `top' = (`hdfe' \ `dof' \ `setup' \ `load' + `partial' \ `ols' \ 0)
    local sum = 0
    forvalues i = 1/6 {
        local sum = `sum' + `top'[`i', 1]
    }
    matrix `top' = `top' \ (`total' - `sum')
    matrix `top' = `top', 100 * `top' / max(`total', 1e-12)
    matrix rownames `top' = "HDFE" "DoF" "Expand factors/Lags" "Partial out" "Solve OLS" "Store alphas" "(Remainder)"
    matrix colnames `top' = "Time" "(% Total)"
    di as text _n `"{bf: Timer results:} Top-level"'
    matlist `top', noblank cspec(& %25s | %10.2fc | %9.1f &) rspec(||&&&&&&|) rowtitle(Step)
    local tp = `load' + `partial'
    matrix `part' = (`load' \ 0 \ `partial' \ 0)
    matrix `part' = `part', 100 * `part' / max(`tp', 1e-12)
    matrix rownames `part' = "Load data" "Standardize/etc" "MAP/LSMR" "(Remainder)"
    matrix colnames `part' = "Time" "(% Total)"
    di as text _n `"{bf: Timer results:} Partial-out"'
    matlist `part', noblank cspec(& %25s | %10.2fc | %9.1f &) rspec(||&&&|) rowtitle(Step)
    return matrix timer = `part'
end

* vce(bootstrap): run the replicates on the fit the plugin retained, in
* chunks of about 50 (dots and Break between chunks; a multiple of the
* replicates the plugin runs at once, __creghdfe_bs_par, so that no thread
* idles at a chunk's end), then summarize them into the
* __creghdfe_bs_* matrices. Returns r(rngstate), the random-number state
* before the first draw. draws(stata) draws bsample's uniforms with Stata's
* generator and the plugin maps them as bsample does; draws(fast) takes a
* master seed from Stata's generator (docs/development/BOOTSTRAP_PLAN.md, 5.8).
program define _creghdfe_bootstrap, rclass
    version 14.1
    syntax, REPS(integer) DRAWS(name) DOTS(integer) LEVEL(real) MSE(integer) ///
        [SEED(string) THREADS(string) COEF(integer 1) FE(integer 0) FEVARS(string) ///
        NFE(integer 1) NREPS(integer 0) CI(integer 0) DSE(integer 0) BSDRAWS(string) ///
        SAVEFILE(string) EVERY(string) STRATAOBS]
    local targets = cond(`fe', " fe", "") + cond(`dse', " dse", "") + cond(`ci', " ci", "")
    local Kb = scalar(__creghdfe_K_keep) + 1
    _ctools_hdfe_boot run, cmd(creghdfe) reps(`reps') draws(`draws') dots(`dots') kb(`Kb') ///
        seed(`"`seed'"') threads(`"`threads'"') targets(`"`targets'"') bsdraws(`bsdraws') ///
        savefile(`"`savefile'"') every(`every') `strataobs'
    return local rngstate `"`r(rngstate)'"'
    matrix __creghdfe_bs_V = J(`Kb', `Kb', .)
    foreach m in b se bias z0 {
        matrix __creghdfe_bs_`m' = J(1, `Kb', .)
    }
    foreach m in ci_normal ci_percentile ci_bc {
        matrix __creghdfe_bs_`m' = J(2, `Kb', .)
    }
    matrix __creghdfe_bs_counts = J(1, 5, .)
    if `fe' matrix __creghdfe_bs_fe_table = J(`nfe', 5, .)
    local finish = cond(`coef', " coef", "") + cond(`fe', " fe", "") + ///
        cond(`nreps', " nreps", "") + cond(`ci', " ci", "") + cond(`dse', " dse", "") + ///
        cond(`mse', " mse", "")
    plugin call ctools_plugin `fevars', ///
        "creghdfe `threads' boot_finish level=`level'`finish'"
end

* wildbootstrap(): the wild bootstrap on the fit the plugin retained. The
* weights come from Stata's generator as boottest draws them (enumerated in
* the plugin when all 2^C Rademacher vectors fit in reps), passed in chunks
* of 50 replicates. Results: matrices __creghdfe_wb_p, _ci, _cistat; r(test)
* (tested regressor positions), r(reps), r(enumerated), r(C), r(path),
* r(rngstate) (before the first draw).
program define _creghdfe_wild, rclass
    version 14.1
    syntax, REPS(integer) TYPE(name) NULL(integer) PTYPE(integer) ENUMERATE(integer) ///
        PARTIALLED(integer) CI(integer) LEVEL(real) DOTS(integer) ///
        [NAMES(string) COEFS(string) SEED(string) THREADS(string) VCEVAR(string) ///
        BOOTVAR(string) KEEP(integer 0) PROJECTED(integer 0)]
    local K : word count `names'
    local test
    if `"`coefs'"' == "" {
        forvalues k = 1/`K' {
            if scalar(__creghdfe_collinear_`k') != 1 & scalar(__creghdfe_fe_collinear_`k') != 1 ///
                local test `test' `k'
        }
    }
    else {
        foreach c of local coefs {
            local p : list posof `"`c'"' in names
            if !`p' {
                di as error "creghdfe: wildbootstrap(coefs()): `c' is not a regressor"
                exit 198
            }
            if scalar(__creghdfe_collinear_`p') == 1 | scalar(__creghdfe_fe_collinear_`p') == 1 {
                di as error "creghdfe: wildbootstrap(coefs()): `c' was omitted"
                exit 198
            }
            local test `test' `p'
        }
    }
    local T : word count `test'
    if `T' == 0 {
        di as error "creghdfe: wildbootstrap() has no coefficient to test"
        exit 198
    }
    matrix __creghdfe_wb_test = J(1, `T', .)
    forvalues j = 1/`T' {
        matrix __creghdfe_wb_test[1, `j'] = `: word `j' of `test''
    }
    if `"`seed'"' != "" set seed `seed'
    return local rngstate `"`c(rngstate)'"'
    local flags = cond("`vcevar'" != "", " vce", "") + cond("`bootvar'" != "", " boot", "") + ///
        cond(!`null', " nonull", "") + cond(`partialled', " partialled", "") + ///
        cond(`projected', " projected", "")
    plugin call ctools_plugin `vcevar' `bootvar', "creghdfe `threads' wild_setup`flags'"
    local C = scalar(__creghdfe_wb_C)
    local path = scalar(__creghdfe_wb_path)
    * Any failure or Break while the weights are drawn and passed leaves
    * neither their variables nor the Mata draws behind
    capture noisily {
        _ctools_hdfe_boot wilddraw, clusters(`C') reps(`reps') type(`type') ///
            `=cond(`enumerate', "", "noenumerate")'
        local B = r(reps)
        local enum = r(enumerate)
        if `dots' _dots 0, title(Wild bootstrap replications) reps(`B')
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
                    "creghdfe `threads' wild_chunk first=`first' count=`n' total=`B' enum"
            }
            else {
                _ctools_hdfe_boot wildweights, first(`=`first' + 1') count(`n') variables(`owned_weights')
                local wvars `r(varlist)'
                capture noisily plugin call ctools_plugin `wvars', ///
                    "creghdfe `threads' wild_chunk first=`first' count=`n' total=`B'"
                local rc = _rc
                if `rc' {
                    capture drop `owned_weights'
                    _ctools_hdfe_boot wildclear
                    exit `rc'
                }
            }
            if `dots' {
                forvalues r = 1/`n' {
                    _dots `=`first' + `r'' 0
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
    if `dots' di as text ""
    _ctools_hdfe_boot wildclear
    matrix __creghdfe_wb_p = J(1, `T', .)
    matrix __creghdfe_wb_ci = J(`T', 2, .)
    matrix __creghdfe_wb_cistat = J(`T', 5, 0)
    plugin call ctools_plugin, "creghdfe `threads' wild_finish level=`level' ptype=`ptype'`=cond(`ci', " ci", "")'`=cond(`keep', " keep", "")'"
    return local test `test'
    return scalar reps = `B'
    return scalar enumerated = `enum'
    return scalar C = `C'
    return scalar path = `path'
end

* wildbootstrap(): its e() results (plan 4.7), with the prefix wb_
program define _creghdfe_wild_post, eclass
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
    matrix `p' = __creghdfe_wb_p
    matrix colnames `t' = `coefs'
    matrix colnames `p' = `coefs'
    ereturn matrix wb_t = `t'
    ereturn matrix wb_p = `p'
    if `ci' {
        matrix `cim' = __creghdfe_wb_ci
        matrix `stm' = __creghdfe_wb_cistat
        matrix rownames `cim' = `coefs'
        matrix colnames `cim' = ll ul
        matrix rownames `stm' = `coefs'
        matrix colnames `stm' = open_lower open_upper disjoint empty failed
        ereturn matrix wb_ci = `cim'
        ereturn matrix wb_cistatus = `stm'
    }
    ereturn scalar wb_reps = `reps'
    ereturn scalar wb_misreps = scalar(__creghdfe_wb_misreps)
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

* wildbootstrap(): the table after the coefficient table (plan 4.7)
program define _creghdfe_wild_table
    version 14.1
    tempname t p ci st
    matrix `t' = e(wb_t)
    matrix `p' = e(wb_p)
    local hasci = 0
    capture confirm matrix e(wb_ci)
    if !_rc {
        local hasci = 1
        matrix `ci' = e(wb_ci)
        matrix `st' = e(wb_cistatus)
    }
    local kind = cond("`e(wb_bootclust)'" != "", "Wild cluster bootstrap", "Wild bootstrap")
    local how = cond("`e(wb_method)'" == "wcr", "restricted (WCR)", "unrestricted (WCU)")
    local wt `e(wb_weighttype)'
    if inlist("`wt'", "rademacher", "mammen", "webb") local wt = strproper("`wt'")
    local who = cond("`e(wb_bootclust)'" != "", "(`e(wb_bootclust)')", "(observations)")
    local pl = cond("`e(wb_ptype)'" == "lower", "P<t", cond("`e(wb_ptype)'" == "upper", "P>t", "P>|t|"))
    di as text ""
    di as text "`kind': `how', `wt' weights"
    di as text "Replications = " as result strtrim(string(e(wb_reps), "%12.0fc")) as text ///
        _col(40) "Bootstrap clusters = " as result strtrim(string(e(wb_N_clust), "%12.0fc")) ///
        as text " `who'"
    di as text "{hline 13}{c TT}{hline 64}"
    di as text %12s abbrev("`e(depvar)'", 12) " {c |}" _col(20) "Coef." _col(34) "t" _col(39) "`pl'" ///
        _col(50) "[`=strofreal(e(wb_level))'% Conf. Interval]"
    di as text "{hline 13}{c +}{hline 64}"
    local names : colnames `p'
    local j 0
    local open
    local disjoint
    local empty
    local trialfail
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
        di as text %12s abbrev("`c'", 12) " {c |}" as result _col(15) %10.0g _b[`c'] ///
            _col(27) %8.2f `t'[1, `j'] _col(38) %6.4f `p'[1, `j'] _col(49) %10.0g `lo' ///
            _col(62) %10.0g `hi'
    }
    di as text "{hline 13}{c BT}{hline 64}"
    local pdesc = cond("`e(wb_ptype)'" == "symmetric", "symmetric", cond("`e(wb_ptype)'" == "equaltail", "equal-tailed", "`e(wb_ptype)'-tailed"))
    di as text "P-values are `pdesc' bootstrap-t" _c
    if `hasci' di as text "; intervals invert the bootstrap test." _c
    di as text ""
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
    if e(wb_misreps) > 0 & e(wb_misreps) < . {
        di as text "Note: " e(wb_misreps) " replicates had no defined statistic and were left out."
    }
end

* vce(bootstrap): the e() results Stata's vce(bootstrap) posts, so that
* estat bootstrap works, and the Wald chi2 of the kept regressors
program define _creghdfe_bootstrap_post, eclass
    version 14.1
    syntax, VMB(name) BSB(name) BSSE(name) BSBIAS(name) BSZ0(name) CINORMAL(name) ///
        CIPERCENTILE(name) CIBC(name) LEVEL(real) DFM(integer) COMMAND(string asis) ///
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
    matrix `reps' = J(1, colsof(e(b)), scalar(__creghdfe_bs_N_reps))
    matrix colnames `reps' = `names'
    matrix rownames `reps' = y1
    matrix `counts' = __creghdfe_bs_counts
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
    ereturn scalar N_reps = scalar(__creghdfe_bs_N_reps)
    ereturn scalar N_misreps = scalar(__creghdfe_bs_N_misreps)
    ereturn scalar level = `level'
    ereturn scalar bs_version = 3
    ereturn scalar k_exp = 0
    ereturn scalar k_eexp = colsof(e(b))
    ereturn scalar k_extra = 0
    ereturn scalar k_eq = 1
    mata: _creghdfe_bs_wald()
    ereturn scalar chi2 = `__chi2'
    ereturn scalar p = chi2tail(`dfm', `__chi2')
    ereturn local chi2type Wald
    local i 0
    foreach name of local names {
        local ++i
        ereturn local exp`i' "_b[`name']"
    }
    ereturn local command `"`command'"'
    ereturn local cmdname creghdfe
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
program define _creghdfe_fe_post, eclass
    version 14.1
    syntax, SEVARS(string) ABSVARS(string) [NREPVARS(string) DSE(string) ///
        LBVARS(string) UBVARS(string) LEVEL(real 95) NINT(integer -1) ///
        RNGSTATE(string) DRAWS(name) CLUSTER(string) STRATA(string)]
    * the normalization concerns the absorbed intercepts (a slope is
    * identified within its level)
    if `nint' < 0 local nint : word count `absvars'
    tempname table
    matrix `table' = __creghdfe_bs_fe_table
    matrix colnames `table' = levels never always min_reps median_reps
    matrix rownames `table' = `absvars'
    ereturn matrix fe_bs_table = `table'
    ereturn scalar fe_N_reps = scalar(__creghdfe_bs_fe_N_reps)
    ereturn scalar fe_split_reps = scalar(__creghdfe_bs_fe_split_reps)
    ereturn scalar fe_norm_warning = (`nint' >= 2)
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
program define _creghdfe_fe_warning
    version 14.1
    local G : word count `e(absvars)'
    di as text ""
    di as text "{p 0 9 2}warning: with `G' absorbed effects (`e(absvars)'), individual fixed effects are " ///
        "identified only up to a normalization within each connected set of levels. The bootstrap " ///
        "standard errors in `e(fe_se)' depend on that normalization and can be large where the sets " ///
        "are thinly connected. " _c
    if "`e(fe_dse)'" != "" di as text "The summed fixed effect d does not depend on it; its " ///
        "standard errors are in `e(fe_dse)'.{p_end}"
    else di as text "The summed fixed effect d does not depend on it: specify " ///
        "fevce(bootstrap, dse({it:newvar})) for its standard errors.{p_end}"
    if e(fe_split_reps) > 0 & e(fe_split_reps) < . {
        di as text "{p 9 9 2}In " e(fe_split_reps) " of " e(fe_N_reps) " replicates a connected set " ///
            "split; levels in a split set were treated as not estimable in those replicates (see " ///
            "e(fe_bs_table)).{p_end}"
    }
end

* Release the fit the plugin retained for vce(bootstrap), after a failure
program define _creghdfe_bs_release
    version 14.1
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    plugin call ctools_plugin, "creghdfe boot_release"
end

capture mata: mata drop _creghdfe_bs_wald()
mata:
/* Wald chi2 of the regressors estimated (not omitted, not _cons) under
 * e(V), in the caller's local __chi2 (missing when V is singular there) */
void _creghdfe_bs_wald()
{
    real rowvector b, keep
    real matrix V, Vi
    string colvector names
    real scalar chi2

    b = st_matrix("e(b)")
    V = st_matrix("e(V)")
    names = st_matrixcolstripe("e(b)")[., 2]
    keep = selectindex((diagonal(V)' :!= 0) :& (names' :!= "_cons"))
    chi2 = .
    if (cols(keep)) {
        Vi = invsym(V[keep, keep])
        if (rank(Vi) == cols(keep)) chi2 = b[keep] * Vi * b[keep]'
    }
    st_local("__chi2", strofreal(chi2, "%21.17g"))
}
end


* group(), individual() and aggregation(), checked as reghdfe checks them:
* s(group_id) and s(individual_id) (variable names, or empty) and
* s(aggregation) (mean or sum)
program define _creghdfe_teams_opt, sclass
    version 14.1
    syntax [, group_id(string) individual_id(string) aggregation(string)]
    sreturn clear
    local opts
    if `"`group_id'"' != "" local opts `"`opts' group(`group_id')"'
    if `"`individual_id'"' != "" local opts `"`opts' individual(`individual_id')"'
    local 0 `", `opts'"'
    syntax [, Group(varname numeric) Individual(varname numeric)]
    if "`individual'" != "" & "`group'" == "" {
        di as error "cannot set the individual() identifiers without the group() identifiers"
        exit 198
    }
    if `"`aggregation'"' == "" local aggregation mean
    if inlist(`"`aggregation'"', "average", "avg") local aggregation mean
    if !inlist(`"`aggregation'"', "mean", "sum") {
        di as error "assertion is false"
        exit 9
    }
    sreturn local group_id `group'
    sreturn local individual_id `individual'
    sreturn local aggregation `aggregation'
end

* group(): reghdfe's ValidateGroups. The variables (and the weight) must be
* constant within each group, and a group and an individual must identify
* the observations of touse. itouse is created as touse, and touse is then
* set to the first row of each group (in group and individual order).
program define _creghdfe_groups, sortpreserve
    version 14.1
    syntax, TOUSE(name) ITOUSE(name) GROUP(varname) [INDIVIDUAL(string) VARS(string) ///
        WEIGHTVAR(string) WEIGHTEXP(string asis)]
    sort `touse' `group' `individual'
    foreach var of local vars {
        capture by `touse' `group': assert `var' == `var'[1] if `touse'
        if _rc {
            di as error "variable `var' is not constant within `group'"
            exit 498
        }
    }
    if "`weightvar'" != "" {
        capture by `touse' `group': assert `weightvar' == `weightvar'[1] if `touse'
        if _rc {
            local weightexp = strtrim(subinstr(`"`weightexp'"', "=", "", 1))
            di as error `"variable `weightexp' is not constant within `group'"'
            exit 498
        }
    }
    if "`individual'" != "" {
        capture by `touse' `group' `individual': assert _n == 1 if `touse'
        if _rc {
            di as error "identifiers for group (`group') and individual (`individual') do not uniquely identify the observations'"
            exit 459
        }
    }
    quietly generate byte `itouse' = `touse'
    quietly by `touse' `group': replace `touse' = (_n == 1) & `itouse'
end
