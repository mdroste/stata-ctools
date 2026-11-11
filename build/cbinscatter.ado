*! version 1.2.6 20261111 github.com/mdroste/stata-ctools

program define cbinscatter, eclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    local __cbin_args `"`0'"'
    * binscatter's abbreviations: n() is nquantiles(), m() msymbols(). Its
    * legacy options nbins(), create_xq, x_q(), symbols(), unique() and
    * linetype(noline) are taken as binscatter takes them.
    syntax varlist(min=2 fv) [aw fw pw iw] [if] [in], ///
        [Nquantiles(integer 20)] ///
        [Controls(varlist fv ts)] ///
        [Absorb(varlist)] ///
        [BY(varname)] ///
        [LINEtype(string)] ///
        [DISCrete] ///
        [GENxq(name)] ///
        [SAVEdata(string)] ///
        [REPORTreg] ///
        [NOgraph] ///
        [Verbose] ///
        [METHod(string)] ///
        [THReads(integer 0)] ///
        /* Graph options */ ///
        [TITLE(string asis)] ///
        [YTItle(string asis)] ///
        [XTItle(string asis)] ///
        [LEGEND(string asis)] ///
        [COLors(string)] ///
        [Msymbols(string)] ///
        [MLAbels(string asis)] ///
        /* binscatter options */ ///
        [MEDians noAddmean rd(numlist ascending) xq(varname numeric)] ///
        [MColors(string) LColors(string) savegraph(string)] ///
        [randvar(varname numeric) randcut(real 1) randn(integer -1)] ///
        [nbins(integer 20) create_xq x_q(varname numeric) symbols(string) unique(string)] ///
        [*]

    * Several y variables, as binscatter y1 y2 ... x: one binned scatter per
    * y over a common sample, drawn in one graph
    if `: word count `varlist'' > 2 {
        _cbinscatter_multi `__cbin_args'
        exit
    }

    * =========================================================================
    * UPFRONT VALIDATION - check all options before any data manipulation
    * This ensures errors don't leave the dataset in a bad state
    * =========================================================================

    * Own only the frame created by this invocation
    tempname graph_frame

    * binscatter's legacy options, mapped as binscatter maps them
    local __xlast : word `: word count `varlist'' of `varlist'
    _cbinscatter_legacy, xvar(`__xlast') nquantiles(`nquantiles') nbins(`nbins') ///
        genxq(`genxq') `create_xq' xq(`xq') x_q(`x_q') msymbols(`msymbols') ///
        symbols(`symbols') linetype(`linetype') unique(`unique')

    * Validate nquantiles early
    if `nquantiles' < 0 {
        di as error "nquantiles() invalid"
        exit 3300
    }
    if `nquantiles' < 2 {
        di as error "nquantiles() must be greater than or equal to 2"
        exit 198
    }

    * Validate linetype option early
    local linetype_num = 1   // Default to linear
    if "`linetype'" == "" {
        local linetype = "linear"
        local linetype_lower = "linear"
        local linetype_num = 1
    }
    else {
        local linetype_lower = lower("`linetype'")
        if "`linetype_lower'" == "none" | "`linetype_lower'" == "n" {
            local linetype_num = 0
        }
        else if "`linetype_lower'" == "linear" | "`linetype_lower'" == "line" | "`linetype_lower'" == "lfit" {
            local linetype_num = 1
        }
        else if "`linetype_lower'" == "quadratic" | "`linetype_lower'" == "qfit" {
            local linetype_num = 2
        }
        else if "`linetype_lower'" == "cubic" {
            local linetype_num = 3
        }
        else if "`linetype_lower'" == "connect" {
            local linetype_num = 0  // No fit, just connected points
        }
        else {
            di as error "cbinscatter: invalid linetype() option: `linetype'"
            di as error "  valid options: none, linear, qfit, cubic, connect"
            exit 198
        }
    }

    if "`reportreg'" != "" & `linetype_num' == 0 {
        di as error "Cannot specify reportreg when no fit line is being created."
        exit 198
    }
    local connect_opt
    if "`linetype_lower'" == "connect" local connect_opt "connect(l)"

    * Validate method option early (Cattaneo et al. "On Binscatter")
    * 0 = classic (default): residualize both Y and X, bin on residualized X
    * 1 = binsreg: bin on raw X, residualize Y only (conditional means)
    local method_num = 0   // Default to classic
    if "`method'" != "" {
        local method_lower = lower("`method'")
        if "`method_lower'" == "classic" | "`method_lower'" == "c" | "`method_lower'" == "residualize" {
            local method_num = 0
        }
        else if "`method_lower'" == "binsreg" | "`method_lower'" == "b" | "`method_lower'" == "conditional" {
            local method_num = 1
        }
        else {
            di as error "cbinscatter: invalid method() option: `method'"
            di as error "  valid options: classic, binsreg"
            exit 198
        }
    }

    * Validate savedata filename if specified (check directory exists)
    if "`savedata'" != "" {
        * Extract directory path
        local savedir = ""
        local lastslash = max(strpos(reverse("`savedata'"), "/"), strpos(reverse("`savedata'"), "\"))
        if `lastslash' > 0 {
            local savedir = substr("`savedata'", 1, strlen("`savedata'") - `lastslash')
        }
        if "`savedir'" != "" {
            * direxists(): confirm file "dir/." fails for directories on Windows
            mata: st_local("__dirok", strofreal(direxists(st_local("savedir"))))
            if !`__dirok' {
                di as error "cbinscatter: directory does not exist for savedata(): `savedir'"
                exit 601
            }
        }
    }

    * replace and nofastxtile are binscatter's; the rest are twoway options.
    * Without a graph they are checked against twoway, so a typo is an error
    * as it is with the graph.
    _cbinscatter_options, `options'
    local options `"`s(options)'"'
    local savereplace `"`s(replace)'"'
    local fastxtile `"`s(fastxtile)'"'

    * binscatter's options for its own (classic) method. binsreg's method
    * bins raw x per group with regression-adjusted means, which they do not
    * apply to.
    local rand_given = ("`randvar'" != "" | `randcut' != 1 | `randn' != -1)
    if `method_num' == 1 {
        local classic_only
        if "`medians'" != "" local classic_only `classic_only' medians
        if "`addmean'" != "" local classic_only `classic_only' noaddmean
        if "`rd'" != "" local classic_only `classic_only' rd()
        if "`xq'" != "" local classic_only `classic_only' xq()
        if "`randvar'" != "" local classic_only `classic_only' randvar()
        if `randcut' != 1 local classic_only `classic_only' randcut()
        if `randn' != -1 local classic_only `classic_only' randn()
        if "`classic_only'" != "" {
            di as error "cbinscatter: `classic_only' not allowed with method(binsreg)"
            exit 198
        }
    }
    if "`fastxtile'" == "nofastxtile" & `rand_given' {
        di as error "Cannot combine randvar, randcut or randn with nofastxtile"
        exit 198
    }
    if ("`genxq'" != "" & ("`xq'" != "" | "`discrete'" != "")) | ("`xq'" != "" & "`discrete'" != "") {
        di as error "Cannot specify more than one of genxq(), xq(), and discrete simultaneously."
        exit 198
    }
    if "`genxq'" != "" confirm new variable `genxq'
    if "`xq'" != "" {
        * As binscatter checks it, over every observation
        capture assert `xq' == int(`xq') & `xq' > 0
        if _rc {
            di as error "xq() must contain only positive integers."
            exit 198
        }
        if `nquantiles' != 20 {
            di as error "Cannot specify nquantiles in combination with discrete or an xq variable."
            exit 198
        }
        * (once for several y: _cbinscatter_multi sets the global after the first)
        if "`controls'`absorb'" != "" & "$CBINSCATTER_NOTES_SHOWN" == "" di as text "warning: xq() is specified in combination with controls() or absorb(). note that binning takes places after residualization, so the xq variable should contain bins of the residuals."
    }
    * binscatter takes quantile bins from a random subsample only when it
    * forms quantile bins (not with xq() or discrete); fastxtile's checks
    local rand_mode = 0
    if `rand_given' & "`xq'`discrete'" == "" {
        if `randn' != -1 {
            if `randcut' != 1 {
                di as error "cannot specify both randcut() and randn()"
                exit 198
            }
            if `randn' < 1 {
                di as error "randn() must be a positive integer"
                exit 198
            }
        }
        else if "`randvar'" == "" & !(`randcut' > 0 & `randcut' < 1) & `randcut' != 1 {
            di as error "if randcut() is specified without randvar(), a uniform r.v. will be generated and randcut() must be in (0,1)"
            exit 198
        }
        local rand_mode = cond("`randvar'" != "", 1, 2)
    }
    if `"`savegraph'"' != "" & "`nograph'" != "" {
        di as error "cbinscatter: savegraph() requires the graph (nograph specified)"
        exit 198
    }
    * savegraph() without replace refuses an existing file, as binscatter does
    if `"`savegraph'"' != "" & "`savereplace'" == "" {
        if regexm(`"`savegraph'"', "\.[a-zA-Z0-9]+$") confirm new file `"`savegraph'"'
        else confirm new file `"`savegraph'.gph"'
    }
    * savedata() without an extension writes binscatter's <name>.csv and
    * <name>.do (method(classic)) besides the bins in <name>.dta; a name with
    * an extension is the .dta file alone
    local save_csv = 0
    if "`savedata'" != "" {
        local __base = subinstr(`"`savedata'"', "\", "/", .)
        local __cut = strpos(reverse(`"`__base'"'), "/")
        if `__cut' local __base = substr(`"`__base'"', strlen(`"`__base'"') - `__cut' + 2, .)
        local save_csv = !strpos(`"`__base'"', ".") & `method_num' == 0
    }
    * As in binscatter, savedata() does not overwrite a file without replace
    if "`savedata'" != "" & "`savereplace'" == "" {
        local __savefile `"`savedata'"'
        * save adds .dta to a file name without an extension (r() is left
        * alone: the if qualifier is not marked yet)
        if !strpos(`"`__base'"', ".") {
            local __savefile `"`__savefile'.dta"'
        }
        if `save_csv' {
            confirm new file `"`savedata'.csv"'
            confirm new file `"`savedata'.do"'
        }
        capture confirm new file `"`__savefile'"'
        if _rc == 602 {
            di as error `"file `__savefile' already exists; specify replace"'
            exit 602
        }
    }
    if "`nograph'" != "" {
        local check_gropts `"`options'"'
        if `"`title'"' != "" local check_gropts `"`check_gropts' title(`title')"'
        if `"`ytitle'"' != "" local check_gropts `"`check_gropts' ytitle(`ytitle')"'
        if `"`xtitle'"' != "" local check_gropts `"`check_gropts' xtitle(`xtitle')"'
        if `"`legend'"' != "" local check_gropts `"`check_gropts' legend(`legend')"'
        _cbinscatter_check_gropts, `check_gropts'
    }

    * =========================================================================
    * END UPFRONT VALIDATION - now proceed with data processing
    * =========================================================================

    local weight_var ""
    local weight_type = 0

    * Parse variable list
    gettoken depvar xvar : varlist
    local depvar = trim("`depvar'")
    local xvar = trim("`xvar'")

    * Validate that y and x are numeric (not string) - match binscatter rc=109
    capture confirm numeric variable `depvar'
    if _rc {
        di as error "type mismatch"
        exit 109
    }
    capture confirm numeric variable `xvar'
    if _rc {
        di as error "type mismatch"
        exit 109
    }

    local num_absorb : word count `absorb'

    * Count the requested sample before excluding incomplete observations
    marksample touse, novarlist

    * Stata counts the sample and drops incomplete observations only where
    * the ado needs them before the plugin runs (by-groups, discrete bins,
    * the verbose header, the reportreg regressions, xq() bins, the size of
    * a random subsample); otherwise the plugin does both while it loads
    local stata_sample = ("`by'" != "" | "`discrete'" != "" | "`verbose'" != "" | ///
        "`reportreg'" != "" | "`xq'" != "" | `rand_mode' > 0)
    if `stata_sample' {
        qui count if `touse'
        local nobs = r(N)
        if `nobs' == 0 {
            di as error "no observations"
            exit 3301
        }
    }

    local check_weights = 0
    if "`weight'" != "" {
        * A weight that is a numeric variable is passed as it is, and the
        * plugin checks it; other expressions are evaluated and checked here
        local wexp = strtrim(substr(`"`exp'"', 2, .))
        capture unab wvar : `wexp', max(1)
        if !_rc capture confirm numeric variable `wvar'
        if !_rc & "`wvar'" == "`wexp'" {
            local weight_var "`wvar'"
            local check_weights = 1
        }
        else {
            tempvar evaluated_weight
            _ctools_weight, generate(`evaluated_weight') touse(`touse') type(`weight') expression(`"`exp'"')
            local weight_var "`evaluated_weight'"
        }
        if "`weight'" == "aweight" local weight_type = 1
        else if "`weight'" == "fweight" local weight_type = 2
        else if "`weight'" == "pweight" local weight_type = 3
        else if "`weight'" == "iweight" local weight_type = 4
    }
    local controls_original `controls'
    if "`controls'" != "" {
        fvrevar `controls' if `touse'
        local controls `r(varlist)'
    }
    local num_controls : word count `controls'

    * Group IDs must be dense after excluding incomplete observations.
    * absorb() values are labels (numeric or string); the plugin maps each
    * variable to dense levels over the final sample.
    if `stata_sample' {
        markout `touse' `depvar' `xvar' `controls' `weight_var'
        if "`absorb'" != "" markout `touse' `absorb', strok
        if "`by'" != "" markout `touse' `by', strok
        if "`xq'" != "" markout `touse' `xq'
    }

    * randvar()/randcut()/randn(): the quantile cutpoints come from a random
    * subsample of the n observations in the sample, as in binscatter's
    * fastxtile. randn() is a share of n. Without randvar(), binscatter
    * draws runiform() for every observation of the data sorted by the
    * sample and x, so the sample's observation of rank k in x gets draw
    * _N - n + k: those draws, in that order, are passed in observations
    * 1..n, and the plugin gives them to the observations by rank.
    local rand_var
    if `rand_mode' > 0 {
        if "`weight_var'" != "" quietly count if `touse' & `weight_var' > 0 & !missing(`weight_var')
        else quietly count if `touse'
        local rand_n = r(N)
        if `randn' != -1 {
            if `randn' > `rand_n' {
                if "$CBINSCATTER_NOTES_SHOWN" == "" di as text "randn() is larger than the population. using the full population."
                local rand_mode = 0
            }
            else {
                local randcut = `randn' / `rand_n'
                if "`randvar'" != "" {
                    quietly summarize `randvar', meanonly
                    if r(min) < 0 | r(max) > 1 {
                        di as error "with randn(), the randvar specified must be in [0,1] and ought to be uniformly distributed"
                        exit 198
                    }
                }
            }
        }
    }
    if `rand_mode' == 1 {
        tempvar rand_var
        quietly generate byte `rand_var' = `touse' & `randvar' <= `randcut'
    }
    else if `rand_mode' == 2 {
        * The draws are undone if x turns out to be discrete (binscatter then
        * draws nothing)
        local rand_state `"`c(rngstate)'"'
        tempvar rand_draw rand_var
        quietly generate `rand_draw' = runiform()
        if `rand_n' > 0 quietly generate double `rand_var' = `rand_draw'[_n + _N - `rand_n'] in 1/`rand_n'
        else quietly generate double `rand_var' = .
    }
    local num_by_groups = 1
    if "`by'" != "" {
        local by_orig "`by'"
        tempvar by_numeric
        quietly _ctools_group `by' if `touse', generate(`by_numeric') label
        local by "`by_numeric'"
        local num_by_groups = r(N_groups)
        if `num_by_groups' == 0 error 2000
        * binscatter's legend key for each group, byvarname=byvalname: a
        * numeric value as a macro expression evaluates it, or its value
        * label; a string's text (from egen group(), label)
        local by_reps `"`r(reps)'"'
        capture confirm numeric variable `by_orig'
        local by_isnum = !_rc
        local by_vlabel
        if `by_isnum' local by_vlabel : value label `by_orig'
        forval g = 1/`num_by_groups' {
            if `by_isnum' {
                local rep : word `g' of `by_reps'
                local byval = `by_orig'[`rep']
                if "`by_vlabel'" == "" local by_key`g' = `byval'
                else local by_key`g' : label `by_vlabel' `byval'
            }
            else local by_key`g' : label (`by') `g'
        }
    }

    * Load the platform-appropriate ctools plugin if not already loaded
    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    * Set up parameters via Stata scalars
    scalar __cbinscatter_nquantiles = `nquantiles'
    scalar __cbinscatter_linetype = `linetype_num'
    scalar __cbinscatter_compute_se = 0
    scalar __cbinscatter_discrete = ("`discrete'" != "")
    scalar __cbinscatter_has_controls = (`num_controls' > 0)
    scalar __cbinscatter_num_controls = `num_controls'
    scalar __cbinscatter_has_absorb = (`num_absorb' > 0)
    scalar __cbinscatter_num_absorb = `num_absorb'
    scalar __cbinscatter_has_by = ("`by'" != "")
    scalar __cbinscatter_has_weights = (`weight_type' > 0)
    scalar __cbinscatter_weight_type = `weight_type'
    scalar __cbinscatter_verbose = ("`verbose'" != "")
    scalar __cbinscatter_reportreg = ("`reportreg'" != "")
    scalar __cbinscatter_method = `method_num'
    scalar __cbinscatter_maxiter = 10000
    scalar __cbinscatter_tolerance = 1e-10
    * The plugin selects the sample from a marker variable (or, without
    * if/in, from the data alone) rather than from an if condition, which
    * Stata answers one observation at a time
    scalar __cbinscatter_touse_mode = cond(`"`if'`in'"' != "", 1, 2)
    scalar __cbinscatter_check_weights = `check_weights'
    scalar __cbinscatter_addmean = ("`addmean'" == "")
    scalar __cbinscatter_medians = ("`medians'" != "")
    scalar __cbinscatter_has_xq = ("`xq'" != "")
    scalar __cbinscatter_rand_mode = `rand_mode'
    scalar __cbinscatter_randcut = `randcut'
    scalar __cbinscatter_genxq = ("`genxq'" != "")
    local n_rd : word count `rd'
    scalar __cbinscatter_n_rd = `n_rd'
    if `n_rd' > 0 {
        local rd_row : subinstr local rd " " ",", all
        matrix __cbinscatter_rd = (`rd_row')
    }

    * Display info if verbose
    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "cbinscatter: C-Accelerated Binned Scatter"
        di as text "{hline 60}"
        di as text "Y variable: " as result "`depvar'"
        di as text "X variable: " as result "`xvar'"
        if `num_controls' > 0 {
            di as text "Controls:   " as result "`controls'"
        }
        if `num_absorb' > 0 {
            di as text "Absorb:     " as result "`absorb'"
        }
        if "`by'" != "" {
            di as text "By:         " as result "`by'"
        }
        di as text "Bins:       " as result `nquantiles'
        di as text "Obs:        " as result `nobs'
        di as text "{hline 60}"
        di ""
    }

    * Create matrices for C plugin to fill
    local bin_rows = `nquantiles' * `num_by_groups'
    if "`xq'" != "" {
        * One bin per (group, xq() value)
        tempvar bin_group
        quietly _ctools_group `by' `xq' if `touse', generate(`bin_group')
        local bin_rows = max(1, min(r(N_groups), c(max_matdim)))
    }
    else if "`discrete'" != "" {
        * One bin per (group, value): counted exactly for a raw x, bounded by
        * the sample for a residualized x and trimmed after the plugin. Both
        * are capped at the largest matrix; the plugin reports more bins.
        if `method_num' == 0 & `num_controls' + `num_absorb' > 0 {
            local bin_rows = `nobs'
        }
        else {
            tempvar bin_group
            quietly _ctools_group `by' `xvar' if `touse', generate(`bin_group')
            local bin_rows = r(N_groups)
        }
        local bin_rows = max(1, min(`bin_rows', c(max_matdim)))
    }
    else if `bin_rows' > c(max_matdim) {
        * A group has no more bins than observations; above that the bins of
        * every group must fit in one matrix
        local bin_rows = min(`nquantiles', _N) * `num_by_groups'
        if `bin_rows' > c(max_matdim) {
            di as error "cbinscatter: nquantiles() bins for each by() group need more rows than a Stata matrix holds (c(max_matdim) = `c(max_matdim)')"
            _ctools_cleanup cbinscatter
            exit 198
        }
    }
    matrix __cbinscatter_bins = J(`bin_rows', 5, .)
    * With rd(), one fit per group and interval
    if `linetype_num' > 0 {
        matrix __cbinscatter_coefs = J(`num_by_groups' * (`n_rd' + 1), 4, .)
        matrix __cbinscatter_fit_stats = J(`num_by_groups' * (`n_rd' + 1), 2, .)
    }

    * Build varlist for plugin: y x [controls] [absorb] [by] [weight] [xq]
    * [random subsample] [touse] [genxq]
    local plugin_varlist `depvar' `xvar'
    if "`controls'" != "" {
        local plugin_varlist `plugin_varlist' `controls'
    }
    if "`absorb'" != "" {
        local plugin_varlist `plugin_varlist' `absorb'
    }
    if "`by'" != "" {
        local plugin_varlist `plugin_varlist' `by'
    }
    if `weight_type' > 0 {
        local plugin_varlist `plugin_varlist' `weight_var'
    }
    if "`xq'" != "" {
        local plugin_varlist `plugin_varlist' `xq'
    }
    if `rand_mode' > 0 {
        local plugin_varlist `plugin_varlist' `rand_var'
    }
    if `"`if'`in'"' != "" {
        local plugin_varlist `plugin_varlist' `touse'
    }
    * genxq(): fastxtile's storage type for nquantiles() bins
    if "`genxq'" != "" {
        local xq_type = cond(`nquantiles' <= 100, "byte", cond(`nquantiles' <= 32740, "int", "long"))
        quietly generate `xq_type' `genxq' = .
        local plugin_varlist `plugin_varlist' `genxq'
    }

    * Build threads option string
    local threads_code ""
    if `threads' > 0 {
        local threads_code "threads(`threads')"
    }

    * Record start time
    quietly _ctools_clock __ct_lap99

    * Call the C plugin
    capture noisily plugin call ctools_plugin `plugin_varlist', ///
        "cbinscatter `threads_code' compute_bins"

    local plugin_rc = _rc
    if `plugin_rc' {
        * The plugin itself reports an empty sample and invalid weights
        capture local reported = __cbinscatter_reported
        if _rc local reported = 0
        _ctools_cleanup cbinscatter
        if "`genxq'" != "" capture drop `genxq'
        if !`reported' di as error "cbinscatter: plugin error (rc=`plugin_rc')"
        exit `plugin_rc'
    }
    if !`stata_sample' local nobs = __cbinscatter_N_touse
    * An x with no more values than nquantiles() is binned as discrete: as
    * binscatter, no random subsample is drawn and genxq() is not made then
    local discrete_used = __cbinscatter_discrete_used
    if `rand_mode' == 2 & `discrete_used' set rngstate `rand_state'
    if "`genxq'" != "" {
        if `discrete_used' {
            drop `genxq'
            di as text `"note: the x-variable has fewer unique values than the number of bins specified (`nquantiles').  It will therefore be treated as discrete, and genxq() will be ignored"'
        }
        else label variable `genxq' "`nquantiles' quantiles of `xvar'"
    }
    if "`discrete'" != "" {
        local bins_used = __cbinscatter_bins_used
        if `bins_used' >= 1 & `bins_used' < rowsof(__cbinscatter_bins) {
            matrix __cbinscatter_bins = __cbinscatter_bins[1..`bins_used', 1...]
        }
    }

    _ctools_lap elapsed __ct_lap99 `__ct_lap99'

    * Retrieve results from scalars
    local N_used = __cbinscatter_N
    local N_dropped = `nobs' - `N_used'
    local actual_num_groups = __cbinscatter_num_groups

    * reportreg prints the underlying fit regressions, as binscatter does.
    * Compute these only on request; plotted estimates remain the plugin's.
    * They use the plugin's sample (marked out above). summarize takes
    * pweights as the aweights they equal for a mean, and areg, which has no
    * iweights, is replaced by regress on dummies. A failure leaves no
    * scratch results behind.
    if "`reportreg'" != "" {
        capture noisily _cbinscatter_reportreg, depvar(`depvar') xvar(`xvar') ///
            touse(`touse') weight(`weight') weightvar(`weight_var') ///
            controls(`controls') absorb(`absorb') by(`by') method(`method_num') ///
            linetype(`linetype_num') groups(`actual_num_groups') rd(`rd') `addmean'
        if _rc {
            local report_rc = _rc
            _ctools_cleanup cbinscatter
            if "`genxq'" != "" capture drop `genxq'
            exit `report_rc'
        }
    }
    * As binscatter, a side of an rd() break point without observations has
    * no fit line
    else if `linetype_num' > 0 & `n_rd' > 0 {
        forvalues r = 1/`=rowsof(__cbinscatter_coefs)' {
            if missing(__cbinscatter_coefs[`r', 1]) di as error "no observations for one of the fit lines. add 'reportreg' for more info."
        }
    }

    * Display timing if requested
    if "`verbose'" != "" {
        local plugin_overhead = `elapsed' - __cbinscatter_time_total
        if `plugin_overhead' < 0 local plugin_overhead = 0
        di as text ""
        di as text "{hline 55}"
        di as text "cbinscatter timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin internals:"
        di as text "    Data load:              " as result %8.4f __cbinscatter_time_load " sec"
        di as text "    Residualize:            " as result %8.4f __cbinscatter_time_resid " sec"
        di as text "    Bin computation:        " as result %8.4f __cbinscatter_time_bins " sec"
        di as text "    Line fitting:           " as result %8.4f __cbinscatter_time_fit " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f __cbinscatter_time_total " sec"
        di as text "  {hline 53}"
        di as text "  Plugin call overhead:     " as result %8.4f `plugin_overhead' " sec"
        di as text "{hline 55}"
        di as text "    Wall clock total:       " as result %8.4f `elapsed' " sec"
        di as text "{hline 55}"

        * Display thread diagnostics
        capture local __threads_max = __cbinscatter_threads_max
        if _rc == 0 {
            capture local __openmp_enabled = __cbinscatter_openmp_enabled
            if _rc != 0 local __openmp_enabled = 0
            di as text ""
            di as text "  Thread diagnostics:"
            di as text "    OpenMP enabled:         " as result %8.0f `__openmp_enabled'
            di as text "    Max threads available:  " as result %8.0f `__threads_max'
            di as text "{hline 55}"
        }
    }

    * Display summary
    if `N_dropped' > 0 {
        di as text "(dropped " as result `N_dropped' as text " observations due to missing values)"
    }

    * Save bin data if requested. The bins are built in a frame of their own
    * (Stata 16 or later), so the data in memory are neither copied nor
    * cleared; before Stata 16 they are preserved (on disk above
    * c(max_preservemem)) and restored.
    if "`savedata'" != "" {
        if c(stata_version) >= 16 {
            tempname binframe
            frame create `binframe'
            frame `binframe': _cbinscatter_savebins `"`savedata'"'
            frame drop `binframe'
        }
        else {
            preserve
            _cbinscatter_savebins `"`savedata'"'
            restore
        }
        di as text "(bin data saved to " as result "`savedata'" as text ")"
        if `save_csv' {
            local bykeys
            forvalues g = 1/`actual_num_groups' {
                local bykeys `"`bykeys' `"`by_key`g''"'"'
            }
            local save_gropts
            if `"`title'"' != "" local save_gropts `"`save_gropts' title(`title')"'
            if `"`ytitle'"' != "" local save_gropts `"`save_gropts' ytitle(`ytitle')"'
            if `"`xtitle'"' != "" local save_gropts `"`save_gropts' xtitle(`xtitle')"'
            if `"`legend'"' != "" local save_gropts `"`save_gropts' legend(`legend')"'
            local save_gropts `"`save_gropts' `options'"'
            local save_coefs
            if `linetype_num' > 0 local save_coefs coefs(__cbinscatter_coefs)
            local save_discrete
            if "`discrete'" != "" | `discrete_used' local save_discrete discrete
            local save_connect
            if "`linetype_lower'" == "connect" local save_connect connect
            local save_by
            if "`by_orig'" != "" local save_by by(`by_orig')
            local save_xq
            if "`xq'" != "" local save_xq xq
            _cbinscatter_savecsv, file(`"`savedata'"') yvars(`depvar') xvar(`xvar') ///
                ngroups(`actual_num_groups') bins(__cbinscatter_bins) linetype(`linetype_num') ///
                nquantiles(`nquantiles') `save_xq' ///
                `save_coefs' `save_connect' `save_discrete' rd(`rd') `save_by' bykeys(`bykeys') ///
                colors(`colors') mcolors(`mcolors') lcolors(`lcolors') msymbols(`msymbols') ///
                gropts(`save_gropts') `savereplace'
        }
    }

    * Create graph if not suppressed
    local graphcmd
    if "`nograph'" == "" {
        * Default colors (same as binscatter)
        local colorlist "navy maroon forest_green dkorange teal cranberry lavender khaki sienna emidblue emerald brown erose gold bluishgray"
        if "`colors'" == "" local colors "`colorlist'"

        * binscatter's marker and line colors: one series without connect
        * marks its bins in the first color and draws its line in the
        * second; otherwise series g takes color g. mcolors() and lcolors()
        * give the lists themselves (a list too short leaves the default).
        local one_series = (`actual_num_groups' == 1 & "`linetype_lower'" != "connect")
        if `"`mcolors'"' == "" {
            if `one_series' local mcolors : word 1 of `colors'
            else local mcolors `colors'
        }
        if `"`lcolors'"' == "" {
            if `one_series' local lcolors : word 2 of `colors'
            else local lcolors `colors'
        }

        * Get bin data from matrix
        tempname binmat
        matrix `binmat' = __cbinscatter_bins
        local nrows = rowsof(`binmat')
        * rd(): one fit per interval between the break points, drawn from
        * the first bin to the first break point, between break points, and
        * from the last to the last bin, with dashed lines at the breaks
        local nint = `n_rd' + 1
        local xlines
        if `n_rd' > 0 local xlines xline(`rd', lpattern(dash) lcolor(gs8))

        * Build scatter commands using scatteri (like binscatter)
        local scatters ""
        local fits ""
        local legend_labels ""
        local legend_order ""
        forval g = 1/`actual_num_groups' {
            local coords ""
            local grp_x_min = .
            local grp_x_max = .
            * scatteri labels each point with text that follows its coordinates
            local marker_label ""
            if `"`mlabels'"' != "" {
                if `actual_num_groups' == 1 {
                    gettoken marker_label : mlabels
                    local marker_label `""`marker_label'""'
                }
                else {
                    local marker_label : word `g' of `mlabels'
                    if `"`marker_label'"' != "" local marker_label `""`marker_label'""'
                }
            }
            forval r = 1/`nrows' {
                if `actual_num_groups' > 1 & `binmat'[`r', 1] != `g' continue
                local yval = `binmat'[`r', 4]
                local xval = `binmat'[`r', 3]
                if `yval' != . & `xval' != . {
                    local coords `"`coords' `yval' `xval' `marker_label'"'
                    * binscatter draws the line from the first bin's x to the
                    * last's (not the extremes, which differ for xq() bins)
                    if `grp_x_min' == . local grp_x_min = `xval'
                    local grp_x_max = `xval'
                }
            }

            local thismcolor : word `g' of `mcolors'
            local thislcolor : word `g' of `lcolors'
            local thissymbol ""
            local symbol : word `g' of `msymbols'
            if "`symbol'" != "" local thissymbol "msymbol(`symbol')"

            local scatters `"`scatters' (scatteri `coords', mcolor(`thismcolor') lcolor(`thislcolor') `thissymbol' `connect_opt')"'
            if `actual_num_groups' > 1 {
                local legend_labels `"`legend_labels' label(`g' `by_orig'=`by_key`g'')"'
                local legend_order "`legend_order' `g'"
            }

            * Fit lines (one per rd() interval), skipping those without a fit
            if `linetype_num' > 0 & `"`coords'"' != "" {
                forval r = 1/`nint' {
                    local row = (`g' - 1) * `nint' + `r'
                    local coef_const = __cbinscatter_coefs[`row', 1]
                    local coef_linear = __cbinscatter_coefs[`row', 2]
                    if missing(`coef_const', `coef_linear') continue
                    local f "`coef_linear'*x+`coef_const'"
                    if `linetype_num' >= 2 {
                        local coef_quad = __cbinscatter_coefs[`row', 3]
                        if missing(`coef_quad') continue
                        local f "`coef_quad'*x^2+`f'"
                    }
                    if `linetype_num' == 3 {
                        local coef_cubic = __cbinscatter_coefs[`row', 4]
                        if missing(`coef_cubic') continue
                        local f "`coef_cubic'*x^3+`f'"
                    }
                    local lo = `grp_x_min'
                    local hi = `grp_x_max'
                    if `r' > 1 local lo : word `=`r' - 1' of `rd'
                    if `r' < `nint' local hi : word `r' of `rd'
                    local lo = `lo'
                    local hi = `hi'
                    local fits `"`fits' (function `f', range(`lo' `hi') lcolor(`thislcolor'))"'
                }
            }
        }

        * Legend off for one group; otherwise the default legend, with any
        * legend() options given merged after it
        if `actual_num_groups' == 1 local legend_opt "legend(off)"
        else local legend_opt `"legend(`legend_labels' order(`legend_order'))"'

        * As binscatter does, the defaults come first and the title and
        * legend options given follow them as typed, so suboptions, several
        * lines and empty titles reach twoway unchanged
        local user_gropts
        if `"`title'"' != "" local user_gropts `"`user_gropts' title(`title')"'
        if `"`ytitle'"' != "" local user_gropts `"`user_gropts' ytitle(`ytitle')"'
        if `"`xtitle'"' != "" local user_gropts `"`user_gropts' xtitle(`xtitle')"'
        if `"`legend'"' != "" local user_gropts `"`user_gropts' legend(`legend')"'

        local graphcmd `"twoway `scatters' `fits', graphregion(fcolor(white)) `xlines' ytitle(`depvar') xtitle(`xvar') `legend_opt' `user_gropts' `options'"'

        * twoway's speed depends on the data in memory, even for scatteri:
        * draw from an empty frame (Stata 16 or later) or from one
        * observation of the preserved data
        if c(stata_version) >= 16 {
            local _orig_frame = c(frame)
            capture frame drop `graph_frame'
            frame create `graph_frame'
            frame change `graph_frame'
            capture noisily `graphcmd'
            local _graph_rc = _rc
            * Always clean up frame, even on error
            frame change `_orig_frame'
            capture frame drop `graph_frame'
        }
        else {
            preserve
            qui keep in 1/1
            capture noisily `graphcmd'
            local _graph_rc = _rc
            restore
        }
        if `_graph_rc' {
            _ctools_cleanup cbinscatter
            di as error "cbinscatter: graph generation failed (rc=`_graph_rc')"
            exit `_graph_rc'
        }

        * savegraph(): graph save for .gph (or no extension), else export
        if `"`savegraph'"' != "" {
            local graphextension
            if regexm(`"`savegraph'"', "\.[a-zA-Z0-9]+$") local graphextension = regexs(0)
            if inlist(`"`graphextension'"', ".gph", "") graph save `"`savegraph'"', `savereplace'
            else graph export `"`savegraph'"', `savereplace'
        }
    }

    * binscatter's e(sample), e(y1_coefs) (each fit's coefficients, x terms
    * first, rows by#:rd#) and, for a numeric by(), e(byvalues)
    markout `touse' `depvar' `xvar' `controls' `weight_var'
    if "`absorb'" != "" markout `touse' `absorb', strok
    if "`by'" != "" markout `touse' `by'
    if "`xq'" != "" markout `touse' `xq'
    tempname ycoefs byvalues
    if `linetype_num' > 0 {
        _cbinscatter_ycoefs __cbinscatter_coefs `ycoefs', linetype(`linetype_num') ///
            ngroups(`actual_num_groups') nrd(`n_rd') xvar(`xvar') `=cond("`by'" != "", "by", "")'
    }
    if "`by_isnum'" == "1" {
        matrix `byvalues' = J(`actual_num_groups', 1, .)
        local by_rows
        forval g = 1/`actual_num_groups' {
            local rep : word `g' of `by_reps'
            matrix `byvalues'[`g', 1] = `by_orig'[`rep']
            local by_rows `by_rows' by`g'
        }
        matrix rownames `byvalues' = `by_rows'
        matrix colnames `byvalues' = `by_orig'
    }

    * Post e() results
    ereturn post, esample(`touse')
    ereturn scalar N = `N_used'
    ereturn scalar N_dropped = `N_dropped'
    ereturn scalar nquantiles = `nquantiles'
    ereturn scalar num_groups = `actual_num_groups'
    ereturn scalar discrete = ("`discrete'" != "" | `discrete_used')
    ereturn matrix bindata = __cbinscatter_bins
    if `linetype_num' > 0 {
        ereturn matrix coefs = __cbinscatter_coefs
        ereturn matrix fit_stats = __cbinscatter_fit_stats
    }
    if "`reportreg'" != "" {
        ereturn matrix reportreg_b = __cbinscatter_reportreg_b
    }
    ereturn local depvar "`depvar'"
    ereturn local xvar "`xvar'"
    if "`controls'" != "" {
        ereturn local controls "`controls_original'"
    }
    if "`absorb'" != "" {
        ereturn local absorb "`absorb'"
    }
    if "`by'" != "" {
        ereturn local by "`by_orig'"
    }
    ereturn local linetype "`linetype'"
    if `n_rd' > 0 {
        * binscatter's e(rdintervals): each interval's bounds (gt, lt_eq)
        tempname rdintervals
        local rd_col : subinstr local rd " " "\", all
        matrix `rdintervals' = (. \ `rd_col'), (`rd_col' \ .)
        local rd_names
        forvalues r = 1/`=`n_rd' + 1' {
            local rd_names `rd_names' rd`r'
        }
        matrix rownames `rdintervals' = `rd_names'
        matrix colnames `rdintervals' = gt lt_eq
        ereturn matrix rdintervals = `rdintervals'
        ereturn local rd "`rd'"
    }
    if `"`graphcmd'"' != "" ereturn local graphcmd `"`graphcmd'"'
    if `linetype_num' > 0 ereturn matrix y1_coefs = `ycoefs'
    if "`by_isnum'" == "1" ereturn matrix byvalues = `byvalues'
    ereturn local cmd "cbinscatter"
    ereturn local cmdline "cbinscatter `0'"

    _ctools_cleanup cbinscatter
end

* binscatter's e(y#_coefs) from cbinscatter's coefficients (columns: constant,
* x, x^2, x^3; one row per group and rd() interval): the x terms from the
* highest power down and _cons, rows by#:rd# (by# and rd# only with by() and
* rd(), and = without rd()), as binscatter names them.
program define _cbinscatter_ycoefs
    version 14.1
    syntax anything(name=mats), linetype(integer) ngroups(integer) nrd(integer) ///
        xvar(string) [by]
    gettoken src mats : mats
    gettoken dst : mats
    tempname row out
    local nint = `nrd' + 1
    forvalues g = 1/`ngroups' {
        forvalues r = 1/`nint' {
            local i = (`g' - 1) * `nint' + `r'
            if `linetype' == 1 matrix `row' = (`src'[`i', 2], `src'[`i', 1])
            else if `linetype' == 2 matrix `row' = (`src'[`i', 3], `src'[`i', 2], `src'[`i', 1])
            else matrix `row' = (`src'[`i', 4], `src'[`i', 3], `src'[`i', 2], `src'[`i', 1])
            if "`by'" != "" matrix roweq `row' = "by`g'"
            if `nrd' matrix rownames `row' = "rd`r'"
            else matrix rownames `row' = "="
            matrix `out' = nullmat(`out') \ `row'
        }
    }
    if `linetype' == 1 matrix colnames `out' = "`xvar'" "_cons"
    else if `linetype' == 2 matrix colnames `out' = "`xvar'^2" "`xvar'" "_cons"
    else matrix colnames `out' = "`xvar'^3" "`xvar'^2" "`xvar'" "_cons"
    matrix `dst' = `out'
end

* reportreg: the fit regressions binscatter prints, on cbinscatter's sample.
* Leaves their coefficients, one row per group, in __cbinscatter_reportreg_b.
program define _cbinscatter_reportreg
    version 14.1
    syntax, depvar(varname) xvar(varname) touse(varname) method(integer) ///
        linetype(integer) groups(integer) [weight(string) weightvar(string) ///
        controls(string) absorb(string) by(string) rd(numlist) noaddmean]
    local ry `depvar'
    local rx `xvar'
    local wt
    local wt_mean
    if "`weight'" != "" {
        local wt [`weight'=`weightvar']
        local wt_mean `wt'
        if "`weight'" == "pweight" local wt_mean [aweight=`weightvar']
    }
    * areg retains singleton FE groups, matching binscatter's residualization.
    * Encode labels for string/negative FE values; additional FEs are dummies,
    * and all of them are with iweights, which areg does not take.
    local report_absorb
    local report_fe
    foreach fe of local absorb {
        tempvar report_group
        quietly _ctools_group `fe' if `touse', generate(`report_group')
        if "`report_absorb'" == "" & "`weight'" != "iweight" local report_absorb `report_group'
        else local report_fe `report_fe' i.`report_group'
    }
    if `method' == 0 & "`controls'`absorb'" != "" {
        foreach role in ry rx {
            tempvar resid
            if "`report_absorb'" == "" {
                quietly regress ``role'' `controls' `report_fe' `wt' if `touse'
                quietly predict double `resid' if `touse', residuals
            }
            else {
                quietly areg ``role'' `controls' `report_fe' `wt' if `touse', absorb(`report_absorb')
                quietly predict double `resid' if `touse', residuals
            }
            if "`addmean'" == "" {
                quietly summarize ``role'' `wt_mean' if `touse', meanonly
                quietly replace `resid' = `resid' + r(mean) if `touse'
            }
            local `role' `resid'
        }
    }
    local fitvars `rx'
    if `linetype' >= 2 {
        tempvar quadratic
        quietly gen double `quadratic' = `rx'^2 if `touse'
        local fitvars `quadratic' `fitvars'
    }
    if `linetype' == 3 {
        tempvar cubic
        quietly gen double `cubic' = `rx'^3 if `touse'
        local fitvars `cubic' `fitvars'
    }
    tempname coefs
    local nrd : word count `rd'
    forvalues g=1/`groups' {
        local group_if
        if "`by'" != "" local group_if & `by' == `g'
        * rd(): one regression on each side of each break point, as in
        * binscatter (x <= the first, between two, above the last)
        forvalues r = 1/`=`nrd' + 1' {
            local rd_if
            local rd_text
            if `nrd' > 0 {
                local lo
                if `r' > 1 local lo : word `=`r' - 1' of `rd'
                local hi : word `r' of `rd'
                if `r' == 1 {
                    local rd_if & `rx' <= `hi'
                    local rd_text ", RD: `xvar'<=`hi'"
                }
                else if `r' <= `nrd' {
                    local rd_if & `rx' > `lo' & `rx' <= `hi'
                    local rd_text ", RD: `xvar'>`lo' & `xvar'<=`hi'"
                }
                else {
                    local rd_if & `rx' > `lo'
                    local rd_text ", RD: `xvar'>`lo'"
                }
            }
            di as text "Fit regression: `depvar' on `xvar' (group `g'`rd_text')"
            if `method' == 1 & "`report_absorb'" != "" {
                areg `ry' `fitvars' `controls' `report_fe' `wt' if `touse' `group_if', absorb(`report_absorb')
            }
            else {
                local fitcontrols
                if `method' == 1 local fitcontrols `controls' `report_fe'
                capture noisily regress `ry' `fitvars' `fitcontrols' `wt' if `touse' `group_if' `rd_if'
                * A side without observations has no line (binscatter
                * goes on); other errors stop
                if _rc == 2000 & `nrd' > 0 {
                    tempname empty
                    matrix `empty' = J(1, `linetype' + 1, .)
                    matrix colnames `empty' = `fitvars' _cons
                    matrix `coefs' = nullmat(`coefs') \ `empty'
                    continue
                }
                if _rc exit _rc
            }
            matrix `coefs' = nullmat(`coefs') \ e(b)
        }
    }
    * Name the x terms x, x^2 and x^3 rather than their temporary variables
    local names : colfullnames `coefs'
    if "`rx'" != "`xvar'" local names : subinstr local names "`rx'" "`xvar'", all
    if `linetype' >= 2 local names : subinstr local names "`quadratic'" "`xvar'^2", all
    if `linetype' == 3 local names : subinstr local names "`cubic'" "`xvar'^3", all
    matrix colnames `coefs' = `names'
    matrix __cbinscatter_reportreg_b = `coefs'
end

* Return binscatter's replace (s(replace), which lets savedata() and
* savegraph() overwrite their files) and nofastxtile (s(fastxtile); quantiles
* are always exact, so it only refuses the random-subsample options), and
* the remaining (twoway) options in s(options).
program define _cbinscatter_options, sclass
    version 14.1
    syntax [, REPLACE nofastxtile *]
    sreturn local options `"`options'"'
    sreturn local replace "`replace'"
    sreturn local fastxtile "`fastxtile'"
end

* binscatter 7.02's legacy options, with its notes: nbins() is nquantiles(),
* create_xq is genxq(q_<x>), x_q() is xq(), symbols() is msymbols(),
* linetype(noline) is none, and unique() is ignored. Where binscatter stops
* with r(0) after its message (both the new and the legacy option given),
* this exits with r(198). Sets the mapped locals in the caller.
program define _cbinscatter_legacy
    version 14.1
    syntax, xvar(string) nquantiles(integer) nbins(integer) [genxq(string) create_xq ///
        xq(string) x_q(string) msymbols(string) symbols(string) linetype(string) unique(string)]
    if `nbins' != 20 {
        if `nquantiles' != 20 {
            di as error "Cannot specify both nquantiles() and nbins(): both are the same option, nbins is supported only for backward compatibility."
            exit 198
        }
        di as text "NOTE: legacy option nbins() has been renamed nquantiles(), and is supported only for backward compatibility."
        c_local nquantiles `nbins'
    }
    if "`create_xq'" != "" {
        if "`genxq'" != "" {
            di as error "Cannot specify both genxq() and create_xq: both are the same option, create_xq is supported only for backward compatibility."
            exit 198
        }
        di as text "NOTE: legacy option create_xq has been renamed genxq(), and is supported only for backward compatibility."
        c_local genxq q_`xvar'
    }
    if "`x_q'" != "" {
        if "`xq'" != "" {
            di as error "Cannot specify both xq() and x_q(): both are the same option, x_q() is supported only for backward compatibility."
            exit 198
        }
        di as text "NOTE: legacy option x_q() has been renamed xq(), and is supported only for backward compatibility."
        c_local xq `x_q'
    }
    if `"`symbols'"' != "" {
        if `"`msymbols'"' != "" {
            di as error "Cannot specify both msymbols() and symbols(): both are the same option, symbols() is supported only for backward compatibility."
            exit 198
        }
        di as text "NOTE: legacy option symbols() has been renamed msymbols(), and is supported only for backward compatibility."
        c_local msymbols `symbols'
    }
    if "`linetype'" == "noline" {
        di as text "NOTE: legacy line type 'noline' has been renamed 'none', and is supported only for backward compatibility."
        c_local linetype none
    }
    if `"`unique'"' != "" {
        di as text "NOTE: unique() is no longer a recognized option, and will be ignored. binscatter now considers the x-variable discrete if it has fewer unique values than nquantiles()"
    }
end

* Check twoway options without drawing a graph: a temporary graph of one
* point in an empty frame (Stata 16 or later), dropped at once. name() and
* saving() are left out, so nothing is stored or written.
program define _cbinscatter_check_gropts
    version 14.1
    syntax [, NAME(string asis) SAVing(string asis) *]
    if `"`options'"' == "" exit
    tempname check
    if c(stata_version) >= 16 {
        tempname frame
        frame create `frame'
        capture noisily frame `frame': twoway scatteri 0 0, nodraw `options' name(`check')
        local rc = _rc
        frame drop `frame'
    }
    else {
        capture noisily twoway scatteri 0 0, nodraw `options' name(`check')
        local rc = _rc
    }
    capture graph drop `check'
    if `rc' {
        di as error "cbinscatter: invalid graph option(s)"
        exit `rc'
    }
end

* Replace the current dataset with the bins in matrix __cbinscatter_bins and
* save them.
program define _cbinscatter_savebins
    version 14.1
    args savedata
    quietly {
        clear
        svmat __cbinscatter_bins, names(col)
        rename c1 by_group
        rename c2 bin_id
        rename c3 x_mean
        rename c4 y_mean
        rename c5 n_obs
        * Drop empty rows
        drop if missing(bin_id)
        save `"`savedata'"', replace
    }
end

* binscatter's savedata(): <file>.csv holds the bin coordinates (the x of
* each by() group, then each y for each group, one row per bin up to
* binscatter's number of bins) and <file>.do reads it, labels the variables
* and redraws the graph with binscatter's twoway command. bins() and coefs()
* name one matrix per y in cbinscatter's layout (e(bindata), e(coefs));
* bykeys() holds each group's value or value label, as compound-quoted words.
program define _cbinscatter_savecsv
    version 14.1
    syntax, file(string) yvars(string) xvar(string) ngroups(integer) bins(string) ///
        linetype(integer) nquantiles(integer) [coefs(string) connect discrete xq rd(numlist) by(string) ///
        bykeys(string asis) colors(string) mcolors(string) lcolors(string) ///
        msymbols(string) gropts(string asis) replace]
    local ynum : word count `yvars'
    local bynum = `ngroups'
    local x_var `xvar'
    local y_vars `yvars'
    forvalues g = 1/`bynum' {
        gettoken key bykeys : bykeys
        local byvalname`g' `"`key'"'
    }
    * binscatter's scatter points: x and y of each group side by side, the
    * groups padded with missing values to the most bins of one
    local nq = 0
    forvalues j = 1/`ynum' {
        tempname pts`j'
        local m : word `j' of `bins'
        mata: _cbinscatter_pts("`m'", `bynum', "`pts`j''", "`discrete'" != "")
        if `j' == 1 local nq = r_cbin_nq
    }
    capture scalar drop r_cbin_nq
    * binscatter's number of bins: the distinct x values (discrete) or xq()
    * values; for quantile bins the bins formed, or with by() the
    * nquantiles() requested
    if "`by'" != "" & "`discrete'`xq'" == "" local nq = `nquantiles'
    local nq = max(`nq', rowsof(`pts1'))

    * The CSV, as binscatter writes it
    tempname fh
    file open `fh' using `"`file'.csv"', write text `replace'
    forvalues row = 0/`nq' {
        forvalues counter_by = 1/`bynum' {
            if `row' == 0 {
                if "`by'" != "" local bynlabel _by`counter_by'
                file write `fh' "`x_var'`bynlabel',"
            }
            else {
                if `row' <= `=rowsof(`pts1')' file write `fh' (`pts1'[`row', `counter_by'*2-1]) ","
                else file write `fh' ".,"
            }
        }
        local counter_depvar = 0
        foreach depvar of local y_vars {
            local ++counter_depvar
            forvalues counter_by = 1/`bynum' {
                if `row' == 0 {
                    if "`by'" != "" local bynlabel _by`counter_by'
                    file write `fh' "`depvar'`bynlabel'"
                }
                else {
                    if `row' <= `=rowsof(`pts`counter_depvar'')' file write `fh' (`pts`counter_depvar''[`row', `counter_by'*2])
                    else file write `fh' "."
                }
                if !(`counter_depvar' == `ynum' & `counter_by' == `bynum') file write `fh' ","
            }
        }
        file write `fh' _n
    }
    file close `fh'
    di as text `"(file `file'.csv written containing saved data)"'

    * binscatter's twoway command for the saved variables
    if `"`colors'"' == "" local colors navy maroon forest_green dkorange teal cranberry lavender ///
        khaki sienna emidblue emerald brown erose gold bluishgray
    if `"`mcolors'"' == "" {
        if (`ynum' == 1 & `bynum' == 1 & "`connect'" == "") local mcolors `: word 1 of `colors''
        else local mcolors `colors'
    }
    if `"`lcolors'"' == "" {
        if (`ynum' == 1 & `bynum' == 1 & "`connect'" == "") local lcolors `: word 2 of `colors''
        else local lcolors `colors'
    }
    local connect_opt
    if "`connect'" != "" local connect_opt "c(l)"
    local symbol_prefix
    local symbol_suffix
    if "`msymbols'" != "" {
        local symbol_prefix "msymbol("
        local symbol_suffix ")"
    }
    local c = 0
    local counter_series = 0
    local savedata_scatters
    local legend_labels
    local order
    forvalues counter_by = 1/`bynum' {
        local xind = `counter_by'*2-1
        local yind = `counter_by'*2
        local counter_depvar = 0
        foreach depvar of local y_vars {
            local ++counter_depvar
            local ++c
            local xval = `pts`counter_depvar''[1, `xind']
            local yval = `pts`counter_depvar''[1, `yind']
            if missing(`xval', `yval') continue
            local ++counter_series
            if ("`by'" == "") local savedata_scatters `savedata_scatters' (scatter `depvar' `x_var'
            else local savedata_scatters `savedata_scatters' (scatter `depvar'_by`counter_by' `x_var'_by`counter_by'
            local scatter_options `connect_opt' mcolor(`: word `c' of `mcolors'') lcolor(`: word `c' of `lcolors'') `symbol_prefix'`: word `c' of `msymbols''`symbol_suffix'
            local savedata_scatters `savedata_scatters', `scatter_options')
            if "`by'" == "" {
                if (`ynum' == 1) local legend_labels off
                else local legend_labels `legend_labels' lab(`counter_series' `depvar')
            }
            else {
                if (`ynum' == 1) local legend_labels `legend_labels' lab(`counter_series' `by'=`byvalname`counter_by'')
                else local legend_labels `legend_labels' lab(`counter_series' `depvar': `by'=`byvalname`counter_by'')
            }
            if ("`by'" != "" | `ynum' > 1) local order `order' `counter_series'
        }
    }
    * Fit lines, from the first to the last bin of the group (and between
    * the rd() break points)
    local fits
    if `linetype' > 0 & "`coefs'" != "" {
        local c = 0
        local rdnum = wordcount("`rd'") + 1
        forvalues counter_by = 1/`bynum' {
            local xind = `counter_by'*2-1
            local row0 = `rdnum' * (`counter_by' - 1)
            local counter_depvar = 0
            foreach depvar of local y_vars {
                local ++counter_depvar
                local ++c
                local cm : word `counter_depvar' of `coefs'
                local lb = `pts`counter_depvar''[1, `xind']
                local ub = .
                forvalues r = `=rowsof(`pts`counter_depvar'')'(-1)1 {
                    if `ub' == . local ub = `pts`counter_depvar''[`r', `xind']
                }
                local bounds `lb' `rd' `ub'
                forvalues counter_rd = 1/`rdnum' {
                    local i = `row0' + `counter_rd'
                    local coef_cubic = 0
                    local coef_quad = 0
                    local coef_lin = `cm'[`i', 2]
                    local coef_cons = `cm'[`i', 1]
                    if `linetype' >= 2 local coef_quad = `cm'[`i', 3]
                    if `linetype' == 3 local coef_cubic = `cm'[`i', 4]
                    if !missing(`coef_cubic', `coef_quad', `coef_lin', `coef_cons') {
                        local leftbound : word `counter_rd' of `bounds'
                        local rightbound : word `=`counter_rd'+1' of `bounds'
                        local f `coef_quad'*x^2+`coef_lin'*x+`coef_cons'
                        if `linetype' == 3 local f `coef_cubic'*x^3+`f'
                        local fits `fits' (function `f', range(`leftbound' `rightbound') lcolor(`: word `c' of `lcolors''))
                    }
                }
            }
        }
    }
    local xlines
    if "`rd'" != "" local xlines xline(`rd', lpattern(dash) lcolor(gs8))
    if (`ynum' == 1) local ytitle `y_vars'
    else if (`ynum' == 2) local ytitle : subinstr local y_vars " " " and "
    else local ytitle : subinstr local y_vars " " "; ", all
    mata: st_local("gropts", strtrim(st_local("gropts")))
    local savedata_graphcmd twoway `savedata_scatters' `fits', graphregion(fcolor(white)) `xlines' xtitle(`x_var') ytitle(`ytitle') legend(`legend_labels' order(`order')) `gropts'

    * The do-file, as binscatter writes it
    file open `fh' using `"`file'.do"', write text `replace'
    file write `fh' `"insheet using `file'.csv"' _n _n
    if "`by'" != "" {
        foreach var in `x_var' `y_vars' {
            forvalues counter_by = 1/`bynum' {
                file write `fh' `"label variable `var'_by`counter_by' "`var'; `by'==`byvalname`counter_by''""' _n
            }
        }
        file write `fh' _n
    }
    file write `fh' `"`savedata_graphcmd'"' _n
    file close `fh'
    di as text `"(file `file'.do written containing commands to process saved data)"'
end

capture mata: mata drop _cbinscatter_pts()
mata:
// binscatter's scatter points of one y from cbinscatter's bins (columns
// group, bin, x, y, n): x and y of group g in columns 2g-1 and 2g, padded
// with missing values. r_cbin_nq is binscatter's number of bins: the
// distinct bins overall, or the distinct x values for discrete bins.
void _cbinscatter_pts(string scalar src, real scalar ng, string scalar dst,
    real scalar discrete)
{
    real matrix B, P, G
    real scalar g, n, nr

    B = st_matrix(src)
    B = select(B, B[., 2] :< .)
    nr = 0
    for (g = 1; g <= ng; g++) {
        n = ng > 1 ? sum(B[., 1] :== g) : rows(B)
        if (n > nr) nr = n
    }
    P = J(max((nr, 1)), 2 * ng, .)
    for (g = 1; g <= ng; g++) {
        G = ng > 1 ? select(B, B[., 1] :== g) : B
        if (rows(G)) P[|1, 2*g-1 \ rows(G), 2*g|] = G[., (3, 4)]
    }
    st_matrix(dst, P)
    if (!rows(B)) st_numscalar("r_cbin_nq", 0)
    else st_numscalar("r_cbin_nq", rows(uniqrows(B[., discrete ? 3 : 2])))
}
end

* Several y variables (binscatter y1 ... yk x). As in binscatter, the sample is
* common to all of them, so the bins of x are the same; each y gets its own
* binned means and fit. Each y runs through cbinscatter on that sample
* without a graph, and the series are drawn together: for each by-group,
* one per y, colors in that order, the legend naming the y (and the group),
* and the y titles joined as binscatter joins them.
program define _cbinscatter_multi, eclass
    version 14.1
    syntax varlist(min=3 numeric) [aw fw pw iw] [if] [in], ///
        [Nquantiles(integer 20) Controls(varlist fv ts) Absorb(varlist) BY(varname) ///
         LINEtype(string) DISCrete GENxq(name) SAVEdata(string) REPLACE ///
         REPORTreg NOgraph Verbose METHod(passthru) THReads(passthru) ///
         TITLE(string asis) YTItle(string asis) XTItle(string asis) LEGEND(string asis) ///
         COLors(string) Msymbols(string) MLAbels(string asis) ///
         MEDians noAddmean rd(numlist ascending) xq(varname numeric) ///
         MColors(string) LColors(string) savegraph(string) nofastxtile ///
         randvar(varname numeric) randcut(real 1) randn(integer -1) ///
         nbins(integer 20) create_xq x_q(varname numeric) symbols(string) unique(string) *]

    local nvars : word count `varlist'
    local xvar : word `nvars' of `varlist'
    local yvars : list varlist - xvar
    local ynum : word count `yvars'

    * binscatter's legacy options, once for all the y
    _cbinscatter_legacy, xvar(`xvar') nquantiles(`nquantiles') nbins(`nbins') ///
        genxq(`genxq') `create_xq' xq(`xq') x_q(`x_q') msymbols(`msymbols') ///
        symbols(`symbols') linetype(`linetype') unique(`unique')

    * The common sample: the qualifiers, every y and x, and what each run
    * marks out itself (controls, absorbed and by variables, weights)
    marksample touse
    if "`controls'" != "" {
        fvrevar `controls' if `touse'
        markout `touse' `r(varlist)'
    }
    if "`absorb'" != "" markout `touse' `absorb', strok
    if "`by'" != "" markout `touse' `by', strok
    if "`xq'" != "" markout `touse' `xq'
    local wgt
    if "`weight'" != "" {
        local wgt [`weight'`exp']
        tempvar wcheck
        quietly generate double `wcheck' = `=substr(`"`exp'"', 2, .)' if `touse'
        markout `touse' `wcheck'
    }
    quietly count if `touse'
    if r(N) == 0 {
        di as error "no observations"
        exit 3301
    }

    * savedata() holds every y's bins: as binscatter refuses an existing
    * file without replace. Without an extension, binscatter's .csv and .do
    * are written too (method(classic)).
    local save_csv = 0
    if `"`savedata'"' != "" {
        local __base = subinstr(`"`savedata'"', "\", "/", .)
        local __cut = strpos(reverse(`"`__base'"'), "/")
        if `__cut' local __base = substr(`"`__base'"', strlen(`"`__base'"') - `__cut' + 2, .)
        local save_csv = !strpos(`"`__base'"', ".") & !regexm(lower(`"`method'"'), "\((binsreg|b|conditional)\)$")
    }
    if `"`savedata'"' != "" & "`replace'" == "" {
        local __savefile `"`savedata'"'
        if !strpos(`"`__base'"', ".") local __savefile `"`__savefile'.dta"'
        if `save_csv' {
            confirm new file `"`savedata'.csv"'
            confirm new file `"`savedata'.do"'
        }
        capture confirm new file `"`__savefile'"'
        if _rc == 602 {
            di as error `"file `__savefile' already exists; specify replace"'
            exit 602
        }
    }

    * savegraph() without replace refuses an existing file, as binscatter does
    if `"`savegraph'"' != "" {
        if "`nograph'" != "" {
            di as error "cbinscatter: savegraph() requires the graph (nograph specified)"
            exit 198
        }
        if "`replace'" == "" {
            if regexm(`"`savegraph'"', "\.[a-zA-Z0-9]+$") confirm new file `"`savegraph'"'
            else confirm new file `"`savegraph'.gph"'
        }
    }

    local copts nquantiles(`nquantiles') `discrete' `reportreg' `verbose' `method' `threads'
    if "`linetype'" != "" local copts `copts' linetype(`linetype')
    local copts `copts' `medians' `addmean' `fastxtile'
    if "`controls'" != "" local copts `copts' controls(`controls')
    if "`absorb'" != "" local copts `copts' absorb(`absorb')
    if "`by'" != "" local copts `copts' by(`by')
    if "`rd'" != "" local copts `copts' rd(`rd')
    if "`xq'" != "" local copts `copts' xq(`xq')
    if "`randvar'" != "" local copts `copts' randvar(`randvar')
    if `randcut' != 1 local copts `copts' randcut(`randcut')
    if `randn' != -1 local copts `copts' randn(`randn')
    * binscatter bins x once for all the y: genxq() comes from the first
    * run, and every run draws the same random subsample (the generator
    * starts each run where it started the first)
    local redraw = (`randcut' != 1 | `randn' != -1) & "`randvar'" == ""
    if `redraw' local rand_state `"`c(rngstate)'"'

    local k = 0
    foreach y of local yvars {
        local ++k
        if `k' > 1 & `redraw' set rngstate `rand_state'
        local first_opts
        if `k' == 1 & "`genxq'" != "" local first_opts genxq(`genxq')
        * binscatter shows its notes once, not once per y
        if `k' == 2 global CBINSCATTER_NOTES_SHOWN 1
        capture noisily cbinscatter `y' `xvar' `wgt' if `touse', nograph `copts' `first_opts'
        local run_rc = _rc
        if `k' == `ynum' | `run_rc' global CBINSCATTER_NOTES_SHOWN
        if `run_rc' {
            if `k' > 1 & "`genxq'" != "" capture drop `genxq'
            exit `run_rc'
        }
        tempname bins`k' coefs`k'
        matrix `bins`k'' = e(bindata)
        local has_fit = 0
        capture confirm matrix e(coefs)
        if !_rc {
            matrix `coefs`k'' = e(coefs)
            local has_fit = 1
        }
        tempname ycoefs`k'
        local has_ycoefs = 0
        capture confirm matrix e(y1_coefs)
        if !_rc {
            matrix `ycoefs`k'' = e(y1_coefs)
            local has_ycoefs = 1
        }
        if `k' == 1 {
            tempname byvalues
            capture confirm matrix e(byvalues)
            local has_byvalues = !_rc
            if `has_byvalues' matrix `byvalues' = e(byvalues)
            local discrete_used = e(discrete)
        }
        local N = e(N)
        local ngroups = e(num_groups)
        local nq = e(nquantiles)
        local ltype `e(linetype)'
        if "`rd'" != "" {
            tempname rdintervals
            matrix `rdintervals' = e(rdintervals)
        }
    }

    * by() legend keys, as for one y: byvar=value (or its value label)
    if "`by'" != "" {
        tempvar bygroup
        quietly _ctools_group `by' if `touse', generate(`bygroup') label
        local by_reps `"`r(reps)'"'
        capture confirm numeric variable `by'
        local by_isnum = !_rc
        local by_vlabel
        if `by_isnum' local by_vlabel : value label `by'
        forval g = 1/`ngroups' {
            if `by_isnum' {
                local rep : word `g' of `by_reps'
                local byval = `by'[`rep']
                if "`by_vlabel'" == "" local by_key`g' = `byval'
                else local by_key`g' : label `by_vlabel' `byval'
            }
            else local by_key`g' : label (`bygroup') `g'
        }
    }

    if `"`savedata'"' != "" {
        tempname allbins
        forval j = 1/`ynum' {
            tempname bj
            matrix `bj' = J(rowsof(`bins`j''), 1, `j'), `bins`j''
            if `j' == 1 matrix `allbins' = `bj'
            else matrix `allbins' = `allbins' \ `bj'
        }
        preserve
        quietly {
            clear
            svmat `allbins'
            rename (`allbins'1 `allbins'2 `allbins'3 `allbins'4 `allbins'5 `allbins'6) ///
                (y_index by_group bin_id x_mean y_mean n_obs)
            drop if missing(bin_id)
            generate str32 depvar = ""
            forval j = 1/`ynum' {
                replace depvar = "`: word `j' of `yvars''" if y_index == `j'
            }
            order depvar, after(y_index)
            save `"`savedata'"', replace
        }
        restore
        di as text "(bin data saved to " as result `"`savedata'"' as text ")"
        if `save_csv' {
            local bykeys
            local coefs_list
            local bins_list
            forval g = 1/`ngroups' {
                local bykeys `"`bykeys' `"`by_key`g''"'"'
            }
            forval j = 1/`ynum' {
                local bins_list `bins_list' `bins`j''
                if `has_fit' local coefs_list `coefs_list' `coefs`j''
            }
            local lt = lower("`ltype'")
            local lt_num = cond(inlist("`lt'", "none", "n", "connect"), 0, ///
                cond(inlist("`lt'", "quadratic", "qfit"), 2, cond("`lt'" == "cubic", 3, 1)))
            local save_coefs
            if `has_fit' & `lt_num' > 0 local save_coefs coefs(`coefs_list')
            local save_discrete
            if `discrete_used' local save_discrete discrete
            local save_connect
            if "`lt'" == "connect" local save_connect connect
            local save_by
            if "`by'" != "" local save_by by(`by')
            local save_gropts
            if `"`title'"' != "" local save_gropts `"`save_gropts' title(`title')"'
            if `"`ytitle'"' != "" local save_gropts `"`save_gropts' ytitle(`ytitle')"'
            if `"`xtitle'"' != "" local save_gropts `"`save_gropts' xtitle(`xtitle')"'
            if `"`legend'"' != "" local save_gropts `"`save_gropts' legend(`legend')"'
            local save_gropts `"`save_gropts' `options'"'
            local save_xq
            if "`xq'" != "" local save_xq xq
            _cbinscatter_savecsv, file(`"`savedata'"') yvars(`yvars') xvar(`xvar') ///
                ngroups(`ngroups') bins(`bins_list') linetype(`lt_num') ///
                nquantiles(`nq') `save_xq' ///
                `save_coefs' `save_connect' `save_discrete' rd(`rd') `save_by' bykeys(`bykeys') ///
                colors(`colors') mcolors(`mcolors') lcolors(`lcolors') msymbols(`msymbols') ///
                gropts(`save_gropts') `replace'
        }
    }

    if "`nograph'" == "" {
        local colorlist "navy maroon forest_green dkorange teal cranberry lavender khaki sienna emidblue emerald brown erose gold bluishgray"
        if "`colors'" == "" local colors "`colorlist'"
        * With several y every series takes its own color, as in binscatter
        if `"`mcolors'"' == "" local mcolors `colors'
        if `"`lcolors'"' == "" local lcolors `colors'
        local connect_opt
        if "`ltype'" == "connect" local connect_opt "connect(l)"
        local n_rd : word count `rd'
        local nint = `n_rd' + 1
        local xlines
        if `n_rd' > 0 local xlines xline(`rd', lpattern(dash) lcolor(gs8))
        local scatters
        local fits
        local legend_labels
        local legend_order
        local c = 0
        local series = 0
        forval g = 1/`ngroups' {
            forval j = 1/`ynum' {
                local ++c
                local yname : word `j' of `yvars'
                local coords
                local xmin = .
                local xmax = .
                forval r = 1/`=rowsof(`bins`j'')' {
                    if `bins`j''[`r', 1] != `g' continue
                    local yv = `bins`j''[`r', 4]
                    local xv = `bins`j''[`r', 3]
                    if missing(`yv', `xv') continue
                    local coords `coords' `yv' `xv'
                    * The line runs from the first bin's x to the last's
                    if missing(`xmin') local xmin = `xv'
                    local xmax = `xv'
                }
                if "`coords'" == "" continue
                local ++series
                local mc : word `c' of `mcolors'
                local lc : word `c' of `lcolors'
                local sym : word `c' of `msymbols'
                local symopt
                if "`sym'" != "" local symopt msymbol(`sym')
                local scatters `"`scatters' (scatteri `coords', mcolor(`mc') lcolor(`lc') `symopt' `connect_opt')"'
                if "`by'" == "" local legend_labels `"`legend_labels' label(`series' `yname')"'
                else local legend_labels `"`legend_labels' label(`series' `yname': `by'=`by_key`g'')"'
                local legend_order `legend_order' `series'
                if `has_fit' & "`ltype'" != "connect" & "`ltype'" != "none" {
                    forval r = 1/`nint' {
                        local row = (`g' - 1) * `nint' + `r'
                        local b0 = `coefs`j''[`row', 1]
                        local b1 = `coefs`j''[`row', 2]
                        if missing(`b0', `b1') continue
                        local f "`b1'*x+`b0'"
                        if colsof(`coefs`j'') >= 3 {
                            local b2 = `coefs`j''[`row', 3]
                            if !missing(`b2') local f "`b2'*x^2+`f'"
                        }
                        if colsof(`coefs`j'') >= 4 {
                            local b3 = `coefs`j''[`row', 4]
                            if !missing(`b3') local f "`b3'*x^3+`f'"
                        }
                        local lo = `xmin'
                        local hi = `xmax'
                        if `r' > 1 local lo : word `=`r' - 1' of `rd'
                        if `r' < `nint' local hi : word `r' of `rd'
                        local lo = `lo'
                        local hi = `hi'
                        local fits `"`fits' (function `f', range(`lo' `hi') lcolor(`lc'))"'
                    }
                }
            }
        }
        if `ynum' == 2 local ytitle_def : subinstr local yvars " " " and "
        else local ytitle_def : subinstr local yvars " " "; ", all
        local user_gropts
        if `"`title'"' != "" local user_gropts `"`user_gropts' title(`title')"'
        if `"`ytitle'"' != "" local user_gropts `"`user_gropts' ytitle(`ytitle')"'
        if `"`xtitle'"' != "" local user_gropts `"`user_gropts' xtitle(`xtitle')"'
        if `"`legend'"' != "" local user_gropts `"`user_gropts' legend(`legend')"'
        local graphcmd `"twoway `scatters' `fits', graphregion(fcolor(white)) `xlines' ytitle(`ytitle_def') xtitle(`xvar') legend(`legend_labels' order(`legend_order')) `user_gropts' `options'"'
        if c(stata_version) >= 16 {
            tempname gframe
            local oframe = c(frame)
            frame create `gframe'
            frame change `gframe'
            capture noisily `graphcmd'
            local grc = _rc
            frame change `oframe'
            frame drop `gframe'
        }
        else {
            preserve
            quietly keep in 1/1
            capture noisily `graphcmd'
            local grc = _rc
            restore
        }
        if `grc' {
            di as error "cbinscatter: graph generation failed (rc=`grc')"
            exit `grc'
        }
        if `"`savegraph'"' != "" {
            local graphextension
            if regexm(`"`savegraph'"', "\.[a-zA-Z0-9]+$") local graphextension = regexs(0)
            if inlist(`"`graphextension'"', ".gph", "") graph save `"`savegraph'"', `replace'
            else graph export `"`savegraph'"', `replace'
        }
    }

    * binscatter's e(sample), e(y#_coefs) and e(byvalues)
    ereturn post, esample(`touse')
    ereturn scalar N = `N'
    ereturn scalar nquantiles = `nq'
    ereturn scalar num_groups = `ngroups'
    ereturn scalar num_depvars = `ynum'
    ereturn scalar discrete = `discrete_used'
    forval j = 1/`ynum' {
        ereturn matrix bindata_y`j' = `bins`j''
        if `has_fit' ereturn matrix coefs_y`j' = `coefs`j''
    }
    if `has_ycoefs' {
        forval j = `ynum'(-1)1 {
            ereturn matrix y`j'_coefs = `ycoefs`j''
        }
    }
    if `has_byvalues' ereturn matrix byvalues = `byvalues'
    ereturn local depvar "`yvars'"
    ereturn local xvar "`xvar'"
    if "`controls'" != "" ereturn local controls "`controls'"
    if "`absorb'" != "" ereturn local absorb "`absorb'"
    if "`by'" != "" ereturn local by "`by'"
    ereturn local linetype "`ltype'"
    if "`rd'" != "" {
        ereturn matrix rdintervals = `rdintervals'
        ereturn local rd "`rd'"
    }
    if "`nograph'" == "" ereturn local graphcmd `"`graphcmd'"'
    ereturn local cmd "cbinscatter"
    ereturn local cmdline `"cbinscatter `0'"'
end
