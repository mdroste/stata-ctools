*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! C-accelerated reshape between wide and long form (ctools)
*! original: Stata's reshape (reshape.ado 4.7.1, October 4, 2023)
*! original: see help creshape for details

program define creshape, rclass
    version 14.1

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    u_mi_not_mi_set creshape other

    gettoken sub : 0, parse(" ,")
    if "`sub'" == "long" | "`sub'" == "wide" {
        gettoken sub 0 : 0, parse(" ,")
        _creshape_impl `sub' `0'
        return add
        exit
    }

    * Everything else declares or reports the settings (creshape i, j, xij,
    * xi, favor, clear, query) or lists problem observations (creshape
    * error). The settings are characteristics of the dataset that reshape
    * and creshape share, so reshape handles these subcommands.
    reshape `0'
end

* creshape long and creshape wide: reshape's basic syntax (stubs and options)
* or its advanced one (no arguments; the settings are already declared).
program define _creshape_impl, rclass
    version 14.1
    gettoken dir 0 : 0
    syntax [anything(name=stubs id="stub names")] [, I(string) J(string) ATwl(string) ///
        String FAVor(string) Verbose THReads(integer 0)]

    if `threads' < 0 {
        di as error "threads() must be nonnegative"
        exit 198
    }
    local timing = ("`verbose'" != "")
    if `timing' quietly _ctools_clock __ct_start
    if c(stata_version) >= 18 {
        _alias_not_allowed
    }

    * =========================================================================
    * Basic syntax: declare the settings, as reshape does
    * =========================================================================
    if `"`stubs'"' != "" {
        if `"`i'"' == "" {
            di as error "option {bf:i()} required"
            _creshape_picture
            exit 198
        }
        unab i : `i'
        gettoken jvar jvals : j
        if "`jvar'" == "" local jvar "_j"
        local jvals `jvals'
        if `"`favor'"' != "" {
            if "`favor'" == "space" local favor "memory"
            if !inlist("`favor'", "memory", "speed") {
                di as error "{p}invalid {bf:favor} method;{break}option {bf:favor()} requires argument {bf:speed} or {bf:memory}{p_end}"
                exit 198
            }
        }
        if "`dir'" == "wide" {
            * the stubs and j are variables of the long data
            capture unab stubs : `stubs'
            capture unab jvar : `jvar'
            if _rc {
                if "`jvar'" == "_j" {
                    di as error "option {bf:j()} required"
                    _creshape_picture
                    exit 198
                }
                di as error "variable {bf:`jvar'} not found"
                di as error "{p 4 4 2}Data are already wide.{p_end}"
                exit 111
            }
        }
        else {
            capture confirm new variable `jvar'
            if _rc {
                if _rc == 110 {
                    di as error "variable {bf:`jvar'} already exists"
                    di as error "{p 4 4 2}Data are already long.{p_end}"
                    exit 110
                }
                confirm new variable `jvar'
            }
        }
        if "`string'" != "" local string ", string"
        if `"`atwl'"' != "" local atwl `", atwl(`atwl')"'
        reshape clear
        reshape i `i'
        reshape j `jvar' `jvals' `string'
        reshape xij `stubs' `atwl'
        if "`favor'" != "" char _dta[ReS_favor] "`favor'"
    }
    else if `"`i'`j'`atwl'`string'`favor'"' != "" {
        di as error "options i(), j(), atwl(), string and favor() require the names of the variables to reshape"
        _creshape_picture
        exit 198
    }

    * =========================================================================
    * The settings
    * =========================================================================
    local __creshape_i : char _dta[ReS_i]
    local __creshape_j : char _dta[ReS_j]
    local __creshape_jv : char _dta[ReS_jv]
    local __creshape_xij : char _dta[ReS_Xij]
    local __creshape_xi : char _dta[Res_Xi]
    local __creshape_atwl : char _dta[ReS_atwl]
    local isstr : char _dta[ReS_str]
    local ver : char _dta[ReS_ver]

    * Settings that are incomplete, or were declared in the syntax of Stata
    * 10 and earlier (reshape groups, vars, cons), are left to reshape.
    if "`__creshape_i'" == "" | "`__creshape_j'" == "" | "`__creshape_xij'" == "" | "`ver'" == "v.1" {
        reshape `dir'
        exit
    }
    if "`isstr'" == "" local isstr 0

    if "`dir'" == "wide" {
        capture confirm variable `__creshape_j', exact
        if _rc {
            di as text "(data are already wide)"
            exit
        }
        confirm variable `__creshape_j' `__creshape_i' `__creshape_xi'
    }
    else {
        confirm variable `__creshape_i' `__creshape_xi'
        capture confirm new variable `__creshape_j'
        if _rc {
            di as text "(data are already long)"
            exit
        }
    }
    local xilist = ("`__creshape_xi'" != "")

    * =========================================================================
    * The plugin stages the data and checks the request; nothing has changed
    * in Stata until it returns
    * =========================================================================
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    unab __creshape_vars : *
    _ctools_strw `__creshape_vars'
    local threadopt
    if `threads' > 0 local threadopt "threads(`threads')"

    local N0 = _N
    local k0 = c(k)
    local changed0 = c(changed)
    local fn0 `"$S_FN"'
    local fndate0 `"$S_FNDATE"'
    local dropped 0
    local saved 0

    if `timing' {
        quietly _ctools_clock __ct_now
        local __t_parse = `__ct_now' - `__ct_start'
        local __ct_start = `__ct_now'
    }

    * Break is held from here: between the plan and the store the plugin has
    * the only copy of the observations.
    local fallback 0
    nobreak {
        capture noisily plugin call ctools_plugin `__creshape_vars', ///
            "creshape `threadopt' plan `dir' string=`isstr' report=1 xilist=`xilist' maxvar=`c(maxvar)' verbose=`timing'"
        local rc = _rc
        if !`rc' & "`__creshape_fallback'" != "" local fallback 1
        if `rc' | `fallback' capture plugin call ctools_plugin, "creshape clear"
        if `timing' {
            quietly _ctools_clock __ct_now
            local __t_plan = `__ct_now' - `__ct_start'
            local __ct_start = `__ct_now'
        }

        if !`rc' & !`fallback' {
            * The description of the result: Stata saves the dataset without
            * its observations, the plugin rewrites that description, Stata
            * loads it.
            capture noisily {
                tempfile __creshape_template __creshape_schema
                if _N quietly drop in 1/l
                local dropped 1
                * (CTOOLS_creshape_fault: validation hook, a failure at a
                * named stage after the data have started to change)
                if "$CTOOLS_creshape_fault" == "template" error 1
                quietly save `"`__creshape_template'"', emptyok orphans
                local saved 1
                plugin call ctools_plugin, "creshape schema"
                if "`__creshape_rc'" != "0" exit `__creshape_rc'
                if "$CTOOLS_creshape_fault" == "schema" error 1
                quietly use `"`__creshape_schema'"', clear
                if `__creshape_N' quietly set obs `__creshape_N'
                if "$CTOOLS_creshape_fault" == "store" error 1
                if `timing' {
                    quietly _ctools_clock __ct_now
                    local __t_schema = `__ct_now' - `__ct_start'
                    local __ct_start = `__ct_now'
                }
                plugin call ctools_plugin *, "creshape `threadopt' store"
                if "`__creshape_rc'" != "0" exit `__creshape_rc'
            }
            local rc = _rc
            if `rc' {
                * Put the data back: the template is the dataset as it was,
                * and the plugin still holds its observations.
                if `saved' capture quietly use `"`__creshape_template'"', clear
                if `dropped' & `N0' capture quietly set obs `N0'
                capture noisily plugin call ctools_plugin *, "creshape restore"
                if _rc di as error "creshape: the original data could not be restored"
                global S_FN `"`fn0'"'
                global S_FNDATE `"`fndate0'"'
                mata: st_updata(`changed0')
            }
            else {
                * as reshape leaves them: the data in memory are no longer a
                * file's
                global S_FN
                global S_FNDATE
            }
        }
    }
    if `rc' {
        _ctools_cleanup creshape
        exit `rc'
    }
    if `fallback' {
        * strL variables, which a plugin cannot write: reshape does it
        if `timing' di as text "creshape: used reshape (`__creshape_fallback')"
        reshape `dir'
        exit
    }

    return scalar N = _N
    return scalar k = c(k)
    return scalar N_before = `N0'
    return scalar k_before = `k0'
    foreach stat in load plan levels sort groups schema store {
        capture local __time_`stat' = _creshape_time_`stat'
        if _rc local __time_`stat' = .
        return scalar time_`stat' = `__time_`stat''
    }
    capture local __threads = _creshape_threads_max
    if _rc local __threads = .
    return scalar threads_max = `__threads'
    return local cmd "creshape"

    if `timing' {
        quietly _ctools_clock __ct_now
        local __t_store = `__ct_now' - `__ct_start'
        local __t_total = `__t_parse' + `__t_plan' + `__t_schema' + `__t_store'
        local __t_c = `__time_load' + `__time_plan' + `__time_schema' + `__time_store'
        di as text ""
        di as text "{hline 55}"
        di as text "creshape timing breakdown:"
        di as text "{hline 55}"
        di as text "  C plugin:"
        di as text "    Copy data from Stata:   " as result %8.4f `__time_load' " sec"
        di as text "    Order and check:        " as result %8.4f `__time_plan' " sec"
        if "`dir'" == "wide" {
            di as text "      Values of j:          " as result %8.4f `__time_levels' " sec"
            di as text "      Order by i and j:     " as result %8.4f `__time_sort' " sec"
            di as text "      Groups and cells:     " as result %8.4f `__time_groups' " sec"
        }
        di as text "    Describe the result:    " as result %8.4f `__time_schema' " sec"
        di as text "    Write data to Stata:    " as result %8.4f `__time_store' " sec"
        di as text "  {hline 53}"
        di as text "    C plugin total:         " as result %8.4f `__t_c' " sec"
        di as text "  Stata:"
        di as text "    Settings and checks:    " as result %8.4f `__t_parse' " sec"
        di as text "    Template, use, set obs: " as result %8.4f (`__t_total' - `__t_parse' - `__t_c') " sec"
        di as text "  {hline 53}"
        di as text "    Wall clock total:       " as result %8.4f `__t_total' " sec"
        di as text "    Threads:                " as result %8.0f `__threads'
        di as text "{hline 55}"
    }
    _ctools_cleanup creshape
end

* reshape's picture of the two forms
program define _creshape_picture
    di as error ""
    di as error _col(9) " {it:long}" _skip(32) "{it:wide}"
    di as error _col(9) "{c TLC}{hline 15}{c TRC}" _skip(19) "{c TLC}{hline 18}{c TRC}"
    di as error _col(9) "{c |} {it:i   j}   a   b {c |}" _skip(19) "{c |} {it:i}   a1 a2  b1 b2 {c |}"
    di as error _col(9) "{c |}{hline 15}{c |} <-- {bf:creshape} ---> {c |}{hline 18}{c |}"
    di as error _col(9) "{c |} 1   1   1   2 {c |}" _skip(19) "{c |} 1   1   3   2  4 {c |}"
    di as error _col(9) "{c |} 1   2   3   4 {c |}" _skip(19) "{c |} 2   5   7   6  8 {c |}"
    di as error _col(9) "{c |} 2   1   5   6 {c |}" _skip(19) "{c BLC}{hline 18}{c BRC}"
    di as error _col(9) "{c |} 2   2   7   8 {c |}"
    di as error _col(9) "{c BLC}{hline 15}{c BRC}"
    di as error ""
    di as error _col(9) "{it:long} to {it:wide}: {bf:creshape wide a b, i(}{it:i}{bf:) j(}{it:j}{bf:)}    ({it:j} existing variable)"
    di as error _col(9) "{it:wide} to {it:long}: {bf:creshape long a b, i(}{it:i}{bf:) j(}{it:j}{bf:)}    ({it:j} new variable)"
end
