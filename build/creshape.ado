*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
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
