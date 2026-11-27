*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! C-accelerated collapse: make a dataset of summary statistics
*! original: Stata's collapse (collapse.ado 6.2.2, June 28, 2020)
*! original: see help ccollapse for details

program define ccollapse
    version 14.1

    * Under version 4 or earlier, collapse runs its old command
    if _caller() < 5 {
        local c = string(_caller())
        version `c': collapse `0'
        exit
    }

    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }

    quietly _ctools_clock __ct_start

    * clist, qualifiers, weights and options are collapse's. frame() is in
    * its syntax, undocumented and unused; it is accepted here the same way.
    syntax [anything(name=clist id=clist equalok)] [if] [in] [aw fw iw pw] ///
        [, BY(varlist) CW FAST CLABEL frame(name) Verbose THReads(integer 0)]
    if !inrange(`threads', 0, 256) {
        di as error "threads() must be between 0 and 256"
        exit 198
    }
    local wgt
    if "`weight'" != "" local wgt `"[`weight'`exp']"'
    local command `"`clist' `if' `in' `wgt'"'
    local clist : subinstr local clist "[]" "", all

    * ---------------------------------------------------------------------
    * The clist: [(stat)] varlist  and  [(stat)] target = varname
    * ---------------------------------------------------------------------
    local stat mean
    local n 0
    local targets
    local sources
    while `"`clist'"' != "" {
        local done `n'
        _ccollapse_stat stat clist : "`stat'" `"`clist'"'
        gettoken first after : clist, parse(" =.") bind
        gettoken second : after, parse(" =.") bind
        if "`second'" == "=" {
            _ccollapse_target lhs rhs clist : `"`clist'"'
            tsrevar `rhs'
            local ++n
            local src`n' `r(varlist)'
            local new`n' `lhs'
            local lab`n' `rhs'
            local stat`n' `stat'
            local fmt`n' : format `src`n''
        }
        else {
            _ccollapse_varlist vars clist : `"`clist'"'
            tsrevar `vars'
            local actual `r(varlist)'
            foreach name of local vars {
                gettoken var actual : actual
                local ++n
                local src`n' `var'
                local lab`n' `name'
                local stat`n' `stat'
                if "`name'" == "`var'" {
                    local new`n' `var'
                    local fmt`n' : format `var'
                }
                else {
                    * a time-series operator: L2.x becomes the target L2x
                    local norm : tsnorm `name', varname
                    local dot = strpos("`norm'", ".")
                    local base = substr("`norm'", `dot' + 1, .)
                    local new`n' = substr("`norm'", 1, `dot' - 1) + "`base'"
                    local fmt`n' : format `base'
                }
            }
        }
        forvalues i = `=`done' + 1'/`n' {
            local targets `targets' `new`i''
            local sources `sources' `src`i''
        }
    }
    if `n' == 0 error 198

    local shared : list by & targets
    if "`shared'" != "" {
        di as error "`: word 1 of `shared'' may not be both target and by()"
        exit 198
    }

    * strL values cannot cross the plugin interface in both directions: with
    * a strL by() or source variable collapse itself does the work
    local strl 0
    foreach v in `by' `sources' {
        if "`: type `v''" == "strL" local strl 1
    }
    if `strl' {
        local options
        if "`by'" != "" local options `options' by(`by')
        if "`frame'" != "" local options `options' frame(`frame')
        if "`verbose'" != "" di as text "ccollapse: used collapse for data with strL variables"
        collapse `command', `options' `cw' `fast' `clabel'
        exit
    }

    * ---------------------------------------------------------------------
    * Checks, in collapse's order: the sample comes before the targets
    * ---------------------------------------------------------------------
    local problem
    local dups : list dups targets
    if "`dups'" != "" {
        forvalues i = 1/`n' {
            forvalues j = `=`i' + 1'/`n' {
                if "`new`i''" == "`new`j''" & "`problem'" == "" {
                    local problem conflict
                    local pi `i'
                    local pj `j'
                }
            }
        }
    }
    local count_note 0
    forvalues i = 1/`n' {
        if "`problem'" != "" continue, break
        local s `stat`i''
        if "`s'" == "count" & "`weight'" == "aweight" local count_note 1
        if ("`weight'" == "pweight" & inlist("`s'", "sd", "semean", "sebinomial", "sepoisson")) | ///
           ("`weight'" == "iweight" & inlist("`s'", "semean", "sebinomial", "sepoisson")) |      ///
           ("`weight'" == "aweight" & inlist("`s'", "sebinomial", "sepoisson")) {
            local problem weight
            local pi `i'
        }
        else if substr("`: type `src`i'''", 1, 3) == "str" & ///
            !inlist("`s'", "first", "last", "firstnm", "lastnm") {
            local problem string
            local pi `i'
        }
    }
    if "`problem'" != "" {
        local needed : list uniq sources
        if "`cw'" != "" local needed cw(`needed')
        else local needed
        _ccollapse_sample `if' `in' `wgt', `needed'
        if "`problem'" == "conflict" {
            di as error "error:" _n _col(8) `"`new`pi'' = (`stat`pi'') `lab`pi''"' ///
                _n _col(8) `"`new`pj'' = (`stat`pj'') `lab`pj''"' _n "name conflict"
            exit 198
        }
        if "`problem'" == "weight" {
            di as error "`stat`pi'' not allowed with `weight's"
            exit 135
        }
        di as error "string variable {bf:`src`pi''} not allowed with statistic {bf:`stat`pi''}"
        exit 109
    }
    if _N == 0 error 2000

    * ---------------------------------------------------------------------
    * Storage types of the targets, and the plugin's list of them
    * ---------------------------------------------------------------------
    local srclist : list uniq sources
    * The plugin reads the target list from the local ccollapse_targets: a
    * long clist would not fit in the plugin's argument string, which Stata
    * cuts off after a few thousand characters
    local types
    local ccollapse_targets
    local sep
    forvalues i = 1/`n' {
        local s `stat`i''
        local ty : type `src`i''
        if substr("`ty'", 1, 3) != "str" {
            if inlist("`s'", "sum", "rawsum", "percent") local ty double
            else if "`s'" == "count" {
                if inlist("`weight'", "fweight", "iweight", "pweight") local ty double
                else local ty long
            }
            else if !inlist("`s'", "max", "min", "first", "last", "firstnm", "lastnm") {
                if !inlist("`ty'", "double", "long") local ty `c(type)'
                else local ty double
            }
        }
        local types `types' `ty'
        if substr("`s'", 1, 2) == "p " local s p`=real(substr("`s'", 3, .))'
        local ccollapse_targets `ccollapse_targets'`sep'`: list posof "`src`i''" in srclist':`s'
        local sep ,
    }
    local tlen : length local ccollapse_targets

    * A weight that is one numeric variable is read as it is; an expression
    * is evaluated once. The sample marker replaces an if condition: Stata
    * answers an if one observation at a time, the plugin reads a marker in
    * parallel. An in range goes to the plugin as it is.
    local wvar
    local wtype none
    if "`weight'" != "" {
        local wtype `weight'
        local wexp = trim(substr(trim(`"`exp'"'), 2, .))
        local direct 0
        if `: word count `wexp'' == 1 {
            capture confirm numeric variable `wexp', exact
            if !_rc & !`: list wexp in by' & !`: list wexp in srclist' local direct 1
        }
        if `direct' local wvar `wexp'
        else {
            tempvar wvar
            quietly generate double `wvar' `exp' `if' `in'
        }
    }
    local touse
    local markeropt
    local range `in'
    if `"`if'"' != "" {
        tempvar touse
        mark `touse' `if' `in'
        local markeropt marker
        local range
    }
    local threadopt
    if `threads' > 0 local threadopt threads(`threads')
    local nby : word count `by'
    local nsrc : word count `srclist'

    _ctools_load
    * Stata scopes plugin registrations to the calling ado program.
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    capture confirm number 0

    quietly _ctools_lap __time_setup __ct_start `__ct_start'

    * compute: every statistic, in the plugin; the dataset is not changed
    _ctools_strw `by' `srclist' `wvar' `touse'
    capture noisily plugin call ctools_plugin `by' `srclist' `wvar' `touse' `range', ///
        "ccollapse `threadopt' compute nby=`nby' nsrc=`nsrc' wtype=`wtype' tlen=`tlen' `markeropt' `cw' `verbose'"
    local rc = _rc
    * the sample errors are collapse's, with Stata's messages
    if inlist(`rc', 401, 402, 2000) error `rc'
    if `rc' exit `rc'
    local ngroups = scalar(_ccollapse_ngroups)
    quietly _ctools_lap __time_compute __ct_start `__ct_start'

    * store: one observation per group. Nothing before this point changed
    * the data, and Break waits until the collapsed dataset is complete.
    * The by() variables are made anew with everything that describes them
    * (type, format, labels, characteristics and notes): dropping the data
    * and starting from ngroups observations takes a fraction of the time
    * that cutting a large dataset down to its first observations does.
    unab every : *
    local bylist : list every & by
    local nbyvars 0
    local bytypes
    local alias_chars fralias_from fralias_varname
    foreach v of local bylist {
        local ++nbyvars
        local bytypes `bytypes' `: type `v''
        local bfmt`nbyvars' : format `v'
        local bval`nbyvars' : value label `v'
        local bchars`nbyvars' : char `v'[]
        * collapse makes an alias variable an ordinary one: the result has
        * no characteristics that point to the linked frame
        local alias 0
        capture local alias : isalias `v'
        if "`alias'" == "1" {
            local bchars`nbyvars' : list bchars`nbyvars' - alias_chars
        }
        mata: st_local("blab`nbyvars'", st_varlabel("`v'"))
        local j 0
        foreach c of local bchars`nbyvars' {
            local ++j
            mata: st_local("bchar`nbyvars'_`j'", st_global("`v'[`c']"))
        }
    }
    local dlabel : data label
    local dchars : char _dta[]
    local j 0
    foreach c of local dchars {
        local ++j
        mata: st_local("dchar`j'", st_global("_dta[`c']"))
    }
    nobreak {
        quietly drop _all
        quietly set obs `ngroups'
        mata: (void) st_addvar(tokens(st_local("bytypes")), tokens(st_local("bylist")), 1)
        mata: (void) st_addvar(tokens(st_local("types")), tokens(st_local("targets")), 1)
        _ctools_strw `by' `targets'
        capture noisily plugin call ctools_plugin `by' `targets', "ccollapse `threadopt' store"
        local rc = _rc
        local i 0
        foreach v of local bylist {
            local ++i
            format `v' `bfmt`i''
            if "`bval`i''" != "" label values `v' `bval`i''
            mata: st_varlabel("`v'", st_local("blab`i'"))
            * Stata lists characteristics newest first: defining them from
            * the end of the list leaves the list as it was
            local names `bchars`i''
            forvalues j = `: word count `names''(-1)1 {
                local c : word `j' of `names'
                mata: st_global("`v'[`c']", st_local("bchar`i'_`j'"))
            }
        }
        * with by(), collapse keeps the data label and the characteristics
        * and notes of the dataset; without by() it starts from an empty one
        if "`by'" != "" {
            label data `"`dlabel'"'
            forvalues j = `: word count `dchars''(-1)1 {
                local c : word `j' of `dchars'
                mata: st_global("_dta[`c']", st_local("dchar`j'"))
            }
        }
        if `rc' == 0 {
            forvalues i = 1/`n' {
                if "`clabel'" == "" label variable `new`i'' `"(`stat`i'') `lab`i''"'
                else label variable `new`i'' `"`stat`i'' of `lab`i''"'
                format `new`i'' `fmt`i''
            }
            if "`by'" != "" sort `by'
        }
    }
    if `rc' exit `rc'
    global S_FN
    global S_FNDATE
    quietly _ctools_lap __time_store __ct_start `__ct_start'

    if "`verbose'" != "" {
        di as text ""
        di as text "{hline 60}"
        di as text "ccollapse timing breakdown:"
        di as text "{hline 60}"
        di as text "  Stata setup:              " as result %8.4f `__time_setup' " sec"
        di as text "  {hline 58}"
        di as text "  C plugin phases:"
        di as text "    Data load:              " as result %8.4f scalar(_ccollapse_time_load) " sec"
        di as text "    Sample (weights, cw):   " as result %8.4f scalar(_ccollapse_time_sample) " sec"
        di as text "    Groups (`ccollapse_method'):" _col(31) as result %8.4f scalar(_ccollapse_time_group) " sec"
        di as text "    Statistics:             " as result %8.4f scalar(_ccollapse_time_stats) " sec"
        di as text "    Data store:             " as result %8.4f scalar(_ccollapse_time_store) " sec"
        di as text "  {hline 58}"
        di as text "  Compute call:             " as result %8.4f `__time_compute' " sec"
        di as text "  Replace dataset:          " as result %8.4f `__time_store' " sec"
        di as text "{hline 60}"
        di as text "    Wall clock total:       " as result %8.4f `__time_setup' + `__time_compute' + `__time_store' " sec"
        di as text "{hline 60}"
        di as text ""
        di as text "  Observations: " as result %12.0fc scalar(_ccollapse_nobs)
        di as text "  Groups:       " as result %12.0fc scalar(_ccollapse_ngroups)
        di as text "  Threads:      " as result %12.0f scalar(_ccollapse_threads)
        di as text "{hline 60}"
    }
    if `count_note' di as text "(note: {bf:aweight}s not used to compute {bf:count}s)"
end

* A leading (stat) of the clist, if there is one; otherwise the statistic in
* force stays. Names are collapse's: semean, sebinomial and sepoisson may be
* abbreviated to three letters, median is p50, p# runs from p1 to p99.
program define _ccollapse_stat
    args into remainder colon current line
    gettoken word rest : line, parse("() ") match(paren)
    if "`paren'" == "" {
        c_local `into' "`current'"
        c_local `remainder' `"`line'"'
        exit
    }
    if `: word count `word'' == 1 {
        local word = trim(`"`word'"')
        local found
        if inlist("`word'", "mean", "sd", "sum", "rawsum", "count", "max", "min") | ///
           inlist("`word'", "iqr", "first", "firstnm", "last", "lastnm", "percent") local found `word'
        else if "`word'" == "median" local found "p 50"
        else if strlen("`word'") >= 3 {
            foreach name in semean sebinomial sepoisson {
                if "`word'" == substr("`name'", 1, strlen("`word'")) local found `name'
            }
        }
        if "`found'" == "" & substr("`word'", 1, 1) == "p" {
            * as collapse reports it, without the p
            local word = substr("`word'", 2, .)
            capture confirm integer number `word'
            if !_rc {
                if 1 <= `word' & `word' <= 99 local found "p `word'"
            }
        }
        if "`found'" != "" {
            * a statistic is followed by at least one variable
            if trim(`"`rest'"') == "" {
                local 0
                syntax varlist(ts)
            }
            c_local `into' "`found'"
            c_local `remainder' `"`rest'"'
            exit
        }
    }
    di as error "(`word') invalid statistic"
    exit 198
end

* target = varname at the start of the clist
program define _ccollapse_target
    args lhs_into rhs_into remainder colon line
    gettoken lhs line : line, parse(" =")
    gettoken equal line : line, parse(" =")
    gettoken 0 line : line, parse(" ")
    confirm name `lhs'
    syntax varname(ts)
    c_local `lhs_into' `lhs'
    c_local `rhs_into' `varlist'
    c_local `remainder' `"`line'"'
end

* The varlist at the start of the clist: it ends before the next (stat), and
* before a name that is followed by an equal sign
program define _ccollapse_varlist
    args into remainder colon line
    local taken
    while 1 {
        gettoken word after : line, parse(" =")
        local next = substr(trim(`"`after'"'), 1, 1)
        if "`next'" == "=" {
            c_local `remainder' `"`line'"'
            continue, break
        }
        local taken `taken' `word'
        local line `"`after'"'
        if "`next'" == "(" | "`next'" == "" {
            c_local `remainder' `"`after'"'
            continue, break
        }
    }
    local 0 `"`taken'"'
    syntax varlist(ts)
    c_local `into' `varlist'
end

* collapse's sample, marked the way collapse marks it: stops with Stata's
* error for a bad weight (r(401), r(402)) or an empty sample (r(2000))
program define _ccollapse_sample
    syntax [if] [in] [aw fw iw pw] [, CW(varlist)]
    tempvar touse
    if "`weight'" != "" {
        tempvar w
        quietly generate double `w' `exp' `if' `in'
        mark `touse' `if' `in' [`weight'=`w']
    }
    else mark `touse' `if' `in'
    if "`cw'" != "" markout `touse' `cw'
    quietly count if `touse'
    if r(N) == 0 error 2000
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
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
