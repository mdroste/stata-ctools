*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! C-accelerated distinct value counts
*! original: distinct 1.4.3 (September 3, 2023) by Gary Longton and Nicholas J. Cox
*! original: see help cdistinct for details
program define cdistinct, rclass byable(recall)
    version 14.1
    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    * Options, defaults and messages below are those of distinct 1.4.3.
    * width() and a bare sort are in its syntax but not in its help file.
    syntax [varlist] [if] [in] [, MISSing Abbrev(integer 32) Joint ///
        MINimum(integer 1) MAXimum(integer -1) Width(integer -1) ///
        SORT SORT2(string asis) Verbose THReads(integer 0)]
    if !inrange(`threads', 0, 256) exit 198
    if `maximum' == -1 local maximum .
    if `minimum' > `maximum' {
        local swap `minimum'
        local minimum `maximum'
        local maximum `swap'
        di as txt "min(`maximum') max(`minimum') interpreted as min(`minimum') max(`maximum')"
    }
    * Column widths: the names as abbreviated, and counts of at least 11
    if "`joint'" == "" {
        local vnl = 1
        foreach v of local varlist {
            local vnl = max(`vnl', length(abbrev("`v'", `abbrev')))
        }
        local w1 = `vnl' + 2
    }
    else local w1 = 11
    if `width' > 0 local w2 = max(11, `width')
    else {
        local w2 = 11
        if _N > 1e7 local w2 = 5 + floor(log10(_N))
    }
    * sort() is read only without joint, as in distinct
    local tosort 0
    foreach opt in total distinct alpha descending {
        local `opt' 0
    }
    if "`joint'" == "" {
        if "`sort'" != "" local sort distinct
        if "`sort2'" != "" local sort `sort' `sort2'
        foreach word of local sort {
            local word = lower("`word'")
            local length = length("`word'")
            if "`word'" == "d" {
                di as err "sort() ambiguous: distinct or descending?"
                exit 198
            }
            if "`word'" == substr("total", 1, `length') {
                local total 1
            }
            else if "`word'" == substr("distinct", 1, `length') {
                local distinct 1
            }
            else if "`word'" == substr("alpha", 1, `length') {
                local alpha 1
            }
            else if "`word'" == substr("descending", 1, `length') {
                local descending 1
            }
            else {
                di as err "sort() option invalid?"
                exit 198
            }
        }
        local tosort = `alpha' + `total' + `distinct'
        if `tosort' > 1 | (`tosort' == 0 & `descending') {
            di as err "sort request invalid"
            exit 198
        }
    }
    marksample touse, novarlist
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')
    * strL values the plugin cannot read (binary, or over 2,045 bytes) are
    * counted through stand-in ranks, "" getting a missing one
    local stubs
    foreach v of local varlist {
        tempvar stub
        local stubs `stubs' `stub'
    }
    _ctools_strlkey `varlist', generate(`stubs') missing
    local keylist `__ctools_keys'
    local sets `varlist'
    if "`joint'" != "" local sets joint
    local ctools_count_N 0
    local ctools_count_distinct 0
    * Rows that pass minimum() and maximum(); the counts of the last variable
    * are returned whether or not its row is shown.
    local rows 0
    local index 0
    foreach v of local sets {
        local keys : word `=`index'+1' of `keylist'
        if "`joint'" != "" local keys `keylist'
        local nk : word count `keys'
        _ctools_strw `keys' `touse'
        if _N plugin call ctools_plugin `keys' `touse', ///
            "cdistinct nkeys=`nk' `missing' `verbose' `threadopt'"
        local ++index
        local N = `ctools_count_N'
        local ndistinct = `ctools_count_distinct'
        if "`joint'" != "" & `N' == 0 error 2000
        if `ndistinct' >= `minimum' & `ndistinct' <= `maximum' {
            local ++rows
            if "`joint'" != "" local rowname`rows' "(jointly)"
            else local rowname`rows' = abbrev("`v'", `abbrev')
            local rowN`rows' `N'
            local rowD`rows' `ndistinct'
            local names `names' `rowname`rows''
            local totals `totals' `N'
            local distincts `distincts' `ndistinct'
        }
    }
    if `rows' == 0 {
        noi display as txt "note: no distinct values satisfy specification"
        return scalar N = 0
        return scalar ndistinct = 0
        exit 0
    }
    * Row order. distinct leaves rows that tie on the sort key in an order
    * that changes from run to run (Mata's order() breaks ties at random);
    * here they keep the order of varlist.
    local order
    if `tosort' {
        local updn = cond(`descending', -1, 1)
        if `alpha' {
            mata: st_local("order", invtokens(strofreal(order(           ///
                (tokens(st_local("names"))', strofreal(1::`rows', "%010.0f")), ///
                (`updn', 2))')))
        }
        else {
            local keys `totals'
            if `distinct' local keys `distincts'
            mata: st_local("order", invtokens(strofreal(order(           ///
                (strtoreal(tokens(st_local("keys")))', (1::`rows')),     ///
                (`updn', 2))', "%18.0g")))
        }
    }
    * The table, drawn as distinct draws it. Timing output is already above.
    di
    tempname mytab
    .`mytab' = ._tab.new, col(3) lmargin(0)
    .`mytab'.width `w1' | `w2' `w2'
    .`mytab'.sep, top
    .`mytab'.titles " " "total" "distinct"
    .`mytab'.sep
    forvalues i = 1/`rows' {
        local r `i'
        if `tosort' local r : word `i' of `order'
        .`mytab'.strcolor "green" "yellow" "yellow"
        .`mytab'.row "`rowname`r''" "`rowN`r''" "`rowD`r''"
    }
    .`mytab'.sep, bottom
    return scalar N = `N'
    return scalar ndistinct = `ndistinct'
end

* ctools-helpers:begin (generated by scripts/sync_release.py from src/ado; edit those files)
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

* _ctools_strlkey — numeric stand-ins for strL keys the plugin cannot read
* 
* Usage:  _ctools_strlkey varlist, generate(namelist) [missing]
* _ctools_strlkey varlist, check
* 
* The plugin reads strL text of at most 2,045 bytes and no binary values.
* Each strL variable of varlist that holds another value anywhere in the
* data gets a stand-in: the corresponding name of generate() becomes a long
* holding the rank of the value among the variable's distinct values, in
* Stata's sort order. Grouping by the ranks, ordering by them and counting
* their distinct values give what the strL values give. With missing, ""
* (a missing string) gets a missing rank, so that commands that leave out
* missing keys still leave it out; otherwise "" ranks first, as it sorts.
* 
* Sets local __ctools_keys in the caller: varlist with each such variable
* replaced by its stand-in. Variables the plugin can read stay as they are.
* The order of the observations, the sort marker and r() are left alone.
* check only sets local __ctools_strl_bad in the caller: the variables of
* varlist that would get a stand-in.

program define _ctools_strlkey
    version 14.1
    syntax varlist, [GENerate(namelist) MISSing CHECK]
    if "`check'" == "" & `: word count `generate'' != `: word count `varlist'' exit 198

    local strl
    foreach v of local varlist {
        if "`: type `v''" == "strL" local strl `strl' `v'
    }
    local strl : list uniq strl
    if "`strl'" == "" | _N == 0 {
        c_local __ctools_keys `varlist'
        c_local __ctools_strl_bad
        exit
    }

    * strlcheck reads lengths and binary flags only. If the plugin cannot
    * answer, every strL key gets a stand-in, which is exact anyway.
    local bad `strl'
    capture _ctools_load
    if !_rc {
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if !_rc | _rc == 110 {
            local ctools_strl_bad ""
            capture plugin call ctools_plugin `strl', "strlcheck"
            if !_rc & `: word count `ctools_strl_bad'' == `: word count `strl'' {
                local bad
                local i = 0
                foreach v of local strl {
                    local ++i
                    if `: word `i' of `ctools_strl_bad'' local bad `bad' `v'
                }
            }
        }
    }

    if "`check'" != "" {
        c_local __ctools_strl_bad `bad'
        exit
    }

    local keys
    local made
    local i = 0
    foreach v of local varlist {
        local ++i
        if !`: list v in bad' {
            local keys `keys' `v'
            continue
        }
        * A variable named twice shares one stand-in
        local at : list posof "`v'" in made
        if `at' {
            local keys `keys' `: word `at' of `madeid''
            continue
        }
        local id : word `i' of `generate'
        _ctools_strlkey_rank `v' `id' `missing'
        local made `made' `v'
        local madeid `madeid' `id'
        local keys `keys' `id'
    }
    c_local __ctools_keys `keys'
end

* Rank the values of v in Stata's sort order (sortpreserve puts the
* observations, and the sort marker, back; on Break as well).
program define _ctools_strlkey_rank, sortpreserve
    args v id missing
    sort `v'
    quietly by `v': generate long `id' = _n == 1
    quietly replace `id' = sum(`id')
    if "`missing'" != "" quietly replace `id' = . if `v' == ""
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
