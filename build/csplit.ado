*! version 1.4.0 20261127 github.com/mdroste/stata-ctools
*! C-accelerated string splitting (ctools)
*! original: Stata's split (split.ado 2.1.0, November 2, 2021), written by Nicholas J. Cox
*! original: see help csplit for details
program define csplit, rclass
    version 14.1
    syntax varname(string) [if] [in], [Generate(string) noTrim Parse(string asis) ///
        DESTRING force float Ignore(string asis) percent Limit(numlist int >0 max=1) ///
        Verbose THReads(integer 0)]
    * The if is evaluated now, as split does, before any helper replaces
    * the r() or e() results it may name.
    if `"`if'"' != "" {
        tempvar iftouse
        mark `iftouse' `if' `in'
        local if "if `iftouse'"
    }
    if "`limit'" == "" local limit 0
    if `threads' < 0 {
        di as error "threads() must be nonnegative"
        exit 198
    }
    if _N > 2147483647 exit 920
    * destring's options and ignore() are checked as split checks them
    if "`destring'" == "" {
        foreach option in force float ignore percent {
            if `"``option''"' != "" {
                di as error "may not specify `option' without specifying destring"
                exit 198
            }
        }
    }
    else if `"`ignore'"' != "" {
        if ustrpos(`"`ignore'"', ",") & !ustrpos(`"`ignore'"', `"""') {
            di as error `"ignore string must be enclosed with quotes to ignore commas or use options"'
            exit 198
        }
        _parse_ignore `ignore'
        sreturn clear
    }
    if "`generate'" == "" local generate `varlist'
    confirm name `generate'
    if `: word count `generate'' != 1 {
        di as error "invalid stub `generate'"
        exit 198
    }
    local trimcode = "`trim'" == ""
    * parse("") is the default, as in split
    if `"`parse'"' == "" | `"`parse'"' == `""""' {
        if "`trim'" != "" {
            di as error "parsing on spaces incompatible with notrim"
            exit 198
        }
        local parse `"" ""'
        local trimcode = 2
    }
    local nparse : word count `parse'
    if !`nparse' {
        di as error "parse() is empty"
        exit 198
    }
    tokenize `"`parse'"'
    forvalues i = 1/`nparse' {
        local csplit_delim`i' `"``i''"'
        if `"`csplit_delim`i''"' == "" {
            di as error "parse() contains an empty delimiter"
            exit 198
        }
    }
    local cap = cond(`limit' == 0, 2046, min(`limit',2046))
    local threadopt
    if `threads' local threadopt threads(`threads')
    local selectopt
    if `"`if'"' == "" local selectopt skipif
    local fillall = (`"`if'`in'"' == "")
    if `fillall' local selectopt `selectopt' fillall
    local timing = "`verbose'" != ""
    * SPI cannot write long strings. Native split retains full strL semantics.
    if "`: type `varlist''" == "strL" {
        local limopt
        if `limit' local limopt limit(`limit')
        local ignopt
        if `"`ignore'"' != "" local ignopt `"ignore(`ignore')"'
        local parseopt
        if `trimcode' != 2 local parseopt `"parse(`parse')"'
        split `varlist' `if' `in', generate(`generate') `trim' `parseopt' ///
            `destring' `force' `float' `ignopt' `percent' `limopt'
        return add
        if `timing' di as text "csplit: used native split for strL input"
        exit
    }
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    _ctools_strw `varlist'
    capture noisily {
        plugin call ctools_plugin `varlist' `if' `in', ///
            "csplit `threadopt' scan `nparse' `cap' `trimcode' `timing' `selectopt'"
        local k : word count `csplit_widths'
        local newvars
        forvalues j = 1/`k' {
            local newvars `newvars' `generate'`j'
        }
        confirm new variable `newvars'
        local staged
        foreach width of local csplit_widths {
            tempvar part
            mata: (void) st_addvar("str`width'", "`part'", `fillall')
            local staged `staged' `part'
        }
        plugin call ctools_plugin `staged', "csplit `threadopt' write"
        rename (`staged') (`newvars')
    }
    local rc = _rc
    capture plugin call ctools_plugin, "csplit clear"
    if `rc' exit `rc'
    * split's report: the new variables, listed by ds, then destring's
    * messages for each. destring runs on the final names, so its
    * characteristic destring_cmd names them (split's list ends in a space).
    return local varlist "`newvars' "
    return local nvars `k'
    return scalar k_new = `k'
    local s = cond(`k' > 1, "s", "")
    if "`destring'" != "" di as result "variable`s' born as string: "
    else di as result "variable`s' created as string: "
    if `k' ds `newvars'
    if "`destring'" != "" {
        local ignopt
        if `"`ignore'"' != "" local ignopt `"ignore(`ignore')"'
        local destring_vars "`newvars' "
        destring `destring_vars', replace `force' `float' `ignopt' `percent'
    }
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
