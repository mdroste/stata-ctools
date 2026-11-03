*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
*! C-accelerated string splitting (ctools)
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
    if "`destring'" == "" & `"`force'`float'`ignore'`percent'"' != "" {
        di as error "destring options require destring"
        exit 198
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
        if "`destring'" != "" {
            local ignopt
            if `"`ignore'"' != "" local ignopt `"ignore(`ignore')"'
            quietly destring `staged', replace `force' `float' `ignopt' `percent'
        }
        rename (`staged') (`newvars')
        * destring recorded the command with the temporary names; split
        * runs it on the final ones (its list ends in a space)
        if "`destring'" != "" {
            foreach v of local newvars {
                local cmd : char `v'[destring_cmd]
                local cmd : subinstr local cmd "`staged'" "`newvars' "
                char `v'[destring_cmd] `"`macval(cmd)'"'
            }
        }
    }
    local rc = _rc
    capture plugin call ctools_plugin, "csplit clear"
    if `rc' exit `rc'
    return local varlist `newvars'
    return local nvars `k'
    return scalar k_new = `k'
    di as text "variables created: " as result "`newvars'"
end
