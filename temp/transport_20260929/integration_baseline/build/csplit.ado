*! C-accelerated string splitting (ctools)
program define csplit, rclass
    version 14.1
    syntax varname(string) [if] [in], [Generate(string) noTrim Parse(string asis) ///
        DESTRING force float Ignore(string asis) percent Limit(numlist int >0 max=1) ///
        Verbose THReads(integer 0)]
    if "`limit'" == "" local limit 0
    if `threads' < 0 exit 198
    if _N > 2147483647 exit 920
    if "`destring'" == "" & `"`force'`float'`ignore'`percent'"' != "" {
        di as error "destring options require destring"
        exit 198
    }
    if "`generate'" == "" local generate `varlist'
    confirm name `generate'
    if `: word count `generate'' != 1 exit 198
    marksample touse, strok
    quietly count if `touse'
    if r(N) == 0 error 2000
    local trimcode = "`trim'" == ""
    if `"`parse'"' == "" | `"`parse'"' == `"""' {
        if "`trim'" != "" exit 198
        local parse `"" ""'
        local trimcode = 2
    }
    local nparse : word count `parse'
    if !`nparse' exit 198
    tokenize `"`parse'"'
    forvalues i = 1/`nparse' {
        local csplit_delim`i' `"``i''"'
        if `"`csplit_delim`i''"' == "" exit 198
    }
    local cap = cond(`limit' == 0, 2046, min(`limit',2046))
    local threadopt
    if `threads' local threadopt threads(`threads')
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
        plugin call ctools_plugin `varlist' if `touse', ///
            "csplit `threadopt' scan `nparse' `cap' `trimcode' `timing'"
        local k : word count `csplit_widths'
        local newvars
        forvalues j = 1/`k' {
            local newvars `newvars' `generate'`j'
        }
        confirm new variable `newvars'
        local staged
        foreach width of local csplit_widths {
            tempvar part
            quietly gen str`width' `part' = ""
            local staged `staged' `part'
        }
        plugin call ctools_plugin `staged', "csplit `threadopt' write"
        if "`destring'" != "" {
            local ignopt
            if `"`ignore'"' != "" local ignopt `"ignore(`ignore')"'
            quietly destring `staged', replace `force' `float' `ignopt' `percent'
        }
        rename (`staged') (`newvars')
    }
    local rc = _rc
    capture plugin call ctools_plugin, "csplit clear"
    if `rc' exit `rc'
    return local varlist `newvars'
    return local nvars `k'
    return scalar k_new = `k'
    di as text "variables created: " as result "`newvars'"
end
