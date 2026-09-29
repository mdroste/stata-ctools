*! C-accelerated linear interpolation (ctools)
program define cipolate, byable(onecall)
    version 14.1
    syntax varlist(numeric min=2 max=2) [if] [in], Generate(name) ///
        [BY(varlist) Epolate Verbose THReads(integer 0)]
    if _by() {
        if "`by'" != "" exit 190
        local by "`_byvars'"
    }
    confirm new variable `generate'
    if `threads' < 0 exit 198
    if _N > 2147483647 exit 920
    tokenize `varlist'
    local y `1'
    local x `2'
    * Alias copies allow x/y also to appear among the by variables.
    tempvar touse result yy xx
    marksample touse, novarlist
    quietly replace `touse' = 0 if missing(`x')
    quietly gen double `result' = .
    quietly count if `touse'
    if r(N) {
        local sourcey `y'
        local sourcex `x'
        if `: list y in by' {
            quietly clonevar `yy' = `y'
            local sourcey `yy'
        }
        if `: list x in by' {
            quietly clonevar `xx' = `x'
            local sourcex `xx'
        }
        _ctools_load
        capture program ctools_plugin, plugin using("`__ctools_plugin'")
        if _rc != 0 & _rc != 110 exit 601
        local nby : word count `by'
        local ext = "`epolate'" != ""
        local timing = "`verbose'" != ""
        local threadopt
        if `threads' local threadopt threads(`threads')
        _ctools_strw `sourcey' `sourcex' `by'
        plugin call ctools_plugin `sourcey' `sourcex' `by' `result' if `touse', ///
            "cipolate `threadopt' `nby' `ext' `timing'"
    }
    label variable `result' "Interpolation of `y' on `x'"
    rename `result' `generate'
    quietly count if missing(`generate')
    if r(N) di as text "(" r(N) " missing values generated)"
end
