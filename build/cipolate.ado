*! version 1.2.1 20261029 github.com/mdroste/stata-ctools
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
    tempvar result yy xx touse
    local sample `"`if' `in'"'
    local longby = 0
    foreach v of local by {
        if "`: type `v''" == "strL" local longby = 1
    }
    local run = _N > 0
    if `longby' {
        * Do not load a long by-key from a missing-x row: the SPI cannot read
        * arbitrary strL values, and native ipolate excludes that row first.
        marksample touse, novarlist
        markout `touse' `x'
        local sample "if `touse'"
        quietly count if `touse'
        local run = r(N) > 0
    }
    * The plugin already leaves missing x values missing. Pass the original
    * selection through once instead of materializing and rescanning touse.
    * With no restriction it writes every result, including missing values.
    if `"`if'`in'"' == "" & !`longby' mata: (void) st_addvar("double", st_local("result"), 1)
    else quietly gen double `result' = .
    local cipolate_nmissing = _N
    if `run' {
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
        local selectopt
        if `"`if'"' == "" & !`longby' local selectopt skipif
        _ctools_strw `sourcey' `sourcex' `by'
        plugin call ctools_plugin `sourcey' `sourcex' `by' `result' `sample', ///
            "cipolate `threadopt' `nby' `ext' `timing' `selectopt'"
    }
    label variable `result' "Interpolation of `y' on `x'"
    rename `result' `generate'
    if `cipolate_nmissing' di as text "(`cipolate_nmissing' missing values generated)"
end
