*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
* The absorbed degrees of freedom after creghdfe, named in e(footnote) (for
* the displays of bootstrap, jackknife and other prefixes) and shown by
* creghdfe itself: reghdfe's footnote table from e(dof_table), with ? where
* the count of redundant coefficients is a lower bound and * for effects
* nested within a cluster variable. Nothing is shown without absorb().
program define _creghdfe_footnote
    version 14.1
    syntax [, width(integer 13)]
    if `"`e(absvars)'"' == "_cons" | `"`e(absvars)'"' == "" exit
    tempname table mytab
    matrix `table' = e(dof_table)
    local vars : rownames `table'
    foreach var of local vars {
        local width = max(`width', strlen("`var'"))
    }
    local rows = rowsof(`table')
    di as text _n "Absorbed degrees of freedom:"
    .`mytab' = ._tab.new, col(5) lmargin(0)
    .`mytab'.width `width' | 12 12 14 1 |
    .`mytab'.pad . 1 1 1 0
    .`mytab'.numfmt . %9.0g %9.0g %9.0g .
    .`mytab'.numcolor . text text result .
    .`mytab'.sep, top
    .`mytab'.titles "Absorbed FE" "Categories" " - Redundant" "  = Num. Coefs" ""
    .`mytab'.sep, middle
    local explain_exact = 0
    local explain_nested = 0
    forvalues i = 1/`rows' {
        local var : word `i' of `vars'
        local var = subinstr("`var'", "1.", "", .)
        local note " "
        if `table'[`i', 4] == 1 {
            local note "?"
            local explain_exact = 1
        }
        if `table'[`i', 5] == 1 {
            local note "*"
            local explain_nested = 1
        }
        .`mytab'.row "`var'" `=`table'[`i', 1]' `=`table'[`i', 2]' `=`table'[`i', 3]' "`note'"
    }
    .`mytab'.sep, bottom
    if `explain_exact' di as text "? = number of redundant parameters may be higher"
    if `explain_nested' di as text "* = FE nested within cluster; treated as redundant for DoF computation"
end
