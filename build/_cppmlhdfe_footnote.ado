*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
* Footnote of a cppmlhdfe fit, named in e(footnote): the absorbed degrees of
* freedom table e(dof_table), formatted as reghdfe's footnote (which
* ppmlhdfe shows); nothing without absorbed effects.
program define _cppmlhdfe_footnote
    version 14.1
    syntax [, width(int 13)]
    if "`e(cmd)'" != "cppmlhdfe" error 301
    if `"`e(absvars)'"' == "_cons" exit

    tempname table
    matrix `table' = e(dof_table)
    local rows = rowsof(`table')
    * the effects' names as typed (the matrix row names would carry c.
    * prefixes on interactions)
    local vars `"`e(extended_absvars)'"'
    if `: word count `vars'' != `rows' local vars : rownames `table'
    foreach var of local vars {
        local var = subinstr("`var'", "1.", "", .)
        if strlen("`var'") > `width' local width = strlen("`var'")
    }

    di as text _n "Absorbed degrees of freedom:"
    tempname mytab
    .`mytab' = ._tab.new, col(5) lmargin(0)
    .`mytab'.width    `width'  | 12   12   14   1 |
    .`mytab'.pad      .          1    1    1   0
    .`mytab'.numfmt   .      %9.0g %9.0g %9.0g  .
    .`mytab'.numcolor .       text text result  .
    .`mytab'.sep, top

    local explain_exact 0
    local explain_nested 0
    .`mytab'.titles "Absorbed FE" "Categories" " - Redundant" "  = Num. Coefs" ""
    .`mytab'.sep, middle
    forvalues i = 1/`rows' {
        local var : word `i' of `vars'
        local var = subinstr("`var'", "1.", "", .)
        local note " "
        if `table'[`i', 4] == 1 {
            local note "?"
            local explain_exact 1
        }
        if `table'[`i', 5] == 1 {
            local note "*"
            local explain_nested 1
        }
        .`mytab'.row "`var'" `=`table'[`i', 1]' `=`table'[`i', 2]' `=`table'[`i', 3]' "`note'"
    }
    .`mytab'.sep, bottom
    if `explain_exact' di as text "? = number of redundant parameters may be higher"
    if `explain_nested' di as text `"* = FE nested within cluster; treated as redundant for DoF computation"'
end
