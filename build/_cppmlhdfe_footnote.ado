*! version 1.2.4 20261103 github.com/mdroste/stata-ctools
* Footnote of a cppmlhdfe fit, named in e(footnote): the absorbed degrees of
* freedom table e(dof_table).
program define _cppmlhdfe_footnote
    version 14.1
    if "`e(cmd)'" != "cppmlhdfe" error 301
    di as text ""
    di as text "Absorbed degrees of freedom:"
    matrix list e(dof_table), noheader format(%9.0g)
end
