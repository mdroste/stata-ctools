*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_clear_sortedby
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    * Stata drops the sort marker when a sort variable changes: change one
    * value with replace, then put the original back through Mata (no copy
    * of the variable is needed)
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_sb_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_sb_hold)
    }
    else {
        mata: __ctools_sb_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_sb_hold)
    }
    mata: mata drop __ctools_sb_hold
end
