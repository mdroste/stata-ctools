*! version 1.1.0 20261016 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_assert
    args codes n1 n2 n3 n4 n5
    local failed = 0
    forvalues c = 1/5 {
        if !`: list c in codes' & `n`c'' > 0 local failed = 1
    }
    if !`failed' exit
    local n : word count `codes'
    local comma = cond(`n' > 2, ",", "")
    di as err "{p 0 4 2}"
    di as err "after {bf:merge}, not all observations"
    local i = 0
    foreach c of local codes {
        local ++i
        if `c' == 1 local text "from master"
        else if `c' == 2 local text "from using"
        else if `c' == 3 local text "matched"
        else if `c' == 4 local text "matched and updated"
        else local text "matched and conflicting"
        if `i' == 1 di as err "`text'`comma'"
        else if `i' < `n' di as err "`text'`comma'"
        else di as err "or `text'"
    }
    di as err "{p_end}"
    di as err "(merged result left in memory)"
    global CTOOLS_cmerge_assert_failed 1
    exit 9
end
