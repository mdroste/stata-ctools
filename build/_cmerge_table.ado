*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_table
    args isupdate m1 m2 m3 m4 m5 mergevar
    if "`mergevar'" != "" {
        forvalues c = 1/5 {
            local v`c' "(`mergevar'==`c')"
        }
    }
    di
    di as txt _col(5) "Result" _col(33) "Number of obs"
    di as txt _col(5) "{hline 41}"
    di as txt _col(5) "Not matched" _col(30) as res %16.0fc (`m1'+`m2')
    if (`m1' | `m2') {
        di as txt _col(9) "from master" _col(30) as res %16.0fc `m1' as txt "  `v1'"
        di as txt _col(9) "from using" _col(30) as res %16.0fc `m2' as txt "  `v2'"
        di
    }
    if "`isupdate'" == "" {
        di as txt _col(5) "Matched" _col(30) as res %16.0fc `m3' as txt "  `v3'"
    }
    else {
        if (`m1' == 0 & `m2' == 0) di
        di as txt _col(5) "Matched" _col(30) as res %16.0fc (`m3'+`m4'+`m5')
        di as txt _col(9) "not updated" _col(30) as res %16.0fc `m3' as txt "  `v3'"
        di as txt _col(9) "missing updated" _col(30) as res %16.0fc `m4' as txt "  `v4'"
        di as txt _col(9) "nonmissing conflict" _col(30) as res %16.0fc `m5' as txt "  `v5'"
    }
    di as txt _col(5) "{hline 41}"
end
