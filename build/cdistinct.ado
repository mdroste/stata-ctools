*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! C-accelerated distinct value counts
program define cdistinct, rclass byable(recall)
    version 14.1
    syntax [varlist] [if] [in], [MISSing Abbrev(integer -1) Joint ///
        MINimum(integer 0) MAXimum(integer -1) Verbose THReads(integer 0)]
    if !inrange(`threads', 0, 256) exit 198
    if `maximum' == -1 local maximum .
    if `minimum' > `maximum' {
        local swap = `minimum'
        local minimum = `maximum'
        local maximum = `swap'
        di as text "min(`maximum') max(`minimum') interpreted as min(`minimum') max(`maximum')"
    }
    marksample touse, novarlist
    if `abbrev' == -1 {
        foreach v of local varlist {
            local abbrev = max(`abbrev', length("`v'"))
        }
    }
    local abbrev = max(5, `abbrev')
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')
    local sets `varlist'
    if "`joint'" != "" local sets joint
    local ctools_count_N 0
    local ctools_count_distinct 0
    local index 0
    foreach v of local sets {
        local keys `v'
        if "`joint'" != "" local keys `varlist'
        local nk : word count `keys'
        _ctools_strw `keys' `touse'
        if _N plugin call ctools_plugin `keys' `touse', ///
            "cdistinct nkeys=`nk' `missing' `verbose' `threadopt'"
        local ++index
        local count_n`index' `ctools_count_N'
        local count_d`index' `ctools_count_distinct'
    }
    * Keep timing output above the table, so its rows stay aligned.
    if "`joint'" != "" {
        di as text _newline "        Observations"
        di as text "      total   distinct"
    }
    else {
        local rule = `abbrev' + 2
        local divider = `abbrev' + 3
        di as text _newline _col(`divider') "{c |}        Observations"
        di as text _col(`divider') "{c |}      total   distinct"
        di as text "{hline `rule'}{c +}{hline 22}"
    }
    local index 0
    foreach v of local sets {
        local ++index
        local display_n `count_n`index''
        local display_d `count_d`index''
        if inrange(`display_d', `minimum', `maximum') {
            if "`joint'" != "" {
                di as result %11.0g `display_n' "  " %9.0g `display_d'
            }
            else {
                di as text " " %`abbrev's abbrev("`v'", `abbrev') " {c |}  " ///
                    as result %9.0g `display_n' "  " %9.0g `display_d'
            }
        }
    }
    return scalar N = `ctools_count_N'
    return scalar ndistinct = `ctools_count_distinct'
end
