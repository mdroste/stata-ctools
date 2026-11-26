*! version 1.3.2 20261124 github.com/mdroste/stata-ctools
*! C-accelerated distinct value counts
*! original: distinct 1.4.3 (September 3, 2023) by Gary Longton and Nicholas J. Cox
*! original: see help cdistinct for details
program define cdistinct, rclass byable(recall)
    version 14.1
    * Check observation limit (Stata plugin API limitation)
    if _N > 2147483647 {
        di as error "ctools does not support datasets exceeding 2^31 (2.147 billion) observations"
        di as error "This is a limitation of Stata's plugin API"
        exit 920
    }
    * Options, defaults and messages below are those of distinct 1.4.3.
    * width() and a bare sort are in its syntax but not in its help file.
    syntax [varlist] [if] [in] [, MISSing Abbrev(integer 32) Joint ///
        MINimum(integer 1) MAXimum(integer -1) Width(integer -1) ///
        SORT SORT2(string asis) Verbose THReads(integer 0)]
    if !inrange(`threads', 0, 256) exit 198
    if `maximum' == -1 local maximum .
    if `minimum' > `maximum' {
        local swap `minimum'
        local minimum `maximum'
        local maximum `swap'
        di as txt "min(`maximum') max(`minimum') interpreted as min(`minimum') max(`maximum')"
    }
    * Column widths: the names as abbreviated, and counts of at least 11
    if "`joint'" == "" {
        local vnl = 1
        foreach v of local varlist {
            local vnl = max(`vnl', length(abbrev("`v'", `abbrev')))
        }
        local w1 = `vnl' + 2
    }
    else local w1 = 11
    if `width' > 0 local w2 = max(11, `width')
    else {
        local w2 = 11
        if _N > 1e7 local w2 = 5 + floor(log10(_N))
    }
    * sort() is read only without joint, as in distinct
    local tosort 0
    foreach opt in total distinct alpha descending {
        local `opt' 0
    }
    if "`joint'" == "" {
        if "`sort'" != "" local sort distinct
        if "`sort2'" != "" local sort `sort' `sort2'
        foreach word of local sort {
            local word = lower("`word'")
            local length = length("`word'")
            if "`word'" == "d" {
                di as err "sort() ambiguous: distinct or descending?"
                exit 198
            }
            if "`word'" == substr("total", 1, `length') {
                local total 1
            }
            else if "`word'" == substr("distinct", 1, `length') {
                local distinct 1
            }
            else if "`word'" == substr("alpha", 1, `length') {
                local alpha 1
            }
            else if "`word'" == substr("descending", 1, `length') {
                local descending 1
            }
            else {
                di as err "sort() option invalid?"
                exit 198
            }
        }
        local tosort = `alpha' + `total' + `distinct'
        if `tosort' > 1 | (`tosort' == 0 & `descending') {
            di as err "sort request invalid"
            exit 198
        }
    }
    marksample touse, novarlist
    _ctools_load
    capture program ctools_plugin, plugin using("`__ctools_plugin'")
    if _rc != 0 & _rc != 110 exit 601
    local threadopt
    if `threads' local threadopt threads(`threads')
    * strL values the plugin cannot read (binary, or over 2,045 bytes) are
    * counted through stand-in ranks, "" getting a missing one
    local stubs
    foreach v of local varlist {
        tempvar stub
        local stubs `stubs' `stub'
    }
    _ctools_strlkey `varlist', generate(`stubs') missing
    local keylist `__ctools_keys'
    local sets `varlist'
    if "`joint'" != "" local sets joint
    local ctools_count_N 0
    local ctools_count_distinct 0
    * Rows that pass minimum() and maximum(); the counts of the last variable
    * are returned whether or not its row is shown.
    local rows 0
    local index 0
    foreach v of local sets {
        local keys : word `=`index'+1' of `keylist'
        if "`joint'" != "" local keys `keylist'
        local nk : word count `keys'
        _ctools_strw `keys' `touse'
        if _N plugin call ctools_plugin `keys' `touse', ///
            "cdistinct nkeys=`nk' `missing' `verbose' `threadopt'"
        local ++index
        local N = `ctools_count_N'
        local ndistinct = `ctools_count_distinct'
        if "`joint'" != "" & `N' == 0 error 2000
        if `ndistinct' >= `minimum' & `ndistinct' <= `maximum' {
            local ++rows
            if "`joint'" != "" local rowname`rows' "(jointly)"
            else local rowname`rows' = abbrev("`v'", `abbrev')
            local rowN`rows' `N'
            local rowD`rows' `ndistinct'
            local names `names' `rowname`rows''
            local totals `totals' `N'
            local distincts `distincts' `ndistinct'
        }
    }
    if `rows' == 0 {
        noi display as txt "note: no distinct values satisfy specification"
        return scalar N = 0
        return scalar ndistinct = 0
        exit 0
    }
    * Row order. distinct leaves rows that tie on the sort key in an order
    * that changes from run to run (Mata's order() breaks ties at random);
    * here they keep the order of varlist.
    local order
    if `tosort' {
        local updn = cond(`descending', -1, 1)
        if `alpha' {
            mata: st_local("order", invtokens(strofreal(order(           ///
                (tokens(st_local("names"))', strofreal(1::`rows', "%010.0f")), ///
                (`updn', 2))')))
        }
        else {
            local keys `totals'
            if `distinct' local keys `distincts'
            mata: st_local("order", invtokens(strofreal(order(           ///
                (strtoreal(tokens(st_local("keys")))', (1::`rows')),     ///
                (`updn', 2))', "%18.0g")))
        }
    }
    * The table, drawn as distinct draws it. Timing output is already above.
    di
    tempname mytab
    .`mytab' = ._tab.new, col(3) lmargin(0)
    .`mytab'.width `w1' | `w2' `w2'
    .`mytab'.sep, top
    .`mytab'.titles " " "total" "distinct"
    .`mytab'.sep
    forvalues i = 1/`rows' {
        local r `i'
        if `tosort' local r : word `i' of `order'
        .`mytab'.strcolor "green" "yellow" "yellow"
        .`mytab'.row "`rowname`r''" "`rowN`r''" "`rowD`r''"
    }
    .`mytab'.sep, bottom
    return scalar N = `N'
    return scalar ndistinct = `ndistinct'
end
