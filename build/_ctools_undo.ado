*! version 1.2.3 20261030 github.com/mdroste/stata-ctools
*! Exact undo of a ctools wrapper's changes without a dataset snapshot.
* The wrappers that use it only add observations and variables at the end of
* the dataset, add value labels and _dta characteristics or notes, change
* variable labels they name, and lose the sort marker; the plugin restores any
* values it overwrote first. The record lives in globals named <prefix>*.
*
*   _ctools_undo record <prefix>
*       before the first change: _N, c(k), the sort marker, c(changed), the
*       value-label names and the _dta characteristics (notes included)
*   _ctools_undo varlabels <prefix> <varlist>
*       add the variable labels of varlist to the record
*   _ctools_undo rollback <prefix>
*       drop what was added and put the recorded state back (call under
*       nobreak, after the plugin has restored overwritten values)
*   _ctools_undo clear <prefix>
*       drop the record after success
program define _ctools_undo
    version 14.1
    gettoken sub 0 : 0
    gettoken p 0 : 0
    if "`sub'" == "record" _ctools_undo_record `p'
    else if "`sub'" == "varlabels" _ctools_undo_varlabels `p' `0'
    else if "`sub'" == "rollback" _ctools_undo_rollback `p'
    else if "`sub'" == "clear" macro drop `p'*
    else exit 198
end

program define _ctools_undo_record
    args p
    macro drop `p'*
    global `p'N = _N
    global `p'k = c(k)
    global `p'sortedby : sortedby
    global `p'changed = c(changed)
    quietly label dir
    global `p'labels `r(names)'
    local names : char _dta[]
    global `p'dn `names'
    local i = 0
    foreach c of local names {
        local ++i
        mata: st_global("`p'd`i'", st_global("_dta[`c']"))
    }
    global `p'nl = 0
    global `p'active 1
end

program define _ctools_undo_varlabels
    gettoken p 0 : 0
    if "${`p'active}" == "" exit
    local i = ${`p'nl}
    foreach v of local 0 {
        local ++i
        global `p'lv`i' `v'
        mata: st_global("`p'll`i'", st_varlabel("`v'"))
    }
    global `p'nl = `i'
end

program define _ctools_undo_rollback
    args p
    if "${`p'active}" == "" {
        macro drop `p'*
        exit
    }
    local N0 = ${`p'N}
    local k0 = ${`p'k}
    if _N > `N0' quietly drop in `=`N0'+1'/l
    if c(k) > `k0' {
        mata: st_local("added", invtokens(st_varname((`k0'+1)..st_nvar())))
        drop `added'
    }
    quietly label dir
    local now `r(names)'
    local before ${`p'labels}
    local added : list now - before
    if "`added'" != "" label drop `added'
    local now : char _dta[]
    local old ${`p'dn}
    foreach c of local now {
        if !`: list c in old' char _dta[`c']
    }
    local i = 0
    foreach c of local old {
        local ++i
        mata: st_global("_dta[`c']", st_global("`p'd`i'"))
    }
    forvalues i = 1/${`p'nl} {
        mata: st_varlabel(st_global("`p'lv`i'"), st_global("`p'll`i'"))
    }
    local sortedby ${`p'sortedby}
    if "`: sortedby'" != "`sortedby'" {
        * The observations are back in their order: a stable sort only sets
        * the marker again
        if "`sortedby'" != "" sort `sortedby', stable
        else _ctools_undo_unsort
    }
    mata: st_updata(${`p'changed})
    macro drop `p'*
end

* Drop the sort marker without moving data: change one value of the first sort
* variable with replace, which Stata notices, then put it back through Mata.
program define _ctools_undo_unsort
    local marked : sortedby
    if "`marked'" == "" | _N == 0 exit
    local v : word 1 of `marked'
    capture confirm string variable `v'
    if _rc {
        mata: __ctools_undo_hold = st_data(1, "`v'")
        quietly replace `v' = cond(`v' == 0, 1, 0) in 1
        mata: st_store(1, "`v'", __ctools_undo_hold)
    }
    else {
        mata: __ctools_undo_hold = st_sdata(1, "`v'")
        quietly replace `v' = cond(`v' == "", "x", "") in 1
        mata: st_sstore(1, "`v'", __ctools_undo_hold)
    }
    mata: mata drop __ctools_undo_hold
end
