/*******************************************************************************
 * rollback_helpers.do - failed-call rollback checks shared by the suites
 *
 * By default a ctools wrapper that fails or is interrupted leaves the data as
 * they were without writing a snapshot of the dataset; option snapshot uses
 * preserve instead and option fast undoes nothing. Each wrapper has a
 * validation hook, the global CTOOLS_<command>_fault, that fails the call at a
 * named stage after it has started changing the data (with r(1), as Break
 * would, or inside the plugin with the memory error it simulates).
 *
 * _ctools_rb_case runs one such case in the default, snapshot and fast modes.
 * In the default and snapshot modes it compares the data afterwards with a
 * copy saved before: values (datasignature and cf), variable order, storage
 * types, formats, variable and value labels, value-label definitions, notes
 * and other characteristics, the data label, the sort marker and c(changed).
 * It also checks that no undo record (globals CTOOLS_u_*), plugin handle or
 * frame is left behind. With fast it checks that the failure is reported and
 * that the data did change (so the fault came after the first change).
 * Characteristics are compared by name: Stata's own restore does not keep
 * their internal order.
 ******************************************************************************/

capture program drop _ctools_rb_state
program define _ctools_rb_state, rclass
    quietly datasignature
    local s `"N=`=_N' k=`=c(k)' sig=`r(datasignature)' sortedby=`: sortedby' changed=`c(changed)'"'
    local s `"`macval(s)' label=`: data label'"'
    local chars : char _dta[]
    local chars : list sort chars
    foreach c of local chars {
        local value : char _dta[`c']
        local s `"`macval(s)' | _dta[`c']=`macval(value)'"'
    }
    if c(k) > 0 {
        foreach v of varlist _all {
            local s `"`macval(s)' || `v' `: type `v'' `: format `v'' [`: value label `v''] [`: variable label `v'']"'
            local chars : char `v'[]
            local chars : list sort chars
            foreach c of local chars {
                local value : char `v'[`c']
                local s `"`macval(s)' | `v'[`c']=`macval(value)'"'
            }
        }
    }
    quietly label dir
    local names `r(names)'
    local s `"`macval(s)' || labels: `names'"'
    if "`names'" != "" {
        tempfile defs
        tempname fh
        quietly label save `names' using `"`defs'"', replace
        file open `fh' using `"`defs'"', read text
        file read `fh' line
        while r(eof) == 0 {
            local s `"`macval(s)' | `macval(line)'"'
            file read `fh' line
        }
        file close `fh'
    }
    return local state `"`macval(s)'"'
end

* Run cmd (whose options end the command) with the fault hook global(name) set
* to fault, in the default, snapshot and fast modes; it must fail with rc.
* Default and snapshot must leave the data in memory as they were.
capture program drop _ctools_rb_case
program define _ctools_rb_case
    syntax , cmd(string asis) global(name) fault(string) testname(string) ///
        [rc(integer 1) modes(string)]
    if "`modes'" == "" local modes "default snapshot fast"
    tempfile copy
    quietly save `copy'
    _ctools_rb_state
    local before `"`r(state)'"'
    local frames ""
    if c(stata_version) >= 16 {
        quietly frames dir
        local frames `"`r(frames)'"'
    }
    foreach mode of local modes {
        quietly use `copy', clear
        local option = cond("`mode'" == "default", "", "`mode'")
        global `global' `fault'
        capture `cmd' `option'
        local got = _rc
        global `global'
        local name "`testname' [`mode']"
        if `got' != `rc' {
            test_fail "`name'" "expected r(`rc'), got r(`got')"
            continue
        }
        _ctools_rb_state
        local after `"`r(state)'"'
        mata: st_local("same", strofreal(st_local("before") == st_local("after")))
        capture cf _all using `copy', all
        local cf_rc = _rc
        * fast undoes nothing, so the data must differ: otherwise the fault
        * came before the first change and the other modes prove nothing
        if "`mode'" == "fast" {
            if `same' & !`cf_rc' test_fail "`name'" "the fault came before any change"
            else test_pass "`name'"
            continue
        }
        local left : all globals "CTOOLS_u_*"
        foreach g in CTOOLS_cwinsor_plugin CTOOLS_crangejoin_plugin {
            if "${`g'}" != "" local left `left' `g'
        }
        local frames_now ""
        if c(stata_version) >= 16 {
            quietly frames dir
            local frames_now `"`r(frames)'"'
        }
        if !`same' | `cf_rc' {
            test_fail "`name'" "data changed (state same=`same', cf r(`cf_rc'))"
        }
        else if "`left'" != "" | `"`frames_now'"' != `"`frames'"' {
            test_fail "`name'" "undo record, plugin handle or frame left behind (`left')"
        }
        else {
            test_pass "`name'"
        }
    }
    quietly use `copy', clear
end
