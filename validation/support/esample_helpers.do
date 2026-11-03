/*******************************************************************************
 * esample_helpers.do - estimation-sample parity checks shared by the suites
 *
 * Stata keeps e(sample) in a hidden column that native commands which move or
 * replace observations (sort, merge, rangejoin) carry along or reset. A ctools
 * wrapper that rewrites the visible columns in place has to give the same
 * result. _ctools_esample_case runs a setup program (which loads the data and
 * leaves estimation results m1 active, with m2 also stored), then the native
 * command and the ctools command, and compares row by row:
 *   - e(sample) right after the command
 *   - predict ... if e(sample)
 *   - e(sample) after estimates restore m2 and estimates restore m1
 *   - e(functions), e(_estimates_name), e(cmd) and c(changed)
 * Option sortby() compares the results in that order, for commands whose
 * output order is not part of the parity contract.
 ******************************************************************************/

capture program drop _ctools_esample_case
program define _ctools_esample_case
    syntax, setup(name) native(string asis) ctools(string asis) testname(string) ///
        [sortby(string)]
    tempfile out_native out_ctools
    foreach side in native ctools {
        `setup'
        capture noisily ``side''
        local rc_`side' = _rc
        local f : e(functions)
        local st_`side' `"rc=`rc_`side'' functions=[`f'] name=[`e(_estimates_name)'] cmd=[`e(cmd)'] changed=`c(changed)'"'
        quietly gen byte __es_after = e(sample)
        quietly predict double __es_p if e(sample)
        quietly estimates restore m2
        quietly estimates restore m1
        quietly gen byte __es_restored = e(sample)
        if "`sortby'" != "" sort `sortby', stable
        quietly save `out_`side''
    }
    if `"`st_native'"' != `"`st_ctools'"' {
        test_fail "`testname'" `"native `st_native' vs ctools `st_ctools'"'
        exit
    }
    quietly use `out_native', clear
    capture cf _all using `out_ctools', all
    if _rc test_fail "`testname'" "e(sample), predictions or restored samples differ from native"
    else test_pass "`testname'"
end
