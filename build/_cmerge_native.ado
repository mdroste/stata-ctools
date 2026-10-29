*! version 1.2.1 20261029 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_native, rclass
    args cmdline
    local 0 `"`cmdline'"'
    gettoken mtype 0 : 0, parse(" ")
    syntax [anything] using/ [, Keep(string) ASSert(string) GENerate(name) ///
        NOGENerate KEEPUSing(string) SORTED FORCE NOREPort Verbose NOLabel ///
        NONotes UPDATE REPLACE PRESERVE_order(integer 0) THReads(integer 0) ///
        SNAPshot FAST]
    local opts `nogenerate' `sorted' `force' `noreport' `nolabel' `nonotes' `update' `replace'
    if "`keep'" != "" local opts `opts' keep(`keep')
    if "`assert'" != "" local opts `opts' assert(`assert')
    if "`generate'" != "" local opts `opts' generate(`generate')
    if "`keepusing'" != "" local opts `opts' keepusing(`keepusing')
    di as txt "(note: strL variables are involved; cmerge runs {bf:merge})"
    local caller = cond("$CTOOLS_cmerge_caller" != "", "$CTOOLS_cmerge_caller", "`c(stata_version)'")
    * merge can fail after sorting the master, so by default this path keeps
    * a snapshot: a frame in memory, which cmerge restores after a failure
    * (preserve before Stata 16). Options snapshot (cmerge preserved already)
    * and fast need none here.
    local preserved = 0
    if "`snapshot'`fast'" == "" {
        if c(stata_version) >= 16 frame copy `c(frame)' $CTOOLS_cmerge_snapframe
        else {
            preserve
            local preserved = 1
        }
    }
    capture noisily version `caller': merge `mtype' `anything' using `"`using'"', `opts'
    local rc = _rc
    if `rc' {
        * merge exits 9 only when assert() fails; its result stays in memory
        if `rc' == 9 {
            global CTOOLS_cmerge_assert_failed 1
            if `preserved' restore, not
        }
        exit `rc'
    }
    * Validation hook: fail as Break would once the result is in place
    if "$CTOOLS_cmerge_fault" == "post" exit 1
    local mvar = cond("`generate'" != "", "`generate'", "_merge")
    return scalar N = _N
    if "`nogenerate'" == "" {
        forvalues c = 1/5 {
            quietly count if `mvar' == `c'
            return scalar N_`c' = r(N)
        }
    }
    return local using `"`using'"'
    return local keyvars "`anything'"
    if `preserved' restore, not
end
