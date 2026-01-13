/*******************************************************************************
 * benchmark_helpers.do
 *
 * Benchmark helper programs for ctools validation test suite
 * Each program runs both ctools and Stata native commands, compares results,
 * and reports [PASS] or [FAIL] with detailed error descriptions
 ******************************************************************************/

/*******************************************************************************
 * benchmark_sort - Compare csort vs sort, stable
 *
 * Syntax: benchmark_sort varlist [, testname(string) algorithm(string)]
 ******************************************************************************/
capture program drop benchmark_sort
program define benchmark_sort
    syntax varlist [, testname(string) ALGorithm(string)]

    if "`testname'" == "" local testname "sort `varlist'"
    if "`algorithm'" != "" local testname "`testname' [alg=`algorithm']"

    preserve

    * Save original data
    tempfile original
    quietly save `original'

    * Run Stata sort, stable
    sort `varlist', stable
    tempfile stata_sorted
    quietly save `stata_sorted'

    * Restore and run csort
    use `original', clear
    if "`algorithm'" != "" {
        capture csort `varlist', algorithm(`algorithm')
    }
    else {
        capture csort `varlist'
    }
    local csort_rc = _rc

    if `csort_rc' != 0 {
        restore
        test_fail "`testname'" "csort returned error `csort_rc'"
        exit
    }

    * Compare datasets
    capture quietly cf _all using `stata_sorted'
    local cf_rc = _rc

    restore

    if `cf_rc' == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "sorted data differs from Stata sort"
    }
end

/*******************************************************************************
 * benchmark_merge - Compare cmerge vs merge
 *
 * Syntax: benchmark_merge mergetype keyvar using filename [, options]
 ******************************************************************************/
capture program drop benchmark_merge
program define benchmark_merge
    syntax anything using/, [keep(string) GENerate(name) NOGENerate KEEPUSing(string) ///
        testname(string) SORTed ASSert(string) UPDATE REPLACE FORCE NOLabel NONotes NOREPort]

    * Parse merge type and key vars
    gettoken mergetype keyvars : anything

    if "`testname'" == "" local testname "merge `mergetype' `keyvars'"

    * Build option string
    local opts ""
    if "`keep'" != "" local opts "`opts' keep(`keep')"
    if "`generate'" != "" local opts "`opts' generate(`generate')"
    if "`nogenerate'" != "" local opts "`opts' nogenerate"
    if "`keepusing'" != "" local opts "`opts' keepusing(`keepusing')"
    if "`sorted'" != "" local opts "`opts' sorted"
    if "`assert'" != "" local opts "`opts' assert(`assert')"
    if "`update'" != "" local opts "`opts' update"
    if "`replace'" != "" local opts "`opts' replace"
    if "`force'" != "" local opts "`opts' force"
    if "`nolabel'" != "" local opts "`opts' nolabel"
    if "`nonotes'" != "" local opts "`opts' nonotes"

    preserve

    * Save original data
    tempfile original
    quietly save `original'

    * Run Stata merge
    capture merge `mergetype' `keyvars' using `using', `opts'
    local stata_rc = _rc

    if `stata_rc' != 0 {
        restore
        test_fail "`testname'" "Stata merge returned error `stata_rc'"
        exit
    }

    * Get _merge counts from Stata
    local stata_m1 = 0
    local stata_m2 = 0
    local stata_m3 = 0
    capture confirm variable _merge
    if _rc == 0 {
        quietly count if _merge == 1
        local stata_m1 = r(N)
        quietly count if _merge == 2
        local stata_m2 = r(N)
        quietly count if _merge == 3
        local stata_m3 = r(N)
    }

    * Sort by all variables for comparison
    qui ds
    local allvars `r(varlist)'
    sort `allvars'
    tempfile stata_merged
    quietly save `stata_merged'

    * Restore and run cmerge
    use `original', clear
    capture cmerge `mergetype' `keyvars' using `using', `opts' noreport
    local cmerge_rc = _rc

    if `cmerge_rc' != 0 {
        restore
        test_fail "`testname'" "cmerge returned error `cmerge_rc'"
        exit
    }

    * Get _merge counts from cmerge
    local cmerge_m1 = 0
    local cmerge_m2 = 0
    local cmerge_m3 = 0
    capture confirm variable _merge
    if _rc == 0 {
        quietly count if _merge == 1
        local cmerge_m1 = r(N)
        quietly count if _merge == 2
        local cmerge_m2 = r(N)
        quietly count if _merge == 3
        local cmerge_m3 = r(N)
    }

    * Check counts match
    local counts_match = (`stata_m1' == `cmerge_m1') & (`stata_m2' == `cmerge_m2') & (`stata_m3' == `cmerge_m3')

    * Sort and compare datasets
    qui ds
    local allvars `r(varlist)'
    sort `allvars'
    capture quietly cf _all using `stata_merged'
    local cf_rc = _rc

    restore

    if `counts_match' & `cf_rc' == 0 {
        test_pass "`testname'"
    }
    else if !`counts_match' {
        test_fail "`testname'" "_merge counts differ: stata(m1=`stata_m1',m2=`stata_m2',m3=`stata_m3') vs cmerge(m1=`cmerge_m1',m2=`cmerge_m2',m3=`cmerge_m3')"
    }
    else {
        test_fail "`testname'" "merged data differs from Stata merge"
    }
end

/*******************************************************************************
 * benchmark_reghdfe - Compare creghdfe vs reghdfe
 *
 * Syntax: benchmark_reghdfe depvar indepvars [weight], absorb(varlist) [options]
 ******************************************************************************/
capture program drop benchmark_reghdfe
program define benchmark_reghdfe
    syntax varlist(min=2 fv) [aw fw pw] [if] [in], Absorb(varlist) [vce(string) testname(string) tol(real 1e-6)]

    gettoken depvar indepvars : varlist

    if "`testname'" == "" local testname "reghdfe `depvar' `indepvars', absorb(`absorb')"

    * Build weight string
    local wtexp ""
    if "`weight'" != "" local wtexp "[`weight'`exp']"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run reghdfe
    capture quietly reghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    local reghdfe_rc = _rc

    if `reghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "reghdfe returned error `reghdfe_rc'"
        exit
    }

    matrix reghdfe_b = e(b)
    local reghdfe_N = e(N)

    * Run creghdfe
    capture quietly creghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    local creghdfe_rc = _rc

    if `creghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "creghdfe returned error `creghdfe_rc'"
        exit
    }

    matrix creghdfe_b = e(b)
    local creghdfe_N = e(N)

    restore

    * Compare N
    if `reghdfe_N' != `creghdfe_N' {
        test_fail "`testname'" "N differs: reghdfe=`reghdfe_N' vs creghdfe=`creghdfe_N'"
        exit
    }

    * Compare coefficients
    tempname diff
    matrix `diff' = reghdfe_b - creghdfe_b
    local cols = colsof(`diff')
    local maxdiff = 0
    forvalues j = 1/`cols' {
        local d = abs(`diff'[1, `j'])
        if `d' > `maxdiff' local maxdiff = `d'
    }

    if `maxdiff' < `tol' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "max coef diff = `maxdiff' (tol=`tol')"
    }
end

/*******************************************************************************
 * benchmark_qreg - Compare cqreg vs qreg
 *
 * Syntax: benchmark_qreg depvar indepvars [, quantile(real) options]
 ******************************************************************************/
capture program drop benchmark_qreg
program define benchmark_qreg
    syntax varlist(min=2 fv) [if] [in], [Quantile(real 0.5) vce(string) testname(string) tol(real 1e-6)]

    gettoken depvar indepvars : varlist

    if "`testname'" == "" local testname "qreg `depvar' `indepvars', q(`quantile')"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run qreg
    capture quietly qreg `depvar' `indepvars' `if' `in', quantile(`quantile') `vceopt'
    local qreg_rc = _rc

    if `qreg_rc' != 0 {
        restore
        test_fail "`testname'" "qreg returned error `qreg_rc'"
        exit
    }

    matrix qreg_b = e(b)
    local qreg_N = e(N)

    * Run cqreg
    capture quietly cqreg `depvar' `indepvars' `if' `in', quantile(`quantile') `vceopt'
    local cqreg_rc = _rc

    if `cqreg_rc' != 0 {
        restore
        test_fail "`testname'" "cqreg returned error `cqreg_rc'"
        exit
    }

    matrix cqreg_b = e(b)
    local cqreg_N = e(N)

    restore

    * Compare N
    if `qreg_N' != `cqreg_N' {
        test_fail "`testname'" "N differs: qreg=`qreg_N' vs cqreg=`cqreg_N'"
        exit
    }

    * Compare coefficients
    tempname diff
    matrix `diff' = qreg_b - cqreg_b
    local cols = colsof(`diff')
    local maxdiff = 0
    forvalues j = 1/`cols' {
        local d = abs(`diff'[1, `j'])
        if `d' > `maxdiff' local maxdiff = `d'
    }

    if `maxdiff' < `tol' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "max coef diff = `maxdiff' (tol=`tol')"
    }
end

/*******************************************************************************
 * benchmark_import - Compare cimport vs import delimited
 *
 * Syntax: benchmark_import using filename [, options]
 ******************************************************************************/
capture program drop benchmark_import
program define benchmark_import
    syntax using/, [DELIMiters(string) VARNames(string) CASE(string) testname(string)]

    if "`testname'" == "" local testname "import `using'"

    * Build option string
    local opts "clear"
    if "`delimiters'" != "" local opts "`opts' delimiters(`delimiters')"
    if "`varnames'" != "" local opts "`opts' varnames(`varnames')"
    if "`case'" != "" local opts "`opts' case(`case')"

    preserve

    * Run Stata import
    capture import delimited `using', `opts'
    local stata_rc = _rc

    if `stata_rc' != 0 {
        restore
        test_fail "`testname'" "import delimited returned error `stata_rc'"
        exit
    }

    local stata_n = _N
    local stata_k = c(k)
    tempfile stata_import
    quietly save `stata_import'

    * Run cimport
    capture cimport delimited `using', `opts'
    local cimport_rc = _rc

    if `cimport_rc' != 0 {
        restore
        test_fail "`testname'" "cimport returned error `cimport_rc'"
        exit
    }

    local cimport_n = _N
    local cimport_k = c(k)

    restore

    * Compare dimensions
    if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "dimensions differ: stata(N=`stata_n',K=`stata_k') vs cimport(N=`cimport_n',K=`cimport_k')"
    }
end

/*******************************************************************************
 * benchmark_export - Compare cexport vs export delimited
 *
 * Syntax: benchmark_export [varlist] using filename [, options]
 ******************************************************************************/
capture program drop benchmark_export
program define benchmark_export
    syntax [varlist] using/, [DELIMiter(string) NOVARNames QUOTE NOLabel replace testname(string)]

    if "`testname'" == "" local testname "export `using'"

    * Build option string
    local opts "replace"
    if "`delimiter'" != "" local opts "`opts' delimiter(`delimiter')"
    if "`novarnames'" != "" local opts "`opts' novarnames"
    if "`quote'" != "" local opts "`opts' quote"
    if "`nolabel'" != "" local opts "`opts' nolabel"

    preserve

    * Generate temp filenames
    tempfile stata_export cexport_export
    local stata_csv "`stata_export'.csv"
    local cexport_csv "`cexport_export'.csv"

    * Run Stata export
    if "`varlist'" != "" {
        capture export delimited `varlist' using "`stata_csv'", `opts'
    }
    else {
        capture export delimited using "`stata_csv'", `opts'
    }
    local stata_rc = _rc

    if `stata_rc' != 0 {
        restore
        test_fail "`testname'" "export delimited returned error `stata_rc'"
        exit
    }

    * Run cexport
    if "`varlist'" != "" {
        capture cexport delimited `varlist' using "`cexport_csv'", `opts'
    }
    else {
        capture cexport delimited using "`cexport_csv'", `opts'
    }
    local cexport_rc = _rc

    if `cexport_rc' != 0 {
        restore
        test_fail "`testname'" "cexport returned error `cexport_rc'"
        exit
    }

    * Reimport and compare
    import delimited using "`stata_csv'", clear
    local stata_n = _N
    local stata_k = c(k)

    import delimited using "`cexport_csv'", clear
    local cexport_n = _N
    local cexport_k = c(k)

    * Cleanup temp files
    capture erase "`stata_csv'"
    capture erase "`cexport_csv'"

    restore

    * Compare dimensions
    if `stata_n' == `cexport_n' & `stata_k' == `cexport_k' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "dimensions differ: stata(N=`stata_n',K=`stata_k') vs cexport(N=`cexport_n',K=`cexport_k')"
    }
end

/*******************************************************************************
 * benchmark_ivreghdfe - Compare civreghdfe vs ivreghdfe
 *
 * Syntax: benchmark_ivreghdfe spec, absorb(varlist) [options]
 ******************************************************************************/
capture program drop benchmark_ivreghdfe
program define benchmark_ivreghdfe
    syntax anything(name=spec), Absorb(varlist) [vce(string) testname(string) tol(real 1e-6)]

    if "`testname'" == "" local testname "ivreghdfe `spec', absorb(`absorb')"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run ivreghdfe
    capture quietly ivreghdfe `spec', absorb(`absorb') `vceopt'
    local ivreghdfe_rc = _rc

    if `ivreghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "ivreghdfe returned error `ivreghdfe_rc'"
        exit
    }

    matrix ivreghdfe_b = e(b)
    matrix ivreghdfe_V = e(V)
    local ivreghdfe_N = e(N)

    * Run civreghdfe
    capture quietly civreghdfe `spec', absorb(`absorb') `vceopt'
    local civreghdfe_rc = _rc

    if `civreghdfe_rc' != 0 {
        restore
        test_fail "`testname'" "civreghdfe returned error `civreghdfe_rc'"
        exit
    }

    matrix civreghdfe_b = e(b)
    matrix civreghdfe_V = e(V)
    local civreghdfe_N = e(N)

    restore

    * Compare N
    if `ivreghdfe_N' != `civreghdfe_N' {
        test_fail "`testname'" "N differs: ivreghdfe=`ivreghdfe_N' vs civreghdfe=`civreghdfe_N'"
        exit
    }

    * Compare coefficients
    tempname diff_b
    matrix `diff_b' = ivreghdfe_b - civreghdfe_b
    local cols = colsof(`diff_b')
    local maxdiff_b = 0
    forvalues j = 1/`cols' {
        local d = abs(`diff_b'[1, `j'])
        if `d' > `maxdiff_b' local maxdiff_b = `d'
    }

    * Compare VCE using relative tolerance
    tempname diff_V
    matrix `diff_V' = ivreghdfe_V - civreghdfe_V
    local rows = rowsof(`diff_V')
    local cols = colsof(`diff_V')
    local maxreldiff_V = 0
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local diff = abs(`diff_V'[`i', `j'])
            local base = abs(ivreghdfe_V[`i', `j'])
            if `base' > 1e-10 {
                local reldiff = `diff' / `base'
            }
            else {
                local reldiff = `diff'
            }
            if `reldiff' > `maxreldiff_V' local maxreldiff_V = `reldiff'
        }
    }

    if `maxdiff_b' < `tol' & `maxreldiff_V' < `tol' {
        test_pass "`testname'"
    }
    else if `maxdiff_b' >= `tol' {
        test_fail "`testname'" "max coef diff = `maxdiff_b' (tol=`tol')"
    }
    else {
        test_fail "`testname'" "max VCE relative diff = `maxreldiff_V' (tol=`tol')"
    }
end

di as text "Benchmark helper programs loaded"
