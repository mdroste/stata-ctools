/*******************************************************************************
 * benchmark_helpers.do
 *
 * Benchmark helper programs for ctools validation test suite
 * Each program runs both ctools and Stata native commands, compares results,
 * and reports PASS/FAIL in a concise format.
 *
 * Usage:
 *   benchmark_sort varlist [, options]
 *   benchmark_merge mergetype keyvar using filename [, options]
 *   benchmark_reghdfe depvar indepvars, absorb(varlist) [options]
 *   benchmark_qreg depvar indepvars [, options]
 *   benchmark_import using filename [, options]
 *   benchmark_export using filename [, options]
 ******************************************************************************/

/*******************************************************************************
 * benchmark_sort
 *
 * Compare csort vs sort, stable
 *
 * Syntax: benchmark_sort varlist [, testname(string)]
 *
 * Example:
 *   sysuse auto, clear
 *   benchmark_sort price
 *   benchmark_sort foreign rep78, testname("multi-var sort")
 ******************************************************************************/
capture program drop benchmark_sort
program define benchmark_sort
    syntax varlist [, testname(string)]

    if "`testname'" == "" local testname "sort `varlist'"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

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
    csort `varlist'

    * Compare datasets
    capture quietly cf _all using `stata_sorted'
    local rc = _rc

    restore

    if `rc' == 0 {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "[PASS] `testname'"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname'"
    }
end

/*******************************************************************************
 * benchmark_merge
 *
 * Compare cmerge vs merge
 * Validates: (1) _merge counts match, (2) datasets are identical via cf _all
 *
 * Syntax: benchmark_merge mergetype keyvar using filename [, keep(string)
 *         generate(string) nogenerate keepusing(varlist) testname(string)]
 *
 * Example:
 *   benchmark_merge 1:1 id using mydata.dta
 *   benchmark_merge m:1 foreign using lookup.dta, keep(match)
 ******************************************************************************/
capture program drop benchmark_merge
program define benchmark_merge
    syntax anything using/, [keep(string) GENerate(name) NOGENerate KEEPUSing(string) testname(string) SORTed]

    * Parse merge type and key vars
    gettoken mergetype keyvars : anything

    if "`testname'" == "" local testname "merge `mergetype' `keyvars'"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Build option string
    local opts ""
    if "`keep'" != "" local opts "`opts' keep(`keep')"
    if "`generate'" != "" local opts "`opts' generate(`generate')"
    if "`nogenerate'" != "" local opts "`opts' nogenerate"
    if "`keepusing'" != "" local opts "`opts' keepusing(`keepusing')"
    if "`sorted'" != "" local opts "`opts' sorted"

    preserve

    * Save original data
    tempfile original
    quietly save `original'

    * Run Stata merge
    merge `mergetype' `keyvars' using `using', `opts'

    * Get _merge counts from Stata
    capture confirm variable _merge
    if _rc == 0 {
        quietly count if _merge == 1
        local stata_m1 = r(N)
        quietly count if _merge == 2
        local stata_m2 = r(N)
        quietly count if _merge == 3
        local stata_m3 = r(N)
    }
    else {
        local stata_m1 = .
        local stata_m2 = .
        local stata_m3 = .
    }

    * Save unsorted for first comparison
    tempfile stata_merged_unsorted
    quietly save `stata_merged_unsorted'

    * Sort by ALL variables for deterministic comparison
    * (keyvars alone may leave order undefined within key groups)
    qui ds
    local allvars `r(varlist)'
    sort `allvars'
    tempfile stata_merged_sorted
    quietly save `stata_merged_sorted'

    * Restore and run cmerge
    use `original', clear
    cmerge `mergetype' `keyvars' using `using', `opts'

    * Get _merge counts from cmerge
    capture confirm variable _merge
    if _rc == 0 {
        quietly count if _merge == 1
        local cmerge_m1 = r(N)
        quietly count if _merge == 2
        local cmerge_m2 = r(N)
        quietly count if _merge == 3
        local cmerge_m3 = r(N)
    }
    else {
        local cmerge_m1 = .
        local cmerge_m2 = .
        local cmerge_m3 = .
    }

    * Check 1: _merge counts match
    local counts_match = 1
    if `stata_m1' != `cmerge_m1' | `stata_m2' != `cmerge_m2' | `stata_m3' != `cmerge_m3' {
        local counts_match = 0
    }

    * Check 2a: cf _all BEFORE sorting (tests if order is identical too)
    capture quietly cf _all using `stata_merged_unsorted'
    local cf_unsorted_rc = _rc

    * Check 2b: cf _all AFTER sorting by ALL variables (tests content only)
    qui ds
    local allvars `r(varlist)'
    sort `allvars'
    capture quietly cf _all using `stata_merged_sorted'
    local cf_sorted_rc = _rc

    restore

    * Report results
    if `counts_match' & `cf_sorted_rc' == 0 {
        global TESTS_PASSED = $TESTS_PASSED + 1
        if `cf_unsorted_rc' == 0 {
            di as result "[PASS] `testname'"
        }
        else {
            di as result "[PASS] `testname' (sort order differs)"
        }
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        if !`counts_match' {
            di as error "[FAIL] `testname' (_merge counts: m1:`stata_m1'!=`cmerge_m1' m2:`stata_m2'!=`cmerge_m2' m3:`stata_m3'!=`cmerge_m3')"
        }
        else {
            di as error "[FAIL] `testname' (cf _all: datasets differ even after sorting)"
        }
    }
end

/*******************************************************************************
 * benchmark_reghdfe
 *
 * Compare creghdfe vs reghdfe
 *
 * Syntax: benchmark_reghdfe depvar indepvars [weight], absorb(varlist)
 *         [vce(string) testname(string) tol(real)]
 *
 * Example:
 *   benchmark_reghdfe price mpg weight, absorb(foreign)
 *   benchmark_reghdfe y x1 x2 [aw=w], absorb(firm year) vce(cluster firm)
 ******************************************************************************/
capture program drop benchmark_reghdfe
program define benchmark_reghdfe
    syntax varlist(min=2 fv) [aw fw pw] [if] [in], Absorb(varlist) [vce(string) testname(string) tol(real 1e-6)]

    gettoken depvar indepvars : varlist

    if "`testname'" == "" local testname "reghdfe `depvar' ... , absorb(`absorb')"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Build weight string
    local wtexp ""
    if "`weight'" != "" local wtexp "[`weight'`exp']"

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run reghdfe
    quietly reghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    matrix reghdfe_b = e(b)
    local reghdfe_N = e(N)

    * Run creghdfe
    quietly creghdfe `depvar' `indepvars' `wtexp' `if' `in', absorb(`absorb') `vceopt'
    matrix creghdfe_b = e(b)
    local creghdfe_N = e(N)

    restore

    * Compare N
    if `reghdfe_N' != `creghdfe_N' {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (N: `reghdfe_N' != `creghdfe_N')"
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
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "[PASS] `testname'"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (max coef diff: " %9.2e `maxdiff' ")"
    }
end

/*******************************************************************************
 * benchmark_qreg
 *
 * Compare cqreg vs qreg
 *
 * Syntax: benchmark_qreg depvar indepvars [, quantile(real) vce(string)
 *         testname(string) tol(real)]
 *
 * Example:
 *   benchmark_qreg price mpg weight
 *   benchmark_qreg price mpg weight, quantile(0.25)
 ******************************************************************************/
capture program drop benchmark_qreg
program define benchmark_qreg
    syntax varlist(min=2 fv) [if] [in], [Quantile(real 0.5) vce(string) testname(string) tol(real 1e-6)]

    gettoken depvar indepvars : varlist

    if "`testname'" == "" local testname "qreg `depvar' ... , q(`quantile')"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Build vce option
    local vceopt ""
    if "`vce'" != "" local vceopt "vce(`vce')"

    preserve

    * Run qreg
    quietly qreg `depvar' `indepvars' `if' `in', quantile(`quantile') `vceopt'
    matrix qreg_b = e(b)
    local qreg_N = e(N)

    * Run cqreg
    quietly cqreg `depvar' `indepvars' `if' `in', quantile(`quantile') `vceopt'
    matrix cqreg_b = e(b)
    local cqreg_N = e(N)

    restore

    * Compare N
    if `qreg_N' != `cqreg_N' {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (N: `qreg_N' != `cqreg_N')"
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
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "[PASS] `testname'"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (max coef diff: " %9.2e `maxdiff' ")"
    }
end

/*******************************************************************************
 * benchmark_import
 *
 * Compare cimport delimited vs import delimited
 *
 * Syntax: benchmark_import using filename [, delimiters(string) varnames(string)
 *         case(string) rowrange(string) testname(string)]
 *
 * Example:
 *   benchmark_import using "data.csv"
 *   benchmark_import using "data.tsv", delimiters(tab)
 ******************************************************************************/
capture program drop benchmark_import
program define benchmark_import
    syntax using/, [DELIMiters(string) VARNames(string) CASE(string) ROWRange(string) testname(string) clear]

    if "`testname'" == "" local testname "import `using'"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Build option string
    local opts "clear"
    if "`delimiters'" != "" local opts "`opts' delimiters(`delimiters')"
    if "`varnames'" != "" local opts "`opts' varnames(`varnames')"
    if "`case'" != "" local opts "`opts' case(`case')"
    if "`rowrange'" != "" local opts "`opts' rowrange(`rowrange')"

    preserve

    * Run Stata import
    import delimited `using', `opts'
    local stata_n = _N
    local stata_k = c(k)
    tempfile stata_import
    quietly save `stata_import'

    * Run cimport
    cimport delimited `using', `opts'
    local cimport_n = _N
    local cimport_k = c(k)

    restore

    * Compare dimensions
    if `stata_n' == `cimport_n' & `stata_k' == `cimport_k' {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "[PASS] `testname' (N=`stata_n', K=`stata_k')"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (N: `stata_n'!=`cimport_n', K: `stata_k'!=`cimport_k')"
    }
end

/*******************************************************************************
 * benchmark_export
 *
 * Compare cexport delimited vs export delimited
 *
 * Syntax: benchmark_export [varlist] using filename [, delimiter(string)
 *         novarnames quote testname(string)]
 *
 * Example:
 *   benchmark_export using "output.csv", replace
 *   benchmark_export price mpg using "output.csv", replace
 ******************************************************************************/
capture program drop benchmark_export
program define benchmark_export
    syntax [varlist] using/, [DELIMiter(string) NOVARNames QUOTE replace testname(string)]

    if "`testname'" == "" local testname "export `using'"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Build option string
    local opts "replace"
    if "`delimiter'" != "" local opts "`opts' delimiter(`delimiter')"
    if "`novarnames'" != "" local opts "`opts' novarnames"
    if "`quote'" != "" local opts "`opts' quote"

    preserve

    * Generate temp filenames
    tempfile stata_export cexport_export
    local stata_csv = subinstr("`stata_export'", ".dta", ".csv", .)
    local cexport_csv = subinstr("`cexport_export'", ".dta", ".csv", .)

    * Run Stata export
    if "`varlist'" != "" {
        export delimited `varlist' using "`stata_csv'", `opts'
    }
    else {
        export delimited using "`stata_csv'", `opts'
    }

    * Run cexport
    if "`varlist'" != "" {
        cexport delimited `varlist' using "`cexport_csv'", `opts'
    }
    else {
        cexport delimited using "`cexport_csv'", `opts'
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
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "[PASS] `testname' (N=`stata_n', K=`cexport_k')"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (N: `stata_n'!=`cexport_n', K: `stata_k'!=`cexport_k')"
    }
end

/*******************************************************************************
 * benchmark_binscatter (for cbinscatter validation)
 *
 * Compare cbinscatter vs binscatter (if installed)
 *
 * Syntax: benchmark_binscatter yvar xvar [, controls(varlist) absorb(varlist)
 *         by(varname) nquantiles(int) testname(string)]
 ******************************************************************************/
capture program drop benchmark_binscatter
program define benchmark_binscatter
    syntax varlist(min=2 max=2) [if] [in] [aw fw pw], [controls(varlist) Absorb(varlist) BY(varname) NQuantiles(integer 20) testname(string) nograph]

    gettoken yvar xvar : varlist

    if "`testname'" == "" local testname "binscatter `yvar' `xvar'"

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Build option string
    local opts "nquantiles(`nquantiles') nograph"
    if "`controls'" != "" local opts "`opts' controls(`controls')"
    if "`absorb'" != "" local opts "`opts' absorb(`absorb')"
    if "`by'" != "" local opts "`opts' by(`by')"

    * Build weight string
    local wtexp ""
    if "`weight'" != "" local wtexp "[`weight'`exp']"

    preserve

    * Run cbinscatter
    quietly cbinscatter `yvar' `xvar' `wtexp' `if' `in', `opts'
    local cbin_N = e(N)
    local cbin_nq = e(nquantiles)
    matrix cbin_data = e(bindata)

    restore

    * Verify cbinscatter produced valid results
    if `cbin_N' > 0 & `cbin_nq' == `nquantiles' {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "[PASS] `testname' (N=`cbin_N', bins=`cbin_nq')"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "[FAIL] `testname' (N=`cbin_N', bins=`cbin_nq')"
    }
end

di as text "Benchmark helper programs loaded"
