/*******************************************************************************
 * validate_setup.do
 *
 * Setup and helper programs for ctools validation test suite
 * This file is included by all individual validation tests
 ******************************************************************************/

version 14.0
clear all
set more off
set trace off

* Load ctools from build directory
adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

* Global counters
global TESTS_PASSED = 0
global TESTS_FAILED = 0
global TESTS_TOTAL = 0

* Helper program to compare scalars with tolerance
capture program drop assert_scalar_equal
program define assert_scalar_equal
    args name val1 val2 tol testname

    if "`tol'" == "" local tol = 1e-6

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    local diff = abs(`val1' - `val2')
    * Handle exact match case (for integers like N)
    if `diff' == 0 | `diff' < `tol' {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (`name' diff = " %12.2e `diff' ")"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "  FAIL: `testname' (`name': `val1' vs `val2', diff = " %12.2e `diff' ")"
    }
end

* Helper program to compare matrices
capture program drop assert_matrix_equal
program define assert_matrix_equal
    args mat1 mat2 tol testname

    if "`tol'" == "" local tol = 1e-6

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    tempname diff
    matrix `diff' = `mat1' - `mat2'
    local maxdiff = 0
    local rows = rowsof(`diff')
    local cols = colsof(`diff')
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local d = abs(`diff'[`i', `j'])
            if `d' > `maxdiff' local maxdiff = `d'
        }
    }

    if `maxdiff' < `tol' {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (max diff = " %12.2e `maxdiff' ")"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "  FAIL: `testname' (max diff = " %12.2e `maxdiff' ")"
    }
end

* Helper program to compare variables
capture program drop assert_var_equal
program define assert_var_equal
    args var1 var2 tol testname

    if "`tol'" == "" local tol = 1e-10

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    quietly count if abs(`var1' - `var2') > `tol'
    local ndiff = r(N)

    if `ndiff' == 0 {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (all values match)"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "  FAIL: `testname' (`ndiff' values differ)"
    }
end

* Helper program to compare string variables
capture program drop assert_strvar_equal
program define assert_strvar_equal
    args var1 var2 testname

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    quietly count if `var1' != `var2'
    local ndiff = r(N)

    if `ndiff' == 0 {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (all values match)"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "  FAIL: `testname' (`ndiff' values differ)"
    }
end

* Helper to check if two datasets are identical (with fallback to sorted comparison)
capture program drop assert_data_equal
program define assert_data_equal
    args file1 file2 testname

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    preserve
    capture quietly {
        use "`file1'", clear
        cf _all using "`file2'"
    }
    local rc = _rc
    restore

    if `rc' == 0 {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (datasets identical)"
    }
    else {
        * Try sorting both datasets by all variables before comparing
        preserve
        capture quietly {
            use "`file1'", clear
            ds
            local allvars `r(varlist)'
            sort `allvars'
            tempfile sorted1
            save `sorted1', replace

            use "`file2'", clear
            sort `allvars'
            cf _all using `sorted1'
        }
        local rc2 = _rc
        restore

        if `rc2' == 0 {
            global TESTS_PASSED = $TESTS_PASSED + 1
            di as result "  PASS: `testname' (datasets identical after sorting)"
        }
        else {
            global TESTS_FAILED = $TESTS_FAILED + 1
            di as error "  FAIL: `testname' (datasets differ even after sorting)"
        }
    }
end

* Helper to compare datasets after sorting (for merges where row order may differ)
capture program drop assert_data_equal_sorted
program define assert_data_equal_sorted
    args file1 file2 sortkey testname

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    preserve
    quietly {
        * Load and sort first file
        use "`file1'", clear
        sort `sortkey', stable
        tempfile sorted1
        save `sorted1'

        * Load and sort second file
        use "`file2'", clear
        sort `sortkey', stable
        tempfile sorted2
        save `sorted2'

        * Compare
        use `sorted1', clear
        cf _all using `sorted2'
    }
    local rc = _rc
    restore

    if `rc' == 0 {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (datasets identical after sorting)"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "  FAIL: `testname' (datasets differ even after sorting)"
    }
end

* Helper to compare merge counts (for merges where row order may differ)
capture program drop assert_merge_counts_equal
program define assert_merge_counts_equal
    args file1 file2 testname

    global TESTS_TOTAL = $TESTS_TOTAL + 1

    preserve
    quietly {
        use "`file1'", clear
        local n1 = _N
        count if _merge == 1
        local m1_1 = r(N)
        count if _merge == 2
        local m1_2 = r(N)
        count if _merge == 3
        local m1_3 = r(N)

        use "`file2'", clear
        local n2 = _N
        count if _merge == 1
        local m2_1 = r(N)
        count if _merge == 2
        local m2_2 = r(N)
        count if _merge == 3
        local m2_3 = r(N)
    }
    restore

    local pass = (`n1' == `n2') & (`m1_1' == `m2_1') & (`m1_2' == `m2_2') & (`m1_3' == `m2_3')

    if `pass' {
        global TESTS_PASSED = $TESTS_PASSED + 1
        di as result "  PASS: `testname' (N=`n1', m1=`m1_1', m2=`m1_2', m3=`m1_3')"
    }
    else {
        global TESTS_FAILED = $TESTS_FAILED + 1
        di as error "  FAIL: `testname' counts differ"
        di as error "    Stata: N=`n1', m1=`m1_1', m2=`m1_2', m3=`m1_3'"
        di as error "    cmerge: N=`n2', m1=`m2_1', m2=`m2_2', m3=`m2_3'"
    }
end

* Helper program to print section header
capture program drop print_section
program define print_section
    args title
    di as text ""
    di as text "{hline 70}"
    di as text "`title'"
    di as text "{hline 70}"
end

* Helper program to print test summary
capture program drop print_summary
program define print_summary
    args component
    di as text ""
    di as text "{hline 70}"
    di as text "SUMMARY: `component'"
    di as text "{hline 70}"
    di as text "Tests passed: " as result $TESTS_PASSED
    di as text "Tests failed: " as result $TESTS_FAILED
    di as text "Total tests:  " as result $TESTS_TOTAL
    if $TESTS_FAILED > 0 {
        di as error "VALIDATION FAILED"
    }
    else {
        di as result "ALL TESTS PASSED"
    }
    di as text "{hline 70}"
end

* Create temp directory for test files
capture mkdir "validation/temp"

di as text ""
di as text "ctools validation framework loaded"
di as text ""
