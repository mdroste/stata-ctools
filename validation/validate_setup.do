/*******************************************************************************
 * validate_setup.do
 *
 * Core validation framework for ctools test suite
 * Provides helper programs and consistent PASS/FAIL output formatting
 *
 * Output format:
 *   [PASS] test description
 *   [FAIL] test description: detailed error message
 *
 * Failed tests are collected and displayed in a summary report at the end.
 ******************************************************************************/

version 14.0
clear all
set more off
set trace off

* Add build directory with ctools to adopath (use ++ for highest priority)
* Works from both project root and validation directory
capture adopath ++ "build"
capture adopath ++ "../build"

* Global counters
global TESTS_PASSED = 0
global TESTS_FAILED = 0
global TESTS_TOTAL = 0

* Default tolerance
global DEFAULT_TOL = 1e-7

* Initialize failure tracking (up to 100 failures)
global FAILURE_COUNT = 0
forvalues i = 1/100 {
    global FAILURE_`i' = ""
}

/*******************************************************************************
 * test_pass - Record a passing test
 ******************************************************************************/
capture program drop test_pass
program define test_pass
    args testname
    global TESTS_PASSED = $TESTS_PASSED + 1
    global TESTS_TOTAL = $TESTS_TOTAL + 1
    di as result "[PASS] `testname'"
end

/*******************************************************************************
 * test_fail - Record a failing test with description
 ******************************************************************************/
capture program drop test_fail
program define test_fail
    args testname reason
    global TESTS_FAILED = $TESTS_FAILED + 1
    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Store failure details
    global FAILURE_COUNT = $FAILURE_COUNT + 1
    local idx = $FAILURE_COUNT
    if `idx' <= 100 {
        global FAILURE_`idx' = "`testname': `reason'"
    }

    * Print failure immediately (with emphasis)
    di as error ""
    di as error ">>> [FAIL] `testname'"
    di as error ">>>        `reason'"
    di as error ""
end

/*******************************************************************************
 * print_failures - Print all recorded failures
 ******************************************************************************/
capture program drop print_failures
program define print_failures
    if $FAILURE_COUNT == 0 {
        exit
    }

    di as error ""
    di as error "{hline 70}"
    di as error "FAILED TESTS DETAIL:"
    di as error "{hline 70}"

    local max_to_show = min($FAILURE_COUNT, 100)
    forvalues i = 1/`max_to_show' {
        di as error "  `i'. ${FAILURE_`i'}"
    }

    if $FAILURE_COUNT > 100 {
        local more = $FAILURE_COUNT - 100
        di as error ""
        di as error "  ... and `more' more failures not shown"
    }

    di as error "{hline 70}"
    di as error ""
end

/*******************************************************************************
 * assert_scalar_equal - Compare two scalar values
 ******************************************************************************/
capture program drop assert_scalar_equal
program define assert_scalar_equal
    args val1 val2 tol testname

    if "`tol'" == "" local tol = $DEFAULT_TOL

    local diff = abs(`val1' - `val2')
    if `diff' == 0 | `diff' < `tol' {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "expected `val2', got `val1' (diff=`diff')"
    }
end

/*******************************************************************************
 * assert_matrix_equal - Compare two matrices element-wise
 ******************************************************************************/
capture program drop assert_matrix_equal
program define assert_matrix_equal
    args mat1 mat2 tol testname

    if "`tol'" == "" local tol = $DEFAULT_TOL

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
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "max element diff=`maxdiff'"
    }
end

/*******************************************************************************
 * assert_var_equal - Compare two numeric variables
 ******************************************************************************/
capture program drop assert_var_equal
program define assert_var_equal
    args var1 var2 tol testname

    if "`tol'" == "" local tol = 1e-10

    quietly count if abs(`var1' - `var2') > `tol'
    local ndiff = r(N)

    if `ndiff' == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "`ndiff' values differ"
    }
end

/*******************************************************************************
 * assert_strvar_equal - Compare two string variables
 ******************************************************************************/
capture program drop assert_strvar_equal
program define assert_strvar_equal
    args var1 var2 testname

    quietly count if `var1' != `var2'
    local ndiff = r(N)

    if `ndiff' == 0 {
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "`ndiff' values differ"
    }
end

/*******************************************************************************
 * assert_data_equal - Compare two datasets (with optional sorting)
 ******************************************************************************/
capture program drop assert_data_equal
program define assert_data_equal
    args file1 file2 testname

    preserve

    * Try direct comparison first
    capture quietly {
        use "`file1'", clear
        cf _all using "`file2'"
    }
    local rc = _rc

    if `rc' == 0 {
        restore
        test_pass "`testname'"
        exit
    }

    * Try sorted comparison
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
        test_pass "`testname'"
    }
    else {
        test_fail "`testname'" "datasets differ"
    }
end

/*******************************************************************************
 * print_section - Print section header
 ******************************************************************************/
capture program drop print_section
program define print_section
    args title
    di as text ""
    di as text "{hline 70}"
    di as text "`title'"
    di as text "{hline 70}"
end

/*******************************************************************************
 * print_summary - Print final test summary with failure details
 ******************************************************************************/
capture program drop print_summary
program define print_summary
    args component

    * Print detailed failure report first
    print_failures

    di as text ""
    di as text "{hline 70}"
    di as text "SUMMARY: `component'"
    di as text "{hline 70}"
    di as text "Tests passed: " as result $TESTS_PASSED as text " / " as result $TESTS_TOTAL
    di as text "Tests failed: " as result $TESTS_FAILED
    di as text "{hline 70}"
    if $TESTS_FAILED > 0 {
        di as error "VALIDATION FAILED - See failed tests above"
    }
    else {
        di as result "ALL TESTS PASSED"
    }
    di as text "{hline 70}"
end

/*******************************************************************************
 * reset_counters - Reset test counters (for running multiple test suites)
 ******************************************************************************/
capture program drop reset_counters
program define reset_counters
    global TESTS_PASSED = 0
    global TESTS_FAILED = 0
    global TESTS_TOTAL = 0
    global FAILURE_COUNT = 0
    forvalues i = 1/100 {
        global FAILURE_`i' = ""
    }
end

* Load benchmark helper programs (works from project root or validation dir)
capture quietly do "validation/benchmark_helpers.do"
if _rc != 0 {
    quietly do "benchmark_helpers.do"
}
