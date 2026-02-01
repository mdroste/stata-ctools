/*******************************************************************************
 * validate_setup.do
 *
 * Core validation framework for ctools test suite
 * Minimal output: only summary table and failures are displayed
 *
 * TOLERANCE POLICY:
 *   All comparisons use SIGNIFICANT FIGURES, not absolute differences.
 *   sigfigs = -log10(|a - b| / max(|a|, |b|))
 *   Default threshold is 7 significant figures ($DEFAULT_SIGFIGS).
 *   This ensures scale-invariant precision: 1e6 vs 1e6+1 and 0.001 vs 0.001001
 *   both require ~3 sigfigs of agreement.
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

* Default significant figures threshold - do NOT change this value except by human decision.
* This ensures numerical precision of ctools implementations.
* A value of 7 means results must agree to at least 7 significant figures.
global DEFAULT_SIGFIGS = 7

* Legacy alias for backwards compatibility
global DEFAULT_TOL = 1e-7

* Initialize failure tracking (up to 100 failures)
global FAILURE_COUNT = 0
forvalues i = 1/100 {
    global FAILURE_`i' = ""
}

/*******************************************************************************
 * test_pass - Record a passing test (silent)
 ******************************************************************************/
capture program drop test_pass
program define test_pass
    args testname
    global TESTS_PASSED = $TESTS_PASSED + 1
    global TESTS_TOTAL = $TESTS_TOTAL + 1
    * Silent - no output
end

/*******************************************************************************
 * test_fail - Record a failing test with description (silent, stores for later)
 ******************************************************************************/
capture program drop test_fail
program define test_fail
    args testname reason
    global TESTS_FAILED = $TESTS_FAILED + 1
    global TESTS_TOTAL = $TESTS_TOTAL + 1

    * Store failure details (silently)
    global FAILURE_COUNT = $FAILURE_COUNT + 1
    local idx = $FAILURE_COUNT
    if `idx' <= 100 {
        global FAILURE_`idx' = "`testname': `reason'"
    }
    * No immediate output - failures shown in summary
end

/*******************************************************************************
 * print_failures - Print all recorded failures (called by validate_all.do)
 ******************************************************************************/
capture program drop print_failures
program define print_failures
    if $FAILURE_COUNT == 0 {
        exit
    }

    local max_to_show = min($FAILURE_COUNT, 100)
    forvalues i = 1/`max_to_show' {
        di as error "  ${FAILURE_`i'}"
    }

    if $FAILURE_COUNT > 100 {
        local more = $FAILURE_COUNT - 100
        di as error ""
        di as error "  ... and `more' more failures not shown"
    }
end

/*******************************************************************************
 * print_section - No-op for backwards compatibility (silent)
 ******************************************************************************/
capture program drop print_section
program define print_section
    args title
    * Silent - no output
end

/*******************************************************************************
 * print_summary - No-op for backwards compatibility (silent)
 * Summary is now printed only by validate_all.do
 ******************************************************************************/
capture program drop print_summary
program define print_summary
    args component
    di as text ""
    di as text "{hline 45}"
    di as text "  `component' validation results"
    di as text "{hline 45}"
    di as text "  Passed: " as result "$TESTS_PASSED"
    di as text "  Failed: " _continue
    if $TESTS_FAILED == 0 {
        di as result "$TESTS_FAILED"
    }
    else {
        di as error "$TESTS_FAILED"
    }
    di as text "  Total:  " as result "$TESTS_TOTAL"
    di as text "{hline 45}"
    if $TESTS_FAILED == 0 {
        di as result "  ALL TESTS PASSED"
    }
    else {
        di as error "  SOME TESTS FAILED"
        print_failures
    }
    di as text "{hline 45}"
    di as text ""
end

/*******************************************************************************
 * sigfigs - Compute number of significant figures of agreement between two values
 *
 * Formula: -log10(|a - b| / max(|a|, |b|))
 *
 * Returns:
 *   - 15 (double precision limit) if values are identical
 *   - 0 if one value is zero and the other is not
 *   - Number of matching significant figures otherwise
 *
 * Examples:
 *   1.234567 vs 1.234568 → ~6.1 sig figs (differ in 7th digit)
 *   1000000 vs 1000001   → ~6.0 sig figs
 *   0.001234 vs 0.001235 → ~3.1 sig figs
 ******************************************************************************/
capture program drop sigfigs
program define sigfigs, rclass
    args val1 val2

    * Handle missing values
    if missing(`val1') | missing(`val2') {
        return scalar sigfigs = 0
        exit
    }

    * If values are identical (including both zero), return max precision
    if `val1' == `val2' {
        return scalar sigfigs = 15
        exit
    }

    * If both values are effectively zero (below machine precision), they match
    * This handles cases like 0 vs 1e-17 which are numerically equivalent
    local eps = 1e-14
    if (abs(`val1') < `eps') & (abs(`val2') < `eps') {
        return scalar sigfigs = 15
        exit
    }

    * If one is zero and other is not, no significant figures agree
    if (`val1' == 0) | (`val2' == 0) {
        return scalar sigfigs = 0
        exit
    }

    * Compute significant figures of agreement
    local absdiff = abs(`val1' - `val2')
    local maxabs = max(abs(`val1'), abs(`val2'))
    local reldiff = `absdiff' / `maxabs'

    * Handle edge case where reldiff is 0 (shouldn't happen, but be safe)
    if `reldiff' == 0 {
        return scalar sigfigs = 15
        exit
    }

    * sigfigs = -log10(reldiff)
    local sf = -log10(`reldiff')

    * Handle missing result from log10 (extreme values)
    if missing(`sf') {
        return scalar sigfigs = 15
        exit
    }

    * Clamp to [0, 15] range (double precision limits)
    if `sf' < 0 {
        local sf = 0
    }
    if `sf' > 15 {
        local sf = 15
    }

    return scalar sigfigs = `sf'
end

/*******************************************************************************
 * assert_scalar_equal - Compare two scalar values by significant figures
 *
 * Tests whether two values agree to at least N significant figures.
 * Default threshold is 7 significant figures ($DEFAULT_SIGFIGS).
 ******************************************************************************/
capture program drop assert_scalar_equal
program define assert_scalar_equal
    args val1 val2 min_sigfigs testname

    if "`min_sigfigs'" == "" local min_sigfigs = $DEFAULT_SIGFIGS

    * Compute significant figures of agreement
    sigfigs `val1' `val2'
    local sf = r(sigfigs)

    if `sf' >= `min_sigfigs' {
        test_pass "`testname'"
    }
    else {
        local sf_fmt : display %5.1f `sf'
        test_fail "`testname'" "expected `val2', got `val1' (sigfigs=`sf_fmt', need `min_sigfigs')"
    }
end

/*******************************************************************************
 * assert_matrix_equal - Compare two matrices element-wise by significant figures
 *
 * Tests whether all matrix elements agree to at least N significant figures.
 * Reports the minimum sigfigs found across all element comparisons.
 ******************************************************************************/
capture program drop assert_matrix_equal
program define assert_matrix_equal
    args mat1 mat2 min_sigfigs testname

    if "`min_sigfigs'" == "" local min_sigfigs = $DEFAULT_SIGFIGS

    local min_sf = 15
    local min_i = 1
    local min_j = 1
    local rows = rowsof(`mat1')
    local cols = colsof(`mat1')
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local v1 = `mat1'[`i', `j']
            local v2 = `mat2'[`i', `j']
            sigfigs `v1' `v2'
            local sf = r(sigfigs)
            if `sf' < `min_sf' {
                local min_sf = `sf'
                local min_i = `i'
                local min_j = `j'
            }
        }
    }

    if `min_sf' >= `min_sigfigs' {
        test_pass "`testname'"
    }
    else {
        local sf_fmt : display %5.1f `min_sf'
        test_fail "`testname'" "element [`min_i',`min_j'] sigfigs=`sf_fmt', need `min_sigfigs'"
    }
end

/*******************************************************************************
 * matrix_min_sigfigs - Return minimum sigfigs across all matrix elements
 *
 * For use in inline comparisons. Returns r(min_sigfigs), r(min_i), r(min_j).
 ******************************************************************************/
capture program drop matrix_min_sigfigs
program define matrix_min_sigfigs, rclass
    args mat1 mat2

    local min_sf = 15
    local min_i = 1
    local min_j = 1
    local rows = rowsof(`mat1')
    local cols = colsof(`mat1')
    forvalues i = 1/`rows' {
        forvalues j = 1/`cols' {
            local v1 = `mat1'[`i', `j']
            local v2 = `mat2'[`i', `j']
            sigfigs `v1' `v2'
            local sf = r(sigfigs)
            if `sf' < `min_sf' {
                local min_sf = `sf'
                local min_i = `i'
                local min_j = `j'
            }
        }
    }

    return scalar min_sigfigs = `min_sf'
    return scalar min_i = `min_i'
    return scalar min_j = `min_j'
end

/*******************************************************************************
 * assert_var_equal - Compare two numeric variables by significant figures
 *
 * Tests whether all observations agree to at least N significant figures.
 ******************************************************************************/
capture program drop assert_var_equal
program define assert_var_equal
    args var1 var2 min_sigfigs testname

    if "`min_sigfigs'" == "" local min_sigfigs = $DEFAULT_SIGFIGS

    * Generate sigfigs for each observation
    tempvar sf
    quietly {
        * Handle identical values (including both zero)
        gen double `sf' = 15 if `var1' == `var2'
        * Handle one zero, other non-zero
        replace `sf' = 0 if (`sf' == .) & ((`var1' == 0) | (`var2' == 0))
        * Compute sigfigs for remaining cases
        replace `sf' = -log10(abs(`var1' - `var2') / max(abs(`var1'), abs(`var2'))) if `sf' == .
        * Clamp to valid range
        replace `sf' = 0 if `sf' < 0
        replace `sf' = 15 if `sf' > 15
    }

    * Count observations that fail the threshold
    quietly count if `sf' < `min_sigfigs'
    local ndiff = r(N)

    if `ndiff' == 0 {
        test_pass "`testname'"
    }
    else {
        quietly summarize `sf'
        local min_sf : display %5.1f r(min)
        test_fail "`testname'" "`ndiff' values below `min_sigfigs' sigfigs (min=`min_sf')"
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
