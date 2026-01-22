*! validation_utils.do - Shared utilities for ctools validation tests
*! Include this file at the start of each validation test

* Setup adopath
adopath ++ "build"

* Initialize globals
global validation_errors = 0
global validation_tests = 0
global validation_cmd = ""

* Test reporting program
capture program drop report_test
program define report_test
    args test_name passed
    global validation_tests = $validation_tests + 1
    if `passed' {
        di as text "  [PASS] `test_name'"
    }
    else {
        di as error "  [FAIL] `test_name'"
        global validation_errors = $validation_errors + 1
    }
end

* Run a test and report result - simplifies common pattern
capture program drop run_test
program define run_test
    gettoken test_name cmd_to_run : 0
    capture `cmd_to_run'
    local passed = (_rc == 0)
    report_test "`test_name'" `passed'
end

* Run a test expecting a specific error code
capture program drop run_test_error
program define run_test_error
    gettoken test_name rest : 0
    gettoken expected_rc cmd_to_run : rest
    capture `cmd_to_run'
    local passed = (_rc == `expected_rc')
    report_test "`test_name'" `passed'
end

* Print validation header
capture program drop validation_header
program define validation_header
    args cmd_name
    global validation_cmd = "`cmd_name'"
    di ""
    di as text "{hline 70}"
    di as text "VALIDATION SUITE: `cmd_name'"
    di as text "{hline 70}"
    di ""
end

* Print validation summary
capture program drop validation_summary
program define validation_summary
    di ""
    di as text "{hline 70}"
    di as text "VALIDATION SUMMARY: $validation_cmd"
    di as text "{hline 70}"
    di as text "Tests run: $validation_tests"
    di as text "Errors:    $validation_errors"
    if $validation_errors == 0 {
        di as result "ALL TESTS PASSED"
    }
    else {
        di as error "$validation_errors TEST(S) FAILED"
    }
    di as text "{hline 70}"
    di ""
end
