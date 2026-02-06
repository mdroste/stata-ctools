/*******************************************************************************
 * validate_all.do
 *
 * Master test runner for ctools validation suite
 * Minimal output: only summary table and failures
 *
 * Usage: runstata validation/validate_all.do
 ******************************************************************************/


quietly {
	
clear all 
set more off

* Track overall results
local total_passed = 0
local total_failed = 0
local total_tests = 0

* Track all failures across all components (with component prefix)
local all_failure_count = 0

* Create temp directory
capture mkdir "temp"

* Define list of validation scripts and their names
local scripts "csort cmerge cimport cexport creghdfe cqreg civreghdfe cdecode cencode cdestring csample cbsample cbinscatter cpsmatch"
local script_names "csort cmerge cimport cexport creghdfe cqreg civreghdfe cdecode cencode cdestring csample cbsample cbinscatter cpsmatch"
local num_scripts : word count `scripts'

* Count total validations (number of scripts)
noi di as text ""
noi di as text "Running all validation scripts..."
noi di as text ""

* Initialize storage for per-script results
forvalues i = 1/`num_scripts' {
    local passed_`i' = 0
    local failed_`i' = 0
    local total_`i' = 0
}

* Run each validation script
local script_idx = 1
foreach script of local scripts {
    noi di as text "  Testing `script'..."
    * Run the script silently
    capture quietly do "validation/validate_`script'.do"
    local rc_`script_idx' = _rc

    * Capture results
    local passed_`script_idx' = $TESTS_PASSED
    local failed_`script_idx' = $TESTS_FAILED
    local total_`script_idx' = $TESTS_TOTAL

    * Accumulate totals
    local total_passed = `total_passed' + $TESTS_PASSED
    local total_failed = `total_failed' + $TESTS_FAILED
    local total_tests = `total_tests' + $TESTS_TOTAL

    * Collect failures with component prefix
    local name : word `script_idx' of `script_names'
    if $FAILURE_COUNT > 0 {
        forvalues i = 1/$FAILURE_COUNT {
            if `i' <= 100 {
                local all_failure_count = `all_failure_count' + 1
                if `all_failure_count' <= 100 {
                    local all_failure_`all_failure_count' = "[`name'] ${FAILURE_`i'}"
                }
            }
        }
    }

    local script_idx = `script_idx' + 1
}

* Print summary table
noi di as text "{hline 55}"
noi di as text "  Command       Passed   Failed   Total    Status"
noi di as text "{hline 55}"

local script_idx = 1
foreach name of local script_names {
    local p = `passed_`script_idx''
    local f = `failed_`script_idx''
    local t = `total_`script_idx''

    * Pad name to 12 chars
    local padded_name = substr("`name'" + "            ", 1, 12)

    if `f' == 0 {
        noi di as text "  `padded_name'" %8.0f `p' %8.0f `f' %8.0f `t' "    " as result "PASS"
    }
    else {
        noi di as text "  `padded_name'" %8.0f `p' %8.0f `f' %8.0f `t' "    " as error "FAIL"
    }
    local script_idx = `script_idx' + 1
}

noi di as text "{hline 55}"
noi di as text "  TOTAL       " %8.0f `total_passed' %8.0f `total_failed' %8.0f `total_tests'
noi di as text "{hline 55}"
noi di as text ""

* Print failures if any
if `total_failed' > 0 {
    noi di as error "Failed Tests:"
    noi di as error "{hline 70}"

    local max_to_show = min(`all_failure_count', 100)
    forvalues i = 1/`max_to_show' {
        noi di as error "  `all_failure_`i''"
    }

    if `all_failure_count' > 100 {
        local more = `all_failure_count' - 100
        noi di as error ""
        noi di as error "  ... and `more' more failures not shown"
    }
    noi di as error "{hline 70}"
    noi di as text ""
}

* Final status
if `total_failed' == 0 {
    noi di as result "ALL TESTS PASSED"
}
else {
    noi di as error "VALIDATION FAILED"
}
noi di as text ""

} /* end quietly */

* Return error code if any tests failed
if `total_failed' > 0 {
    exit 1
}
