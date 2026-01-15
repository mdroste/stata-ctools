/*******************************************************************************
 * validate_all.do
 *
 * Master test runner for ctools validation suite
 * Runs all individual validation tests and provides summary
 *
 * Usage:
 *   stata -b do validate_all.do
 *
 * Individual tests can also be run separately:
 *   stata -b do validate_csort.do
 *   stata -b do validate_cmerge.do
 *   stata -b do validate_cimport.do
 *   stata -b do validate_cexport.do
 *   stata -b do validate_creghdfe.do
 *   stata -b do validate_cqreg.do
 *   stata -b do validate_civreghdfe.do
 ******************************************************************************/

version 14.0
clear all
set more off
set trace off

quietly {

noi di as text ""
noi di as text "######################################################################"
noi di as text "#                                                                    #"
noi di as text "#               CTOOLS VALIDATION TEST SUITE                         #"
noi di as text "#                                                                    #"
noi di as text "######################################################################"
noi di as text ""
noi di as text "Running comprehensive validation tests for all ctools commands."
noi di as text "Each command is compared against its native Stata equivalent."
noi di as text ""
noi di as text "Date: " c(current_date) " " c(current_time)
noi di as text "Stata version: " c(stata_version)
noi di as text "Platform: " c(os) " " c(machine_type)
noi di as text ""

* Track overall results
local total_passed = 0
local total_failed = 0
local total_tests = 0

* Track all failures across all components
local all_failures = ""

* Create temp directory
capture mkdir "temp"

/*******************************************************************************
 * RUN CSORT VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: csort validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_csort.do"
local csort_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local csort_passed = $TESTS_PASSED
local csort_failed = $TESTS_FAILED
local csort_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [csort] ${FAILURE_`i'}"
        }
    }
}

if `csort_rc' == 0 & `csort_failed' == 0 {
    noi di as result "csort: ALL TESTS PASSED (`csort_passed'/`csort_total')"
}
else {
    noi di as error "csort: FAILED (`csort_failed' of `csort_total' tests failed)"
}

/*******************************************************************************
 * RUN CMERGE VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: cmerge validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_cmerge.do"
local cmerge_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cmerge_passed = $TESTS_PASSED
local cmerge_failed = $TESTS_FAILED
local cmerge_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [cmerge] ${FAILURE_`i'}"
        }
    }
}

if `cmerge_rc' == 0 & `cmerge_failed' == 0 {
    noi di as result "cmerge: ALL TESTS PASSED (`cmerge_passed'/`cmerge_total')"
}
else {
    noi di as error "cmerge: FAILED (`cmerge_failed' of `cmerge_total' tests failed)"
}

/*******************************************************************************
 * RUN CIMPORT VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: cimport validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_cimport.do"
local cimport_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cimport_passed = $TESTS_PASSED
local cimport_failed = $TESTS_FAILED
local cimport_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [cimport] ${FAILURE_`i'}"
        }
    }
}

if `cimport_rc' == 0 & `cimport_failed' == 0 {
    noi di as result "cimport: ALL TESTS PASSED (`cimport_passed'/`cimport_total')"
}
else {
    noi di as error "cimport: FAILED (`cimport_failed' of `cimport_total' tests failed)"
}

/*******************************************************************************
 * RUN CEXPORT VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: cexport validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_cexport.do"
local cexport_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cexport_passed = $TESTS_PASSED
local cexport_failed = $TESTS_FAILED
local cexport_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [cexport] ${FAILURE_`i'}"
        }
    }
}

if `cexport_rc' == 0 & `cexport_failed' == 0 {
    noi di as result "cexport: ALL TESTS PASSED (`cexport_passed'/`cexport_total')"
}
else {
    noi di as error "cexport: FAILED (`cexport_failed' of `cexport_total' tests failed)"
}

/*******************************************************************************
 * RUN CREGHDFE VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: creghdfe validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_creghdfe.do"
local creghdfe_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local creghdfe_passed = $TESTS_PASSED
local creghdfe_failed = $TESTS_FAILED
local creghdfe_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [creghdfe] ${FAILURE_`i'}"
        }
    }
}

if `creghdfe_rc' == 0 & `creghdfe_failed' == 0 {
    noi di as result "creghdfe: ALL TESTS PASSED (`creghdfe_passed'/`creghdfe_total')"
}
else {
    noi di as error "creghdfe: FAILED (`creghdfe_failed' of `creghdfe_total' tests failed)"
}

/*******************************************************************************
 * RUN CQREG VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: cqreg validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_cqreg.do"
local cqreg_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cqreg_passed = $TESTS_PASSED
local cqreg_failed = $TESTS_FAILED
local cqreg_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [cqreg] ${FAILURE_`i'}"
        }
    }
}

if `cqreg_rc' == 0 & `cqreg_failed' == 0 {
    noi di as result "cqreg: ALL TESTS PASSED (`cqreg_passed'/`cqreg_total')"
}
else {
    noi di as error "cqreg: FAILED (`cqreg_failed' of `cqreg_total' tests failed)"
}

/*******************************************************************************
 * RUN CIVREGHDFE VALIDATION
 ******************************************************************************/
noi di as text ""
noi di as text "======================================================================"
noi di as text "  Running: civreghdfe validation tests"
noi di as text "======================================================================"

capture noisily do "validation/validate_civreghdfe.do"
local civreghdfe_rc = _rc

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local civreghdfe_passed = $TESTS_PASSED
local civreghdfe_failed = $TESTS_FAILED
local civreghdfe_total = $TESTS_TOTAL

* Collect failures
if $FAILURE_COUNT > 0 {
    forvalues i = 1/$FAILURE_COUNT {
        if `i' <= 100 {
            local all_failures = `"`all_failures'"' + char(10) + "  [civreghdfe] ${FAILURE_`i'}"
        }
    }
}

if `civreghdfe_rc' == 0 & `civreghdfe_failed' == 0 {
    noi di as result "civreghdfe: ALL TESTS PASSED (`civreghdfe_passed'/`civreghdfe_total')"
}
else {
    noi di as error "civreghdfe: FAILED (`civreghdfe_failed' of `civreghdfe_total' tests failed)"
}

/*******************************************************************************
 * FINAL SUMMARY
 ******************************************************************************/
noi di as text ""
noi di as text "######################################################################"
noi di as text "#                     VALIDATION SUMMARY                             #"
noi di as text "######################################################################"
noi di as text ""

* Show all failures if any
if `total_failed' > 0 {
    noi di as error "{hline 70}"
    noi di as error "ALL FAILED TESTS:"
    noi di as error "{hline 70}"
    noi di as error `"`all_failures'"'
    noi di as error "{hline 70}"
    noi di as text ""
}

noi di as text "Component Results:"
noi di as text "{hline 55}"
noi di as text "  Command       Passed   Failed   Total    Status"
noi di as text "{hline 55}"

* csort
if `csort_failed' == 0 {
    noi di as text "  csort      " %8.0f `csort_passed' %8.0f `csort_failed' %8.0f `csort_total' "    " as result "PASS"
}
else {
    noi di as text "  csort      " %8.0f `csort_passed' %8.0f `csort_failed' %8.0f `csort_total' "    " as error "FAIL"
}

* cmerge
if `cmerge_failed' == 0 {
    noi di as text "  cmerge     " %8.0f `cmerge_passed' %8.0f `cmerge_failed' %8.0f `cmerge_total' "    " as result "PASS"
}
else {
    noi di as text "  cmerge     " %8.0f `cmerge_passed' %8.0f `cmerge_failed' %8.0f `cmerge_total' "    " as error "FAIL"
}

* cimport
if `cimport_failed' == 0 {
    noi di as text "  cimport    " %8.0f `cimport_passed' %8.0f `cimport_failed' %8.0f `cimport_total' "    " as result "PASS"
}
else {
    noi di as text "  cimport    " %8.0f `cimport_passed' %8.0f `cimport_failed' %8.0f `cimport_total' "    " as error "FAIL"
}

* cexport
if `cexport_failed' == 0 {
    noi di as text "  cexport    " %8.0f `cexport_passed' %8.0f `cexport_failed' %8.0f `cexport_total' "    " as result "PASS"
}
else {
    noi di as text "  cexport    " %8.0f `cexport_passed' %8.0f `cexport_failed' %8.0f `cexport_total' "    " as error "FAIL"
}

* creghdfe
if `creghdfe_failed' == 0 {
    noi di as text "  creghdfe   " %8.0f `creghdfe_passed' %8.0f `creghdfe_failed' %8.0f `creghdfe_total' "    " as result "PASS"
}
else {
    noi di as text "  creghdfe   " %8.0f `creghdfe_passed' %8.0f `creghdfe_failed' %8.0f `creghdfe_total' "    " as error "FAIL"
}

* cqreg
if `cqreg_failed' == 0 {
    noi di as text "  cqreg      " %8.0f `cqreg_passed' %8.0f `cqreg_failed' %8.0f `cqreg_total' "    " as result "PASS"
}
else {
    noi di as text "  cqreg      " %8.0f `cqreg_passed' %8.0f `cqreg_failed' %8.0f `cqreg_total' "    " as error "FAIL"
}

* civreghdfe
if `civreghdfe_failed' == 0 {
    noi di as text "  civreghdfe " %8.0f `civreghdfe_passed' %8.0f `civreghdfe_failed' %8.0f `civreghdfe_total' "    " as result "PASS"
}
else {
    noi di as text "  civreghdfe " %8.0f `civreghdfe_passed' %8.0f `civreghdfe_failed' %8.0f `civreghdfe_total' "    " as error "FAIL"
}

noi di as text "{hline 55}"
noi di as text "  TOTAL      " %8.0f `total_passed' %8.0f `total_failed' %8.0f `total_tests'
noi di as text "{hline 55}"
noi di as text ""

* Overall result
if `total_failed' == 0 {
    noi di as result "######################################################################"
    noi di as result "#                                                                    #"
    noi di as result "#                   ALL VALIDATION TESTS PASSED                      #"
    noi di as result "#                                                                    #"
    noi di as result "######################################################################"
}
else {
    noi di as error "######################################################################"
    noi di as error "#                                                                    #"
    noi di as error "#                   SOME VALIDATION TESTS FAILED                     #"
    noi di as error "#                                                                    #"
    noi di as error "######################################################################"
    noi di as error ""
    noi di as error "See 'ALL FAILED TESTS' section above for details."
}

noi di as text ""
noi di as text "Completed: " c(current_date) " " c(current_time)
noi di as text ""

}

* Return error code if any tests failed
if `total_failed' > 0 {
    exit 1
}
