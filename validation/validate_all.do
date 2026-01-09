/*******************************************************************************
 * validate_all.do
 *
 * Master test runner for ctools validation suite
 * Runs all individual validation tests and provides summary
 *
 * Usage:
 *   stata -b do validation/validate_all.do
 *
 * Individual tests can also be run separately:
 *   stata -b do validation/validate_csort.do
 *   stata -b do validation/validate_cmerge.do
 *   etc.
 ******************************************************************************/

version 14.0
clear all
set more off
set trace off

di as text ""
di as text "######################################################################"
di as text "#                                                                    #"
di as text "#               CTOOLS VALIDATION TEST SUITE                         #"
di as text "#                                                                    #"
di as text "######################################################################"
di as text ""
di as text "Running comprehensive validation tests for all ctools commands."
di as text "Each command is compared against its native Stata equivalent."
di as text ""
di as text "Date: " c(current_date) " " c(current_time)
di as text "Stata version: " c(stata_version)
di as text "Platform: " c(os) " " c(machine_type)
di as text ""

* Track overall results
local total_passed = 0
local total_failed = 0
local total_tests = 0

* Create temp directory
capture mkdir "validation/temp"

/*******************************************************************************
 * RUN CSORT VALIDATION
 ******************************************************************************/
di as text ""
di as text "======================================================================"
di as text "  Running: csort validation tests"
di as text "======================================================================"

capture noisily do "validation/validate_csort.do"
local csort_rc = _rc

if `csort_rc' == 0 {
    di as result "csort: ALL TESTS PASSED"
}
else {
    di as error "csort: SOME TESTS FAILED"
}

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local csort_passed = $TESTS_PASSED
local csort_failed = $TESTS_FAILED
local csort_total = $TESTS_TOTAL

/*******************************************************************************
 * RUN CMERGE VALIDATION
 ******************************************************************************/
di as text ""
di as text "======================================================================"
di as text "  Running: cmerge validation tests"
di as text "======================================================================"

capture noisily do "validation/validate_cmerge.do"
local cmerge_rc = _rc

if `cmerge_rc' == 0 {
    di as result "cmerge: ALL TESTS PASSED"
}
else {
    di as error "cmerge: SOME TESTS FAILED"
}

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cmerge_passed = $TESTS_PASSED
local cmerge_failed = $TESTS_FAILED
local cmerge_total = $TESTS_TOTAL

/*******************************************************************************
 * RUN CIMPORT VALIDATION
 ******************************************************************************/
di as text ""
di as text "======================================================================"
di as text "  Running: cimport validation tests"
di as text "======================================================================"

capture noisily do "validation/validate_cimport.do"
local cimport_rc = _rc

if `cimport_rc' == 0 {
    di as result "cimport: ALL TESTS PASSED"
}
else {
    di as error "cimport: SOME TESTS FAILED"
}

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cimport_passed = $TESTS_PASSED
local cimport_failed = $TESTS_FAILED
local cimport_total = $TESTS_TOTAL

/*******************************************************************************
 * RUN CEXPORT VALIDATION
 ******************************************************************************/
di as text ""
di as text "======================================================================"
di as text "  Running: cexport validation tests"
di as text "======================================================================"

capture noisily do "validation/validate_cexport.do"
local cexport_rc = _rc

if `cexport_rc' == 0 {
    di as result "cexport: ALL TESTS PASSED"
}
else {
    di as error "cexport: SOME TESTS FAILED"
}

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cexport_passed = $TESTS_PASSED
local cexport_failed = $TESTS_FAILED
local cexport_total = $TESTS_TOTAL

/*******************************************************************************
 * RUN CREGHDFE VALIDATION
 ******************************************************************************/
di as text ""
di as text "======================================================================"
di as text "  Running: creghdfe validation tests"
di as text "======================================================================"

capture noisily do "validation/validate_creghdfe.do"
local creghdfe_rc = _rc

if `creghdfe_rc' == 0 {
    di as result "creghdfe: ALL TESTS PASSED"
}
else {
    di as error "creghdfe: SOME TESTS FAILED"
}

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local creghdfe_passed = $TESTS_PASSED
local creghdfe_failed = $TESTS_FAILED
local creghdfe_total = $TESTS_TOTAL

/*******************************************************************************
 * RUN CQREG VALIDATION
 ******************************************************************************/
di as text ""
di as text "======================================================================"
di as text "  Running: cqreg validation tests"
di as text "======================================================================"

capture noisily do "validation/validate_cqreg.do"
local cqreg_rc = _rc

if `cqreg_rc' == 0 {
    di as result "cqreg: ALL TESTS PASSED"
}
else {
    di as error "cqreg: SOME TESTS FAILED"
}

local total_passed = `total_passed' + $TESTS_PASSED
local total_failed = `total_failed' + $TESTS_FAILED
local total_tests = `total_tests' + $TESTS_TOTAL

local cqreg_passed = $TESTS_PASSED
local cqreg_failed = $TESTS_FAILED
local cqreg_total = $TESTS_TOTAL

/*******************************************************************************
 * FINAL SUMMARY
 ******************************************************************************/
di as text ""
di as text "######################################################################"
di as text "#                     VALIDATION SUMMARY                             #"
di as text "######################################################################"
di as text ""
di as text "Component Results:"
di as text "{hline 50}"
di as text "  Command      Passed   Failed   Total    Status"
di as text "{hline 50}"

* csort
if `csort_failed' == 0 {
    di as text "  csort    " %8.0f `csort_passed' %8.0f `csort_failed' %8.0f `csort_total' "    " as result "PASS"
}
else {
    di as text "  csort    " %8.0f `csort_passed' %8.0f `csort_failed' %8.0f `csort_total' "    " as error "FAIL"
}

* cmerge
if `cmerge_failed' == 0 {
    di as text "  cmerge   " %8.0f `cmerge_passed' %8.0f `cmerge_failed' %8.0f `cmerge_total' "    " as result "PASS"
}
else {
    di as text "  cmerge   " %8.0f `cmerge_passed' %8.0f `cmerge_failed' %8.0f `cmerge_total' "    " as error "FAIL"
}

* cimport
if `cimport_failed' == 0 {
    di as text "  cimport  " %8.0f `cimport_passed' %8.0f `cimport_failed' %8.0f `cimport_total' "    " as result "PASS"
}
else {
    di as text "  cimport  " %8.0f `cimport_passed' %8.0f `cimport_failed' %8.0f `cimport_total' "    " as error "FAIL"
}

* cexport
if `cexport_failed' == 0 {
    di as text "  cexport  " %8.0f `cexport_passed' %8.0f `cexport_failed' %8.0f `cexport_total' "    " as result "PASS"
}
else {
    di as text "  cexport  " %8.0f `cexport_passed' %8.0f `cexport_failed' %8.0f `cexport_total' "    " as error "FAIL"
}

* creghdfe
if `creghdfe_failed' == 0 {
    di as text "  creghdfe " %8.0f `creghdfe_passed' %8.0f `creghdfe_failed' %8.0f `creghdfe_total' "    " as result "PASS"
}
else {
    di as text "  creghdfe " %8.0f `creghdfe_passed' %8.0f `creghdfe_failed' %8.0f `creghdfe_total' "    " as error "FAIL"
}

* cqreg
if `cqreg_failed' == 0 {
    di as text "  cqreg    " %8.0f `cqreg_passed' %8.0f `cqreg_failed' %8.0f `cqreg_total' "    " as result "PASS"
}
else {
    di as text "  cqreg    " %8.0f `cqreg_passed' %8.0f `cqreg_failed' %8.0f `cqreg_total' "    " as error "FAIL"
}

di as text "{hline 50}"
di as text "  TOTAL    " %8.0f `total_passed' %8.0f `total_failed' %8.0f `total_tests'
di as text "{hline 50}"
di as text ""

* Overall result
if `total_failed' == 0 {
    di as result "######################################################################"
    di as result "#                                                                    #"
    di as result "#                   ALL VALIDATION TESTS PASSED                      #"
    di as result "#                                                                    #"
    di as result "######################################################################"
}
else {
    di as error "######################################################################"
    di as error "#                                                                    #"
    di as error "#                   SOME VALIDATION TESTS FAILED                     #"
    di as error "#                                                                    #"
    di as error "######################################################################"
}

di as text ""
di as text "Completed: " c(current_date) " " c(current_time)
di as text ""

* Return error code if any tests failed
if `total_failed' > 0 {
    exit 1
}
