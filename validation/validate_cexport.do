/*******************************************************************************
 * validate_cexport.do
 *
 * Comprehensive validation tests for cexport vs export delimited
 * Tests CSV export functionality across various scenarios
 ******************************************************************************/

do "validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CEXPORT VALIDATION TEST SUITE"
di as text "======================================================================"

* Create temp directory for test files
capture mkdir "temp"

/*******************************************************************************
 * Basic export tests
 ******************************************************************************/
print_section "Basic export tests"

sysuse auto, clear
benchmark_export using "temp/test.csv", testname("auto dataset")

sysuse census, clear
benchmark_export using "temp/test.csv", testname("census dataset")

sysuse voter, clear
benchmark_export using "temp/test.csv", testname("voter dataset")

sysuse sp500, clear
benchmark_export using "temp/test.csv", testname("sp500 dataset")

sysuse uslifeexp, clear
benchmark_export using "temp/test.csv", testname("uslifeexp dataset")

/*******************************************************************************
 * Selected variables
 ******************************************************************************/
print_section "Selected variables"

sysuse auto, clear
benchmark_export make price mpg using "temp/test.csv", testname("three variables")

sysuse auto, clear
benchmark_export make using "temp/test.csv", testname("single string var")

sysuse auto, clear
benchmark_export price using "temp/test.csv", testname("single numeric var")

/*******************************************************************************
 * Delimiter options
 ******************************************************************************/
print_section "Delimiter options"

sysuse auto, clear
keep make price mpg weight
benchmark_export using "temp/test.csv", delimiter(";") testname("semicolon delimiter")

/*******************************************************************************
 * Other options
 ******************************************************************************/
print_section "Other options"

sysuse auto, clear
keep in 1/5
keep make price
benchmark_export using "temp/test.csv", novarnames testname("novarnames")

sysuse auto, clear
keep in 1/5
keep make price
benchmark_export using "temp/test.csv", quote testname("quote")

/*******************************************************************************
 * nolabel option
 ******************************************************************************/
print_section "nolabel option"

sysuse auto, clear
keep in 1/5
keep make price foreign

* Test with labels (default)
export delimited using "temp/stata_label.csv", replace
cexport delimited using "temp/cexport_label.csv", replace

import delimited using "temp/stata_label.csv", clear stringcols(3)
local stata_foreign1 = foreign[1]

import delimited using "temp/cexport_label.csv", clear stringcols(3)
local cexport_foreign1 = foreign[1]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if "`stata_foreign1'" == "`cexport_foreign1'" {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] value labels exported (default)"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] value labels exported (default)"
}

* Test with nolabel
sysuse auto, clear
keep in 1/5
keep make price foreign

export delimited using "temp/stata_nolabel.csv", nolabel replace
cexport delimited using "temp/cexport_nolabel.csv", nolabel replace

import delimited using "temp/stata_nolabel.csv", clear
local stata_foreign1 = foreign[1]

import delimited using "temp/cexport_nolabel.csv", clear
local cexport_foreign1 = foreign[1]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_foreign1' == `cexport_foreign1' {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] nolabel option"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] nolabel option"
}

/*******************************************************************************
 * Tab delimiter
 ******************************************************************************/
print_section "Tab delimiter"

sysuse auto, clear
keep in 1/10
keep make price mpg

export delimited using "temp/stata_tab.tsv", delimiter(tab) replace
cexport delimited using "temp/cexport_tab.tsv", delimiter(tab) replace

import delimited using "temp/stata_tab.tsv", delimiters(tab) clear
local stata_n = _N
local stata_k = c(k)

import delimited using "temp/cexport_tab.tsv", delimiters(tab) clear
local cexport_n = _N
local cexport_k = c(k)

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_n' == `cexport_n' & `stata_k' == `cexport_k' {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] tab delimiter"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] tab delimiter (stata: N=`stata_n' K=`stata_k', cexport: N=`cexport_n' K=`cexport_k')"
}

/*******************************************************************************
 * if/in conditions
 ******************************************************************************/
print_section "if/in conditions"

* Export with if condition
sysuse auto, clear
export delimited using "temp/stata_if.csv" if foreign == 1, replace
cexport delimited using "temp/cexport_if.csv" if foreign == 1, replace

import delimited using "temp/stata_if.csv", clear
local stata_n = _N

import delimited using "temp/cexport_if.csv", clear
local cexport_n = _N

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_n' == `cexport_n' & `stata_n' == 22 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] if condition (foreign==1)"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] if condition"
}

* Export with in range
sysuse auto, clear
export delimited using "temp/stata_in.csv" in 10/20, replace
cexport delimited using "temp/cexport_in.csv" in 10/20, replace

import delimited using "temp/stata_in.csv", clear
local stata_n = _N

import delimited using "temp/cexport_in.csv", clear
local cexport_n = _N

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_n' == `cexport_n' & `stata_n' == 11 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] in range (10/20)"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] in range"
}

/*******************************************************************************
 * Special data types
 ******************************************************************************/
print_section "Special data types"

* Negative numbers
clear
set obs 5
gen id = _n
gen neg_int = -1 * _n * 100
gen neg_float = -1 * _n * 1.5

cexport delimited using "temp/cexport_neg.csv", replace
import delimited using "temp/cexport_neg.csv", clear
local neg1 = neg_int[3]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `neg1' == -300 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] negative numbers"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] negative numbers"
}

* Missing values
clear
set obs 5
gen id = _n
gen value = _n * 10
replace value = . in 2
replace value = . in 4

export delimited using "temp/stata_miss.csv", replace
cexport delimited using "temp/cexport_miss.csv", replace

import delimited using "temp/stata_miss.csv", clear
quietly count if missing(value)
local stata_miss = r(N)

import delimited using "temp/cexport_miss.csv", clear
quietly count if missing(value)
local cexport_miss = r(N)

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `stata_miss' == `cexport_miss' & `stata_miss' == 2 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] missing values"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] missing values"
}

* Numeric precision
clear
set obs 5
gen double precise = _n * 1.123456789012345

export delimited using "temp/stata_prec.csv", replace
cexport delimited using "temp/cexport_prec.csv", replace

import delimited using "temp/stata_prec.csv", clear
local stata_p1 = precise[1]

import delimited using "temp/cexport_prec.csv", clear
local cexport_p1 = precise[1]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if abs(`stata_p1' - `cexport_p1') < 1e-10 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] numeric precision"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] numeric precision"
}

/*******************************************************************************
 * Large dataset
 ******************************************************************************/
print_section "Large dataset"

clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform() * 1000
gen str20 label = "item" + string(runiformint(1, 500))

benchmark_export using "temp/test.csv", testname("10K rows")

* Verify statistics match
export delimited using "temp/stata_large.csv", replace
cexport delimited using "temp/cexport_large.csv", replace

import delimited using "temp/stata_large.csv", clear
quietly summarize value
local stata_mean = r(mean)

import delimited using "temp/cexport_large.csv", clear
quietly summarize value
local cexport_mean = r(mean)

global TESTS_TOTAL = $TESTS_TOTAL + 1
if abs(`stata_mean' - `cexport_mean') < 0.01 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] large dataset statistics"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] large dataset statistics"
}

/*******************************************************************************
 * Round-trip integrity
 ******************************************************************************/
print_section "Round-trip integrity"

sysuse census, clear
local orig_n = _N
quietly summarize pop
local orig_mean = r(mean)

cexport delimited using "temp/census_rt.csv", replace
import delimited using "temp/census_rt.csv", clear
local reimport_n = _N
quietly summarize pop
local reimport_mean = r(mean)

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `orig_n' == `reimport_n' & abs(`orig_mean' - `reimport_mean') < 1 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] round-trip integrity"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] round-trip integrity"
}

/*******************************************************************************
 * Panel data
 ******************************************************************************/
print_section "Panel data"

webuse nlswork, clear
keep in 1/1000
benchmark_export using "temp/test.csv", testname("nlswork panel")

/*******************************************************************************
 * Cleanup temp files
 ******************************************************************************/
local files : dir "temp" files "*.csv"
foreach f of local files {
    capture erase "temp/`f'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cexport"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
