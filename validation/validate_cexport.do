/*******************************************************************************
 * validate_cexport.do
 *
 * Comprehensive validation tests for cexport vs export delimited
 * Tests CSV export functionality across various scenarios
 ******************************************************************************/

do "validation/validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CEXPORT VALIDATION TEST SUITE"
di as text "======================================================================"

* Create temp directory for test files
capture mkdir "validation/temp"

/*******************************************************************************
 * TEST 1: Basic CSV export (auto dataset)
 ******************************************************************************/
print_section "Test 1: Basic CSV export"

sysuse auto, clear

* Stata export
export delimited using "validation/temp/stata_auto.csv", replace

* ctools export
cexport delimited using "validation/temp/cexport_auto.csv", replace

* Reimport and compare
import delimited using "validation/temp/stata_auto.csv", clear
tempfile stata_reimport
quietly save `stata_reimport'

import delimited using "validation/temp/cexport_auto.csv", clear
tempfile cexport_reimport
quietly save `cexport_reimport'

use `stata_reimport', clear
local stata_n = _N
local stata_k = c(k)

use `cexport_reimport', clear
local cexport_n = _N
local cexport_k = c(k)

if `stata_n' == `cexport_n' & `stata_k' == `cexport_k' {
    di as result "  PASS: Basic export dimensions match (N=`stata_n', K=`cexport_k')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Dimensions differ (stata: N=`stata_n' K=`stata_k', cexport: N=`cexport_n' K=`cexport_k')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 2: Semicolon-delimited export
 ******************************************************************************/
print_section "Test 2: Semicolon-delimited export"

sysuse auto, clear
keep make price mpg weight

* Stata export with semicolon delimiter
export delimited using "validation/temp/stata_semi.csv", delimiter(";") replace

* ctools export with semicolon delimiter
cexport delimited using "validation/temp/cexport_semi.csv", delimiter(";") replace

* Reimport and compare
import delimited using "validation/temp/stata_semi.csv", delimiters(";") clear
tempfile stata_reimport
quietly save `stata_reimport'

import delimited using "validation/temp/cexport_semi.csv", delimiters(";") clear
tempfile cexport_reimport
quietly save `cexport_reimport'

use `stata_reimport', clear
local stata_n = _N

use `cexport_reimport', clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    di as result "  PASS: Semicolon-delimited export (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Row counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 3: Export selected variables
 ******************************************************************************/
print_section "Test 3: Export selected variables"

sysuse auto, clear

* Stata export (selected vars)
export delimited make price mpg using "validation/temp/stata_select.csv", replace

* ctools export (selected vars)
cexport delimited make price mpg using "validation/temp/cexport_select.csv", replace

* Reimport and compare
import delimited using "validation/temp/stata_select.csv", clear
local stata_k = c(k)
tempfile stata_reimport
quietly save `stata_reimport'

import delimited using "validation/temp/cexport_select.csv", clear
local cexport_k = c(k)
tempfile cexport_reimport
quietly save `cexport_reimport'

if `stata_k' == `cexport_k' & `stata_k' == 3 {
    di as result "  PASS: Selected variables exported (K=`stata_k')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Variable count incorrect (stata=`stata_k', cexport=`cexport_k', expected=3)"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 4: Export without header (novarnames)
 ******************************************************************************/
print_section "Test 4: Export without header (novarnames)"

sysuse auto, clear
keep in 1/5
keep make price

* Stata export
export delimited using "validation/temp/stata_nohead.csv", novarnames replace

* ctools export
cexport delimited using "validation/temp/cexport_nohead.csv", novarnames replace

* Read first line of each file
file open stata_f using "validation/temp/stata_nohead.csv", read
file read stata_f stata_line1
file close stata_f

file open cexport_f using "validation/temp/cexport_nohead.csv", read
file read cexport_f cexport_line1
file close cexport_f

* First line should not contain "make" or "price" (the header names)
local stata_has_header = (strpos("`stata_line1'", "make") > 0 | strpos("`stata_line1'", "price") > 0)
local cexport_has_header = (strpos("`cexport_line1'", "make") > 0 | strpos("`cexport_line1'", "price") > 0)

if !`stata_has_header' & !`cexport_has_header' {
    di as result "  PASS: novarnames option works correctly"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Header still present with novarnames"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 5: Export with quote option
 ******************************************************************************/
print_section "Test 5: Export with quote option"

sysuse auto, clear
keep in 1/5
keep make price

* Stata export with quote
export delimited using "validation/temp/stata_quote.csv", quote replace

* ctools export with quote
cexport delimited using "validation/temp/cexport_quote.csv", quote replace

* Check that strings are quoted
file open cexport_f using "validation/temp/cexport_quote.csv", read
file read cexport_f line1  // header
file read cexport_f line2  // first data row
file close cexport_f

* The make field should be quoted (starts with ")
if strpos("`line2'", `"""') > 0 {
    di as result "  PASS: quote option adds quotes to strings"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Strings not quoted (line2=`line2')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 6: Export with if condition
 ******************************************************************************/
print_section "Test 6: Export with if condition"

sysuse auto, clear

* Stata export
export delimited using "validation/temp/stata_if.csv" if foreign == 1, replace

* ctools export
cexport delimited using "validation/temp/cexport_if.csv" if foreign == 1, replace

* Reimport and compare counts
import delimited using "validation/temp/stata_if.csv", clear
local stata_n = _N

import delimited using "validation/temp/cexport_if.csv", clear
local cexport_n = _N

* Expected: 22 foreign cars
if `stata_n' == `cexport_n' & `stata_n' == 22 {
    di as result "  PASS: Export with if condition (N=`stata_n' foreign cars)"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Row counts differ or incorrect (stata=`stata_n', cexport=`cexport_n', expected=22)"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 7: Export with in range
 ******************************************************************************/
print_section "Test 7: Export with in range"

sysuse auto, clear

* Stata export
export delimited using "validation/temp/stata_in.csv" in 10/20, replace

* ctools export
cexport delimited using "validation/temp/cexport_in.csv" in 10/20, replace

* Reimport and compare counts
import delimited using "validation/temp/stata_in.csv", clear
local stata_n = _N

import delimited using "validation/temp/cexport_in.csv", clear
local cexport_n = _N

* Expected: 11 rows (10 through 20 inclusive)
if `stata_n' == `cexport_n' & `stata_n' == 11 {
    di as result "  PASS: Export with in range (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Row counts differ (stata=`stata_n', cexport=`cexport_n', expected=11)"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 8: Export numeric precision
 ******************************************************************************/
print_section "Test 8: Export numeric precision"

clear
set obs 5
gen double precise = _n * 1.123456789012345
gen float less_precise = _n * 1.12345
gen int integer_val = _n * 100

* Stata export
export delimited using "validation/temp/stata_numeric.csv", replace

* ctools export
cexport delimited using "validation/temp/cexport_numeric.csv", replace

* Reimport and compare
import delimited using "validation/temp/stata_numeric.csv", clear
local stata_precise1 = precise[1]

import delimited using "validation/temp/cexport_numeric.csv", clear
local cexport_precise1 = precise[1]

if abs(`stata_precise1' - `cexport_precise1') < 1e-10 {
    di as result "  PASS: Numeric precision maintained"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Numeric precision differs (stata=`stata_precise1', cexport=`cexport_precise1')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 9: Export with missing values
 ******************************************************************************/
print_section "Test 9: Export with missing values"

clear
set obs 5
gen id = _n
gen value = _n * 10
replace value = . in 2
replace value = . in 4

* Stata export
export delimited using "validation/temp/stata_missing.csv", replace

* ctools export
cexport delimited using "validation/temp/cexport_missing.csv", replace

* Reimport and check missing values preserved
import delimited using "validation/temp/stata_missing.csv", clear
quietly count if missing(value)
local stata_miss = r(N)

import delimited using "validation/temp/cexport_missing.csv", clear
quietly count if missing(value)
local cexport_miss = r(N)

if `stata_miss' == `cexport_miss' & `stata_miss' == 2 {
    di as result "  PASS: Missing values exported correctly (`stata_miss' missing)"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Missing values differ (stata=`stata_miss', cexport=`cexport_miss', expected=2)"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 10: Export large dataset
 ******************************************************************************/
print_section "Test 10: Export large dataset (10000 rows)"

clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform() * 1000
gen str20 label = "item" + string(runiformint(1, 500))

* Stata export
export delimited using "validation/temp/stata_large.csv", replace

* ctools export
cexport delimited using "validation/temp/cexport_large.csv", replace

* Reimport and compare
import delimited using "validation/temp/stata_large.csv", clear
quietly summarize value
local stata_n = _N
local stata_mean = r(mean)

import delimited using "validation/temp/cexport_large.csv", clear
quietly summarize value
local cexport_n = _N
local cexport_mean = r(mean)

if `stata_n' == `cexport_n' & abs(`stata_mean' - `cexport_mean') < 0.01 {
    di as result "  PASS: Large dataset export (N=`stata_n', mean value=`stata_mean')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Large dataset export differs"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 11: Round-trip integrity (census dataset)
 ******************************************************************************/
print_section "Test 11: Round-trip integrity (census)"

sysuse census, clear
local orig_n = _N
quietly summarize pop
local orig_pop_mean = r(mean)
local orig_pop_sd = r(sd)

* ctools export
cexport delimited using "validation/temp/cexport_census.csv", replace

* Reimport
import delimited using "validation/temp/cexport_census.csv", clear
local reimport_n = _N
quietly summarize pop
local reimport_pop_mean = r(mean)
local reimport_pop_sd = r(sd)

local n_match = (`orig_n' == `reimport_n')
local mean_match = (abs(`orig_pop_mean' - `reimport_pop_mean') < 1)
local sd_match = (abs(`orig_pop_sd' - `reimport_pop_sd') < 1)

if `n_match' & `mean_match' & `sd_match' {
    di as result "  PASS: Round-trip integrity maintained"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Round-trip integrity lost"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 12: String variables with special characters
 ******************************************************************************/
print_section "Test 12: String variables with special characters"

clear
set obs 3
gen str50 name = ""
replace name = "Normal" in 1
replace name = "With,Comma" in 2
replace name = "With Space" in 3
gen value = _n * 10

* ctools export
cexport delimited using "validation/temp/cexport_special.csv", replace

* Reimport
import delimited using "validation/temp/cexport_special.csv", clear

local name2 = name[2]
if strpos("`name2'", ",") > 0 | "`name2'" == "With,Comma" {
    di as result "  PASS: Special characters in strings preserved"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Special characters not preserved (got: `name2')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 13: Export nlswork panel data
 ******************************************************************************/
print_section "Test 13: Export panel data (nlswork)"

webuse nlswork, clear
keep in 1/1000

* Stata export
export delimited using "validation/temp/stata_nlswork.csv", replace

* ctools export
cexport delimited using "validation/temp/cexport_nlswork.csv", replace

* Reimport and compare
import delimited using "validation/temp/stata_nlswork.csv", clear
local stata_n = _N

import delimited using "validation/temp/cexport_nlswork.csv", clear
local cexport_n = _N

if `stata_n' == `cexport_n' {
    di as result "  PASS: Panel data export (N=`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Panel data export differs (stata=`stata_n', cexport=`cexport_n')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 14: Export with negative numbers
 ******************************************************************************/
print_section "Test 14: Export with negative numbers"

clear
set obs 5
gen id = _n
gen neg_int = -1 * _n * 100
gen neg_float = -1 * _n * 1.5

* ctools export
cexport delimited using "validation/temp/cexport_negative.csv", replace

* Reimport
import delimited using "validation/temp/cexport_negative.csv", clear
local neg1 = neg_int[3]

if `neg1' == -300 {
    di as result "  PASS: Negative numbers exported correctly"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Negative number incorrect (got `neg1', expected -300)"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * Cleanup temp files
 ******************************************************************************/
local files : dir "validation/temp" files "*.csv"
foreach f of local files {
    capture erase "validation/temp/`f'"
}
local files : dir "validation/temp" files "*.tsv"
foreach f of local files {
    capture erase "validation/temp/`f'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cexport"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
