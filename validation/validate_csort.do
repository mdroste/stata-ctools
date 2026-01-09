/*******************************************************************************
 * validate_csort.do
 *
 * Comprehensive validation tests for csort vs native sort
 * Tests sorting functionality across various data types and scenarios
 *
 * Note: csort uses inherently stable radix sort, so we compare against
 * Stata's "sort, stable" for exact matches.
 ******************************************************************************/

do "validation/validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CSORT VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * First, test if csort works at all
 ******************************************************************************/
print_section "Preliminary: Check if csort plugin works"

sysuse auto, clear
capture noisily csort price
local csort_works = (_rc == 0)

if !`csort_works' {
    di as error "  csort plugin is not functioning (rc=`=_rc')"
    di as error "  Skipping all csort tests"
    di as text ""
    di as text "{hline 70}"
    di as text "SUMMARY: csort"
    di as text "{hline 70}"
    di as text "Tests passed: 0"
    di as text "Tests failed: 1"
    di as text "Total tests:  1"
    di as error "VALIDATION FAILED - csort plugin not working"
    di as text "{hline 70}"
    global TESTS_PASSED = 0
    global TESTS_FAILED = 1
    global TESTS_TOTAL = 1
    exit 1
}

di as result "  csort plugin is functional"
global TESTS_PASSED = $TESTS_PASSED + 1
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 1: Single numeric variable sort (auto dataset)
 ******************************************************************************/
print_section "Test 1: Single numeric variable sort"

sysuse auto, clear
preserve
sort price, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort price
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Single numeric sort (price)"

/*******************************************************************************
 * TEST 2: Single string variable sort
 ******************************************************************************/
print_section "Test 2: Single string variable sort"

sysuse auto, clear
preserve
sort make, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort make
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Single string sort (make)"

/*******************************************************************************
 * TEST 3: Multiple variable sort
 ******************************************************************************/
print_section "Test 3: Multiple variable sort"

sysuse auto, clear
preserve
sort foreign rep78 price, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort foreign rep78 price
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Multi-variable sort (foreign rep78 price)"

/*******************************************************************************
 * TEST 4: Sort with missing values
 ******************************************************************************/
print_section "Test 4: Sort with missing values"

sysuse auto, clear
* rep78 has missing values
preserve
sort rep78, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort rep78
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Sort with missing values (rep78)"

/*******************************************************************************
 * TEST 5: Sort with duplicates
 ******************************************************************************/
print_section "Test 5: Sort with duplicates"

sysuse auto, clear
* foreign has many duplicates (0 and 1)
preserve
sort foreign, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort foreign
tempfile csort_sorted
quietly save `csort_sorted'
restore

* Check that all foreign=0 come before foreign=1
use `csort_sorted', clear
gen _order = _n
quietly summarize _order if foreign == 0
local max_domestic = r(max)
quietly summarize _order if foreign == 1
local min_foreign = r(min)

if `max_domestic' < `min_foreign' {
    di as result "  PASS: Duplicates sorted correctly (domestic before foreign)"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Duplicates not sorted correctly"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 6: Large synthetic dataset
 ******************************************************************************/
print_section "Test 6: Large synthetic dataset (10,000 obs)"

clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform()

* Sort by group and value
preserve
sort group value, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort group value
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Large dataset sort (10K obs)"

/*******************************************************************************
 * TEST 7: Negative numbers
 ******************************************************************************/
print_section "Test 7: Negative numbers"

clear
set obs 1000
gen x = runiform() * 200 - 100  // Range -100 to 100

preserve
sort x, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort x
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Negative numbers sort"

/*******************************************************************************
 * TEST 8: Census dataset
 ******************************************************************************/
print_section "Test 8: Census dataset"

sysuse census, clear

preserve
sort state, stable
tempfile stata_sorted
quietly save `stata_sorted'
restore

preserve
csort state
tempfile csort_sorted
quietly save `csort_sorted'
restore

assert_data_equal `stata_sorted' `csort_sorted' "Census state sort"

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "csort"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
