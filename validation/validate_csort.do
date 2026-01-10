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
 * Preliminary: Check if csort plugin works
 ******************************************************************************/
print_section "Preliminary: Check if csort plugin works"

sysuse auto, clear
capture noisily csort price
local csort_works = (_rc == 0)

if !`csort_works' {
    di as error "  csort plugin is not functioning (rc=`=_rc')"
    di as error "  Skipping all csort tests"
    global TESTS_PASSED = 0
    global TESTS_FAILED = 1
    global TESTS_TOTAL = 1
    print_summary "csort"
    exit 1
}

di as result "  csort plugin is functional"
global TESTS_PASSED = $TESTS_PASSED + 1
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * Auto dataset tests
 ******************************************************************************/
print_section "Auto dataset"

sysuse auto, clear

* Single numeric variable
benchmark_sort price, testname("single numeric (price)")

sysuse auto, clear
benchmark_sort weight, testname("single numeric (weight)")

sysuse auto, clear
benchmark_sort mpg, testname("single numeric (mpg)")

* Single string variable
sysuse auto, clear
benchmark_sort make, testname("single string (make)")

* Categorical variable with duplicates
sysuse auto, clear
benchmark_sort foreign, testname("binary with duplicates (foreign)")

* Variable with missing values
sysuse auto, clear
benchmark_sort rep78, testname("with missing values (rep78)")

* Multi-variable sorts
sysuse auto, clear
benchmark_sort foreign price, testname("two vars (foreign price)")

sysuse auto, clear
benchmark_sort foreign rep78 price, testname("three vars (foreign rep78 price)")

sysuse auto, clear
benchmark_sort make price, testname("string + numeric (make price)")

/*******************************************************************************
 * Census dataset tests
 ******************************************************************************/
print_section "Census dataset"

sysuse census, clear
benchmark_sort state, testname("state name")

sysuse census, clear
benchmark_sort region, testname("region code")

sysuse census, clear
benchmark_sort pop, testname("population")

sysuse census, clear
benchmark_sort region state, testname("region state")

sysuse census, clear
benchmark_sort region pop, testname("region pop")

/*******************************************************************************
 * Synthetic data tests
 ******************************************************************************/
print_section "Synthetic data"

* Large dataset
clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform()
gen str20 label = "item" + string(runiformint(1, 500))

benchmark_sort value, testname("10K random floats")

clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform()
gen str20 label = "item" + string(runiformint(1, 500))

benchmark_sort group, testname("10K random ints")

clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform()
gen str20 label = "item" + string(runiformint(1, 500))

benchmark_sort label, testname("10K random strings")

clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform()
gen str20 label = "item" + string(runiformint(1, 500))

benchmark_sort group value, testname("10K group + value")

* Negative numbers
clear
set obs 1000
gen x = runiform() * 200 - 100

benchmark_sort x, testname("negative numbers")

* All same values (degenerate case)
clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same values")

* Already sorted
clear
set obs 1000
gen x = _n

benchmark_sort x, testname("already sorted")

* Reverse sorted
clear
set obs 1000
gen x = 1000 - _n

benchmark_sort x, testname("reverse sorted")

* Integer edge cases
clear
set obs 100
gen long x = runiformint(-2147483647, 2147483646)

benchmark_sort x, testname("extreme integers")

* Float precision
clear
set obs 100
gen double x = runiform() * 1e-10

benchmark_sort x, testname("small floats")

/*******************************************************************************
 * Edge cases
 ******************************************************************************/
print_section "Edge cases"

* Empty string handling
clear
set obs 100
gen str10 x = ""
replace x = "test" if _n > 50

benchmark_sort x, testname("empty strings")

* Very long strings
clear
set obs 100
gen str200 x = "a" * runiformint(1, 200)

benchmark_sort x, testname("variable length strings")

* Single observation
clear
set obs 1
gen x = 42

benchmark_sort x, testname("single observation")

* Two observations
clear
set obs 2
gen x = 2 - _n

benchmark_sort x, testname("two observations")

/*******************************************************************************
 * Larger scale tests
 ******************************************************************************/
print_section "Larger scale tests"

clear
set seed 54321
set obs 50000
gen id = _n
gen group = runiformint(1, 500)
gen value = rnormal()

benchmark_sort value, testname("50K random normal")

clear
set seed 54321
set obs 50000
gen id = _n
gen group = runiformint(1, 500)
gen value = rnormal()

benchmark_sort group id, testname("50K group + id")

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "csort"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
