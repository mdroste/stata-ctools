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
 * Algorithm-specific tests: MSD Radix Sort
 ******************************************************************************/
print_section "MSD Radix Sort Algorithm"

* Numeric sorting
sysuse auto, clear
benchmark_sort price, testname("numeric") algorithm(msd)

sysuse auto, clear
benchmark_sort mpg weight, testname("multi-numeric") algorithm(msd)

* String sorting (MSD is optimized for strings)
sysuse auto, clear
benchmark_sort make, testname("string (make)") algorithm(msd)

sysuse census, clear
benchmark_sort state, testname("string (state)") algorithm(msd)

* Mixed types
sysuse auto, clear
benchmark_sort foreign make, testname("mixed types") algorithm(msd)

* With missing values
sysuse auto, clear
benchmark_sort rep78, testname("with missings") algorithm(msd)

* Larger dataset
clear
set seed 12345
set obs 10000
gen str20 label = "item" + string(runiformint(1, 500))
gen value = rnormal()

benchmark_sort label, testname("10K strings") algorithm(msd)

clear
set seed 12345
set obs 10000
gen str20 label = "item" + string(runiformint(1, 500))
gen value = rnormal()

benchmark_sort value, testname("10K numeric") algorithm(msd)

/*******************************************************************************
 * Algorithm-specific tests: Timsort
 ******************************************************************************/
print_section "Timsort Algorithm"

* Numeric sorting
sysuse auto, clear
benchmark_sort price, testname("numeric") algorithm(timsort)

sysuse auto, clear
benchmark_sort mpg weight, testname("multi-numeric") algorithm(timsort)

* String sorting
sysuse auto, clear
benchmark_sort make, testname("string (make)") algorithm(timsort)

sysuse census, clear
benchmark_sort state, testname("string (state)") algorithm(timsort)

* Mixed types
sysuse auto, clear
benchmark_sort foreign make, testname("mixed types") algorithm(timsort)

* With missing values
sysuse auto, clear
benchmark_sort rep78, testname("with missings") algorithm(timsort)

* Already sorted (Timsort should excel here)
clear
set obs 5000
gen x = _n

benchmark_sort x, testname("already sorted") algorithm(timsort)

* Reverse sorted
clear
set obs 5000
gen x = 5000 - _n

benchmark_sort x, testname("reverse sorted") algorithm(timsort)

* Partially sorted (panel-like data)
clear
set seed 99999
set obs 10000
gen panel_id = mod(_n - 1, 100) + 1
gen time = ceil(_n / 100)
gen value = rnormal()
sort panel_id  /* Already sorted by panel_id */

benchmark_sort panel_id time, testname("panel data") algorithm(timsort)

* Larger dataset
clear
set seed 12345
set obs 10000
gen value = rnormal()

benchmark_sort value, testname("10K numeric") algorithm(timsort)

/*******************************************************************************
 * Edge cases for all algorithms
 ******************************************************************************/
print_section "Edge cases - all algorithms"

* Single observation
clear
set obs 1
gen x = 42

benchmark_sort x, testname("single obs") algorithm(lsd)

clear
set obs 1
gen x = 42

benchmark_sort x, testname("single obs") algorithm(msd)

clear
set obs 1
gen x = 42

benchmark_sort x, testname("single obs") algorithm(timsort)

* All same values
clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same") algorithm(lsd)

clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same") algorithm(msd)

clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same") algorithm(timsort)

* Negative numbers
clear
set obs 1000
gen x = runiform() * 200 - 100

benchmark_sort x, testname("negatives") algorithm(msd)

clear
set obs 1000
gen x = runiform() * 200 - 100

benchmark_sort x, testname("negatives") algorithm(timsort)

/*******************************************************************************
 * Algorithm-specific tests: Sample Sort
 ******************************************************************************/
print_section "Sample Sort Algorithm"

* Numeric sorting
sysuse auto, clear
benchmark_sort price, testname("numeric") algorithm(sample)

sysuse auto, clear
benchmark_sort mpg weight, testname("multi-numeric") algorithm(sample)

* String sorting
sysuse auto, clear
benchmark_sort make, testname("string (make)") algorithm(sample)

sysuse census, clear
benchmark_sort state, testname("string (state)") algorithm(sample)

* Mixed types
sysuse auto, clear
benchmark_sort foreign make, testname("mixed types") algorithm(sample)

* With missing values
sysuse auto, clear
benchmark_sort rep78, testname("with missings") algorithm(sample)

* Larger dataset
clear
set seed 12345
set obs 10000
gen value = rnormal()

benchmark_sort value, testname("10K numeric") algorithm(sample)

/*******************************************************************************
 * Algorithm-specific tests: Counting Sort
 ******************************************************************************/
print_section "Counting Sort Algorithm"

* Integer data
sysuse auto, clear
benchmark_sort foreign, testname("binary (0/1)") algorithm(counting)

sysuse auto, clear
benchmark_sort rep78, testname("small int range") algorithm(counting)

* Year-like data
clear
set obs 5000
set seed 12345
gen year = runiformint(1990, 2023)
benchmark_sort year, testname("year data") algorithm(counting)

* State codes
clear
set obs 5000
set seed 54321
gen state = runiformint(1, 50)
benchmark_sort state, testname("state codes") algorithm(counting)

* Multi-key with integers
sysuse auto, clear
benchmark_sort foreign rep78, testname("multi-int") algorithm(counting)

/*******************************************************************************
 * Algorithm-specific tests: Parallel Merge Sort
 ******************************************************************************/
print_section "Parallel Merge Sort Algorithm"

* Numeric sorting
sysuse auto, clear
benchmark_sort price, testname("numeric") algorithm(merge)

sysuse auto, clear
benchmark_sort mpg weight, testname("multi-numeric") algorithm(merge)

* String sorting
sysuse auto, clear
benchmark_sort make, testname("string (make)") algorithm(merge)

sysuse census, clear
benchmark_sort state, testname("string (state)") algorithm(merge)

* Mixed types
sysuse auto, clear
benchmark_sort foreign make, testname("mixed types") algorithm(merge)

* With missing values
sysuse auto, clear
benchmark_sort rep78, testname("with missings") algorithm(merge)

* Larger dataset
clear
set seed 12345
set obs 10000
gen value = rnormal()

benchmark_sort value, testname("10K numeric") algorithm(merge)

/*******************************************************************************
 * Edge cases for new algorithms
 ******************************************************************************/
print_section "Edge cases - new algorithms"

* Single observation
clear
set obs 1
gen x = 42

benchmark_sort x, testname("single obs") algorithm(sample)

clear
set obs 1
gen x = 42

benchmark_sort x, testname("single obs") algorithm(merge)

* All same values
clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same") algorithm(sample)

clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same") algorithm(merge)

/*******************************************************************************
 * Algorithm-specific tests: IPS4o (In-place Parallel Super Scalar Samplesort)
 ******************************************************************************/
print_section "IPS4o Algorithm"

* Numeric sorting
sysuse auto, clear
benchmark_sort price, testname("numeric") algorithm(ips4o)

sysuse auto, clear
benchmark_sort mpg weight, testname("multi-numeric") algorithm(ips4o)

* String sorting
sysuse auto, clear
benchmark_sort make, testname("string (make)") algorithm(ips4o)

sysuse census, clear
benchmark_sort state, testname("string (state)") algorithm(ips4o)

* Mixed types
sysuse auto, clear
benchmark_sort foreign make, testname("mixed types") algorithm(ips4o)

* With missing values
sysuse auto, clear
benchmark_sort rep78, testname("with missings") algorithm(ips4o)

* Larger dataset
clear
set seed 12345
set obs 10000
gen value = rnormal()

benchmark_sort value, testname("10K numeric") algorithm(ips4o)

clear
set seed 12345
set obs 10000
gen str20 label = "item" + string(runiformint(1, 500))

benchmark_sort label, testname("10K strings") algorithm(ips4o)

* Edge cases for IPS4o
clear
set obs 1
gen x = 42

benchmark_sort x, testname("single obs") algorithm(ips4o)

clear
set obs 1000
gen x = 5

benchmark_sort x, testname("all same") algorithm(ips4o)

clear
set obs 1000
gen x = runiform() * 200 - 100

benchmark_sort x, testname("negatives") algorithm(ips4o)

* Already sorted (edge case)
clear
set obs 5000
gen x = _n

benchmark_sort x, testname("already sorted") algorithm(ips4o)

* Reverse sorted
clear
set obs 5000
gen x = 5000 - _n

benchmark_sort x, testname("reverse sorted") algorithm(ips4o)

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "csort"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
