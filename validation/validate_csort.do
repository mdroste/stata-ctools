/*******************************************************************************
 * validate_csort.do
 *
 * Comprehensive validation tests for csort vs native Stata sort
 * Tests all sort algorithms and data type scenarios
 *
 * Algorithms tested:
 *   - lsd (default): LSD radix sort
 *   - msd: MSD radix sort
 *   - timsort: Timsort
 *   - sample: Sample sort
 *   - counting: Counting sort
 *   - merge: Parallel merge sort
 *   - ips4o: In-place parallel super scalar samplesort
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text ""
noi di as text "======================================================================"
noi di as text "              CSORT VALIDATION TEST SUITE"
noi di as text "======================================================================"

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture noisily csort price
if _rc != 0 {
    noi test_fail "csort plugin load" "plugin returned error `=_rc'"
    noi print_summary "csort"
    exit 1
}
noi test_pass "csort plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Default algorithm (LSD) - comprehensive tests
 ******************************************************************************/
noi print_section "Default Algorithm (LSD Radix Sort)"

* Single numeric variable
sysuse auto, clear
noi benchmark_sort price, testname("numeric: price")

sysuse auto, clear
noi benchmark_sort mpg, testname("numeric: mpg")

sysuse auto, clear
noi benchmark_sort weight, testname("numeric: weight")

* Single string variable
sysuse auto, clear
noi benchmark_sort make, testname("string: make")

* Binary/categorical variable
sysuse auto, clear
noi benchmark_sort foreign, testname("binary: foreign")

* Variable with missing values
sysuse auto, clear
noi benchmark_sort rep78, testname("with missings: rep78")

* Multiple sort keys
sysuse auto, clear
noi benchmark_sort foreign price, testname("two keys: foreign price")

sysuse auto, clear
noi benchmark_sort foreign rep78 price, testname("three keys: foreign rep78 price")

sysuse auto, clear
noi benchmark_sort make price, testname("mixed keys: make price")

/*******************************************************************************
 * SECTION 3: MSD Radix Sort algorithm
 ******************************************************************************/
noi print_section "MSD Radix Sort Algorithm"

sysuse auto, clear
noi benchmark_sort price, testname("numeric: price") algorithm(msd)

sysuse auto, clear
noi benchmark_sort make, testname("string: make") algorithm(msd)

sysuse census, clear
noi benchmark_sort state, testname("string: state") algorithm(msd)

sysuse auto, clear
noi benchmark_sort foreign make, testname("mixed: foreign make") algorithm(msd)

sysuse auto, clear
noi benchmark_sort rep78, testname("with missings") algorithm(msd)

/*******************************************************************************
 * SECTION 4: Timsort algorithm
 ******************************************************************************/
noi print_section "Timsort Algorithm"

sysuse auto, clear
noi benchmark_sort price, testname("numeric: price") algorithm(timsort)

sysuse auto, clear
noi benchmark_sort mpg weight, testname("multi-numeric") algorithm(timsort)

sysuse auto, clear
noi benchmark_sort make, testname("string: make") algorithm(timsort)

sysuse auto, clear
noi benchmark_sort rep78, testname("with missings") algorithm(timsort)

* Already sorted (timsort excels here)
clear
set obs 5000
gen x = _n
noi benchmark_sort x, testname("already sorted") algorithm(timsort)

* Reverse sorted
clear
set obs 5000
gen x = 5000 - _n
noi benchmark_sort x, testname("reverse sorted") algorithm(timsort)

/*******************************************************************************
 * SECTION 5: Sample Sort algorithm
 ******************************************************************************/
noi print_section "Sample Sort Algorithm"

sysuse auto, clear
noi benchmark_sort price, testname("numeric: price") algorithm(sample)

sysuse auto, clear
noi benchmark_sort make, testname("string: make") algorithm(sample)

sysuse census, clear
noi benchmark_sort state, testname("string: state") algorithm(sample)

sysuse auto, clear
noi benchmark_sort foreign make, testname("mixed: foreign make") algorithm(sample)

sysuse auto, clear
noi benchmark_sort rep78, testname("with missings") algorithm(sample)

/*******************************************************************************
 * SECTION 6: Counting Sort algorithm
 ******************************************************************************/
noi print_section "Counting Sort Algorithm"

sysuse auto, clear
noi benchmark_sort foreign, testname("binary (0/1)") algorithm(counting)

sysuse auto, clear
noi benchmark_sort rep78, testname("small int range") algorithm(counting)

clear
set obs 5000
set seed 12345
gen year = runiformint(1990, 2023)
noi benchmark_sort year, testname("year data") algorithm(counting)

sysuse auto, clear
noi benchmark_sort foreign rep78, testname("multi-int") algorithm(counting)

/*******************************************************************************
 * SECTION 7: Parallel Merge Sort algorithm
 ******************************************************************************/
noi print_section "Parallel Merge Sort Algorithm"

sysuse auto, clear
noi benchmark_sort price, testname("numeric: price") algorithm(merge)

sysuse auto, clear
noi benchmark_sort make, testname("string: make") algorithm(merge)

sysuse census, clear
noi benchmark_sort state, testname("string: state") algorithm(merge)

sysuse auto, clear
noi benchmark_sort foreign make, testname("mixed: foreign make") algorithm(merge)

sysuse auto, clear
noi benchmark_sort rep78, testname("with missings") algorithm(merge)

/*******************************************************************************
 * SECTION 8: IPS4o algorithm
 ******************************************************************************/
noi print_section "IPS4o Algorithm"

sysuse auto, clear
noi benchmark_sort price, testname("numeric: price") algorithm(ips4o)

sysuse auto, clear
noi benchmark_sort make, testname("string: make") algorithm(ips4o)

sysuse census, clear
noi benchmark_sort state, testname("string: state") algorithm(ips4o)

sysuse auto, clear
noi benchmark_sort foreign make, testname("mixed: foreign make") algorithm(ips4o)

sysuse auto, clear
noi benchmark_sort rep78, testname("with missings") algorithm(ips4o)

* Already sorted
clear
set obs 5000
gen x = _n
noi benchmark_sort x, testname("already sorted") algorithm(ips4o)

* Reverse sorted
clear
set obs 5000
gen x = 5000 - _n
noi benchmark_sort x, testname("reverse sorted") algorithm(ips4o)

/*******************************************************************************
 * SECTION 9: Edge cases - all algorithms
 ******************************************************************************/
noi print_section "Edge Cases"

* Single observation
foreach alg in lsd msd timsort sample merge ips4o {
    clear
    set obs 1
    gen x = 42
    noi benchmark_sort x, testname("single obs") algorithm(`alg')
}

* Two observations
foreach alg in lsd msd timsort sample merge ips4o {
    clear
    set obs 2
    gen x = 2 - _n
    noi benchmark_sort x, testname("two obs") algorithm(`alg')
}

* All same values
foreach alg in lsd msd timsort sample merge ips4o {
    clear
    set obs 1000
    gen x = 5
    noi benchmark_sort x, testname("all same values") algorithm(`alg')
}

/*******************************************************************************
 * SECTION 10: Numeric edge cases
 ******************************************************************************/
noi print_section "Numeric Edge Cases"

* Negative numbers
clear
set obs 1000
gen x = runiform() * 200 - 100
noi benchmark_sort x, testname("negative numbers")

* Extreme integers
clear
set obs 100
gen long x = runiformint(-2147483647, 2147483646)
noi benchmark_sort x, testname("extreme integers")

* Small floats
clear
set obs 100
gen double x = runiform() * 1e-10
noi benchmark_sort x, testname("small floats (1e-10)")

* Large floats
clear
set obs 100
gen double x = runiform() * 1e10
noi benchmark_sort x, testname("large floats (1e10)")

* Mixed positive/negative
clear
set obs 500
gen x = cond(_n <= 250, -_n, _n - 250)
noi benchmark_sort x, testname("mixed pos/neg")

/*******************************************************************************
 * SECTION 11: String edge cases
 ******************************************************************************/
noi print_section "String Edge Cases"

* Empty strings
clear
set obs 100
gen str10 x = ""
replace x = "test" if _n > 50
noi benchmark_sort x, testname("empty strings")

* Variable length strings
clear
set obs 100
gen str200 x = substr("a" * 200, 1, runiformint(1, 200))
noi benchmark_sort x, testname("variable length strings")

* Unicode/special characters (if supported)
clear
set obs 50
gen str20 x = "item" + string(_n)
replace x = "Item_" + string(_n) if _n > 25
noi benchmark_sort x, testname("mixed case strings")

/*******************************************************************************
 * SECTION 12: Large dataset tests
 ******************************************************************************/
noi print_section "Large Dataset Tests"

* 10K observations
clear
set seed 12345
set obs 10000
gen id = _n
gen group = runiformint(1, 100)
gen value = runiform()
gen str20 label = "item" + string(runiformint(1, 500))

noi benchmark_sort value, testname("10K random floats")

clear
set seed 12345
set obs 10000
gen group = runiformint(1, 100)
noi benchmark_sort group, testname("10K random ints")

clear
set seed 12345
set obs 10000
gen str20 label = "item" + string(runiformint(1, 500))
noi benchmark_sort label, testname("10K random strings")

* 50K observations
clear
set seed 54321
set obs 50000
gen value = rnormal()
noi benchmark_sort value, testname("50K random normal")

clear
set seed 54321
set obs 50000
gen id = _n
gen group = runiformint(1, 500)
noi benchmark_sort group id, testname("50K group + id")

/*******************************************************************************
 * SECTION 13: Census dataset tests
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
noi benchmark_sort state, testname("state name")

sysuse census, clear
noi benchmark_sort region, testname("region code")

sysuse census, clear
noi benchmark_sort pop, testname("population")

sysuse census, clear
noi benchmark_sort region state, testname("region state")

sysuse census, clear
noi benchmark_sort region pop, testname("region pop")

/*******************************************************************************
 * SECTION 14: verbose option
 ******************************************************************************/
noi print_section "Option Tests"

* Test verbose option (just verify it doesn't error)
sysuse auto, clear
capture csort price, verbose
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "csort"

if $TESTS_FAILED > 0 {
    exit 1
}

}
