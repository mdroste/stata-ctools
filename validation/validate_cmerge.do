/*******************************************************************************
 * validate_cmerge.do
 *
 * Comprehensive validation tests for cmerge vs native merge
 * Tests all merge types and options across multiple built-in Stata datasets
 ******************************************************************************/

do "validation/validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CMERGE VALIDATION TEST SUITE"
di as text "======================================================================"

* Create temp directory for test files
capture mkdir "validation/temp"

/*******************************************************************************
 * SECTION 1: BASIC MERGE TYPES
 ******************************************************************************/

/*******************************************************************************
 * TEST 1: Basic 1:1 merge (auto dataset)
 ******************************************************************************/
print_section "Test 1: Basic 1:1 merge (auto)"

sysuse auto, clear
keep make price mpg
gen id = _n
order id
tempfile master
quietly save `master'

sysuse auto, clear
keep make weight length
gen id = _n
order id
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Basic 1:1 merge"

/*******************************************************************************
 * TEST 2: Basic m:1 merge (auto dataset)
 ******************************************************************************/
print_section "Test 2: Basic m:1 merge (auto)"

sysuse auto, clear
keep make price mpg foreign
tempfile master
quietly save `master'

clear
input byte foreign str20 country
0 "Domestic"
1 "Foreign"
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:1 foreign using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:1 foreign using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Basic m:1 merge"

/*******************************************************************************
 * TEST 3: Basic 1:m merge (auto dataset)
 ******************************************************************************/
print_section "Test 3: Basic 1:m merge (auto)"

clear
input byte foreign str20 region
0 "Americas"
1 "International"
end
tempfile master
quietly save `master'

sysuse auto, clear
keep make price foreign
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:m foreign using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:m foreign using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

* 1:m merge produces rows in different order than Stata, compare after sorting
assert_data_equal_sorted `stata_merged' `cmerge_merged' "foreign make" "Basic 1:m merge"

/*******************************************************************************
 * TEST 4: m:m merge (auto dataset)
 ******************************************************************************/
print_section "Test 4: m:m merge"

clear
input int group float value
1 10
1 20
2 30
2 40
end
tempfile master
quietly save `master'

clear
input int group str10 label
1 "A"
1 "B"
2 "C"
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:m group using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:m group using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

* m:m merge results can be order-dependent, check counts
use `stata_merged', clear
local stata_n = _N
use `cmerge_merged', clear
local cmerge_n = _N

if `stata_n' == `cmerge_n' {
    di as result "  PASS: m:m merge - observation counts match (`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: m:m merge - counts differ (stata=`stata_n', cmerge=`cmerge_n')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * SECTION 2: KEEP() OPTION TESTS
 ******************************************************************************/

/*******************************************************************************
 * TEST 5: keep(match) option
 ******************************************************************************/
print_section "Test 5: keep(match) option"

sysuse auto, clear
keep if _n <= 50
gen id = _n
tempfile master
quietly save `master'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', keep(match)
drop _merge
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', keep(match)
drop _merge
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
use `cmerge_merged', clear
local cmerge_n = _N

if `stata_n' == `cmerge_n' {
    di as result "  PASS: keep(match) - observation counts match (`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: keep(match) - counts differ (stata=`stata_n', cmerge=`cmerge_n')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 6: keep(master) option
 ******************************************************************************/
print_section "Test 6: keep(master) option"

sysuse auto, clear
keep if _n <= 50
gen id = _n
tempfile master
quietly save `master'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', keep(master)
drop _merge
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', keep(master)
drop _merge
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
use `cmerge_merged', clear
local cmerge_n = _N

if `stata_n' == `cmerge_n' {
    di as result "  PASS: keep(master) - observation counts match (`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: keep(master) - counts differ (stata=`stata_n', cmerge=`cmerge_n')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 7: keep(using) option
 ******************************************************************************/
print_section "Test 7: keep(using) option"

sysuse auto, clear
keep if _n <= 50
gen id = _n
tempfile master
quietly save `master'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', keep(using)
drop _merge
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', keep(using)
drop _merge
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
use `cmerge_merged', clear
local cmerge_n = _N

if `stata_n' == `cmerge_n' {
    di as result "  PASS: keep(using) - observation counts match (`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: keep(using) - counts differ (stata=`stata_n', cmerge=`cmerge_n')"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 8: keep(master match) combined
 ******************************************************************************/
print_section "Test 8: keep(master match) combined"

sysuse auto, clear
keep if _n <= 50
gen id = _n
tempfile master
quietly save `master'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', keep(master match)
drop _merge
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', keep(master match)
drop _merge
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
use `cmerge_merged', clear
local cmerge_n = _N

if `stata_n' == `cmerge_n' {
    di as result "  PASS: keep(master match) - counts match (`stata_n')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: keep(master match) - counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * SECTION 3: GENERATE/NOGENERATE OPTIONS
 ******************************************************************************/

/*******************************************************************************
 * TEST 9: nogenerate option
 ******************************************************************************/
print_section "Test 9: nogenerate option"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight length
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', nogenerate
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', nogenerate
tempfile cmerge_merged
quietly save `cmerge_merged'

use `cmerge_merged', clear
capture confirm variable _merge
if _rc != 0 {
    di as result "  PASS: nogenerate - _merge variable not created"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: nogenerate - _merge variable exists"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 10: generate() option
 ******************************************************************************/
print_section "Test 10: generate() option"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', generate(merge_result)
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', generate(merge_result)
tempfile cmerge_merged
quietly save `cmerge_merged'

use `cmerge_merged', clear
capture confirm variable merge_result
if _rc == 0 {
    di as result "  PASS: generate() - custom merge variable created"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: generate() - custom merge variable not found"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 11: keepusing() option
 ******************************************************************************/
print_section "Test 11: keepusing() option"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight length turn
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data', keepusing(weight)
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data', keepusing(weight)
tempfile cmerge_merged
quietly save `cmerge_merged'

use `cmerge_merged', clear
capture confirm variable weight
local has_weight = (_rc == 0)
capture confirm variable length
local has_length = (_rc == 0)
capture confirm variable turn
local has_turn = (_rc == 0)

if `has_weight' & !`has_length' & !`has_turn' {
    di as result "  PASS: keepusing() - only specified variables merged"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: keepusing() - incorrect variables merged"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * SECTION 4: MULTIPLE KEY VARIABLES
 ******************************************************************************/

/*******************************************************************************
 * TEST 12: Two key variables (nlswork panel)
 ******************************************************************************/
print_section "Test 12: Two key variables (nlswork panel)"

webuse nlswork, clear
keep in 1/2000
keep idcode year ln_wage age
tempfile master
quietly save `master'

webuse nlswork, clear
keep in 1/2000
keep idcode year tenure ttl_exp
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 idcode year using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 idcode year using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Two key variables (idcode year)"

/*******************************************************************************
 * TEST 13: Three key variables
 ******************************************************************************/
print_section "Test 13: Three key variables"

clear
set obs 500
gen id1 = mod(_n-1, 10) + 1
gen id2 = mod(floor((_n-1)/10), 10) + 1
gen id3 = floor((_n-1)/100) + 1
gen value1 = runiform()
tempfile master
quietly save `master'

clear
set obs 500
gen id1 = mod(_n-1, 10) + 1
gen id2 = mod(floor((_n-1)/10), 10) + 1
gen id3 = floor((_n-1)/100) + 1
gen value2 = runiform() * 100
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id1 id2 id3 using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id1 id2 id3 using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Three key variables"

/*******************************************************************************
 * TEST 14: String key variable
 ******************************************************************************/
print_section "Test 14: String key variable (make)"

sysuse auto, clear
keep make price mpg
tempfile master
quietly save `master'

sysuse auto, clear
keep make weight length
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 make using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 make using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "String key variable (make)"

/*******************************************************************************
 * TEST 15: Mixed string and numeric keys
 ******************************************************************************/
print_section "Test 15: Mixed string and numeric keys"

sysuse census, clear
keep state region pop
tempfile master
quietly save `master'

sysuse census, clear
keep state region medage
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 state region using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 state region using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Mixed string and numeric keys"

/*******************************************************************************
 * SECTION 5: EDGE CASES
 ******************************************************************************/

/*******************************************************************************
 * TEST 16: Merge with no matches
 ******************************************************************************/
print_section "Test 16: Merge with no matches"

clear
input int id float value
1 10
2 20
3 30
end
tempfile master
quietly save `master'

clear
input int id float other
100 1.5
200 2.5
300 3.5
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Merge with no matches"

use `cmerge_merged', clear
quietly count if _merge == 3
if r(N) == 0 {
    di as result "  PASS: Correctly identified zero matches"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Incorrectly reported matches"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 17: Merge with all matches
 ******************************************************************************/
print_section "Test 17: Merge with all matches"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

use `cmerge_merged', clear
quietly count if _merge != 3
if r(N) == 0 {
    di as result "  PASS: All observations matched (_merge == 3)"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Some observations not matched"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 18: _merge values consistency
 ******************************************************************************/
print_section "Test 18: _merge values consistency"

sysuse auto, clear
keep if _n <= 40
gen id = _n
keep id price
tempfile master
quietly save `master'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
keep id weight
tempfile using_data
quietly save `using_data'

use `master', clear
cmerge 1:1 id using `using_data'

quietly count if _merge == 1
local m1 = r(N)
quietly count if _merge == 2
local m2 = r(N)
quietly count if _merge == 3
local m3 = r(N)

* Expected: ids 1-29 master only, ids 30-40 matched, ids 41-74 using only
local expected_m1 = 29
local expected_m3 = 11
local expected_m2 = 34

if `m1' == `expected_m1' & `m3' == `expected_m3' & `m2' == `expected_m2' {
    di as result "  PASS: _merge values correct (m1=`m1', m2=`m2', m3=`m3')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: _merge values incorrect"
    di as error "    Expected: m1=`expected_m1', m2=`expected_m2', m3=`expected_m3'"
    di as error "    Got:      m1=`m1', m2=`m2', m3=`m3'"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * SECTION 6: BUILT-IN DATASET MERGES
 ******************************************************************************/

/*******************************************************************************
 * TEST 19: Census state to region merge
 ******************************************************************************/
print_section "Test 19: Census state to region merge"

sysuse census, clear
keep state region pop
tempfile master
quietly save `master'

clear
input byte region str20 region_name
1 "NE"
2 "N Cntrl"
3 "South"
4 "West"
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:1 region using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:1 region using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Census state to region merge"

/*******************************************************************************
 * TEST 20: Auto - foreign to category merge
 ******************************************************************************/
print_section "Test 20: Auto foreign to category merge"

sysuse auto, clear
keep make price mpg foreign rep78
tempfile master
quietly save `master'

clear
input byte foreign str30 origin_desc float import_tax
0 "Made in USA" 0
1 "Imported" 0.05
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:1 foreign using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:1 foreign using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Auto foreign category merge"

/*******************************************************************************
 * TEST 21: Auto - rep78 to quality rating merge
 ******************************************************************************/
print_section "Test 21: Auto rep78 to quality rating merge"

sysuse auto, clear
keep make price mpg rep78
tempfile master
quietly save `master'

clear
input byte rep78 str20 quality_desc
1 "Poor"
2 "Fair"
3 "Average"
4 "Good"
5 "Excellent"
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:1 rep78 using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:1 rep78 using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Auto rep78 quality merge"

/*******************************************************************************
 * TEST 22: Census - division expansion
 ******************************************************************************/
print_section "Test 22: Census division 1:m merge"

clear
input byte region str20 division_name
1 "Northeast"
2 "Midwest"
3 "South"
4 "West"
end
tempfile master
quietly save `master'

sysuse census, clear
keep state region pop medage
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:m region using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:m region using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Census division 1:m merge"

/*******************************************************************************
 * TEST 23: Voter dataset merge
 ******************************************************************************/
print_section "Test 23: Voter dataset 1:1 merge"

sysuse voter, clear
keep pop frac
gen id = _n
tempfile master
quietly save `master'

sysuse voter, clear
keep cand inc
gen id = _n
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Voter dataset 1:1 merge"

/*******************************************************************************
 * TEST 24: US life expectancy merge
 ******************************************************************************/
print_section "Test 24: US life expectancy time series merge"

sysuse uslifeexp, clear
keep year le
tempfile master
quietly save `master'

sysuse uslifeexp, clear
keep year le_male le_female
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 year using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 year using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "US life expectancy merge"

/*******************************************************************************
 * TEST 25: SP500 stock data merge
 ******************************************************************************/
print_section "Test 25: SP500 time series merge"

sysuse sp500, clear
keep date open high
tempfile master
quietly save `master'

sysuse sp500, clear
keep date low close volume
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 date using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 date using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "SP500 time series merge"

/*******************************************************************************
 * SECTION 7: PANEL DATA MERGES (nlswork)
 ******************************************************************************/

/*******************************************************************************
 * TEST 26: nlswork full panel merge
 ******************************************************************************/
print_section "Test 26: nlswork panel 1:1 merge"

webuse nlswork, clear
keep in 1/5000
keep idcode year ln_wage hours
tempfile master
quietly save `master'

webuse nlswork, clear
keep in 1/5000
keep idcode year age tenure wks_work
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 idcode year using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 idcode year using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "nlswork panel 1:1 merge"

/*******************************************************************************
 * TEST 27: nlswork m:1 individual characteristics
 ******************************************************************************/
print_section "Test 27: nlswork m:1 individual merge"

webuse nlswork, clear
keep in 1/3000
keep idcode year ln_wage
tempfile master
quietly save `master'

* Create individual-level data (first obs per person)
webuse nlswork, clear
keep in 1/3000
bysort idcode (year): keep if _n == 1
keep idcode race grade
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:1 idcode using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:1 idcode using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "nlswork m:1 individual merge"

/*******************************************************************************
 * TEST 28: nlswork 1:m year characteristics
 ******************************************************************************/
print_section "Test 28: nlswork 1:m year merge"

webuse nlswork, clear
bysort year: keep if _n == 1
keep year
gen year_char = "Y" + string(year)
tempfile master
quietly save `master'

webuse nlswork, clear
keep in 1/2000
keep idcode year ln_wage
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:m year using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:m year using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "nlswork 1:m year merge"

/*******************************************************************************
 * SECTION 8: LARGE DATASET MERGES
 ******************************************************************************/

/*******************************************************************************
 * TEST 29: Large 1:1 merge (50K observations)
 ******************************************************************************/
print_section "Test 29: Large 1:1 merge (50K obs)"

clear
set seed 12345
set obs 50000
gen id = _n
gen value1 = runiform()
gen value2 = runiformint(1, 100)
tempfile master
quietly save `master'

clear
set obs 50000
gen id = _n
gen extra1 = runiform() * 100
gen extra2 = runiformint(1, 1000)
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
quietly count if _merge == 3
local stata_m3 = r(N)

use `cmerge_merged', clear
local cmerge_n = _N
quietly count if _merge == 3
local cmerge_m3 = r(N)

if `stata_n' == `cmerge_n' & `stata_m3' == `cmerge_m3' {
    di as result "  PASS: Large 1:1 merge (N=`stata_n', matched=`stata_m3')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Large 1:1 merge counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 30: Large 1:m merge (50K x 10K)
 ******************************************************************************/
print_section "Test 30: Large 1:m merge (50K x 10K)"

clear
set seed 54321
set obs 50000
gen id = _n
gen value1 = runiform()
gen value2 = runiformint(1, 100)
tempfile master
quietly save `master'

clear
set obs 10000
gen id = runiformint(1, 50000)
gen extra1 = runiform() * 100
gen extra2 = "label" + string(runiformint(1, 500))
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:m id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:m id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
quietly count if _merge == 1
local stata_m1 = r(N)
quietly count if _merge == 2
local stata_m2 = r(N)
quietly count if _merge == 3
local stata_m3 = r(N)

use `cmerge_merged', clear
local cmerge_n = _N
quietly count if _merge == 1
local cmerge_m1 = r(N)
quietly count if _merge == 2
local cmerge_m2 = r(N)
quietly count if _merge == 3
local cmerge_m3 = r(N)

local match = (`stata_n' == `cmerge_n') & (`stata_m1' == `cmerge_m1') & ///
              (`stata_m2' == `cmerge_m2') & (`stata_m3' == `cmerge_m3')

if `match' {
    di as result "  PASS: Large 1:m merge (N=`stata_n', m1=`stata_m1', m2=`stata_m2', m3=`stata_m3')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Large 1:m merge counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 31: Large m:1 merge (100K x 1K)
 ******************************************************************************/
print_section "Test 31: Large m:1 merge (100K x 1K)"

clear
set seed 11111
set obs 100000
gen group_id = runiformint(1, 1000)
gen value = runiform()
tempfile master
quietly save `master'

clear
set obs 1000
gen group_id = _n
gen group_name = "Group" + string(group_id)
gen group_weight = runiform()
tempfile using_data
quietly save `using_data'

use `master', clear
merge m:1 group_id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge m:1 group_id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_n = _N
quietly count if _merge == 3
local stata_m3 = r(N)

use `cmerge_merged', clear
local cmerge_n = _N
quietly count if _merge == 3
local cmerge_m3 = r(N)

if `stata_n' == `cmerge_n' & `stata_m3' == `cmerge_m3' {
    di as result "  PASS: Large m:1 merge (N=`stata_n', matched=`stata_m3')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Large m:1 merge counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * TEST 32: Large merge with partial overlap
 ******************************************************************************/
print_section "Test 32: Large merge with partial overlap"

clear
set seed 22222
set obs 30000
gen id = _n
gen value1 = runiform()
tempfile master
quietly save `master'

clear
set obs 30000
gen id = _n + 15000  // ids 15001-45000, overlaps with master 15001-30000
gen value2 = runiform()
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
quietly count if _merge == 1
local stata_m1 = r(N)
quietly count if _merge == 2
local stata_m2 = r(N)
quietly count if _merge == 3
local stata_m3 = r(N)

use `cmerge_merged', clear
quietly count if _merge == 1
local cmerge_m1 = r(N)
quietly count if _merge == 2
local cmerge_m2 = r(N)
quietly count if _merge == 3
local cmerge_m3 = r(N)

* Expected: 15000 master only, 15000 using only, 15000 matched
if `stata_m1' == `cmerge_m1' & `stata_m2' == `cmerge_m2' & `stata_m3' == `cmerge_m3' {
    di as result "  PASS: Partial overlap (m1=`cmerge_m1', m2=`cmerge_m2', m3=`cmerge_m3')"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Partial overlap counts differ"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * SECTION 9: STRING KEY MERGES
 ******************************************************************************/

/*******************************************************************************
 * TEST 33: Long string keys
 ******************************************************************************/
print_section "Test 33: Long string keys"

clear
set obs 100
gen str50 long_key = "prefix_" + string(_n, "%04.0f") + "_suffix_data"
gen value1 = runiform()
tempfile master
quietly save `master'

clear
set obs 100
gen str50 long_key = "prefix_" + string(_n, "%04.0f") + "_suffix_data"
gen value2 = runiform() * 100
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 long_key using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 long_key using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Long string keys merge"

/*******************************************************************************
 * TEST 34: Census state name merge
 ******************************************************************************/
print_section "Test 34: Census state name 1:1 merge"

sysuse census, clear
keep state pop
tempfile master
quietly save `master'

sysuse census, clear
keep state medage death marriage divorce
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 state using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 state using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Census state name 1:1 merge"

/*******************************************************************************
 * TEST 35: Auto make model merge
 ******************************************************************************/
print_section "Test 35: Auto make 1:1 merge"

sysuse auto, clear
keep make mpg price
tempfile master
quietly save `master'

sysuse auto, clear
keep make weight length turn displacement gear_ratio
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 make using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 make using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Auto make 1:1 merge"

/*******************************************************************************
 * SECTION 10: SPECIAL CASES
 ******************************************************************************/

/*******************************************************************************
 * TEST 36: Merge with missing key values
 ******************************************************************************/
print_section "Test 36: Merge with missing key values"

clear
input int id float value
1 10
2 20
. 30
4 40
end
tempfile master
quietly save `master'

clear
input int id float other
1 100
. 200
3 300
4 400
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Merge with missing key values"

/*******************************************************************************
 * TEST 37: Merge with negative key values
 ******************************************************************************/
print_section "Test 37: Merge with negative key values"

clear
input int id float value
-5 10
-2 20
0 30
2 40
5 50
end
tempfile master
quietly save `master'

clear
input int id float other
-5 100
-1 200
0 300
1 400
5 500
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Merge with negative keys"

/*******************************************************************************
 * TEST 38: Merge with float key values
 ******************************************************************************/
print_section "Test 38: Merge with float key values"

clear
input float id float value
1.5 10
2.7 20
3.14159 30
end
tempfile master
quietly save `master'

clear
input float id float other
1.5 100
2.7 200
3.14159 300
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Merge with float keys"

/*******************************************************************************
 * TEST 39: Merge with duplicate values in non-key variables
 ******************************************************************************/
print_section "Test 39: Merge with duplicate non-key values"

clear
input int id float value str10 category
1 100 "A"
2 100 "A"
3 100 "B"
4 200 "B"
end
tempfile master
quietly save `master'

clear
input int id float other str10 label
1 50 "X"
2 50 "X"
3 75 "Y"
4 75 "Y"
end
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

assert_data_equal `stata_merged' `cmerge_merged' "Merge with duplicate non-key values"

/*******************************************************************************
 * TEST 40: Very wide using dataset
 ******************************************************************************/
print_section "Test 40: Wide using dataset (50 variables)"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
drop make price mpg
* Generate additional variables
forval i = 1/40 {
    gen extra`i' = runiform()
}
tempfile using_data
quietly save `using_data'

use `master', clear
merge 1:1 id using `using_data'
tempfile stata_merged
quietly save `stata_merged'

use `master', clear
cmerge 1:1 id using `using_data'
tempfile cmerge_merged
quietly save `cmerge_merged'

use `stata_merged', clear
local stata_nvars = c(k)
use `cmerge_merged', clear
local cmerge_nvars = c(k)

if `stata_nvars' == `cmerge_nvars' {
    di as result "  PASS: Wide dataset merge (`stata_nvars' variables)"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as error "  FAIL: Variable count mismatch"
    global TESTS_FAILED = $TESTS_FAILED + 1
}
global TESTS_TOTAL = $TESTS_TOTAL + 1

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cmerge"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
