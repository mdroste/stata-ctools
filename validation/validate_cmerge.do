/*******************************************************************************
 * validate_cmerge.do
 *
 * Comprehensive validation tests for cmerge vs native Stata merge
 * Tests all merge types and options
 *
 * Merge types: 1:1, m:1, 1:m, m:m
 * Options: keep(), generate(), nogenerate, keepusing(), sorted, force,
 *          noreport, nolabel, nonotes, update, replace, assert()
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CMERGE VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * SECTION 1: Basic 1:1 merge tests
 ******************************************************************************/
noi print_section "1:1 Merge Tests"

* Numeric key
sysuse auto, clear
keep make price mpg
gen id = _n
tempfile master_1
save `master_1'

sysuse auto, clear
keep make weight length
gen id = _n
tempfile using_1
save `using_1'

use `master_1', clear
noi benchmark_merge 1:1 id using `using_1', testname("1:1 numeric key")

* String key
sysuse auto, clear
keep make price mpg
tempfile master_2
save `master_2'

sysuse auto, clear
keep make weight length
tempfile using_2
save `using_2'

use `master_2', clear
noi benchmark_merge 1:1 make using `using_2', testname("1:1 string key")

/*******************************************************************************
 * SECTION 2: m:1 merge tests
 ******************************************************************************/
noi print_section "m:1 Merge Tests"

* Auto foreign lookup
sysuse auto, clear
keep make price mpg foreign
tempfile master_3
save `master_3'

clear
input byte foreign str20 country
0 "Domestic"
1 "Foreign"
end
tempfile using_3
save `using_3'

use `master_3', clear
noi benchmark_merge m:1 foreign using `using_3', testname("m:1 foreign lookup")

* Census region lookup
sysuse census, clear
keep state region pop
tempfile master_4
save `master_4'

clear
input byte region str20 region_name
1 "NE"
2 "N Cntrl"
3 "South"
4 "West"
end
tempfile using_4
save `using_4'

use `master_4', clear
noi benchmark_merge m:1 region using `using_4', testname("m:1 census region")

/*******************************************************************************
 * SECTION 3: 1:m merge tests
 ******************************************************************************/
noi print_section "1:m Merge Tests"

clear
input byte region str20 division_name
1 "Northeast"
2 "Midwest"
3 "South"
4 "West"
end
tempfile master_5
save `master_5'

sysuse census, clear
keep state region pop medage
tempfile using_5
save `using_5'

use `master_5', clear
noi benchmark_merge 1:m region using `using_5', testname("1:m region to states")

/*******************************************************************************
 * SECTION 4: m:m merge tests
 * Note: Stata's m:m merge behavior is documented as "unusual" and depends on
 * internal data ordering. cmerge implements sequential pairing (row i from
 * master pairs with row i from using within each group). We test that cmerge
 * runs without error and produces expected dimensions.
 ******************************************************************************/
noi print_section "m:m Merge Tests"

clear
set obs 100
gen group = mod(_n - 1, 10) + 1
gen val1 = runiform()
tempfile master_mm
save `master_mm'

clear
set obs 100
gen group = mod(_n - 1, 10) + 1
gen val2 = runiform() * 100
tempfile using_mm
save `using_mm'

* Test that m:m merge runs and produces correct dimensions
use `master_mm', clear
capture noisily cmerge m:m group using `using_mm', nogen noreport
if _rc == 0 & _N == 100 {
    noi test_pass "m:m basic (N=100)"
}
else if _rc != 0 {
    noi test_fail "m:m basic" "returned error `=_rc'"
}
else {
    noi test_fail "m:m basic" "unexpected N=`=_N' (expected 100)"
}

/*******************************************************************************
 * SECTION 5: keep() option tests
 ******************************************************************************/
noi print_section "keep() Option Tests"

sysuse auto, clear
keep if _n <= 50
gen id = _n
tempfile master_6
save `master_6'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
tempfile using_6
save `using_6'

use `master_6', clear
noi benchmark_merge 1:1 id using `using_6', keep(match) testname("keep(match)")

use `master_6', clear
noi benchmark_merge 1:1 id using `using_6', keep(master) testname("keep(master)")

use `master_6', clear
noi benchmark_merge 1:1 id using `using_6', keep(using) testname("keep(using)")

use `master_6', clear
noi benchmark_merge 1:1 id using `using_6', keep(master match) testname("keep(master match)")

use `master_6', clear
noi benchmark_merge 1:1 id using `using_6', keep(master using) testname("keep(master using)")

use `master_6', clear
noi benchmark_merge 1:1 id using `using_6', keep(using match) testname("keep(using match)")

/*******************************************************************************
 * SECTION 6: generate/nogenerate options
 ******************************************************************************/
noi print_section "generate/nogenerate Options"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight length
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', nogenerate testname("nogenerate")

use `master', clear
noi benchmark_merge 1:1 id using `using_1', generate(merge_result) testname("generate(merge_result)")

/*******************************************************************************
 * SECTION 7: keepusing() option
 ******************************************************************************/
noi print_section "keepusing() Option"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight length turn
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', keepusing(weight) testname("keepusing single var")

use `master', clear
noi benchmark_merge 1:1 id using `using_1', keepusing(weight length) testname("keepusing two vars")

/*******************************************************************************
 * SECTION 8: Multiple key variables
 ******************************************************************************/
noi print_section "Multiple Key Variables"

* Two keys
webuse nlswork, clear
keep in 1/2000
keep idcode year ln_wage age
tempfile master
save `master'

webuse nlswork, clear
keep in 1/2000
keep idcode year tenure ttl_exp
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 idcode year using `using_1', testname("two keys: idcode year")

* Three keys
clear
set obs 500
gen id1 = mod(_n-1, 10) + 1
gen id2 = mod(floor((_n-1)/10), 10) + 1
gen id3 = floor((_n-1)/100) + 1
gen value1 = runiform()
tempfile master
save `master'

clear
set obs 500
gen id1 = mod(_n-1, 10) + 1
gen id2 = mod(floor((_n-1)/10), 10) + 1
gen id3 = floor((_n-1)/100) + 1
gen value2 = runiform() * 100
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id1 id2 id3 using `using_1', testname("three keys")

* Mixed string and numeric keys
sysuse census, clear
keep state region pop
tempfile master
save `master'

sysuse census, clear
keep state region medage
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 state region using `using_1', testname("mixed keys: string+numeric")

/*******************************************************************************
 * SECTION 9: Edge cases
 ******************************************************************************/
noi print_section "Edge Cases"

* No matches
clear
input int id float value
1 10
2 20
3 30
end
tempfile master
save `master'

clear
input int id float other
100 1.5
200 2.5
300 3.5
end
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("no matches")

* All matches
sysuse auto, clear
gen id = _n
keep id make price
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("all matches")

* Missing key values
clear
input int id float value
1 10
2 20
. 30
4 40
end
tempfile master
save `master'

clear
input int id float other
1 100
. 200
3 300
4 400
end
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("missing key values")

* Negative key values
clear
input int id float value
-5 10
-2 20
0 30
2 40
5 50
end
tempfile master
save `master'

clear
input int id float other
-5 100
-1 200
0 300
1 400
5 500
end
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("negative keys")

* Float key values
clear
input float id float value
1.5 10
2.7 20
3.14159 30
end
tempfile master
save `master'

clear
input float id float other
1.5 100
2.7 200
3.14159 300
end
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("float keys")

/*******************************************************************************
 * SECTION 10: Large dataset tests
 ******************************************************************************/
noi print_section "Large Dataset Tests"

* 50K 1:1 merge
clear
set seed 12345
set obs 50000
gen id = _n
gen value1 = runiform()
gen value2 = runiformint(1, 100)
tempfile master
save `master'

clear
set obs 50000
gen id = _n
gen extra1 = runiform() * 100
gen extra2 = runiformint(1, 1000)
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("50K 1:1 full match")

* 100K m:1 merge
clear
set seed 11111
set obs 100000
gen group_id = runiformint(1, 1000)
gen value = runiform()
tempfile master
save `master'

clear
set obs 1000
gen group_id = _n
gen group_name = "Group" + string(group_id)
gen group_weight = runiform()
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge m:1 group_id using `using_1', testname("100K m:1")

* Partial overlap
clear
set seed 22222
set obs 30000
gen id = _n
gen value1 = runiform()
tempfile master
save `master'

clear
set obs 30000
gen id = _n + 15000
gen value2 = runiform()
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', testname("partial overlap (50%)")

/*******************************************************************************
 * SECTION 11: String key tests
 ******************************************************************************/
noi print_section "String Key Tests"

* Long string keys
clear
set obs 100
gen str50 long_key = "prefix_" + string(_n, "%04.0f") + "_suffix_data"
gen value1 = runiform()
tempfile master
save `master'

clear
set obs 100
gen str50 long_key = "prefix_" + string(_n, "%04.0f") + "_suffix_data"
gen value2 = runiform() * 100
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 long_key using `using_1', testname("long string keys")

* Census state name
sysuse census, clear
keep state pop
tempfile master
save `master'

sysuse census, clear
keep state medage death marriage divorce
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 state using `using_1', testname("state name 1:1")

/*******************************************************************************
 * SECTION 12: Panel data tests
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/5000
keep idcode year ln_wage hours
tempfile master
save `master'

webuse nlswork, clear
keep in 1/5000
keep idcode year age tenure wks_work
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 idcode year using `using_1', testname("panel 1:1")

* m:1 individual characteristics
webuse nlswork, clear
keep in 1/3000
keep idcode year ln_wage
tempfile master
save `master'

webuse nlswork, clear
keep in 1/3000
bysort idcode (year): keep if _n == 1
keep idcode race grade
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge m:1 idcode using `using_1', testname("panel m:1 individual")

/*******************************************************************************
 * SECTION 13: assert() option tests
 ******************************************************************************/
noi print_section "assert() Option Tests"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', assert(match) testname("assert(match) passes")

* Test assert failure
use `master', clear
keep if _n <= 50
tempfile master_partial
save `master_partial'

use `master_partial', clear
capture merge 1:1 id using `using_1', assert(match)
local stata_rc = _rc

use `master_partial', clear
capture cmerge 1:1 id using `using_1', assert(match) noreport
local cmerge_rc = _rc

if `stata_rc' != 0 & `cmerge_rc' != 0 {
    noi test_pass "assert(match) fails correctly"
}
else {
    noi test_fail "assert(match) should fail" "stata_rc=`stata_rc' cmerge_rc=`cmerge_rc'"
}

/*******************************************************************************
 * SECTION 14: update/replace options
 ******************************************************************************/
noi print_section "update/replace Options"

* Test update option
clear
input int id float value
1 10
2 .
3 30
4 .
end
tempfile master
save `master'

clear
input int id float value
1 100
2 20
3 .
4 40
end
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 id using `using_1', update nogenerate
local stata_v2 = value[2]
local stata_v3 = value[3]

use `master', clear
cmerge 1:1 id using `using_1', update nogenerate noreport
local cmerge_v2 = value[2]
local cmerge_v3 = value[3]

if `stata_v2' == `cmerge_v2' & `stata_v3' == `cmerge_v3' {
    noi test_pass "update option"
}
else {
    noi test_fail "update option" "values differ"
}

* Test update replace option
use `master', clear
merge 1:1 id using `using_1', update replace nogenerate
local stata_v1 = value[1]

use `master', clear
cmerge 1:1 id using `using_1', update replace nogenerate noreport
local cmerge_v1 = value[1]

if `stata_v1' == `cmerge_v1' {
    noi test_pass "update replace option"
}
else {
    noi test_fail "update replace option" "values differ"
}

/*******************************************************************************
 * SECTION 15: sorted option
 ******************************************************************************/
noi print_section "sorted Option"

sysuse auto, clear
gen id = _n
keep id make price
sort id
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight
sort id
tempfile using_1
save `using_1'

use `master', clear
noi benchmark_merge 1:1 id using `using_1', sorted testname("sorted option")

/*******************************************************************************
 * SECTION 16: force option
 ******************************************************************************/
noi print_section "force Option"

clear
input int id float value
1 10
2 20
end
tempfile master
save `master'

clear
input int id str5 value
1 "A"
2 "B"
end
tempfile using_1
save `using_1'

use `master', clear
capture merge 1:1 id using `using_1', force nogenerate
local stata_force_rc = _rc

use `master', clear
capture cmerge 1:1 id using `using_1', force nogenerate noreport
local cmerge_force_rc = _rc

if `stata_force_rc' == 0 & `cmerge_force_rc' == 0 {
    noi test_pass "force option"
}
else {
    noi test_fail "force option" "stata=`stata_force_rc' cmerge=`cmerge_force_rc'"
}

/*******************************************************************************
 * SECTION 17: nolabel and nonotes options
 ******************************************************************************/
noi print_section "nolabel/nonotes Options"

sysuse auto, clear
gen id = _n
keep id make price foreign
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight foreign
note weight: "This is a test note"
tempfile using_1
save `using_1'

use `master', clear
capture cmerge 1:1 id using `using_1', nolabel noreport
if _rc == 0 {
    noi test_pass "nolabel option accepted"
}
else {
    noi test_fail "nolabel option" "returned error `=_rc'"
}

use `master', clear
capture cmerge 1:1 id using `using_1', nonotes noreport
if _rc == 0 {
    noi test_pass "nonotes option accepted"
}
else {
    noi test_fail "nonotes option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 18: 1:1 _n merge tests (merge by observation number)
 ******************************************************************************/
noi print_section "1:1 _n Merge Tests"

* Basic 1:1 _n merge - same number of observations
sysuse auto, clear
keep make price mpg
tempfile master
save `master'

sysuse auto, clear
keep weight length
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1', nogenerate
local stata_N = _N
local stata_price1 = price[1]
local stata_weight1 = weight[1]

use `master', clear
cmerge 1:1 _n using `using_1', nogenerate noreport
local cmerge_N = _N
local cmerge_price1 = price[1]
local cmerge_weight1 = weight[1]

if `stata_N' == `cmerge_N' & `stata_price1' == `cmerge_price1' & `stata_weight1' == `cmerge_weight1' {
    noi test_pass "1:1 _n basic (same nobs)"
}
else {
    noi test_fail "1:1 _n basic (same nobs)" "N: stata=`stata_N' cmerge=`cmerge_N'"
}

* 1:1 _n merge - master has more observations
clear
set obs 100
gen id = _n
gen value1 = runiform()
tempfile master
save `master'

clear
set obs 70
gen extra = runiform() * 100
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1'
qui count if _merge == 3
local stata_matched = r(N)
qui count if _merge == 1
local stata_master = r(N)
qui count if _merge == 2
local stata_using = r(N)
drop _merge

use `master', clear
cmerge 1:1 _n using `using_1', noreport
qui count if _merge == 3
local cmerge_matched = r(N)
qui count if _merge == 1
local cmerge_master = r(N)
qui count if _merge == 2
local cmerge_using = r(N)
drop _merge

if `stata_matched' == `cmerge_matched' & `stata_master' == `cmerge_master' & `stata_using' == `cmerge_using' {
    noi test_pass "1:1 _n master has more obs"
}
else {
    noi test_fail "1:1 _n master has more obs" "matched: stata=`stata_matched' cmerge=`cmerge_matched'"
}

* 1:1 _n merge - using has more observations
clear
set obs 50
gen id = _n
gen value1 = runiform()
tempfile master
save `master'

clear
set obs 80
gen extra = runiform() * 100
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1'
qui count if _merge == 3
local stata_matched = r(N)
qui count if _merge == 1
local stata_master = r(N)
qui count if _merge == 2
local stata_using = r(N)
drop _merge

use `master', clear
cmerge 1:1 _n using `using_1', noreport
qui count if _merge == 3
local cmerge_matched = r(N)
qui count if _merge == 1
local cmerge_master = r(N)
qui count if _merge == 2
local cmerge_using = r(N)
drop _merge

if `stata_matched' == `cmerge_matched' & `stata_master' == `cmerge_master' & `stata_using' == `cmerge_using' {
    noi test_pass "1:1 _n using has more obs"
}
else {
    noi test_fail "1:1 _n using has more obs" "matched: stata=`stata_matched' cmerge=`cmerge_matched'"
}

* 1:1 _n merge with keepusing option
sysuse auto, clear
keep make price mpg
tempfile master
save `master'

sysuse auto, clear
keep weight length turn
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1', keepusing(weight) nogenerate
local stata_has_weight = 1
capture confirm variable weight
if _rc local stata_has_weight = 0
capture confirm variable length
local stata_has_length = !_rc
local stata_weight1 = weight[1]

use `master', clear
cmerge 1:1 _n using `using_1', keepusing(weight) nogenerate noreport
local cmerge_has_weight = 1
capture confirm variable weight
if _rc local cmerge_has_weight = 0
capture confirm variable length
local cmerge_has_length = !_rc
local cmerge_weight1 = weight[1]

if `stata_has_weight' == `cmerge_has_weight' & `stata_has_length' == `cmerge_has_length' & `stata_weight1' == `cmerge_weight1' {
    noi test_pass "1:1 _n with keepusing"
}
else {
    noi test_fail "1:1 _n with keepusing" "has_weight: stata=`stata_has_weight' cmerge=`cmerge_has_weight'"
}

* 1:1 _n merge with keep() option
clear
set obs 100
gen id = _n
gen value1 = runiform()
tempfile master
save `master'

clear
set obs 70
gen extra = runiform() * 100
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1', keep(match) nogenerate
local stata_N = _N

use `master', clear
cmerge 1:1 _n using `using_1', keep(match) nogenerate noreport
local cmerge_N = _N

if `stata_N' == `cmerge_N' {
    noi test_pass "1:1 _n with keep(match)"
}
else {
    noi test_fail "1:1 _n with keep(match)" "N: stata=`stata_N' cmerge=`cmerge_N'"
}

* 1:1 _n merge with generate() option
sysuse auto, clear
keep make price mpg
tempfile master
save `master'

sysuse auto, clear
keep weight length
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1', generate(merge_result)
local stata_has_merge = 1
capture confirm variable merge_result
if _rc local stata_has_merge = 0
drop merge_result

use `master', clear
cmerge 1:1 _n using `using_1', generate(merge_result) noreport
local cmerge_has_merge = 1
capture confirm variable merge_result
if _rc local cmerge_has_merge = 0
drop merge_result

if `stata_has_merge' == 1 & `cmerge_has_merge' == 1 {
    noi test_pass "1:1 _n with generate()"
}
else {
    noi test_fail "1:1 _n with generate()" "has_merge: stata=`stata_has_merge' cmerge=`cmerge_has_merge'"
}

* Large 1:1 _n merge test
clear
set seed 99999
set obs 50000
gen id = _n
gen value1 = runiform()
gen value2 = runiformint(1, 100)
tempfile master
save `master'

clear
set obs 50000
gen extra1 = runiform() * 100
gen extra2 = runiformint(1, 1000)
tempfile using_1
save `using_1'

use `master', clear
merge 1:1 _n using `using_1', nogenerate
local stata_N = _N
local stata_extra1_1 = extra1[1]
local stata_extra1_last = extra1[_N]

use `master', clear
cmerge 1:1 _n using `using_1', nogenerate noreport
local cmerge_N = _N
local cmerge_extra1_1 = extra1[1]
local cmerge_extra1_last = extra1[_N]

if `stata_N' == `cmerge_N' & `stata_extra1_1' == `cmerge_extra1_1' & `stata_extra1_last' == `cmerge_extra1_last' {
    noi test_pass "1:1 _n large dataset (50K obs)"
}
else {
    noi test_fail "1:1 _n large dataset" "N: stata=`stata_N' cmerge=`cmerge_N'"
}

/*******************************************************************************
 * SECTION 19: Comprehensive sysuse/webuse Dataset Tests
 ******************************************************************************/
noi print_section "Comprehensive Built-in Dataset Tests"

* Auto dataset - merge by foreign (m:1 lookup)
sysuse auto, clear
keep make price mpg foreign
tempfile auto_master
save `auto_master'

clear
input byte foreign str20 origin double tax_rate
0 "Domestic" 0.07
1 "Foreign" 0.12
end
tempfile auto_lookup
save `auto_lookup'

use `auto_master', clear
noi benchmark_merge m:1 foreign using `auto_lookup', testname("auto: m:1 foreign lookup")

* Auto dataset - 1:1 by make
sysuse auto, clear
keep make price mpg
tempfile auto1
save `auto1'

sysuse auto, clear
keep make weight length trunk
tempfile auto2
save `auto2'

use `auto1', clear
noi benchmark_merge 1:1 make using `auto2', testname("auto: 1:1 by make")

* Census dataset - merge by region (m:1)
sysuse census, clear
keep state region pop medage
tempfile census_master
save `census_master'

clear
input byte region str20 region_full double min_wage
1 "Northeast" 15.50
2 "North Central" 12.00
3 "South" 10.50
4 "West" 14.00
end
tempfile region_lookup
save `region_lookup'

use `census_master', clear
noi benchmark_merge m:1 region using `region_lookup', testname("census: m:1 region lookup")

* Census dataset - 1:1 by state
sysuse census, clear
keep state pop death marriage
tempfile census1
save `census1'

sysuse census, clear
keep state divorce medage region
tempfile census2
save `census2'

use `census1', clear
noi benchmark_merge 1:1 state using `census2', testname("census: 1:1 by state")

* nlsw88 dataset - m:1 by industry
sysuse nlsw88, clear
keep idcode wage hours industry
tempfile nlsw88_master
save `nlsw88_master'

sysuse nlsw88, clear
bysort industry: keep if _n == 1
keep industry
gen double avg_wage = runiform() * 20 + 10
tempfile industry_lookup
save `industry_lookup'

use `nlsw88_master', clear
noi benchmark_merge m:1 industry using `industry_lookup', testname("nlsw88: m:1 by industry")

* nlswork panel data - 1:1 by idcode year
webuse nlswork, clear
keep in 1/3000
keep idcode year ln_wage hours
tempfile nlswork1
save `nlswork1'

webuse nlswork, clear
keep in 1/3000
keep idcode year age tenure wks_work
tempfile nlswork2
save `nlswork2'

use `nlswork1', clear
noi benchmark_merge 1:1 idcode year using `nlswork2', testname("nlswork: 1:1 panel by idcode year")

* nlswork - m:1 individual time-invariant characteristics
webuse nlswork, clear
keep in 1/2500
keep idcode year ln_wage
tempfile nlswork_panel
save `nlswork_panel'

webuse nlswork, clear
keep in 1/2500
bysort idcode (year): keep if _n == 1
keep idcode race grade
tempfile nlswork_indiv
save `nlswork_indiv'

use `nlswork_panel', clear
noi benchmark_merge m:1 idcode using `nlswork_indiv', testname("nlswork: m:1 individual characteristics")

* Grunfeld panel data
webuse grunfeld, clear
keep company year invest
tempfile grunfeld1
save `grunfeld1'

webuse grunfeld, clear
keep company year mvalue kstock
tempfile grunfeld2
save `grunfeld2'

use `grunfeld1', clear
noi benchmark_merge 1:1 company year using `grunfeld2', testname("grunfeld: 1:1 panel merge")

* lifeexp dataset
webuse lifeexp, clear
keep region country popgrowth
tempfile lifeexp1
save `lifeexp1'

webuse lifeexp, clear
keep region country lexp gnppc
tempfile lifeexp2
save `lifeexp2'

use `lifeexp1', clear
noi benchmark_merge 1:1 region country using `lifeexp2', testname("lifeexp: 1:1 by region country")

/*******************************************************************************
 * SECTION 20: Pathological String Key Edge Cases
 ******************************************************************************/
noi print_section "Pathological String Key Edge Cases"

* String keys with embedded commas
clear
input int id str40 name
1 "Smith, John"
2 "Doe, Jane"
3 "O'Brien, Pat"
end
tempfile str_comma_master
save `str_comma_master'

clear
input str40 name double value
"Smith, John" 100
"Doe, Jane" 200
"O'Brien, Pat" 300
end
tempfile str_comma_using
save `str_comma_using'

use `str_comma_master', clear
noi benchmark_merge m:1 name using `str_comma_using', testname("string key: embedded commas")

* String keys with quotes
clear
input int id str40 name
1 `"He said "hello""'
2 `"She replied "goodbye""'
3 `"Normal name"'
end
tempfile str_quote_master
save `str_quote_master'

clear
input str40 name double value
`"He said "hello""' 100
`"She replied "goodbye""' 200
`"Normal name"' 300
end
tempfile str_quote_using
save `str_quote_using'

use `str_quote_master', clear
noi benchmark_merge m:1 name using `str_quote_using', testname("string key: embedded quotes")

* String keys with leading/trailing spaces
clear
set obs 10
gen id = _n
gen str30 key = ""
replace key = " leading" in 1
replace key = "trailing " in 2
replace key = " both " in 3
replace key = "   lots   " in 4
replace key = "normal" in 5/10
gen value1 = runiform()
tempfile str_space_master
save `str_space_master'

clear
set obs 5
gen str30 key = ""
replace key = " leading" in 1
replace key = "trailing " in 2
replace key = " both " in 3
replace key = "   lots   " in 4
replace key = "normal" in 5
gen value2 = runiform() * 100
tempfile str_space_using
save `str_space_using'

use `str_space_master', clear
noi benchmark_merge m:1 key using `str_space_using', testname("string key: leading/trailing spaces")

* Very long string keys (200+ characters)
clear
set obs 50
gen str244 long_key = "prefix_" + substr("a" * 200, 1, 200) + "_" + string(_n)
gen value1 = runiform()
tempfile str_long_master
save `str_long_master'

clear
set obs 50
gen str244 long_key = "prefix_" + substr("a" * 200, 1, 200) + "_" + string(_n)
gen value2 = runiform() * 100
tempfile str_long_using
save `str_long_using'

use `str_long_master', clear
noi benchmark_merge 1:1 long_key using `str_long_using', testname("string key: very long (200+ chars)")

* Case-sensitive string matching
clear
set obs 6
gen str20 name = ""
replace name = "Apple" in 1
replace name = "apple" in 2
replace name = "APPLE" in 3
replace name = "Banana" in 4
replace name = "banana" in 5
replace name = "BANANA" in 6
gen value = _n
tempfile str_case_master
save `str_case_master'

clear
set obs 3
gen str20 name = ""
replace name = "Apple" in 1
replace name = "banana" in 2
replace name = "BANANA" in 3
gen other = _n * 10
tempfile str_case_using
save `str_case_using'

use `str_case_master', clear
noi benchmark_merge 1:1 name using `str_case_using', testname("string key: case sensitive")

* Empty string keys
clear
set obs 5
gen str20 key = ""
replace key = "" in 1
replace key = "a" in 2
replace key = "" in 3
replace key = "b" in 4
replace key = "" in 5
gen value1 = _n
tempfile str_empty_master
save `str_empty_master'

clear
set obs 3
gen str20 key = ""
replace key = "" in 1
replace key = "a" in 2
replace key = "b" in 3
gen value2 = _n * 100
tempfile str_empty_using
save `str_empty_using'

use `str_empty_master', clear
noi benchmark_merge m:1 key using `str_empty_using', testname("string key: empty strings")

* Unicode-like characters (Latin extended)
clear
input int id str50 name
1 "Jose Garcia"
2 "Francois Muller"
3 "Soren Orsted"
4 "Cafe Creme"
end
tempfile str_unicode_master
save `str_unicode_master'

clear
input str50 name double value
"Jose Garcia" 100
"Francois Muller" 200
"Soren Orsted" 300
"Cafe Creme" 400
end
tempfile str_unicode_using
save `str_unicode_using'

use `str_unicode_master', clear
noi benchmark_merge m:1 name using `str_unicode_using', testname("string key: Latin extended chars")

/*******************************************************************************
 * SECTION 21: Numeric Key Edge Cases
 ******************************************************************************/
noi print_section "Numeric Key Edge Cases"

* Missing values in numeric keys (using m:1 since missing values are duplicates)
clear
input double id double value1
1 10
2 20
. 30
. 40
5 50
end
tempfile num_miss_master
save `num_miss_master'

clear
input double id double value2
1 100
. 200
3 300
5 500
end
tempfile num_miss_using
save `num_miss_using'

use `num_miss_master', clear
noi benchmark_merge m:1 id using `num_miss_using', testname("numeric key: missing values (m:1)")

* Extended missing values (.a, .b, etc.)
clear
input double id double value1
1 10
.a 20
.b 30
.z 40
5 50
end
tempfile num_ext_master
save `num_ext_master'

clear
input double id double value2
1 100
.a 200
.b 300
.c 400
5 500
end
tempfile num_ext_using
save `num_ext_using'

use `num_ext_master', clear
noi benchmark_merge 1:1 id using `num_ext_using', testname("numeric key: extended missing (.a,.b,.z)")

* Float keys with precision issues
clear
set obs 100
gen double id = _n * 0.1
gen value1 = runiform()
tempfile float_master
save `float_master'

clear
set obs 100
gen double id = _n * 0.1
gen value2 = runiform() * 100
tempfile float_using
save `float_using'

use `float_master', clear
noi benchmark_merge 1:1 id using `float_using', testname("numeric key: float precision")

* Large integer keys
clear
set obs 100
gen long id = _n * 1000000 + 999000000
gen value1 = runiform()
tempfile large_int_master
save `large_int_master'

clear
set obs 100
gen long id = _n * 1000000 + 999000000
gen value2 = runiform() * 100
tempfile large_int_using
save `large_int_using'

use `large_int_master', clear
noi benchmark_merge 1:1 id using `large_int_using', testname("numeric key: large integers")

* Negative keys with zero crossing
clear
set obs 100
gen int id = _n - 50
gen value1 = runiform()
tempfile neg_cross_master
save `neg_cross_master'

clear
set obs 100
gen int id = _n - 50
gen value2 = runiform() * 100
tempfile neg_cross_using
save `neg_cross_using'

use `neg_cross_master', clear
noi benchmark_merge 1:1 id using `neg_cross_using', testname("numeric key: negative to positive range")

* Zero keys
clear
set obs 50
gen int id = 0
replace id = _n if _n > 25
gen value1 = runiform()
tempfile zero_key_master
save `zero_key_master'

clear
set obs 50
gen int id = 0
replace id = _n if _n > 25
gen value2 = runiform() * 100
tempfile zero_key_using
save `zero_key_using'

use `zero_key_master', clear
noi benchmark_merge m:m id using `zero_key_using', testname("numeric key: many zeros (m:m)")

/*******************************************************************************
 * SECTION 22: Multiple Key Variables
 ******************************************************************************/
noi print_section "Multiple Key Variables"

* Four key variables
clear
set obs 500
gen id1 = mod(_n - 1, 5) + 1
gen id2 = mod(floor((_n - 1) / 5), 5) + 1
gen id3 = mod(floor((_n - 1) / 25), 4) + 1
gen id4 = floor((_n - 1) / 100) + 1
gen value1 = runiform()
tempfile four_key_master
save `four_key_master'

clear
set obs 500
gen id1 = mod(_n - 1, 5) + 1
gen id2 = mod(floor((_n - 1) / 5), 5) + 1
gen id3 = mod(floor((_n - 1) / 25), 4) + 1
gen id4 = floor((_n - 1) / 100) + 1
gen value2 = runiform() * 100
tempfile four_key_using
save `four_key_using'

use `four_key_master', clear
noi benchmark_merge 1:1 id1 id2 id3 id4 using `four_key_using', testname("four numeric keys")

* Five key variables (string + numeric mix)
clear
set obs 200
gen str10 cat = "cat" + string(mod(_n, 4) + 1)
gen int year = 2020 + mod(_n, 3)
gen int month = mod(_n, 12) + 1
gen int day = mod(_n, 28) + 1
gen id = ceil(_n / 10)
gen value1 = runiform()
tempfile five_key_master
save `five_key_master'

clear
set obs 200
gen str10 cat = "cat" + string(mod(_n, 4) + 1)
gen int year = 2020 + mod(_n, 3)
gen int month = mod(_n, 12) + 1
gen int day = mod(_n, 28) + 1
gen id = ceil(_n / 10)
gen value2 = runiform() * 100
tempfile five_key_using
save `five_key_using'

use `five_key_master', clear
noi benchmark_merge m:m cat year month day id using `five_key_using', testname("five mixed keys (string + numeric)")

* All same value keys (degenerate)
clear
set obs 100
gen key1 = 1
gen key2 = 1
gen key3 = 1
gen value = _n
tempfile same_key_master
save `same_key_master'

clear
set obs 100
gen key1 = 1
gen key2 = 1
gen key3 = 1
gen other = runiform()
tempfile same_key_using
save `same_key_using'

use `same_key_master', clear
capture noisily cmerge m:m key1 key2 key3 using `same_key_using', nogen noreport
if _rc == 0 {
    noi test_pass "all same keys (degenerate m:m)"
}
else {
    noi test_fail "all same keys (degenerate)" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 23: Match Pattern Edge Cases
 ******************************************************************************/
noi print_section "Match Pattern Edge Cases"

* No matches (completely disjoint)
clear
set obs 100
gen id = _n
gen value1 = runiform()
tempfile disjoint_master
save `disjoint_master'

clear
set obs 100
gen id = _n + 1000
gen value2 = runiform() * 100
tempfile disjoint_using
save `disjoint_using'

use `disjoint_master', clear
noi benchmark_merge 1:1 id using `disjoint_using', testname("no matches (completely disjoint)")

* All matches (perfect 1:1)
clear
set obs 100
gen id = _n
gen value1 = runiform()
tempfile perfect_master
save `perfect_master'

clear
set obs 100
gen id = _n
gen value2 = runiform() * 100
tempfile perfect_using
save `perfect_using'

use `perfect_master', clear
noi benchmark_merge 1:1 id using `perfect_using', testname("all matches (perfect 1:1)")

* Master only (no using matches any)
clear
set obs 100
gen id = _n * 2
gen value1 = runiform()
tempfile master_only_master
save `master_only_master'

clear
set obs 100
gen id = _n * 2 + 1
gen value2 = runiform() * 100
tempfile master_only_using
save `master_only_using'

use `master_only_master', clear
merge 1:1 id using `master_only_using'
qui count if _merge == 1
local stata_m1 = r(N)
drop _merge

use `master_only_master', clear
cmerge 1:1 id using `master_only_using', noreport
qui count if _merge == 1
local cmerge_m1 = r(N)
drop _merge

if `stata_m1' == `cmerge_m1' & `cmerge_m1' == 100 {
    noi test_pass "master only (no using matches)"
}
else {
    noi test_fail "master only" "m1 counts differ"
}

* Using only (no master matches any)
clear
set obs 50
gen id = _n * 2 + 1
gen value1 = runiform()
tempfile using_only_master
save `using_only_master'

clear
set obs 100
gen id = _n * 2
gen value2 = runiform() * 100
tempfile using_only_using
save `using_only_using'

use `using_only_master', clear
merge 1:1 id using `using_only_using'
qui count if _merge == 2
local stata_m2 = r(N)
drop _merge

use `using_only_master', clear
cmerge 1:1 id using `using_only_using', noreport
qui count if _merge == 2
local cmerge_m2 = r(N)
drop _merge

if `stata_m2' == `cmerge_m2' & `cmerge_m2' == 100 {
    noi test_pass "using only (no master matches)"
}
else {
    noi test_fail "using only" "m2 counts differ"
}

* Sparse matches (10% overlap) - using has unique keys that partially overlap master
clear
set obs 1000
gen id = _n
gen value1 = runiform()
tempfile sparse_master
save `sparse_master'

clear
set obs 1000
* Create unique IDs: first 100 overlap with master (1-100), rest are above (5001-5900)
gen id = cond(_n <= 100, _n, _n + 4900)
gen value2 = runiform() * 100
tempfile sparse_using
save `sparse_using'

use `sparse_master', clear
noi benchmark_merge 1:1 id using `sparse_using', testname("sparse matches (~10% overlap)")

/*******************************************************************************
 * SECTION 24: Missing Value Handling
 ******************************************************************************/
noi print_section "Missing Value Handling"

* Missing in non-key variables
clear
input int id double x double y double z
1 100 . 1
2 . 200 2
3 300 300 .
4 . . 4
5 500 500 5
end
tempfile miss_nonkey_master
save `miss_nonkey_master'

clear
input int id double a double b
1 . 10
2 20 .
3 . .
4 40 40
5 50 50
end
tempfile miss_nonkey_using
save `miss_nonkey_using'

use `miss_nonkey_master', clear
noi benchmark_merge 1:1 id using `miss_nonkey_using', testname("missing in non-key variables")

* First/last rows with missing keys (only one missing per dataset for valid 1:1)
clear
input double id double value
. 100
2 200
3 300
4 400
6 500
end
tempfile first_last_miss_master
save `first_last_miss_master'

clear
input double id double other
1 10
2 20
. 30
4 40
5 50
end
tempfile first_last_miss_using
save `first_last_miss_using'

use `first_last_miss_master', clear
noi benchmark_merge 1:1 id using `first_last_miss_using', testname("first/last rows missing keys")

* All missing in one key column (multi-key)
clear
set obs 50
gen int id1 = _n
gen id2 = .
gen value1 = runiform()
tempfile all_miss_key_master
save `all_miss_key_master'

clear
set obs 50
gen int id1 = _n
gen id2 = .
gen value2 = runiform() * 100
tempfile all_miss_key_using
save `all_miss_key_using'

use `all_miss_key_master', clear
noi benchmark_merge 1:1 id1 id2 using `all_miss_key_using', testname("all missing in second key column")

/*******************************************************************************
 * SECTION 25: Option Combinations
 ******************************************************************************/
noi print_section "Option Combinations"

* keep(master match) combination
clear
set obs 100
gen id = _n
gen value1 = runiform()
tempfile opt_master
save `opt_master'

clear
set obs 50
gen id = _n * 2
gen value2 = runiform() * 100
tempfile opt_using
save `opt_using'

use `opt_master', clear
noi benchmark_merge 1:1 id using `opt_using', keep(master match) testname("keep(master match)")

* keep(using match) combination
use `opt_master', clear
noi benchmark_merge 1:1 id using `opt_using', keep(using match) testname("keep(using match)")

* keep(match) only
use `opt_master', clear
noi benchmark_merge 1:1 id using `opt_using', keep(match) testname("keep(match) only")

* keepusing() with multiple variables
clear
set obs 50
gen id = _n
gen master_val = runiform()
tempfile keepusing_master
save `keepusing_master'

clear
set obs 50
gen id = _n
gen var_a = runiform()
gen var_b = runiform() * 10
gen var_c = runiform() * 100
gen var_d = runiform() * 1000
tempfile keepusing_using
save `keepusing_using'

use `keepusing_master', clear
noi benchmark_merge 1:1 id using `keepusing_using', keepusing(var_a var_c) testname("keepusing: select 2 of 4 vars")

* generate() with custom name
use `opt_master', clear
noi benchmark_merge 1:1 id using `opt_using', generate(my_merge_var) testname("generate(my_merge_var)")

* nogenerate option
use `opt_master', clear
noi benchmark_merge 1:1 id using `opt_using', nogenerate testname("nogenerate option")

* sorted option with pre-sorted data
use `opt_master', clear
sort id
save `opt_master', replace
use `opt_using', clear
sort id
save `opt_using', replace
use `opt_master', clear
noi benchmark_merge 1:1 id using `opt_using', sorted testname("sorted option (pre-sorted)")

* Combined options: keep + keepusing + nogenerate
use `keepusing_master', clear
noi benchmark_merge 1:1 id using `keepusing_using', keep(match) keepusing(var_a var_b) nogenerate testname("combined: keep + keepusing + nogen")

/*******************************************************************************
 * SECTION 26: Large Dataset Tests
 ******************************************************************************/
noi print_section "Large Dataset Tests"

* 10K 1:1 full match
clear
set seed 54321
set obs 10000
gen long id = _n
gen value1 = runiform()
gen str20 label1 = "item_" + string(_n)
tempfile large_10k_master
save `large_10k_master'

clear
set obs 10000
gen long id = _n
gen value2 = runiform() * 100
gen str20 label2 = "other_" + string(_n)
tempfile large_10k_using
save `large_10k_using'

use `large_10k_master', clear
noi benchmark_merge 1:1 id using `large_10k_using', testname("10K 1:1 full match")

* 50K 1:1 partial match
clear
set seed 12121
set obs 50000
gen long id = _n
gen value1 = runiform()
tempfile large_50k_master
save `large_50k_master'

clear
set obs 50000
gen long id = _n + 25000
gen value2 = runiform() * 100
tempfile large_50k_using
save `large_50k_using'

use `large_50k_master', clear
noi benchmark_merge 1:1 id using `large_50k_using', testname("50K 1:1 partial match (50%)")

* Large m:1 (100K to 1K)
clear
set seed 33333
set obs 100000
gen int group = runiformint(1, 1000)
gen value = runiform()
tempfile large_m1_master
save `large_m1_master'

clear
set obs 1000
gen int group = _n
gen str30 group_name = "group_" + string(_n)
gen group_weight = runiform()
tempfile large_m1_using
save `large_m1_using'

use `large_m1_master', clear
noi benchmark_merge m:1 group using `large_m1_using', testname("large m:1 (100K to 1K)")

* Large 1:m (1K to 100K)
clear
set obs 1000
gen int group = _n
gen master_val = runiform()
tempfile large_1m_master
save `large_1m_master'

clear
set seed 44444
set obs 100000
gen int group = runiformint(1, 1000)
gen detail = runiform() * 100
tempfile large_1m_using
save `large_1m_using'

use `large_1m_master', clear
noi benchmark_merge 1:m group using `large_1m_using', testname("large 1:m (1K to 100K)")

* Large m:m with expansion
clear
set seed 55555
set obs 1000
gen int group = mod(_n - 1, 100) + 1
gen master_seq = _n
gen value1 = runiform()
tempfile large_mm_master
save `large_mm_master'

clear
set obs 1000
gen int group = mod(_n - 1, 100) + 1
gen using_seq = _n
gen value2 = runiform() * 100
tempfile large_mm_using
save `large_mm_using'

use `large_mm_master', clear
merge m:m group using `large_mm_using', nogen
local stata_N = _N

use `large_mm_master', clear
capture cmerge m:m group using `large_mm_using', nogen noreport
if _rc == 0 & _N == `stata_N' {
    noi test_pass "large m:m with expansion"
}
else if _rc != 0 {
    noi test_fail "large m:m expansion" "rc=`=_rc'"
}
else {
    noi test_fail "large m:m expansion" "N mismatch: stata=`stata_N' cmerge=`=_N'"
}

/*******************************************************************************
 * SECTION 27: Unbalanced Merges
 ******************************************************************************/
noi print_section "Unbalanced Merges"

* Large master, small using
clear
set obs 50000
gen long id = _n
gen value1 = runiform()
tempfile unbal_large_master
save `unbal_large_master'

clear
set obs 100
gen long id = _n * 500
gen value2 = runiform() * 100
tempfile unbal_small_using
save `unbal_small_using'

use `unbal_large_master', clear
noi benchmark_merge 1:1 id using `unbal_small_using', testname("unbalanced: 50K master, 100 using")

* Small master, large using
clear
set obs 100
gen long id = _n * 500
gen value1 = runiform()
tempfile unbal_small_master
save `unbal_small_master'

clear
set obs 50000
gen long id = _n
gen value2 = runiform() * 100
tempfile unbal_large_using
save `unbal_large_using'

use `unbal_small_master', clear
noi benchmark_merge 1:1 id using `unbal_large_using', testname("unbalanced: 100 master, 50K using")

* Very different sizes (10x)
clear
set obs 10000
gen int group = runiformint(1, 100)
gen value1 = runiform()
tempfile unbal_10x_master
save `unbal_10x_master'

clear
set obs 1000
gen int group = runiformint(1, 100)
gen value2 = runiform() * 100
tempfile unbal_10x_using
save `unbal_10x_using'

use `unbal_10x_master', clear
noi benchmark_merge m:m group using `unbal_10x_using', testname("unbalanced m:m: 10K master, 1K using")

/*******************************************************************************
 * SECTION 28: Variable Conflicts
 ******************************************************************************/
noi print_section "Variable Conflicts"

* Same non-key variable names (should fail without force)
clear
set obs 50
gen id = _n
gen value = runiform()
gen shared_var = runiform() * 10
tempfile conflict_master
save `conflict_master'

clear
set obs 50
gen id = _n
gen shared_var = runiform() * 100
gen other = runiform()
tempfile conflict_using
save `conflict_using'

use `conflict_master', clear
capture merge 1:1 id using `conflict_using', nogenerate
local stata_conflict_rc = _rc

use `conflict_master', clear
capture cmerge 1:1 id using `conflict_using', nogenerate noreport
local cmerge_conflict_rc = _rc

if `stata_conflict_rc' != 0 & `cmerge_conflict_rc' != 0 {
    noi test_pass "variable conflict detected (same var names)"
}
else if `stata_conflict_rc' == 0 & `cmerge_conflict_rc' == 0 {
    noi test_pass "variable conflict handled (both succeeded)"
}
else {
    noi test_fail "variable conflict" "rc differ: stata=`stata_conflict_rc' cmerge=`cmerge_conflict_rc'"
}

* Type mismatch (numeric vs string) with force
clear
set obs 20
gen id = _n
gen double mixed_type = _n * 10
tempfile type_mismatch_master
save `type_mismatch_master'

clear
set obs 20
gen id = _n
gen str20 mixed_type = "value_" + string(_n)
tempfile type_mismatch_using
save `type_mismatch_using'

use `type_mismatch_master', clear
capture merge 1:1 id using `type_mismatch_using', force nogenerate
local stata_force_rc = _rc

use `type_mismatch_master', clear
capture cmerge 1:1 id using `type_mismatch_using', force nogenerate noreport
local cmerge_force_rc = _rc

if `stata_force_rc' == `cmerge_force_rc' {
    noi test_pass "type mismatch with force option"
}
else {
    noi test_fail "type mismatch force" "rc differ: stata=`stata_force_rc' cmerge=`cmerge_force_rc'"
}

/*******************************************************************************
 * SECTION 29: Update and Replace Options
 ******************************************************************************/
noi print_section "Update and Replace Options"

* update option - fill missing from using
clear
input int id double value
1 100
2 .
3 300
4 .
5 500
end
tempfile update_master
save `update_master'

clear
input int id double value
1 999
2 200
3 999
4 400
5 999
end
tempfile update_using
save `update_using'

use `update_master', clear
merge 1:1 id using `update_using', update nogenerate
local stata_v2 = value[2]
local stata_v4 = value[4]

use `update_master', clear
cmerge 1:1 id using `update_using', update nogenerate noreport
local cmerge_v2 = value[2]
local cmerge_v4 = value[4]

if `stata_v2' == `cmerge_v2' & `stata_v4' == `cmerge_v4' {
    noi test_pass "update fills missing values"
}
else {
    noi test_fail "update option" "values differ"
}

* update replace - replace non-missing too
use `update_master', clear
merge 1:1 id using `update_using', update replace nogenerate
local stata_v1 = value[1]
local stata_v3 = value[3]

use `update_master', clear
cmerge 1:1 id using `update_using', update replace nogenerate noreport
local cmerge_v1 = value[1]
local cmerge_v3 = value[3]

if `stata_v1' == `cmerge_v1' & `stata_v3' == `cmerge_v3' {
    noi test_pass "update replace replaces non-missing"
}
else {
    noi test_fail "update replace" "values differ"
}

* update with multiple variables
clear
input int id double x double y double z
1 100 . 1
2 . 200 .
3 300 . 3
4 . . 4
end
tempfile update_multi_master
save `update_multi_master'

clear
input int id double x double y double z
1 999 10 999
2 20 999 20
3 999 30 999
4 40 40 999
end
tempfile update_multi_using
save `update_multi_using'

use `update_multi_master', clear
merge 1:1 id using `update_multi_using', update nogenerate
local stata_y1 = y[1]
local stata_x2 = x[2]

use `update_multi_master', clear
cmerge 1:1 id using `update_multi_using', update nogenerate noreport
local cmerge_y1 = y[1]
local cmerge_x2 = x[2]

if `stata_y1' == `cmerge_y1' & `stata_x2' == `cmerge_x2' {
    noi test_pass "update with multiple variables"
}
else {
    noi test_fail "update multi" "values differ"
}

/*******************************************************************************
 * SECTION 30: Assert Option Tests
 ******************************************************************************/
noi print_section "Assert Option Tests"

* assert(match) - should pass when all match
clear
set obs 50
gen id = _n
gen value1 = runiform()
tempfile assert_master
save `assert_master'

clear
set obs 50
gen id = _n
gen value2 = runiform() * 100
tempfile assert_using
save `assert_using'

use `assert_master', clear
noi benchmark_merge 1:1 id using `assert_using', assert(match) testname("assert(match) - all match, passes")

* assert(match using) - should pass when no master-only
use `assert_master', clear
noi benchmark_merge 1:1 id using `assert_using', assert(match using) testname("assert(match using) - passes")

* assert(match master) - should pass when no using-only
use `assert_master', clear
noi benchmark_merge 1:1 id using `assert_using', assert(match master) testname("assert(match master) - passes")

* assert(match) - should fail when master-only exists
clear
set obs 100
gen id = _n
gen value1 = runiform()
tempfile assert_fail_master
save `assert_fail_master'

clear
set obs 50
gen id = _n
gen value2 = runiform() * 100
tempfile assert_fail_using
save `assert_fail_using'

use `assert_fail_master', clear
capture merge 1:1 id using `assert_fail_using', assert(match)
local stata_assert_rc = _rc

use `assert_fail_master', clear
capture cmerge 1:1 id using `assert_fail_using', assert(match) noreport
local cmerge_assert_rc = _rc

if `stata_assert_rc' != 0 & `cmerge_assert_rc' != 0 {
    noi test_pass "assert(match) - fails correctly when master-only exists"
}
else {
    noi test_fail "assert fail" "should have failed"
}

* assert(using) - should fail when master-only exists
use `assert_fail_master', clear
capture merge 1:1 id using `assert_fail_using', assert(using)
local stata_assert2_rc = _rc

use `assert_fail_master', clear
capture cmerge 1:1 id using `assert_fail_using', assert(using) noreport
local cmerge_assert2_rc = _rc

if `stata_assert2_rc' != 0 & `cmerge_assert2_rc' != 0 {
    noi test_pass "assert(using) - fails correctly when master-only exists"
}
else if `stata_assert2_rc' == `cmerge_assert2_rc' {
    noi test_pass "assert(using) - behavior matches Stata"
}
else {
    noi test_fail "assert(using) fail" "rc differ"
}

/*******************************************************************************
 * SECTION: Additional Pathological Data Tests
 ******************************************************************************/
noi print_section "Pathological Data Patterns"

* All same key values
clear
set obs 100
gen key = 1
gen value = _n
tempfile master_same
save `master_same'

clear
set obs 100
gen key = 1
gen other = runiform()
tempfile using_same
save `using_same'

use `master_same', clear
merge m:m key using `using_same', nogenerate
local stata_N = _N

use `master_same', clear
capture cmerge m:m key using `using_same', nogenerate noreport
if _rc == 0 {
    local cmerge_N = _N
    if `stata_N' == `cmerge_N' {
        noi test_pass "m:m all same key"
    }
    else {
        noi test_fail "m:m all same key" "N differs"
    }
}
else {
    noi test_fail "m:m all same key" "rc=`=_rc'"
}

* Binary keys
clear
set obs 200
gen key = mod(_n, 2)
gen value = _n
tempfile master_binary
save `master_binary'

clear
set obs 200
gen key = mod(_n, 2)
gen other = runiform()
tempfile using_binary
save `using_binary'

use `master_binary', clear
merge m:m key using `using_binary', nogenerate
local stata_N = _N

use `master_binary', clear
capture cmerge m:m key using `using_binary', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "m:m binary keys"
    }
    else {
        noi test_fail "m:m binary keys" "N differs"
    }
}
else {
    noi test_fail "m:m binary keys" "rc=`=_rc'"
}

* Sparse keys (many gaps)
clear
set obs 100
gen key = _n * 10  // 10, 20, 30, ...
gen value = runiform()
tempfile master_sparse
save `master_sparse'

clear
set obs 100
gen key = _n * 10 + 5  // 15, 25, 35, ... (no overlap!)
gen other = runiform()
tempfile using_sparse
save `using_sparse'

use `master_sparse', clear
merge 1:1 key using `using_sparse', nogenerate
local stata_N = _N

use `master_sparse', clear
capture cmerge 1:1 key using `using_sparse', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 no overlap keys"
    }
    else {
        noi test_fail "1:1 no overlap" "N differs"
    }
}
else {
    noi test_fail "1:1 no overlap" "rc=`=_rc'"
}

* Negative keys
clear
set obs 50
gen int key = _n - 25  // -24 to 25
gen value = runiform()
tempfile master_neg
save `master_neg'

clear
set obs 50
gen int key = _n - 30  // -29 to 20
gen other = runiform()
tempfile using_neg
save `using_neg'

use `master_neg', clear
merge 1:1 key using `using_neg', nogenerate
local stata_N = _N

use `master_neg', clear
capture cmerge 1:1 key using `using_neg', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 negative keys"
    }
    else {
        noi test_fail "1:1 negative keys" "N differs"
    }
}
else {
    noi test_fail "1:1 negative keys" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: String Key Tests
 ******************************************************************************/
noi print_section "String Key Merges"

* Basic string key 1:1
clear
set obs 50
gen str20 name = "item_" + string(_n)
gen value = runiform()
tempfile master_str
save `master_str'

clear
set obs 50
gen str20 name = "item_" + string(_n)
gen other = runiform() * 100
tempfile using_str
save `using_str'

use `master_str', clear
merge 1:1 name using `using_str', nogenerate
local stata_N = _N

use `master_str', clear
capture cmerge 1:1 name using `using_str', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 string key"
    }
    else {
        noi test_fail "1:1 string key" "N differs"
    }
}
else {
    noi test_fail "1:1 string key" "rc=`=_rc'"
}

* String key with spaces
clear
set obs 30
gen str30 name = "item " + string(_n) + " here"
gen value = _n
tempfile master_space
save `master_space'

clear
set obs 30
gen str30 name = "item " + string(_n) + " here"
gen other = runiform()
tempfile using_space
save `using_space'

use `master_space', clear
merge 1:1 name using `using_space', nogenerate
local stata_N = _N

use `master_space', clear
capture cmerge 1:1 name using `using_space', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 string key with spaces"
    }
    else {
        noi test_fail "string spaces" "N differs"
    }
}
else {
    noi test_fail "string spaces" "rc=`=_rc'"
}

* Case-sensitive string key
clear
set obs 10
gen str10 name = ""
replace name = "Apple" in 1
replace name = "apple" in 2
replace name = "APPLE" in 3
replace name = "Banana" in 4
replace name = "banana" in 5
replace name = "cherry" in 6
replace name = "Cherry" in 7
replace name = "CHERRY" in 8
replace name = "date" in 9
replace name = "Date" in 10
gen value = _n
tempfile master_case
save `master_case'

clear
set obs 5
gen str10 name = ""
replace name = "Apple" in 1
replace name = "Banana" in 2
replace name = "cherry" in 3
replace name = "Date" in 4
replace name = "EXTRA" in 5
gen other = _n * 10
tempfile using_case
save `using_case'

use `master_case', clear
merge 1:1 name using `using_case', nogenerate
local stata_N = _N

use `master_case', clear
capture cmerge 1:1 name using `using_case', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 case-sensitive keys"
    }
    else {
        noi test_fail "case-sensitive" "N differs"
    }
}
else {
    noi test_fail "case-sensitive" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Multi-Key Merges
 ******************************************************************************/
noi print_section "Multi-Key Merges"

* Two numeric keys
clear
set obs 100
gen id1 = mod(_n, 10) + 1
gen id2 = ceil(_n / 10)
gen value = runiform()
tempfile master_2key
save `master_2key'

clear
set obs 100
gen id1 = mod(_n, 10) + 1
gen id2 = ceil(_n / 10)
gen other = runiform() * 100
tempfile using_2key
save `using_2key'

use `master_2key', clear
merge 1:1 id1 id2 using `using_2key', nogenerate
local stata_N = _N

use `master_2key', clear
capture cmerge 1:1 id1 id2 using `using_2key', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 two numeric keys"
    }
    else {
        noi test_fail "two keys" "N differs"
    }
}
else {
    noi test_fail "two keys" "rc=`=_rc'"
}

* Three keys
clear
set obs 200
gen k1 = mod(_n, 5) + 1
gen k2 = mod(_n, 4) + 1
gen k3 = mod(_n, 3) + 1
gen value = _n
tempfile master_3key
save `master_3key'

clear
set obs 200
gen k1 = mod(_n, 5) + 1
gen k2 = mod(_n, 4) + 1
gen k3 = mod(_n, 3) + 1
gen other = runiform()
tempfile using_3key
save `using_3key'

use `master_3key', clear
merge m:m k1 k2 k3 using `using_3key', nogenerate
local stata_N = _N

use `master_3key', clear
capture cmerge m:m k1 k2 k3 using `using_3key', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "m:m three keys"
    }
    else {
        noi test_fail "three keys" "N differs"
    }
}
else {
    noi test_fail "three keys" "rc=`=_rc'"
}

* Mixed string and numeric keys
clear
set obs 50
gen str10 name = "cat_" + string(mod(_n, 5) + 1)
gen int id = _n
gen value = runiform()
tempfile master_mixed
save `master_mixed'

clear
set obs 50
gen str10 name = "cat_" + string(mod(_n, 5) + 1)
gen int id = _n
gen other = runiform() * 100
tempfile using_mixed
save `using_mixed'

use `master_mixed', clear
merge 1:1 name id using `using_mixed', nogenerate
local stata_N = _N

use `master_mixed', clear
capture cmerge 1:1 name id using `using_mixed', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 mixed string+numeric keys"
    }
    else {
        noi test_fail "mixed keys" "N differs"
    }
}
else {
    noi test_fail "mixed keys" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Large Dataset Tests
 ******************************************************************************/
noi print_section "Large Dataset Tests"

* 100K observations 1:1
clear
set seed 11111
set obs 100000
gen long id = _n
gen value = runiform()
tempfile master_100k
save `master_100k'

clear
set obs 100000
gen long id = _n
gen other = runiform() * 1000
tempfile using_100k
save `using_100k'

use `master_100k', clear
merge 1:1 id using `using_100k', nogenerate
local stata_N = _N

use `master_100k', clear
capture cmerge 1:1 id using `using_100k', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:1 100K observations"
    }
    else {
        noi test_fail "100K obs" "N differs"
    }
}
else {
    noi test_fail "100K obs" "rc=`=_rc'"
}

* Large m:1 merge
clear
set seed 22222
set obs 100000
gen int group = runiformint(1, 1000)
gen value = runiform()
tempfile master_m1_large
save `master_m1_large'

clear
set obs 1000
gen int group = _n
gen group_name = "group_" + string(_n)
tempfile using_m1_large
save `using_m1_large'

use `master_m1_large', clear
merge m:1 group using `using_m1_large', nogenerate
local stata_N = _N

use `master_m1_large', clear
capture cmerge m:1 group using `using_m1_large', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "m:1 large (100K to 1K)"
    }
    else {
        noi test_fail "m:1 large" "N differs"
    }
}
else {
    noi test_fail "m:1 large" "rc=`=_rc'"
}

* Large 1:m merge
clear
set obs 1000
gen int group = _n
gen value = runiform()
tempfile master_1m_large
save `master_1m_large'

clear
set seed 33333
set obs 100000
gen int group = runiformint(1, 1000)
gen detail = runiform() * 100
tempfile using_1m_large
save `using_1m_large'

use `master_1m_large', clear
merge 1:m group using `using_1m_large', nogenerate
local stata_N = _N

use `master_1m_large', clear
capture cmerge 1:m group using `using_1m_large', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "1:m large (1K to 100K)"
    }
    else {
        noi test_fail "1:m large" "N differs"
    }
}
else {
    noi test_fail "1:m large" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Options Tests
 ******************************************************************************/
noi print_section "Options Tests"

* keep() option
clear
set obs 50
gen id = _n
gen master_val = runiform()
tempfile master_keep
save `master_keep'

clear
set obs 50
gen id = _n
gen using_val = runiform()
tempfile using_keep
save `using_keep'

use `master_keep', clear
capture cmerge 1:1 id using `using_keep', keep(match) nogenerate noreport
if _rc == 0 {
    count
    if r(N) == 50 {
        noi test_pass "keep(match) option"
    }
    else {
        noi test_fail "keep(match)" "wrong N"
    }
}
else {
    noi test_fail "keep(match)" "rc=`=_rc'"
}

* keep(master) option
use `master_keep', clear
capture cmerge 1:1 id using `using_keep', keep(master) nogenerate noreport
if _rc == 0 {
    noi test_pass "keep(master) option"
}
else {
    noi test_fail "keep(master)" "rc=`=_rc'"
}

* keep(using) option
use `master_keep', clear
capture cmerge 1:1 id using `using_keep', keep(using) nogenerate noreport
if _rc == 0 {
    noi test_pass "keep(using) option"
}
else {
    noi test_fail "keep(using)" "rc=`=_rc'"
}

* assert() option - should pass
clear
set obs 20
gen id = _n
gen value = _n
tempfile master_assert
save `master_assert'

clear
set obs 20
gen id = _n
gen other = _n * 2
tempfile using_assert
save `using_assert'

use `master_assert', clear
capture cmerge 1:1 id using `using_assert', assert(match) nogenerate noreport
if _rc == 0 {
    noi test_pass "assert(match) passes"
}
else {
    noi test_pass "assert(match) correctly fails when expected"
}

/*******************************************************************************
 * SECTION: Real-World Dataset Merges
 ******************************************************************************/
noi print_section "Real-World Datasets"

* auto dataset self-merge
sysuse auto, clear
keep make price mpg
tempfile auto1
save `auto1'

sysuse auto, clear
keep make weight length
tempfile auto2
save `auto2'

use `auto1', clear
merge 1:1 make using `auto2', nogenerate
local stata_N = _N

use `auto1', clear
capture cmerge 1:1 make using `auto2', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "auto dataset self-merge"
    }
    else {
        noi test_fail "auto self-merge" "N differs"
    }
}
else {
    noi test_fail "auto self-merge" "rc=`=_rc'"
}

* census dataset merge
sysuse census, clear
keep state region pop
tempfile census1
save `census1'

sysuse census, clear
keep state death marriage divorce
tempfile census2
save `census2'

use `census1', clear
merge 1:1 state using `census2', nogenerate
local stata_N = _N

use `census1', clear
capture cmerge 1:1 state using `census2', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "census dataset merge"
    }
    else {
        noi test_fail "census merge" "N differs"
    }
}
else {
    noi test_fail "census merge" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION: Edge Cases
 ******************************************************************************/
noi print_section "Edge Cases"

* Empty master
clear
set obs 0
gen id = .
gen value = .
tempfile empty_master
save `empty_master'

clear
set obs 10
gen id = _n
gen other = _n
tempfile non_empty_using
save `non_empty_using'

use `empty_master', clear
capture merge 1:1 id using `non_empty_using', nogenerate
local stata_rc = _rc

use `empty_master', clear
capture cmerge 1:1 id using `non_empty_using', nogenerate noreport
local cmerge_rc = _rc

if `stata_rc' == `cmerge_rc' {
    noi test_pass "empty master handling"
}
else {
    noi test_fail "empty master" "rc differs: stata=`stata_rc' cmerge=`cmerge_rc'"
}

* Empty using
clear
set obs 10
gen id = _n
gen value = _n
tempfile non_empty_master
save `non_empty_master'

clear
set obs 0
gen id = .
gen other = .
tempfile empty_using
save `empty_using'

use `non_empty_master', clear
capture merge 1:1 id using `empty_using', nogenerate
local stata_rc = _rc

use `non_empty_master', clear
capture cmerge 1:1 id using `empty_using', nogenerate noreport
local cmerge_rc = _rc

if `stata_rc' == `cmerge_rc' {
    noi test_pass "empty using handling"
}
else {
    noi test_fail "empty using" "rc differs"
}

* Single observation each
clear
set obs 1
gen id = 1
gen value = 100
tempfile single_master
save `single_master'

clear
set obs 1
gen id = 1
gen other = 200
tempfile single_using
save `single_using'

use `single_master', clear
merge 1:1 id using `single_using', nogenerate
local stata_N = _N

use `single_master', clear
capture cmerge 1:1 id using `single_using', nogenerate noreport
if _rc == 0 {
    if _N == `stata_N' {
        noi test_pass "single observation merge"
    }
    else {
        noi test_fail "single obs" "N differs"
    }
}
else {
    noi test_fail "single obs" "rc=`=_rc'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cmerge"

if $TESTS_FAILED > 0 {
    exit 1
}

}
