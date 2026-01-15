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

use `master_mm', clear
noi benchmark_merge m:m group using `using_mm', testname("m:m basic")

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
 * SUMMARY
 ******************************************************************************/
noi print_summary "cmerge"

if $TESTS_FAILED > 0 {
    exit 1
}

}
