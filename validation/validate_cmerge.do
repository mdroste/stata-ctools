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
 * Basic 1:1 merge tests
 ******************************************************************************/
print_section "1:1 merge tests"

* Auto dataset - numeric key
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("1:1 numeric key")

* Auto dataset - string key
sysuse auto, clear
keep make price mpg
tempfile master
quietly save `master'

sysuse auto, clear
keep make weight length
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 make using `using_1', testname("1:1 string key (make)")

/*******************************************************************************
 * m:1 merge tests
 ******************************************************************************/
print_section "m:1 merge tests"

* Auto foreign lookup
sysuse auto, clear
keep make price mpg foreign
tempfile master
quietly save `master'

clear
input byte foreign str20 country
0 "Domestic"
1 "Foreign"
end
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge m:1 foreign using `using_1', testname("m:1 foreign lookup")

* Census region lookup
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge m:1 region using `using_1', testname("m:1 census region")

/*******************************************************************************
 * 1:m merge tests
 ******************************************************************************/
print_section "1:m merge tests"

* Region to states
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:m region using `using_1', testname("1:m region to states")

/*******************************************************************************
 * keep() option tests
 ******************************************************************************/
print_section "keep() option tests"

sysuse auto, clear
keep if _n <= 50
gen id = _n
tempfile master
quietly save `master'

sysuse auto, clear
keep if _n >= 30
gen id = _n + 29
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', keep(match) testname("keep(match)")

use `master', clear
benchmark_merge 1:1 id using `using_1', keep(master) testname("keep(master)")

use `master', clear
benchmark_merge 1:1 id using `using_1', keep(using) testname("keep(using)")

use `master', clear
benchmark_merge 1:1 id using `using_1', keep(master match) testname("keep(master match)")

/*******************************************************************************
 * generate/nogenerate options
 ******************************************************************************/
print_section "generate/nogenerate options"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight length
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', nogenerate testname("nogenerate")

use `master', clear
benchmark_merge 1:1 id using `using_1', generate(merge_result) testname("generate(merge_result)")

/*******************************************************************************
 * keepusing() option
 ******************************************************************************/
print_section "keepusing() option"

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight length turn
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', keepusing(weight) testname("keepusing(weight)")

/*******************************************************************************
 * Multiple key variables
 ******************************************************************************/
print_section "Multiple key variables"

* Two keys
webuse nlswork, clear
keep in 1/2000
keep idcode year ln_wage age
tempfile master
quietly save `master'

webuse nlswork, clear
keep in 1/2000
keep idcode year tenure ttl_exp
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 idcode year using `using_1', testname("two keys (idcode year)")

* Three keys
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id1 id2 id3 using `using_1', testname("three keys")

* Mixed string and numeric keys
sysuse census, clear
keep state region pop
tempfile master
quietly save `master'

sysuse census, clear
keep state region medage
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 state region using `using_1', testname("mixed keys (string + numeric)")

/*******************************************************************************
 * Edge cases
 ******************************************************************************/
print_section "Edge cases"

* No matches
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("no matches")

* All matches
sysuse auto, clear
gen id = _n
keep id make price
tempfile master
quietly save `master'

sysuse auto, clear
gen id = _n
keep id weight
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("all matches")

* Missing key values
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("missing key values")

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
quietly save `master'

clear
input int id float other
-5 100
-1 200
0 300
1 400
5 500
end
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("negative keys")

* Float key values
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("float keys")

/*******************************************************************************
 * Large dataset tests
 ******************************************************************************/
print_section "Large dataset tests"

* 50K 1:1 merge
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("50K 1:1 full match")

* 100K m:1 merge
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge m:1 group_id using `using_1', testname("100K m:1")

* Partial overlap
clear
set seed 22222
set obs 30000
gen id = _n
gen value1 = runiform()
tempfile master
quietly save `master'

clear
set obs 30000
gen id = _n + 15000
gen value2 = runiform()
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("partial overlap (50%)")

/*******************************************************************************
 * String key tests
 ******************************************************************************/
print_section "String key tests"

* Long string keys
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
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 long_key using `using_1', testname("long string keys")

* Census state name
sysuse census, clear
keep state pop
tempfile master
quietly save `master'

sysuse census, clear
keep state medage death marriage divorce
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 state using `using_1', testname("state name 1:1")

/*******************************************************************************
 * Panel data tests
 ******************************************************************************/
print_section "Panel data (nlswork)"

webuse nlswork, clear
keep in 1/5000
keep idcode year ln_wage hours
tempfile master
quietly save `master'

webuse nlswork, clear
keep in 1/5000
keep idcode year age tenure wks_work
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 idcode year using `using_1', testname("panel 1:1")

* m:1 individual characteristics
webuse nlswork, clear
keep in 1/3000
keep idcode year ln_wage
tempfile master
quietly save `master'

webuse nlswork, clear
keep in 1/3000
bysort idcode (year): keep if _n == 1
keep idcode race grade
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge m:1 idcode using `using_1', testname("panel m:1 individual")

/*******************************************************************************
 * Other built-in datasets
 ******************************************************************************/
print_section "Other datasets"

* Voter dataset
sysuse voter, clear
keep pop frac
gen id = _n
tempfile master
quietly save `master'

sysuse voter, clear
keep cand inc
gen id = _n
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 id using `using_1', testname("voter 1:1")

* US life expectancy
sysuse uslifeexp, clear
keep year le
tempfile master
quietly save `master'

sysuse uslifeexp, clear
keep year le_male le_female
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 year using `using_1', testname("uslifeexp 1:1")

* SP500
sysuse sp500, clear
keep date open high
tempfile master
quietly save `master'

sysuse sp500, clear
keep date low close volume
tempfile using_1
quietly save `using_1'

use `master', clear
benchmark_merge 1:1 date using `using_1', testname("sp500 1:1")

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cmerge"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
