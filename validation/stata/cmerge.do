* cmerge: named cases share this file, each do-call has its own local scope.
gettoken __case 0 : 0
macro shift
if "`__case'" == "" local __case "cmerge"

if "`__case'" == "cmerge" {
* validate_cmerge.do
* Comprehensive validation tests for cmerge vs native Stata merge
* Tests all merge types and options
* Merge types: 1:1, m:1, 1:m, m:m
* Options: keep(), generate(), nogenerate, keepusing(), sorted, force,
* noreport, nolabel, nonotes, update, replace, assert()

* Load setup (works from project root or validation dir)
do "validation/support/validate_setup.do"

quietly {

noi di as text "Running validation tests for cmerge..."

* Basic 1:1 merge tests

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
benchmark_merge 1:1 id using `using_1', testname("1:1 numeric key")

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
benchmark_merge 1:1 make using `using_2', testname("1:1 string key")

* SECTION 2: m:1 merge tests
print_section "m:1 Merge Tests"

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
benchmark_merge m:1 foreign using `using_3', testname("m:1 foreign lookup")

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
benchmark_merge m:1 region using `using_4', testname("m:1 census region")

* SECTION 3: 1:m merge tests
print_section "1:m Merge Tests"

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
benchmark_merge 1:m region using `using_5', testname("1:m region to states")

* SECTION 4: m:m merge tests
* An m:m merge pairs the j-th master row of a key group with the j-th using
* row (the shorter group's last row repeats), so its result depends on the
* order of rows within each key group. Stata's merge gives an unsorted master
* or using file a plain -sort- by the keys, which orders tied keys at random
* (it draws on the sortseed state, which every sort advances): two native
* m:m merges of the same unsorted files in one session pair different rows.
* The m:m benchmarks below therefore stably sort both files by the keys
* first, so that native merge pairs rows in file order and cmerge must match
* it exactly. benchmark_mm_unsorted covers unsorted files: cmerge must pair
* their rows in file order, i.e. match native merge of the stably sorted
* files.
print_section "m:m Merge Tests"

capture program drop benchmark_mm_unsorted
program define benchmark_mm_unsorted
    syntax varlist using/, testname(string)

    preserve
    tempfile original sorted_using reference
    quietly save `original'

    quietly use `"`using'"', clear
    sort `varlist', stable
    quietly save `sorted_using'
    quietly use `original', clear
    sort `varlist', stable
    quietly merge m:m `varlist' using `sorted_using'
    quietly ds
    local allvars `r(varlist)'
    sort `allvars'
    quietly save `reference'

    quietly use `original', clear
    local sorted_before : sortedby
    capture quietly cmerge m:m `varlist' using `"`using'"', noreport
    local rc = _rc
    local cf_rc = 0
    if `rc' == 0 {
        quietly ds
        local allvars `r(varlist)'
        sort `allvars'
        capture quietly cf _all using `reference'
        local cf_rc = _rc
    }
    restore

    if "`sorted_before'" != "" {
        test_fail "`testname'" "master is already sorted by `sorted_before'; test needs unsorted files"
    }
    else if `rc' != 0 {
        test_fail "`testname'" "cmerge returned error `rc'"
    }
    else if `cf_rc' != 0 {
        test_fail "`testname'" "merged data differs from Stata merge of the stably sorted files"
    }
    else {
        test_pass "`testname'"
    }
end

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
benchmark_mm_unsorted group using `using_mm', testname("m:m basic: unsorted files pair in file order")

use `using_mm', clear
sort group, stable
save `using_mm', replace
use `master_mm', clear
sort group, stable
benchmark_merge m:m group using `using_mm', testname("m:m basic")

* SECTION 5: keep() option tests
print_section "keep() Option Tests"

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
benchmark_merge 1:1 id using `using_6', keep(match) testname("keep(match)")

use `master_6', clear
benchmark_merge 1:1 id using `using_6', keep(master) testname("keep(master)")

use `master_6', clear
benchmark_merge 1:1 id using `using_6', keep(using) testname("keep(using)")

use `master_6', clear
benchmark_merge 1:1 id using `using_6', keep(master match) testname("keep(master match)")

use `master_6', clear
benchmark_merge 1:1 id using `using_6', keep(master using) testname("keep(master using)")

use `master_6', clear
benchmark_merge 1:1 id using `using_6', keep(using match) testname("keep(using match)")

* SECTION 6: generate/nogenerate options
print_section "generate/nogenerate Options"

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
benchmark_merge 1:1 id using `using_1', nogenerate testname("nogenerate")

use `master', clear
benchmark_merge 1:1 id using `using_1', generate(merge_result) testname("generate(merge_result)")

* SECTION 7: keepusing() option
print_section "keepusing() Option"

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
benchmark_merge 1:1 id using `using_1', keepusing(weight) testname("keepusing single var")

use `master', clear
benchmark_merge 1:1 id using `using_1', keepusing(weight length) testname("keepusing two vars")

* SECTION 8: Multiple key variables
print_section "Multiple Key Variables"

* Two keys
ctools_fixture nlswork, clear
keep in 1/2000
keep idcode year ln_wage age
tempfile master
save `master'

ctools_fixture nlswork, clear
keep in 1/2000
keep idcode year tenure ttl_exp
tempfile using_1
save `using_1'

use `master', clear
benchmark_merge 1:1 idcode year using `using_1', testname("two keys: idcode year")

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
benchmark_merge 1:1 id1 id2 id3 using `using_1', testname("three keys")

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
benchmark_merge 1:1 state region using `using_1', testname("mixed keys: string+numeric")

* SECTION 9: Edge cases
print_section "Edge Cases"

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
benchmark_merge 1:1 id using `using_1', testname("no matches")

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
benchmark_merge 1:1 id using `using_1', testname("negative keys")

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
benchmark_merge 1:1 id using `using_1', testname("float keys")

* SECTION 10: Large dataset tests
print_section "Large Dataset Tests"

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
benchmark_merge 1:1 id using `using_1', testname("50K 1:1 full match")

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
benchmark_merge m:1 group_id using `using_1', testname("100K m:1")

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
benchmark_merge 1:1 id using `using_1', testname("partial overlap (50%)")

* SECTION 11: String key tests
print_section "String Key Tests"

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
benchmark_merge 1:1 long_key using `using_1', testname("long string keys")

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
benchmark_merge 1:1 state using `using_1', testname("state name 1:1")

* SECTION 12: Panel data tests
print_section "Panel Data (nlswork)"

ctools_fixture nlswork, clear
keep in 1/5000
keep idcode year ln_wage hours
tempfile master
save `master'

ctools_fixture nlswork, clear
keep in 1/5000
keep idcode year age tenure wks_work
tempfile using_1
save `using_1'

use `master', clear
benchmark_merge 1:1 idcode year using `using_1', testname("panel 1:1")

* m:1 individual characteristics
ctools_fixture nlswork, clear
keep in 1/3000
keep idcode year ln_wage
tempfile master
save `master'

ctools_fixture nlswork, clear
keep in 1/3000
bysort idcode (year): keep if _n == 1
keep idcode race grade
tempfile using_1
save `using_1'

use `master', clear
benchmark_merge m:1 idcode using `using_1', testname("panel m:1 individual")

* SECTION 13: assert() option tests
print_section "assert() Option Tests"

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
benchmark_merge 1:1 id using `using_1', assert(match) testname("assert(match) passes")

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

if `stata_rc' != 0 & `cmerge_rc' != 0 test_pass "assert(match) fails correctly"
else test_fail "assert(match) should fail" "stata_rc=`stata_rc' cmerge_rc=`cmerge_rc'"

* SECTION 14: update/replace options
print_section "update/replace Options"

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
benchmark_merge 1:1 id using `using_1', update nogenerate testname("update option")

* Test update replace option
use `master', clear
benchmark_merge 1:1 id using `using_1', update replace nogenerate testname("update replace option")

* SECTION 15: sorted option
print_section "sorted Option"

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
benchmark_merge 1:1 id using `using_1', sorted testname("sorted option")

* SECTION 16: force option
print_section "force Option"

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
benchmark_merge 1:1 id using `using_1', force nogenerate testname("force option")

* SECTION 17: nolabel and nonotes options
print_section "nolabel/nonotes Options"

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
benchmark_merge 1:1 id using `using_1', nolabel testname("nolabel option")

use `master', clear
benchmark_merge 1:1 id using `using_1', nonotes testname("nonotes option")

* SECTION 18: 1:1 _n merge tests (merge by observation number)
print_section "1:1 _n Merge Tests"

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
benchmark_merge 1:1 _n using `using_1', nogenerate testname("1:1 _n basic (same nobs)")

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
benchmark_merge 1:1 _n using `using_1', testname("1:1 _n master has more obs")

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
benchmark_merge 1:1 _n using `using_1', testname("1:1 _n using has more obs")

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
benchmark_merge 1:1 _n using `using_1', keepusing(weight) nogenerate testname("1:1 _n with keepusing")

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
benchmark_merge 1:1 _n using `using_1', keep(match) nogenerate testname("1:1 _n with keep(match)")

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
benchmark_merge 1:1 _n using `using_1', generate(merge_result) testname("1:1 _n with generate()")

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
benchmark_merge 1:1 _n using `using_1', nogenerate testname("1:1 _n large dataset (50K obs)")

* SECTION 19: Comprehensive sysuse/webuse Dataset Tests
print_section "Comprehensive Built-in Dataset Tests"

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
benchmark_merge m:1 foreign using `auto_lookup', testname("auto: m:1 foreign lookup")

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
benchmark_merge 1:1 make using `auto2', testname("auto: 1:1 by make")

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
benchmark_merge m:1 region using `region_lookup', testname("census: m:1 region lookup")

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
benchmark_merge 1:1 state using `census2', testname("census: 1:1 by state")

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
benchmark_merge m:1 industry using `industry_lookup', testname("nlsw88: m:1 by industry")

* nlswork panel data - 1:1 by idcode year
ctools_fixture nlswork, clear
keep in 1/3000
keep idcode year ln_wage hours
tempfile nlswork1
save `nlswork1'

ctools_fixture nlswork, clear
keep in 1/3000
keep idcode year age tenure wks_work
tempfile nlswork2
save `nlswork2'

use `nlswork1', clear
benchmark_merge 1:1 idcode year using `nlswork2', testname("nlswork: 1:1 panel by idcode year")

* nlswork - m:1 individual time-invariant characteristics
ctools_fixture nlswork, clear
keep in 1/2500
keep idcode year ln_wage
tempfile nlswork_panel
save `nlswork_panel'

ctools_fixture nlswork, clear
keep in 1/2500
bysort idcode (year): keep if _n == 1
keep idcode race grade
tempfile nlswork_indiv
save `nlswork_indiv'

use `nlswork_panel', clear
benchmark_merge m:1 idcode using `nlswork_indiv', testname("nlswork: m:1 individual characteristics")

* Grunfeld panel data
ctools_fixture grunfeld, clear
keep company year invest
tempfile grunfeld1
save `grunfeld1'

ctools_fixture grunfeld, clear
keep company year mvalue kstock
tempfile grunfeld2
save `grunfeld2'

use `grunfeld1', clear
benchmark_merge 1:1 company year using `grunfeld2', testname("grunfeld: 1:1 panel merge")

* lifeexp dataset
ctools_fixture lifeexp, clear
keep region country popgrowth
tempfile lifeexp1
save `lifeexp1'

ctools_fixture lifeexp, clear
keep region country lexp gnppc
tempfile lifeexp2
save `lifeexp2'

use `lifeexp1', clear
benchmark_merge 1:1 region country using `lifeexp2', testname("lifeexp: 1:1 by region country")

* SECTION 20: Pathological String Key Edge Cases
print_section "Pathological String Key Edge Cases"

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
benchmark_merge m:1 name using `str_comma_using', testname("string key: embedded commas")

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
benchmark_merge m:1 name using `str_quote_using', testname("string key: embedded quotes")

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
benchmark_merge m:1 key using `str_space_using', testname("string key: leading/trailing spaces")

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
benchmark_merge 1:1 long_key using `str_long_using', testname("string key: very long (200+ chars)")

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
benchmark_merge 1:1 name using `str_case_using', testname("string key: case sensitive")

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
benchmark_merge m:1 key using `str_empty_using', testname("string key: empty strings")

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
benchmark_merge m:1 name using `str_unicode_using', testname("string key: Latin extended chars")

* SECTION 21: Numeric Key Edge Cases
print_section "Numeric Key Edge Cases"

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
benchmark_merge m:1 id using `num_miss_using', testname("numeric key: missing values (m:1)")

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
benchmark_merge 1:1 id using `num_ext_using', testname("numeric key: extended missing (.a,.b,.z)")

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
benchmark_merge 1:1 id using `float_using', testname("numeric key: float precision")

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
benchmark_merge 1:1 id using `large_int_using', testname("numeric key: large integers")

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
benchmark_merge 1:1 id using `neg_cross_using', testname("numeric key: negative to positive range")

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
benchmark_merge m:m id using `zero_key_using', testname("numeric key: many zeros (m:m)")

* SECTION 22: Multiple Key Variables
print_section "Multiple Key Variables"

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
benchmark_merge 1:1 id1 id2 id3 id4 using `four_key_using', testname("four numeric keys")

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
benchmark_merge m:m cat year month day id using `five_key_using', testname("five mixed keys (string + numeric)")

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
benchmark_merge m:m key1 key2 key3 using `same_key_using', testname("all same keys (degenerate m:m)")

* SECTION 23: Match Pattern Edge Cases
print_section "Match Pattern Edge Cases"

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
benchmark_merge 1:1 id using `disjoint_using', testname("no matches (completely disjoint)")

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
benchmark_merge 1:1 id using `perfect_using', testname("all matches (perfect 1:1)")

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
benchmark_merge 1:1 id using `master_only_using', testname("master only (no using matches)")

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
benchmark_merge 1:1 id using `using_only_using', testname("using only (no master matches)")

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
benchmark_merge 1:1 id using `sparse_using', testname("sparse matches (~10% overlap)")

* SECTION 24: Missing Value Handling
print_section "Missing Value Handling"

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
benchmark_merge 1:1 id using `miss_nonkey_using', testname("missing in non-key variables")

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
benchmark_merge 1:1 id using `first_last_miss_using', testname("first/last rows missing keys")

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
benchmark_merge 1:1 id1 id2 using `all_miss_key_using', testname("all missing in second key column")

* SECTION 25: Option Combinations
print_section "Option Combinations"

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
benchmark_merge 1:1 id using `opt_using', keep(master match) testname("keep(master match)")

* keep(using match) combination
use `opt_master', clear
benchmark_merge 1:1 id using `opt_using', keep(using match) testname("keep(using match)")

* keep(match) only
use `opt_master', clear
benchmark_merge 1:1 id using `opt_using', keep(match) testname("keep(match) only")

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
benchmark_merge 1:1 id using `keepusing_using', keepusing(var_a var_c) testname("keepusing: select 2 of 4 vars")

* generate() with custom name
use `opt_master', clear
benchmark_merge 1:1 id using `opt_using', generate(my_merge_var) testname("generate(my_merge_var)")

* nogenerate option
use `opt_master', clear
benchmark_merge 1:1 id using `opt_using', nogenerate testname("nogenerate option")

* sorted option with pre-sorted data
use `opt_master', clear
sort id
save `opt_master', replace
use `opt_using', clear
sort id
save `opt_using', replace
use `opt_master', clear
benchmark_merge 1:1 id using `opt_using', sorted testname("sorted option (pre-sorted)")

* Combined options: keep + keepusing + nogenerate
use `keepusing_master', clear
benchmark_merge 1:1 id using `keepusing_using', keep(match) keepusing(var_a var_b) nogenerate testname("combined: keep + keepusing + nogen")

* SECTION 26: Large Dataset Tests
print_section "Large Dataset Tests"

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
benchmark_merge 1:1 id using `large_10k_using', testname("10K 1:1 full match")

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
benchmark_merge 1:1 id using `large_50k_using', testname("50K 1:1 partial match (50%)")

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
benchmark_merge m:1 group using `large_m1_using', testname("large m:1 (100K to 1K)")

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
benchmark_merge 1:m group using `large_1m_using', testname("large 1:m (1K to 100K)")

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
benchmark_mm_unsorted group using `large_mm_using', testname("large m:m with expansion: unsorted files pair in file order")

* Stable sorts make native merge pair in file order (see SECTION 4)
use `large_mm_using', clear
sort group, stable
save `large_mm_using', replace
use `large_mm_master', clear
sort group, stable
benchmark_merge m:m group using `large_mm_using', testname("large m:m with expansion")

* SECTION 27: Unbalanced Merges
print_section "Unbalanced Merges"

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
benchmark_merge 1:1 id using `unbal_small_using', testname("unbalanced: 50K master, 100 using")

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
benchmark_merge 1:1 id using `unbal_large_using', testname("unbalanced: 100 master, 50K using")

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
benchmark_mm_unsorted group using `unbal_10x_using', testname("unbalanced m:m: unsorted files pair in file order")

* Stable sorts make native merge pair in file order (see SECTION 4)
use `unbal_10x_using', clear
sort group, stable
save `unbal_10x_using', replace
use `unbal_10x_master', clear
sort group, stable
benchmark_merge m:m group using `unbal_10x_using', testname("unbalanced m:m: 10K master, 1K using")

* SECTION 28: Variable Conflicts
print_section "Variable Conflicts"

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
    test_pass "variable conflict detected (same var names)"
}
else if `stata_conflict_rc' == 0 & `cmerge_conflict_rc' == 0 {
    test_pass "variable conflict handled (both succeeded)"
}
else {
    test_fail "variable conflict" "rc differ: stata=`stata_conflict_rc' cmerge=`cmerge_conflict_rc'"
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

if `stata_force_rc' == `cmerge_force_rc' test_pass "type mismatch with force option"
else test_fail "type mismatch force" "rc differ: stata=`stata_force_rc' cmerge=`cmerge_force_rc'"

* SECTION 29: Update and Replace Options
print_section "Update and Replace Options"

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
benchmark_merge 1:1 id using `update_using', update nogenerate testname("update fills missing values")

* update replace - replace non-missing too
use `update_master', clear
benchmark_merge 1:1 id using `update_using', update replace nogenerate testname("update replace replaces non-missing")

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
benchmark_merge 1:1 id using `update_multi_using', update nogenerate testname("update multi vars")

* SECTION 30: Assert Option Tests
print_section "Assert Option Tests"

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
benchmark_merge 1:1 id using `assert_using', assert(match) testname("assert(match) - all match, passes")

* assert(match using) - should pass when no master-only
use `assert_master', clear
benchmark_merge 1:1 id using `assert_using', assert(match using) testname("assert(match using) - passes")

* assert(match master) - should pass when no using-only
use `assert_master', clear
benchmark_merge 1:1 id using `assert_using', assert(match master) testname("assert(match master) - passes")

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

if `stata_assert_rc' != 0 & `cmerge_assert_rc' != 0 test_pass "assert(match) - fails correctly when master-only exists"
else test_fail "assert fail" "should have failed"

* assert(using) - should fail when master-only exists
use `assert_fail_master', clear
capture merge 1:1 id using `assert_fail_using', assert(using)
local stata_assert2_rc = _rc

use `assert_fail_master', clear
capture cmerge 1:1 id using `assert_fail_using', assert(using) noreport
local cmerge_assert2_rc = _rc

if `stata_assert2_rc' != 0 & `cmerge_assert2_rc' != 0 {
    test_pass "assert(using) - fails correctly when master-only exists"
}
else if `stata_assert2_rc' == `cmerge_assert2_rc' {
    test_pass "assert(using) - behavior matches Stata"
}
else {
    test_fail "assert(using) fail" "rc differ"
}

* SECTION: Intentional Error Tests
* These tests verify that cmerge returns the same error codes as Stata's merge
* when given invalid inputs or error conditions.
print_section "Intentional Error Tests"

* Using file doesn't exist
sysuse auto, clear
test_error_match, stata_cmd(merge 1:1 make using "nonexistent_file.dta") ctools_cmd(cmerge 1:1 make using "nonexistent_file.dta", noreport) testname("using file doesn't exist")

* Key variable doesn't exist
sysuse auto, clear
tempfile using_file
save `using_file'
test_error_match, stata_cmd(merge 1:1 nonexistent_var using `using_file') ctools_cmd(cmerge 1:1 nonexistent_var using `using_file', noreport) testname("key variable doesn't exist")

* Invalid merge type
sysuse auto, clear
tempfile using_file
save `using_file'
test_error_match, stata_cmd(merge invalid_type make using `using_file') ctools_cmd(cmerge invalid_type make using `using_file', noreport) testname("invalid merge type")

* _merge variable already exists
sysuse auto, clear
gen _merge = 1
tempfile using_file
save `using_file'
test_error_match, stata_cmd(merge 1:1 make using `using_file') ctools_cmd(cmerge 1:1 make using `using_file', noreport) testname("_merge variable already exists")

* Empty master dataset
clear
set obs 0
gen id = .
tempfile using_file
clear
set obs 5
gen id = _n
save `using_file'
clear
set obs 0
gen id = .
capture merge 1:1 id using `using_file'
local stata_rc = _rc
clear
set obs 0
gen id = .
capture cmerge 1:1 id using `using_file', noreport
local cmerge_rc = _rc
if `stata_rc' == `cmerge_rc' test_pass "[error] empty master dataset (rc=`stata_rc')"
else test_fail "[error] empty master dataset" "stata rc=`stata_rc', cmerge rc=`cmerge_rc'"

* SECTION 31: Extended missing keys out of order
* Both sides list . .a .b out of Stata's order, so the sort must order them.
print_section "Extended Missing Keys (unsorted)"

clear
input double id str1 m
.b "x"
.a "y"
1  "z"
.  "w"
end
tempfile ext_master ext_using
save `ext_master'
clear
input double id str1 u
.a "p"
2  "r"
.b "q"
.  "s"
end
save `ext_using'
use `ext_master', clear
benchmark_merge 1:1 id using `ext_using', testname("numeric key: extended missing out of order (1:1)")

clear
input double id str1 m
.b "x"
.a "y"
.b "v"
1  "z"
.a "t"
.  "w"
end
benchmark_merge m:1 id using `ext_using', testname("numeric key: extended missing out of order (m:1)")

* SECTION 32: A failed or interrupted merge leaves the master as it was
* By default cmerge writes no snapshot of the master: it plans the join
* before changing anything, records the schema changes it makes, and the
* plugin keeps a copy of the master observations it overwrites (paths that
* cannot be undone that way keep an in-memory snapshot frame). Option
* snapshot preserves the whole dataset instead; option fast keeps nothing and
* only reports the failure. Each case runs in all three modes. The
* validation hook $CTOOLS_cmerge_fault fails a merge at one stage with r(1),
* as Break would:
* record  after the record, before any change to the master
* schema  after the recasts, new variables, labels, notes and new rows,
* just before the plugin writes
* store   inside the plugin, after half of the output rows were written
* (reported as the memory error it simulates, r(920))
* post    after the plugin wrote the whole merge (a Break pressed during
* the plugin is acted on there)
* keep    after keep() dropped observations
* Join errors (r(459), r(5)) fail in the plugin's plan, before any change.
* In the default and snapshot modes each case compares the data in memory
* afterwards with a copy saved before:
* values (datasignature and cf), variable order, storage types, formats,
* variable and value labels, value-label definitions, notes and other
* characteristics, the sort marker and c(changed). Characteristics are compared
* by name: Stata's own restore does not keep their internal order either.
print_section "Rollback of Failed Merges"

capture program drop _cmerge_state
program define _cmerge_state, rclass
    quietly datasignature
    local s `"N=`=_N' k=`=c(k)' sig=`r(datasignature)' sortedby=`: sortedby' changed=`c(changed)'"'
    local s `"`macval(s)' label=`: data label'"'
    local chars : char _dta[]
    local chars : list sort chars
    foreach c of local chars {
        local value : char _dta[`c']
        local s `"`macval(s)' | _dta[`c']=`macval(value)'"'
    }
    if c(k) > 0 {
        foreach v of varlist _all {
            local s `"`macval(s)' || `v' `: type `v'' `: format `v'' [`: value label `v''] [`: variable label `v'']"'
            local chars : char `v'[]
            local chars : list sort chars
            foreach c of local chars {
                local value : char `v'[`c']
                local s `"`macval(s)' | `v'[`c']=`macval(value)'"'
            }
        }
    }
    quietly label dir
    local names `r(names)'
    local s `"`macval(s)' || labels: `names'"'
    if "`names'" != "" {
        tempfile defs
        tempname fh
        quietly label save `names' using `"`defs'"', replace
        file open `fh' using `"`defs'"', read text
        file read `fh' line
        while r(eof) == 0 {
            local s `"`macval(s)' | `macval(line)'"'
            file read `fh' line
        }
        file close `fh'
    }
    return local state `"`macval(s)'"'
end

* Run a cmerge command (whose options end the command) that must fail with
* return code rc (by default the code of the injected fault), in the default,
* snapshot and fast modes. Default and snapshot must leave the data saved in
* copy, with no record or snapshot frame left behind; fast only reports.
capture program drop _cmerge_rollback_case
program define _cmerge_rollback_case
    syntax , cmd(string asis) testname(string) [fault(string) rc(integer -1)]
    if `rc' == -1 local rc = cond("`fault'" == "store", 920, 1)
    tempfile copy
    quietly save `copy'
    _cmerge_state
    local before `"`r(state)'"'
    quietly frames dir
    local frames `"`r(frames)'"'
    foreach mode in default snapshot fast {
        quietly use `copy', clear
        local option = cond("`mode'" == "default", "", "`mode'")
        global CTOOLS_cmerge_fault `fault'
        capture `cmd' `option'
        local got = _rc
        global CTOOLS_cmerge_fault
        local name "`testname' [`mode']"
        if `got' != `rc' {
            test_fail "`name'" "expected r(`rc'), got r(`got')"
            continue
        }
        if "`mode'" == "fast" {
            test_pass "`name'"
            continue
        }
        _cmerge_state
        local after `"`r(state)'"'
        mata: st_local("same", strofreal(st_local("before") == st_local("after")))
        capture cf _all using `copy', all
        local cf_rc = _rc
        quietly frames dir
        if !`same' | `cf_rc' {
            test_fail "`name'" "master changed (state same=`same', cf r(`cf_rc'))"
        }
        else if "$CTOOLS_cmu_active$CTOOLS_cmerge_plugin$CTOOLS_cmerge_snapframe" != "" | ///
            `"`r(frames)'"' != `"`frames'"' {
            test_fail "`name'" "rollback record or snapshot frame left behind"
        }
        else {
            test_pass "`name'"
        }
    }
    quietly use `copy', clear
end

* Master: an unsorted key with duplicates (grp), a unique key (id), shared
* variables the merge promotes (grp byte->int, x float->double, s str3->str10
* with a left-aligned format), a value-labeled variable, notes and
* characteristics on variables and _dta, a data label and a sort marker.
clear
set seed 4711
set obs 3000
gen long id = 3001 - _n
gen byte grp = mod(_n * 7, 23) + 1
gen float x = cond(mod(_n, 5) == 0, ., _n / 3)
gen str3 s = cond(mod(_n, 4) == 0, "", substr("abcdef", mod(_n, 4) + 1, 3))
format s %-5s
gen byte lab = mod(_n, 3)
label define mlab 0 "zero" 1 "one" 2 "two"
label values lab mlab
gen double w = runiform()
label variable w "master weight"
label data "Master data"
note: master dataset note
note x: master note on x
char _dta[mchar] master characteristic
char w[wchar] w characteristic
sort lab w
tempfile rb_master
save `rb_master'

* Using: one row per grp 1..25 (two not in the master), the shared x and s at
* wider types, a new labeled variable, notes on a shared and a new variable,
* and _dta notes and characteristics.
clear
set obs 25
gen int grp = _n
gen double x = cond(mod(_n, 3) == 0, ., _n * 100 + .25)
gen str10 s = "u" + string(_n)
gen byte newv = mod(_n, 2)
label define ulab 0 "even" 1 "odd"
label values newv ulab
gen str5 newstr = "n" + string(_n)
label data "Using data"
note: using dataset note
note x: using note on x
note newv: using note on newv
char _dta[uchar] using characteristic
tempfile rb_using
save `rb_using'

* 1:1 using on id and a using with duplicate grp values
clear
set obs 3100
gen long id = _n
gen double y = _n * 2
tempfile rb_using_id
save `rb_using_id'
use `rb_using', clear
expand 2 in 1/3
tempfile rb_using_dup
save `rb_using_dup'

* Staged path (the master is not in key order): every stage
foreach f in record schema store post keep {
    use `rb_master', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_using', keep(master match using) noreport) ///
        testname("rollback m:1 staged: fault at `f'")
}
foreach f in schema store post {
    use `rb_master', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_using', update replace noreport) ///
        testname("rollback m:1 update replace: fault at `f'")
    use `rb_master', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge 1:1 id using `rb_using_id', generate(how) noreport) ///
        testname("rollback 1:1 staged: fault at `f'")
    use `rb_master', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge 1:1 _n using `rb_using_id', noreport) ///
        testname("rollback 1:1 _n: fault at `f'")
}

* Join errors are found by the plugin's plan, before any change
use `rb_master', clear
expand 2 in 1/5
sort lab w
_cmerge_rollback_case, rc(459) cmd(cmerge 1:1 id using `rb_using_id', noreport) ///
    testname("rollback: duplicate master keys r(459)")
use `rb_master', clear
_cmerge_rollback_case, rc(459) cmd(cmerge m:1 grp using `rb_using_dup', noreport) ///
    testname("rollback: duplicate using keys r(459)")
use `rb_master', clear
_cmerge_rollback_case, rc(5) cmd(cmerge m:1 grp using `rb_using', sorted noreport) ///
    testname("rollback: master not sorted r(5)")

* Identity output (master in key order, every key matched): the plugin
* rewrites no master observation unless keep() or update needs a copy
use `rb_master', clear
drop if grp > 20
sort grp id
tempfile rb_sorted
save `rb_sorted'
use `rb_using', clear
keep if grp <= 20
tempfile rb_using20
save `rb_using20'
foreach f in schema post {
    use `rb_sorted', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_using20', noreport) ///
        testname("rollback m:1 identity: fault at `f'")
}
foreach f in post keep {
    use `rb_sorted', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_using20', keep(match) update noreport) ///
        testname("rollback m:1 identity keep update: fault at `f'")
}

* Strings too wide to stage: the merge keeps a snapshot instead
use `rb_master', clear
gen str2045 wide = "w" + string(_n)
sort lab w
tempfile rb_wide
save `rb_wide'
foreach f in record schema store post {
    use `rb_wide', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_using', noreport) ///
        testname("rollback m:1 unstaged strings: fault at `f'")
}

* Empty data (merged without the plugin's join): a 1:1 merge that sorts the
* master keeps a snapshot frame; the others are undone from the record
use `rb_using_id', clear
drop if 1
tempfile rb_empty
save `rb_empty'
use `rb_master', clear
_cmerge_rollback_case, fault(post) ///
    cmd(cmerge 1:1 id using `rb_empty', noreport) ///
    testname("rollback: empty using 1:1, fault after the merge")
use `rb_using', clear
drop if 1
tempfile rb_empty_grp
save `rb_empty_grp'
foreach f in record post {
    use `rb_master', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_empty_grp', noreport) ///
        testname("rollback: empty using m:1, fault at `f'")
    use `rb_master', clear
    drop if 1
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_using', noreport) ///
        testname("rollback: empty master m:1, fault at `f'")
}

* strL variables: merge itself runs the merge (it sorts the master before it
* checks the keys), under an in-memory snapshot frame by default
use `rb_master', clear
gen strL memo = "memo " + string(_n)
sort lab w
tempfile rb_strl
save `rb_strl'
_cmerge_rollback_case, rc(459) cmd(cmerge m:1 grp using `rb_using_dup', noreport) ///
    testname("rollback strL: duplicate using keys r(459)")
use `rb_strl', clear
_cmerge_rollback_case, fault(post) cmd(cmerge m:1 grp using `rb_using', noreport) ///
    testname("rollback strL: fault after the merge")

* Many observations: the plugin restores them with all threads
clear
set obs 400000
gen long id = 400001 - _n
gen int grp = mod(_n * 7919, 1000) + 1
gen double x = _n / 7
gen str12 s = "m" + string(mod(_n, 9973))
gen float z = runiform()
sort z
tempfile rb_big
save `rb_big'
clear
set obs 1200
gen int grp = _n
gen double x = _n * 3
gen str12 s = "u" + string(_n)
gen double gv = _n / 9
tempfile rb_big_using
save `rb_big_using'
foreach f in store post keep {
    use `rb_big', clear
    _cmerge_rollback_case, fault(`f') ///
        cmd(cmerge m:1 grp using `rb_big_using', keep(master match) update threads(4) noreport) ///
        testname("rollback 400K m:1: fault at `f'")
}

* After failures, merges succeed and match merge; the record is released
use `rb_master', clear
benchmark_merge m:1 grp using `rb_using', update testname("merge after rollbacks matches merge")
quietly frames dir
local frames `"`r(frames)'"'
foreach source in rb_master rb_strl {
    foreach mode in "" snapshot fast {
        use ``source'', clear
        quietly cmerge m:1 grp using `rb_using', keep(match) update noreport `mode'
        quietly datasignature
        local sig `r(datasignature)'
        if "`mode'" == "" local default_sig `sig'
        quietly frames dir
        if "$CTOOLS_cmu_active$CTOOLS_cmerge_plugin$CTOOLS_cmerge_snapframe" == "" & ///
            `"`r(frames)'"' == `"`frames'"' & "`sig'" == "`default_sig'" {
            test_pass "successful merge releases its record (`source' `mode')"
        }
        else {
            test_fail "successful merge releases its record (`source' `mode')" ///
                "globals or frames left behind, or result differs from the default"
        }
    }
}
use `rb_master', clear
benchmark_merge m:1 grp using `rb_using', keep(match) testname("merge after rollbacks: keep(match)")
use `rb_master', clear
capture cmerge m:1 grp using `rb_using', snapshot fast
if _rc == 198 test_pass "snapshot and fast together are rejected"
else test_fail "snapshot and fast together are rejected" "r(`=_rc')"

* A failed assert() leaves the merged result in memory, as merge does
use `rb_master', clear
capture merge m:1 grp using `rb_using', assert(match)
local stata_rc = _rc
local stata_N = _N
quietly datasignature
local stata_sig `r(datasignature)'
use `rb_master', clear
capture cmerge m:1 grp using `rb_using', assert(match) noreport
local cmerge_rc = _rc
local cmerge_N = _N
capture confirm variable _merge
local has_merge = !_rc
if `stata_rc' == 9 & `cmerge_rc' == 9 & `cmerge_N' == `stata_N' & `has_merge' & ///
    "$CTOOLS_cmu_active" == "" {
    test_pass "assert() failure keeps the merged result"
}
else {
    test_fail "assert() failure keeps the merged result" ///
        "rc `stata_rc'/`cmerge_rc', N `stata_N'/`cmerge_N', _merge `has_merge'"
}

* End of cmerge validation
noi print_summary "cmerge"
}

* Reached only after the complete component script.
global CTOOLS_COMPONENT_COMPLETE "cmerge"
    global CTOOLS_CASE_COMPLETE "cmerge"
    exit
}

di as error "Unknown cmerge case: `__case'"
exit 198
