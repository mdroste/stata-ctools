/*******************************************************************************
 * quick_validate.do - Quick validation for cmerge
 ******************************************************************************/

version 14.0
clear all
set more off

adopath + "build"
cap program drop ctools_plugin
program ctools_plugin, plugin using("build/ctools_mac_arm.plugin")

local passed = 0
local failed = 0

quietly {

noi di "========================================"
noi di "CMERGE QUICK VALIDATION"
noi di "========================================"

/* Test 1: Basic 1:1 merge */
noi di ""
noi di "Test 1: Basic 1:1 merge"

sysuse auto, clear
gen id = _n
tempfile master
save `master'

clear
sysuse auto
gen id = _n
keep id foreign
tempfile using_data
save `using_data'

use `master', clear
drop foreign
merge 1:1 id using `using_data'
tempfile stata_result
save `stata_result'

use `master', clear
drop foreign
cmerge 1:1 id using `using_data'

* Sort both and compare
sort id
tempfile cmerge_result
save `cmerge_result'

use `stata_result', clear
sort id
capture cf _all using `cmerge_result'
if _rc == 0 {
    noi di as result "  PASS"
    local ++passed
}
else {
    noi di as error "  FAIL"
    local ++failed
}

/* Test 2: m:1 merge */
noi di ""
noi di "Test 2: m:1 merge"

sysuse census, clear
keep state region pop
tempfile master
save `master'

clear
input byte region str20 region_name
1 "NE"
2 "N Cntrl"
3 "South"
4 "West"
end
tempfile using_data
save `using_data'

use `master', clear
merge m:1 region using `using_data'
sort state
tempfile stata_result
save `stata_result'

use `master', clear
cmerge m:1 region using `using_data'
sort state
capture cf _all using `stata_result'
if _rc == 0 {
    noi di as result "  PASS"
    local ++passed
}
else {
    noi di as error "  FAIL"
    local ++failed
}

/* Test 3: 1:m merge */
noi di ""
noi di "Test 3: 1:m merge"

sysuse auto, clear
keep foreign make price
bysort foreign: keep if _n <= 5
tempfile master
save `master'

clear
input byte foreign str10 label
0 "Domestic"
1 "Foreign"
end
tempfile using_data
save `using_data'

use `using_data', clear
merge 1:m foreign using `master'
sort foreign make
tempfile stata_result
save `stata_result'

use `using_data', clear
cmerge 1:m foreign using `master'
sort foreign make
capture cf _all using `stata_result'
if _rc == 0 {
    noi di as result "  PASS"
    local ++passed
}
else {
    noi di as error "  FAIL"
    local ++failed
}

/* Test 4: Merge with no matches */
noi di ""
noi di "Test 4: Merge with no matches"

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
tempfile using_data
save `using_data'

use `master', clear
merge 1:1 id using `using_data'
sort id
tempfile stata_result
save `stata_result'

use `master', clear
cmerge 1:1 id using `using_data'
sort id
capture cf _all using `stata_result'
if _rc == 0 {
    noi di as result "  PASS"
    local ++passed
}
else {
    noi di as error "  FAIL"
    local ++failed
}

/* Test 5: m:m merge */
noi di ""
noi di "Test 5: m:m merge"

clear
input int grp float value
1 10
1 20
2 30
2 40
end
tempfile master
save `master'

clear
input int grp str1 label
1 "A"
1 "B"
2 "C"
end
tempfile using_data
save `using_data'

use `master', clear
merge m:m grp using `using_data'
local stata_n = _N

use `master', clear
cmerge m:m grp using `using_data'
local cmerge_n = _N

if `stata_n' == `cmerge_n' {
    noi di as result "  PASS (N=`stata_n')"
    local ++passed
}
else {
    noi di as error "  FAIL (stata=`stata_n', cmerge=`cmerge_n')"
    local ++failed
}

noi di ""
noi di "========================================"
noi di "RESULTS: `passed' passed, `failed' failed"
noi di "========================================"

}

if `failed' > 0 {
    exit 1
}
