/*******************************************************************************
 * validate_csample.do
 *
 * Validation tests for csample command
 * Tests random sampling without replacement functionality
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

* Plugin functionality check
clear
set obs 100
gen id = _n
set seed 12345
local N_before = _N
capture csample 50
if _rc != 0 {
    test_fail "csample plugin load" "plugin returned error `=_rc'"
    exit 1
}
if _N < `N_before' {
    test_pass "csample plugin loads and runs (50% of 100 = `=_N' obs)"
}
else {
    test_fail "csample plugin load" "no observations dropped"
    exit 1
}

/*******************************************************************************
 * SECTION 2: Basic percentage sampling
 ******************************************************************************/
print_section "Basic Percentage Sampling"

* 50% sample
clear
set obs 1000
gen id = _n
set seed 12345
csample 50
local n_after = _N
if `n_after' > 400 & `n_after' < 600 {
    test_pass "50% sample yields ~500 obs (got `n_after')"
}
else {
    test_fail "50% sample yields ~500 obs" "got `n_after'"
}

* 25% sample
clear
set obs 1000
gen id = _n
set seed 12345
csample 25
local n_after = _N
if `n_after' > 200 & `n_after' < 300 {
    test_pass "25% sample yields ~250 obs (got `n_after')"
}
else {
    test_fail "25% sample yields ~250 obs" "got `n_after'"
}

/*******************************************************************************
 * SECTION 3: Fixed count sampling
 ******************************************************************************/
print_section "Fixed Count Sampling"

clear
set obs 1000
gen id = _n
set seed 12345
csample, count(100)
if _N == 100 {
    test_pass "count(100) yields exactly 100 obs"
}
else {
    test_fail "count(100) yields exactly 100 obs" "got `=_N'"
}

clear
set obs 500
gen id = _n
set seed 12345
csample, count(50)
if _N == 50 {
    test_pass "count(50) yields exactly 50 obs"
}
else {
    test_fail "count(50) yields exactly 50 obs" "got `=_N'"
}

/*******************************************************************************
 * SECTION 4: By-group sampling
 ******************************************************************************/
print_section "By-Group Sampling"

* By-group percentage
clear
set obs 1000
gen group = mod(_n, 10)
gen id = _n
set seed 12345
csample 50, by(group)

local all_ok = 1
forvalues g = 0/9 {
    count if group == `g'
    local n_g = r(N)
    if `n_g' < 35 | `n_g' > 65 {
        local all_ok = 0
    }
}
if `all_ok' {
    test_pass "50% by-group sample (~50 per group)"
}
else {
    test_fail "50% by-group sample" "some groups outside [35,65] range"
}

* By-group fixed count
clear
set obs 1000
gen group = mod(_n, 10)
gen id = _n
set seed 12345
csample, count(10) by(group)

local all_ok = 1
forvalues g = 0/9 {
    count if group == `g'
    local n_g = r(N)
    if `n_g' != 10 {
        local all_ok = 0
    }
}
if `all_ok' & _N == 100 {
    test_pass "count(10) by-group yields 10 per group, 100 total"
}
else {
    test_fail "count(10) by-group" "expected 100 total, got `=_N'"
}

/*******************************************************************************
 * SECTION 5: Reproducibility with set seed
 ******************************************************************************/
print_section "Reproducibility"

clear
set obs 1000
gen id = _n
set seed 99999
csample 30
tempfile sample1
save `sample1'
local N1 = _N

clear
set obs 1000
gen id = _n
set seed 99999
csample 30
tempfile sample2
save `sample2'
local N2 = _N

* Check same count
if `N1' == `N2' {
    test_pass "same seed yields same count"
}
else {
    test_fail "same seed yields same count" "N1=`N1' N2=`N2'"
}

* Check same observations
use `sample1', clear
merge 1:1 id using `sample2'
count if _merge != 3
local n_diff = r(N)
if `n_diff' == 0 {
    test_pass "same seed yields identical samples"
}
else {
    test_fail "same seed yields identical samples" "`n_diff' obs differ"
}

/*******************************************************************************
 * SECTION 6: Edge cases
 ******************************************************************************/
print_section "Edge Cases"

* 100% sample
clear
set obs 100
gen id = _n
set seed 12345
csample 100
if _N == 100 {
    test_pass "100% sample keeps all obs"
}
else {
    test_fail "100% sample keeps all obs" "got `=_N'"
}

* 0% sample
clear
set obs 100
gen id = _n
set seed 12345
csample 0
if _N == 0 {
    test_pass "0% sample drops all obs"
}
else {
    test_fail "0% sample drops all obs" "got `=_N'"
}

* Small dataset
clear
set obs 10
gen id = _n
set seed 12345
csample 50
if _N >= 3 & _N <= 7 {
    test_pass "small dataset sample (N=10, 50%)"
}
else {
    test_fail "small dataset sample" "got `=_N', expected 3-7"
}

/*******************************************************************************
 * SECTION 7: if/in conditions
 ******************************************************************************/
print_section "If/In Conditions"

clear
set obs 200
gen id = _n
gen x = runiform()
gen keep_flag = id <= 100

set seed 12345
csample 50 if keep_flag == 1

local n_kept = _N
count if id > 100
local n_wrong = r(N)

if `n_wrong' == 0 & `n_kept' > 35 & `n_kept' < 65 {
    test_pass "if condition respected"
}
else {
    test_fail "if condition respected" "n_wrong=`n_wrong', n_kept=`n_kept'"
}

* End of csample validation
}
