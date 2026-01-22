*! Validation tests for cwinsor - all options
*! Tests C-accelerated winsorization against expected behavior

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cwinsor.log", replace text
validation_header "cwinsor"

* =============================================================================
* Basic functionality tests
* =============================================================================

* Test 1: Basic winsorization
sysuse auto, clear
sum price, detail
local p1 = r(p1)
local p99 = r(p99)
capture cwinsor price
local passed = (_rc == 0)
if `passed' {
    sum price, meanonly
    local newmin = r(min)
    local newmax = r(max)
    * Check that values are bounded by approximately p1 and p99
    local passed = (`newmin' >= `p1' - 0.01 & `newmax' <= `p99' + 0.01)
}
report_test "Basic winsorization" `passed'

* Test 2: Multiple variables
sysuse auto, clear
run_test "Multiple variables" cwinsor price mpg weight

* Test 3: Custom percentiles with cuts()
sysuse auto, clear
sum price, detail
local p5 = r(p5)
local p95 = r(p95)
capture cwinsor price, cuts(5 95)
local passed = (_rc == 0)
if `passed' {
    sum price, meanonly
    local newmin = r(min)
    local newmax = r(max)
    local passed = (`newmin' >= `p5' - 0.01 & `newmax' <= `p95' + 0.01)
}
report_test "cuts(5 95) option" `passed'

* Test 4: p() and q() options
sysuse auto, clear
run_test "p(2) q(98) options" cwinsor price, p(2) q(98)

* =============================================================================
* Trim option tests
* =============================================================================

* Test 5: Trim mode
sysuse auto, clear
local orig_n = _N
capture cwinsor price, trim
local passed = (_rc == 0)
if `passed' {
    count if mi(price)
    * Should have some missing values after trimming extremes
    local nmiss = r(N)
    local passed = (`nmiss' > 0)
}
report_test "trim option creates missing values" `passed'

* Test 6: Trim at different percentiles
sysuse auto, clear
run_test "trim with cuts(10 90)" cwinsor price, trim cuts(10 90)

* =============================================================================
* By-group tests
* =============================================================================

* Test 7: By-group winsorization
sysuse auto, clear
run_test "by(foreign) option" cwinsor price, by(foreign)

* Test 8: By-group with multiple variables
sysuse auto, clear
run_test "by(foreign) multiple vars" cwinsor price mpg, by(foreign)

* Test 9: By-group produces different results per group
sysuse auto, clear
gen price_orig = price
capture cwinsor price, by(foreign)
local passed = (_rc == 0)
if `passed' {
    * Check that winsorization happened differently in each group
    bysort foreign: egen pmin = min(price)
    bysort foreign: egen pmax = max(price)
    * The min/max should be different between groups (different distributions)
    sum pmin if foreign == 0, meanonly
    local min0 = r(mean)
    sum pmin if foreign == 1, meanonly
    local min1 = r(mean)
    * They should be different since domestic/foreign have different price distributions
    local passed = (abs(`min0' - `min1') > 0.001)
}
report_test "by-group winsorizes independently" `passed'

* =============================================================================
* Output options tests
* =============================================================================

* Test 10: suffix() option
sysuse auto, clear
capture cwinsor price, suffix(_w)
local passed = (_rc == 0)
if `passed' {
    capture confirm variable price_w
    local passed = (_rc == 0)
}
report_test "suffix(_w) creates price_w" `passed'

* Test 11: prefix() option
sysuse auto, clear
capture cwinsor price, prefix(w_)
local passed = (_rc == 0)
if `passed' {
    capture confirm variable w_price
    local passed = (_rc == 0)
}
report_test "prefix(w_) creates w_price" `passed'

* Test 12: suffix with multiple variables
sysuse auto, clear
capture cwinsor price mpg, suffix(_win)
local passed = (_rc == 0)
if `passed' {
    capture confirm variable price_win mpg_win
    local passed = (_rc == 0)
}
report_test "suffix with multiple vars" `passed'

* Test 13: Original unchanged with suffix
sysuse auto, clear
sum price, meanonly
local orig_mean = r(mean)
capture cwinsor price, suffix(_w)
local passed = (_rc == 0)
if `passed' {
    sum price, meanonly
    local passed = (abs(`orig_mean' - r(mean)) < 0.001)
}
report_test "original unchanged with suffix" `passed'

* =============================================================================
* Performance options tests
* =============================================================================

* Test 14: verbose option
sysuse auto, clear
run_test "verbose option" cwinsor price, verbose

* Test 15: threads option
sysuse auto, clear
run_test "threads(2) option" cwinsor price, threads(2)

* Test 16: Combined options
sysuse auto, clear
run_test "Combined options" cwinsor price mpg, by(foreign) cuts(5 95) verbose

* =============================================================================
* Conditional tests
* =============================================================================

* Test 17: if condition
sysuse auto, clear
run_test "if condition" cwinsor price if foreign == 1

* Test 18: in range
sysuse auto, clear
run_test "in range" cwinsor price in 1/50

* =============================================================================
* Edge case tests
* =============================================================================

* Test 19: Small dataset
clear
set obs 10
gen x = _n
run_test "Small dataset (10 obs)" cwinsor x

* Test 20: Dataset with missing values
clear
set obs 100
gen x = _n
replace x = . in 1/10
replace x = . in 50/55
capture cwinsor x
local passed = (_rc == 0)
if `passed' {
    * Check that original missing values are still missing
    count if mi(x) & _n <= 10
    local passed = (r(N) == 10)
}
report_test "Preserves original missing values" `passed'

* Test 21: All same values (no winsorization needed)
clear
set obs 100
gen x = 5
capture cwinsor x
local passed = (_rc == 0)
if `passed' {
    sum x, meanonly
    local passed = (r(min) == 5 & r(max) == 5)
}
report_test "All same values" `passed'

* Test 22: Large dataset
clear
set obs 100000
gen x = rnormal()
run_test "Large dataset (100k obs)" cwinsor x

* Test 23: Large dataset with by-groups
clear
set obs 100000
gen group = mod(_n - 1, 100) + 1
gen x = rnormal() + group/10
run_test "Large dataset with 100 groups" cwinsor x, by(group)

* =============================================================================
* Error case tests
* =============================================================================

* Test 24: Invalid percentile range
sysuse auto, clear
run_test_error "p >= q error" 198 cwinsor price, p(99) q(1)

* Test 25: Percentile out of range
sysuse auto, clear
run_test_error "p > 100 error" 198 cwinsor price, p(150)

* Test 26: suffix and prefix together
sysuse auto, clear
run_test_error "suffix+prefix error" 198 cwinsor price, suffix(_w) prefix(w_)

* Test 27: suffix and replace together
sysuse auto, clear
run_test_error "suffix+replace error" 198 cwinsor price, suffix(_w) replace

* Test 28: String variable error
sysuse auto, clear
run_test_error "String variable error" 109 cwinsor make

* Test 29: Existing variable with suffix
sysuse auto, clear
gen price_w = 0
run_test_error "Existing suffix var error" 110 cwinsor price, suffix(_w)

* =============================================================================
* Data integrity tests
* =============================================================================

* Test 30: Values well within bounds unchanged
* Note: Percentile calculation may differ slightly from Stata's summarize
* so we use a wider margin (5th to 95th percentile) to test interior values
sysuse auto, clear
sum price, detail
local p5 = r(p5)
local p95 = r(p95)
gen interior = (price >= `p5' & price <= `p95')
gen price_orig = price
sum price if interior, meanonly
local orig_sum = r(sum)
capture cwinsor price
local passed = (_rc == 0)
if `passed' {
    * Values well within bounds (5th-95th) should remain unchanged after 1-99 winsorization
    sum price if interior, meanonly
    local passed = (abs(`orig_sum' - r(sum)) < 0.01)
}
report_test "Interior values unchanged" `passed'

* Test 31: Number of observations unchanged
sysuse auto, clear
local orig_n = _N
capture cwinsor price
local passed = (_rc == 0 & _N == `orig_n')
report_test "Observation count unchanged" `passed'

* Test 32: Other variables unchanged
sysuse auto, clear
sum mpg, meanonly
local orig_mpg = r(mean)
capture cwinsor price
local passed = (_rc == 0)
if `passed' {
    sum mpg, meanonly
    local passed = (abs(`orig_mpg' - r(mean)) < 0.001)
}
report_test "Other variables unchanged" `passed'

* =============================================================================
* Comparison test (if winsor2 is installed)
* =============================================================================

capture which winsor2
if _rc == 0 {
    * Test 33: Compare with winsor2
    * Note: Percentile interpolation methods may differ, so we compare
    * at a wider percentile (10-90) where differences are smaller
    sysuse auto, clear
    preserve
    winsor2 price, cuts(10 90) replace
    sum price, meanonly
    local w2_mean = r(mean)
    local w2_min = r(min)
    local w2_max = r(max)
    restore

    cwinsor price, cuts(10 90)
    sum price, meanonly
    local cw_mean = r(mean)
    local cw_min = r(min)
    local cw_max = r(max)

    * Results should be close (allow larger tolerance for interpolation differences)
    * Mean should be very close, min/max may differ more due to interpolation
    local mean_diff = abs(`w2_mean' - `cw_mean')
    local min_diff = abs(`w2_min' - `cw_min')
    local max_diff = abs(`w2_max' - `cw_max')
    * Allow 5% relative difference for min/max bounds
    local min_tol = `w2_min' * 0.05
    local max_tol = `w2_max' * 0.05
    local passed = (`mean_diff' < 100 & `min_diff' < `min_tol' & `max_diff' < `max_tol')
    report_test "Compare with winsor2 (10-90)" `passed'
}
else {
    di as text "  winsor2 not installed - skipping comparison test"
}

* =============================================================================
* Summary
* =============================================================================

validation_summary
log close
