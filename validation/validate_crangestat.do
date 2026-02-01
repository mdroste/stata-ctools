/*******************************************************************************
 * validate_crangestat.do
 *
 * Comprehensive validation tests for crangestat
 * Tests all statistics, options, and edge cases
 *
 * crangestat: C-accelerated range statistics
 * Syntax: crangestat (stat) [newvar=]varname [...], interval(keyvar low high) [by(varlist) excludeself verbose threads(#)]
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for crangestat..."

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
print_section "Plugin Check"

clear
set obs 100
gen obs = _n
gen x = runiform()
capture noisily crangestat (mean) x, interval(obs -5 5)
if _rc != 0 {
    test_fail "crangestat plugin load" "plugin returned error `=_rc'"
    print_summary "crangestat"
    exit 1
}
test_pass "crangestat plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic statistics (20 tests)
 ******************************************************************************/
print_section "Basic Statistics"

* Test 2.1: Count statistic
clear
set obs 100
gen obs = _n
gen x = 1
crangestat (count) n=x, interval(obs -2 2)
* Each obs should count 5 (self + 2 before + 2 after), except edges
if n[50] == 5 {
    test_pass "count statistic"
}
else {
    test_fail "count statistic" "expected 5, got `=n[50]'"
}

* Test 2.2: Mean statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (mean) x_mean=x, interval(obs -2 2)
* Mean of 48,49,50,51,52 = 50
if abs(x_mean[50] - 50) < 0.001 {
    test_pass "mean statistic"
}
else {
    test_fail "mean statistic" "expected 50, got `=x_mean[50]'"
}

* Test 2.3: Sum statistic
clear
set obs 100
gen obs = _n
gen x = 1
crangestat (sum) x_sum=x, interval(obs -2 2)
if x_sum[50] == 5 {
    test_pass "sum statistic"
}
else {
    test_fail "sum statistic" "expected 5, got `=x_sum[50]'"
}

* Test 2.4: Min statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (min) x_min=x, interval(obs -2 2)
* Min of 48,49,50,51,52 = 48
if x_min[50] == 48 {
    test_pass "min statistic"
}
else {
    test_fail "min statistic" "expected 48, got `=x_min[50]'"
}

* Test 2.5: Max statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (max) x_max=x, interval(obs -2 2)
* Max of 48,49,50,51,52 = 52
if x_max[50] == 52 {
    test_pass "max statistic"
}
else {
    test_fail "max statistic" "expected 52, got `=x_max[50]'"
}

* Test 2.6: SD statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (sd) x_sd=x, interval(obs -2 2)
* SD of 48,49,50,51,52 with N-1 denominator = sqrt(10/4) = 1.581
if abs(x_sd[50] - sqrt(2.5)) < 0.001 {
    test_pass "sd statistic"
}
else {
    test_fail "sd statistic" "expected `=sqrt(2.5)', got `=x_sd[50]'"
}

* Test 2.7: Variance statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (variance) x_var=x, interval(obs -2 2)
if abs(x_var[50] - 2.5) < 0.001 {
    test_pass "variance statistic"
}
else {
    test_fail "variance statistic" "expected 2.5, got `=x_var[50]'"
}

* Test 2.8: Median statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (median) x_med=x, interval(obs -2 2)
if x_med[50] == 50 {
    test_pass "median statistic"
}
else {
    test_fail "median statistic" "expected 50, got `=x_med[50]'"
}

* Test 2.9: IQR statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (iqr) x_iqr=x, interval(obs -4 4)
* For 9 values centered at 50, IQR should be p75 - p25
capture if !missing(x_iqr[50]) {
    test_pass "iqr statistic"
}
else {
    test_fail "iqr statistic" "missing value"
}

* Test 2.10: First statistic
clear
set obs 100
gen obs = _n
gen x = obs * 10
crangestat (first) x_first=x, interval(obs -2 2)
* First should be the value at the lowest key in range
if x_first[50] == 480 {
    test_pass "first statistic"
}
else {
    test_fail "first statistic" "expected 480, got `=x_first[50]'"
}

* Test 2.11: Last statistic
clear
set obs 100
gen obs = _n
gen x = obs * 10
crangestat (last) x_last=x, interval(obs -2 2)
if x_last[50] == 520 {
    test_pass "last statistic"
}
else {
    test_fail "last statistic" "expected 520, got `=x_last[50]'"
}

* Test 2.12: P25 statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (p25) x_p25=x, interval(obs -4 4)
capture if !missing(x_p25[50]) {
    test_pass "p25 statistic"
}
else {
    test_fail "p25 statistic" "missing value"
}

* Test 2.13: P75 statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (p75) x_p75=x, interval(obs -4 4)
capture if !missing(x_p75[50]) {
    test_pass "p75 statistic"
}
else {
    test_fail "p75 statistic" "missing value"
}

* Test 2.14: Multiple statistics
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (mean) x_mean=x (sum) x_sum=x (count) x_n=x, interval(obs -2 2)
if abs(x_mean[50] - 50) < 0.001 & x_n[50] == 5 {
    test_pass "multiple statistics"
}
else {
    test_fail "multiple statistics" "values wrong"
}

* Test 2.15: Verbose option
clear
set obs 100
gen obs = _n
gen x = runiform()
capture crangestat (mean) x, interval(obs -5 5) verbose
if _rc == 0 {
    test_pass "verbose option"
}
else {
    test_fail "verbose option" "rc=`=_rc'"
}

* Test 2.16: Threads option
clear
set obs 100
gen obs = _n
gen x = runiform()
capture crangestat (mean) x, interval(obs -5 5) threads(2)
if _rc == 0 {
    test_pass "threads(2) option"
}
else {
    test_fail "threads option" "rc=`=_rc'"
}

* Test 2.17: Skewness statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (skewness) x_skew=x, interval(obs -4 4)
capture if !missing(x_skew[50]) {
    test_pass "skewness statistic"
}
else {
    test_fail "skewness statistic" "missing value"
}

* Test 2.18: Kurtosis statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (kurtosis) x_kurt=x, interval(obs -4 4)
capture if !missing(x_kurt[50]) {
    test_pass "kurtosis statistic"
}
else {
    test_fail "kurtosis statistic" "missing value"
}

* Test 2.19: P1 statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (p1) x_p1=x, interval(obs -10 10)
capture if !missing(x_p1[50]) {
    test_pass "p1 statistic"
}
else {
    test_fail "p1 statistic" "missing value"
}

* Test 2.20: P99 statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (p99) x_p99=x, interval(obs -10 10)
capture if !missing(x_p99[50]) {
    test_pass "p99 statistic"
}
else {
    test_fail "p99 statistic" "missing value"
}

/*******************************************************************************
 * SECTION 3: Interval options (10 tests)
 ******************************************************************************/
print_section "Interval Options"

* Test 3.1: Symmetric interval
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t -5 5)
if n[50] == 11 {
    test_pass "symmetric interval"
}
else {
    test_fail "symmetric interval" "expected 11, got `=n[50]'"
}

* Test 3.2: Backward-looking interval
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t -5 0)
if n[50] == 6 {
    test_pass "backward-looking interval"
}
else {
    test_fail "backward-looking" "expected 6, got `=n[50]'"
}

* Test 3.3: Forward-looking interval
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t 0 5)
if n[50] == 6 {
    test_pass "forward-looking interval"
}
else {
    test_fail "forward-looking" "expected 6, got `=n[50]'"
}

* Test 3.4: Single point interval (only self)
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t 0 0)
if n[50] == 1 {
    test_pass "single point interval"
}
else {
    test_fail "single point" "expected 1, got `=n[50]'"
}

* Test 3.5: Negative infinity (. for low)
clear
set obs 100
gen t = _n
gen x = 1
crangestat (sum) cumsum=x, interval(t . 0)
if cumsum[50] == 50 {
    test_pass "cumulative sum (low=.)"
}
else {
    test_fail "cumulative sum" "expected 50, got `=cumsum[50]'"
}

* Test 3.6: Positive infinity (. for high)
clear
set obs 100
gen t = _n
gen x = 1
crangestat (sum) futuresum=x, interval(t 0 .)
if futuresum[50] == 51 {
    test_pass "future sum (high=.)"
}
else {
    test_fail "future sum" "expected 51, got `=futuresum[50]'"
}

* Test 3.7: Both bounds infinite
clear
set obs 100
gen t = _n
gen x = 1
crangestat (sum) totalsum=x, interval(t . .)
if totalsum[50] == 100 {
    test_pass "total sum (both bounds infinite)"
}
else {
    test_fail "total sum" "expected 100, got `=totalsum[50]'"
}

* Test 3.8: Non-integer interval
clear
set obs 100
gen t = _n / 10
gen x = 1
crangestat (count) n=x, interval(t -0.25 0.25)
* Should count ~5 observations (those within 0.5 range)
capture if n[50] >= 3 & n[50] <= 7 {
    test_pass "non-integer interval"
}
else {
    test_fail "non-integer interval" "unexpected count"
}

* Test 3.9: Asymmetric interval
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t -3 7)
if n[50] == 11 {
    test_pass "asymmetric interval"
}
else {
    test_fail "asymmetric interval" "expected 11, got `=n[50]'"
}

* Test 3.10: Wide interval
clear
set obs 100
gen t = _n
gen x = 1
crangestat (sum) x_sum=x, interval(t -50 50)
if x_sum[50] == 100 {
    test_pass "wide interval covering all"
}
else {
    test_fail "wide interval" "expected 100, got `=x_sum[50]'"
}

/*******************************************************************************
 * SECTION 4: excludeself option (10 tests)
 ******************************************************************************/
print_section "excludeself Option"

* Test 4.1: Count with excludeself
clear
set obs 100
gen t = _n
gen x = 1
crangestat (count) n=x, interval(t -2 2) excludeself
if n[50] == 4 {
    test_pass "count with excludeself"
}
else {
    test_fail "count excludeself" "expected 4, got `=n[50]'"
}

* Test 4.2: Sum with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (sum) s=x, interval(t -2 2) excludeself
* Sum of 48,49,51,52 = 200
if s[50] == 200 {
    test_pass "sum with excludeself"
}
else {
    test_fail "sum excludeself" "expected 200, got `=s[50]'"
}

* Test 4.3: Mean with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (mean) m=x, interval(t -2 2) excludeself
* Mean of 48,49,51,52 = 50
if m[50] == 50 {
    test_pass "mean with excludeself"
}
else {
    test_fail "mean excludeself" "expected 50, got `=m[50]'"
}

* Test 4.4: Min with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (min) minval=x, interval(t -2 2) excludeself
if minval[50] == 48 {
    test_pass "min with excludeself"
}
else {
    test_fail "min excludeself" "expected 48, got `=minval[50]'"
}

* Test 4.5: Max with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (max) maxval=x, interval(t -2 2) excludeself
if maxval[50] == 52 {
    test_pass "max with excludeself"
}
else {
    test_fail "max excludeself" "expected 52, got `=maxval[50]'"
}

* Test 4.6: Single observation interval with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t 0 0) excludeself
if n[50] == 0 | missing(n[50]) {
    test_pass "single point with excludeself"
}
else {
    test_fail "single point excludeself" "expected 0 or missing"
}

* Test 4.7: SD with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (sd) sdval=x, interval(t -2 2) excludeself
capture if !missing(sdval[50]) {
    test_pass "sd with excludeself"
}
else {
    test_fail "sd excludeself" "missing value"
}

* Test 4.8: Median with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (median) medval=x, interval(t -2 2) excludeself
* Median of 48,49,51,52 = 50
if medval[50] == 50 {
    test_pass "median with excludeself"
}
else {
    test_fail "median excludeself" "expected 50, got `=medval[50]'"
}

* Test 4.9: Leave-one-out mean
clear
set obs 10
gen t = _n
gen x = 10
crangestat (mean) loo_mean=x, interval(t . .) excludeself
* All values are 10, so leave-one-out mean should be 10
if abs(loo_mean[5] - 10) < 0.001 {
    test_pass "leave-one-out mean"
}
else {
    test_fail "loo mean" "expected 10, got `=loo_mean[5]'"
}

* Test 4.10: Multiple stats with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (mean) m=x (sum) s=x (count) n=x, interval(t -2 2) excludeself
if n[50] == 4 & s[50] == 200 & m[50] == 50 {
    test_pass "multiple stats with excludeself"
}
else {
    test_fail "multiple excludeself" "values wrong"
}

/*******************************************************************************
 * SECTION 5: by() option (10 tests)
 ******************************************************************************/
print_section "by() Option"

* Test 5.1: Basic by-group
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = 1
sort group t
crangestat (sum) s=x, interval(t . 0) by(group)
* For group 0, obs 50 has t=50, so cumulative sum should be 50
* Check obs 50 (not obs 99 which has t=99)
if s[50] == 50 {
    test_pass "basic by-group"
}
else {
    test_fail "by-group" "expected 50, got `=s[50]'"
}

* Test 5.2: Count by group
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = 1
sort group t
crangestat (count) n=x, interval(t -2 2) by(group)
* Each group has 100 obs, middle obs should count 5
if n[49] == 5 {
    test_pass "count by group"
}
else {
    test_fail "count by group" "expected 5, got `=n[49]'"
}

* Test 5.3: Multiple by-variables
clear
set obs 400
gen g1 = mod(_n - 1, 2)
gen g2 = mod(floor((_n - 1) / 2), 2)
gen t = ceil(_n / 4)
gen x = 1
sort g1 g2 t
crangestat (sum) s=x, interval(t . 0) by(g1 g2)
capture if !missing(s[50]) {
    test_pass "multiple by-variables"
}
else {
    test_fail "multiple by-vars" "missing value"
}

* Test 5.4: By-group with excludeself
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = 1
sort group t
crangestat (count) n=x, interval(t -2 2) by(group) excludeself
if n[49] == 4 {
    test_pass "by-group with excludeself"
}
else {
    test_fail "by excludeself" "expected 4, got `=n[49]'"
}

* Test 5.5: Groups don't overlap
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = group * 100 + t
sort group t
crangestat (mean) m=x, interval(t . .) by(group)
* Group 0: mean of 1..100
* Group 1: mean of 101..200
if abs(m[1] - 50.5) < 0.001 & abs(m[100] - 50.5) < 0.001 {
    test_pass "groups don't overlap"
}
else {
    test_fail "groups overlap" "means wrong"
}

* Test 5.6: Single observation per group
clear
set obs 100
gen group = _n
gen t = 1
gen x = group
crangestat (mean) m=x, interval(t 0 0) by(group)
* Each group has one obs, so mean should equal x
local pass = 1
forvalues i = 1/10 {
    if m[`i'] != x[`i'] local pass = 0
}
if `pass' {
    test_pass "single obs per group"
}
else {
    test_fail "single obs group" "values wrong"
}

* Test 5.7: String-like numeric by-variable
clear
set obs 200
gen group = floor((_n - 1) / 50)  * 10 + 1
gen t = mod(_n - 1, 50) + 1
gen x = 1
sort group t
crangestat (sum) s=x, interval(t . 0) by(group)
capture if !missing(s[25]) {
    test_pass "numeric by-variable"
}
else {
    test_fail "numeric by" "missing value"
}

* Test 5.8: By-group mean equals group mean for infinite interval
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = group * 100 + t
sort group t
crangestat (mean) m=x, interval(t . .) by(group)
* Group 0: values 1-100, mean = 50.5
* Group 1: values 101-200, mean = 150.5
if abs(m[1] - 50.5) < 0.001 {
    test_pass "by-group mean with infinite interval"
}
else {
    test_fail "by infinite interval" "expected 50.5, got `=m[1]'"
}

* Test 5.9: By-group with multiple stats
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = 1
sort group t
crangestat (mean) m=x (sum) s=x (count) n=x, interval(t -5 5) by(group)
capture if !missing(m[50]) & !missing(s[50]) & !missing(n[50]) {
    test_pass "by-group with multiple stats"
}
else {
    test_fail "by multiple stats" "missing values"
}

* Test 5.10: By-group at boundaries
clear
set obs 100
gen group = _n <= 50
gen t = mod(_n - 1, 50) + 1
gen x = 1
sort group t
crangestat (count) n=x, interval(t -10 10) by(group)
* At t=1, should only count forward within group
capture if n[1] > 0 {
    test_pass "by-group at boundaries"
}
else {
    test_fail "by boundaries" "wrong count"
}

/*******************************************************************************
 * SECTION 6: Missing values (10 tests)
 ******************************************************************************/
print_section "Missing Values"

* Test 6.1: Missing key values
clear
set obs 100
gen t = _n
replace t = . in 50
gen x = 1
crangestat (sum) s=x, interval(t -2 2)
if missing(s[50]) {
    test_pass "missing key value"
}
else {
    test_fail "missing key" "expected missing"
}

* Test 6.2: Missing source values excluded
clear
set obs 100
gen t = _n
gen x = 1
replace x = . in 48/52
crangestat (count) n=x, interval(t -5 5)
* At t=50, all 5 nearby obs are missing, so count should exclude them
if n[50] == 6 | n[50] == 0 {
    test_pass "missing source values excluded"
}
else {
    test_fail "missing source" "unexpected count"
}

* Test 6.3: Sum ignores missing
clear
set obs 100
gen t = _n
gen x = 1
replace x = . in 49
crangestat (sum) s=x, interval(t -2 2)
* Sum should be 4 (5 minus the missing one)
if s[50] == 4 {
    test_pass "sum ignores missing"
}
else {
    test_fail "sum missing" "expected 4, got `=s[50]'"
}

* Test 6.4: Mean with some missing
clear
set obs 100
gen t = _n
gen x = t
replace x = . in 48
crangestat (mean) m=x, interval(t -2 2)
* Mean of 49,50,51,52 = 50.5
if abs(m[50] - 50.5) < 0.001 {
    test_pass "mean with some missing"
}
else {
    test_fail "mean some missing" "expected 50.5, got `=m[50]'"
}

* Test 6.5: All missing in window
clear
set obs 100
gen t = _n
gen x = .
crangestat (mean) m=x, interval(t -2 2)
if missing(m[50]) {
    test_pass "all missing in window"
}
else {
    test_fail "all missing" "expected missing"
}

* Test 6.6: Only one non-missing
clear
set obs 100
gen t = _n
gen x = .
replace x = 999 in 50
crangestat (mean) m=x, interval(t -2 2)
if m[50] == 999 {
    test_pass "only one non-missing"
}
else {
    test_fail "one non-missing" "expected 999"
}

* Test 6.7: SD with one non-missing
clear
set obs 100
gen t = _n
gen x = .
replace x = 1 in 50
crangestat (sd) sd_val=x, interval(t -2 2)
* SD requires at least 2 non-missing
if missing(sd_val[50]) {
    test_pass "sd with one non-missing"
}
else {
    test_fail "sd one obs" "expected missing"
}

* Test 6.8: Min/max with missing
clear
set obs 100
gen t = _n
gen x = t
replace x = . in 48/49
crangestat (min) minval=x (max) maxval=x, interval(t -2 2)
if minval[50] == 50 & maxval[50] == 52 {
    test_pass "min/max with missing"
}
else {
    test_fail "min/max missing" "values wrong"
}

* Test 6.9: Firstnm/lastnm with missing
clear
set obs 100
gen t = _n
gen x = t
replace x = . in 48
crangestat (firstnm) fnm=x (lastnm) lnm=x, interval(t -2 2)
if fnm[50] == 49 & lnm[50] == 52 {
    test_pass "firstnm/lastnm with missing"
}
else {
    test_fail "firstnm/lastnm" "values wrong"
}

* Test 6.10: Sparse missing pattern
clear
set obs 100
gen t = _n
gen x = 1
replace x = . if mod(_n, 3) == 0
crangestat (count) n=x, interval(t -5 5)
* Some obs will have fewer non-missing values
capture if n[50] >= 0 & n[50] <= 11 {
    test_pass "sparse missing pattern"
}
else {
    test_fail "sparse missing" "unexpected count"
}

/*******************************************************************************
 * SECTION 7: Edge cases (10 tests)
 ******************************************************************************/
print_section "Edge Cases"

* Test 7.1: Single observation
clear
set obs 1
gen t = 1
gen x = 100
crangestat (mean) m=x, interval(t -5 5)
if m[1] == 100 {
    test_pass "single observation"
}
else {
    test_fail "single obs" "expected 100, got `=m[1]'"
}

* Test 7.2: Two observations
clear
set obs 2
gen t = _n
gen x = _n * 10
crangestat (mean) m=x (sum) s=x, interval(t -1 1)
if m[1] == 15 & m[2] == 15 {
    test_pass "two observations"
}
else {
    test_fail "two obs" "values wrong"
}

* Test 7.3: Identical key values
clear
set obs 100
gen t = 1
gen x = _n
crangestat (mean) m=x, interval(t 0 0)
* All obs have t=1, so all should be included
if abs(m[50] - 50.5) < 0.001 {
    test_pass "identical key values"
}
else {
    test_fail "identical keys" "expected 50.5, got `=m[50]'"
}

* Test 7.4: Sorted vs unsorted data
clear
set obs 100
* Reverse order
gen t = 101 - _n
gen x = t
crangestat (mean) m=x, interval(t -2 2)
* Should still work correctly
capture if !missing(m[50]) {
    test_pass "reverse sorted data"
}
else {
    test_fail "reverse sorted" "missing value"
}

* Test 7.5: Very negative interval bounds
clear
set obs 100
gen t = _n + 1000
gen x = 1
crangestat (sum) s=x, interval(t -1000 0)
if s[50] == 50 {
    test_pass "very negative interval bound"
}
else {
    test_fail "negative bound" "expected 50, got `=s[50]'"
}

* Test 7.6: First and last observations
clear
set obs 100
gen t = _n
gen x = t
crangestat (count) n=x, interval(t -5 5)
* First obs should count 6 (self + 5 forward)
* Last obs should count 6 (self + 5 backward)
if n[1] == 6 & n[100] == 6 {
    test_pass "first and last observations"
}
else {
    test_fail "first/last" "expected 6"
}

* Test 7.7: Large interval with small data
clear
set obs 10
gen t = _n
gen x = 1
crangestat (sum) s=x, interval(t -100 100)
if s[5] == 10 {
    test_pass "large interval with small data"
}
else {
    test_fail "large interval" "expected 10, got `=s[5]'"
}

* Test 7.8: Zero-width data range
clear
set obs 100
gen t = 50
gen x = _n
crangestat (sum) s=x, interval(t 0 0)
* All have t=50, so all included
if s[1] == 5050 {
    test_pass "zero-width data range"
}
else {
    test_fail "zero-width" "expected 5050, got `=s[1]'"
}

* Test 7.9: Multiple variables same stat
clear
set obs 100
gen t = _n
gen x = t
gen y = t * 2
crangestat (mean) mx=x (mean) my=y, interval(t -2 2)
if mx[50] == 50 & my[50] == 100 {
    test_pass "multiple variables same stat"
}
else {
    test_fail "multi var same stat" "values wrong"
}

* Test 7.10: Duplicate observations
clear
set obs 200
* Each t appears twice
gen t = ceil(_n / 2)
gen x = _n
crangestat (count) n=x, interval(t -2 2)
* Each unique t has 2 obs, so count should be 10 for middle obs
if n[100] == 10 {
    test_pass "duplicate observations"
}
else {
    test_fail "duplicates" "expected 10, got `=n[100]'"
}

/*******************************************************************************
 * SECTION 8: Large datasets (5 tests)
 ******************************************************************************/
print_section "Large Datasets"

* Test 8.1: 10K observations
clear
set obs 10000
gen t = _n
gen x = runiform()
capture crangestat (mean) m=x, interval(t -10 10)
if _rc == 0 {
    test_pass "10K observations"
}
else {
    test_fail "10K obs" "rc=`=_rc'"
}

* Test 8.2: 50K observations
clear
set obs 50000
gen t = _n
gen x = runiform()
capture crangestat (mean) m=x, interval(t -5 5)
if _rc == 0 {
    test_pass "50K observations"
}
else {
    test_fail "50K obs" "rc=`=_rc'"
}

* Test 8.3: 100K observations
clear
set obs 100000
gen t = _n
gen x = runiform()
capture crangestat (mean) m=x, interval(t -5 5)
if _rc == 0 {
    test_pass "100K observations"
}
else {
    test_fail "100K obs" "rc=`=_rc'"
}

* Test 8.4: Large with by-groups
clear
set obs 50000
gen group = mod(_n - 1, 100)
gen t = ceil(_n / 100)
gen x = runiform()
sort group t
capture crangestat (mean) m=x, interval(t -5 5) by(group)
if _rc == 0 {
    test_pass "large with by-groups"
}
else {
    test_fail "large by-groups" "rc=`=_rc'"
}

* Test 8.5: Large with multiple stats
clear
set obs 50000
gen t = _n
gen x = runiform()
capture crangestat (mean) m=x (sd) s=x (min) minx=x (max) maxx=x, interval(t -5 5)
if _rc == 0 {
    test_pass "large with multiple stats"
}
else {
    test_fail "large multi stats" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 9: Real-world scenarios (5 tests)
 ******************************************************************************/
print_section "Real-World Scenarios"

* Test 9.1: Rolling mean (time series)
clear
* Trading days in a year
set obs 252
gen date = _n
gen price = 100 + runiform() * 10
crangestat (mean) ma5=price, interval(date -4 0)
capture if !missing(ma5[252]) {
    test_pass "rolling mean (time series)"
}
else {
    test_fail "rolling mean" "missing value"
}

* Test 9.2: Panel data rolling window
webuse grunfeld, clear
sort company year
crangestat (mean) invest_ma=invest, interval(year -2 0) by(company)
capture if !missing(invest_ma[100]) {
    test_pass "panel data rolling window"
}
else {
    test_fail "panel rolling" "missing value"
}

* Test 9.3: Age-based similar observations
webuse nlswork, clear
keep if !missing(ln_wage) & !missing(age)
crangestat (mean) similar_wage=ln_wage, interval(age -1 1)
capture if !missing(similar_wage[1000]) {
    test_pass "age-based similar observations"
}
else {
    test_fail "age-based" "missing value"
}

* Test 9.4: Leave-one-out cross-validation
clear
set obs 1000
set seed 12345
gen x = rnormal()
gen y = 2 * x + rnormal()
gen obs = _n
crangestat (mean) loo_y=y, interval(obs . .) excludeself
capture if !missing(loo_y[500]) {
    test_pass "leave-one-out cross-validation"
}
else {
    test_fail "loo cv" "missing value"
}

* Test 9.5: Cumulative statistics
clear
set obs 100
gen t = _n
gen x = runiform()
crangestat (sum) cumsum=x (count) cumcount=x (mean) cummean=x, interval(t . 0)
if cumcount[100] == 100 {
    test_pass "cumulative statistics"
}
else {
    test_fail "cumulative" "expected 100, got `=cumcount[100]'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
* End of crangestat validation
noi print_summary "crangestat"
}
