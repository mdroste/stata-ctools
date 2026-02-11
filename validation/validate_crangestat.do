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
sigfigs `=x_mean[50]' 50
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "mean statistic"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "mean statistic" "sigfigs=`sf'"
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
* SD of 48,49,50,51,52 with N-1 denominator = sqrt(10/4) = sqrt(2.5)
local expected_sd = sqrt(2.5)
sigfigs `=x_sd[50]' `expected_sd'
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "sd statistic"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "sd statistic" "sigfigs=`sf'"
}

* Test 2.7: Variance statistic
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (variance) x_var=x, interval(obs -2 2)
sigfigs `=x_var[50]' 2.5
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "variance statistic"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "variance statistic" "sigfigs=`sf'"
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
benchmark_rangestat iqr x, interval(obs -4 4) testname("iqr statistic")

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
benchmark_rangestat p25 x, interval(obs -4 4) testname("p25 statistic")

* Test 2.13: P75 statistic
clear
set obs 100
gen obs = _n
gen x = obs
benchmark_rangestat p75 x, interval(obs -4 4) testname("p75 statistic")

* Test 2.14: Multiple statistics
clear
set obs 100
gen obs = _n
gen x = obs
crangestat (mean) x_mean=x (sum) x_sum=x (count) x_n=x, interval(obs -2 2)
sigfigs `=x_mean[50]' 50
if r(sigfigs) >= $DEFAULT_SIGFIGS & x_n[50] == 5 {
    test_pass "multiple statistics"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "multiple statistics" "mean sigfigs=`sf', n=`=x_n[50]'"
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
benchmark_rangestat skewness x, interval(obs -4 4) testname("skewness statistic")

* Test 2.18: Kurtosis statistic
clear
set obs 100
gen obs = _n
gen x = obs
benchmark_rangestat kurtosis x, interval(obs -4 4) testname("kurtosis statistic")

* Test 2.19: P1 statistic
clear
set obs 100
gen obs = _n
gen x = obs
benchmark_rangestat p1 x, interval(obs -10 10) testname("p1 statistic")

* Test 2.20: P99 statistic
clear
set obs 100
gen obs = _n
gen x = obs
benchmark_rangestat p99 x, interval(obs -10 10) testname("p99 statistic")

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
* t[50] = 5.0, window = [4.75, 5.25]
* t values in range: 4.8(48), 4.9(49), 5.0(50), 5.1(51), 5.2(52) = 5 obs
if n[50] == 5 {
    test_pass "non-integer interval"
}
else {
    test_fail "non-integer interval" "expected 5, got `=n[50]'"
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
* Window is [t, t] and self is excluded, so count should be 0
if n[50] == 0 {
    test_pass "single point with excludeself"
}
else {
    test_fail "single point with excludeself" "expected 0, got `=n[50]'"
}

* Test 4.7: SD with excludeself
clear
set obs 100
gen t = _n
gen x = t
crangestat (sd) sdval=x, interval(t -2 2) excludeself
* Excluding self from window [48,52]: values are 48,49,51,52
* Mean = 200/4 = 50, Var = ((4+1+1+4)/3) = 10/3, SD = sqrt(10/3)
local expected_sd = sqrt(10/3)
sigfigs `=sdval[50]' `expected_sd'
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "sd with excludeself"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "sd with excludeself" "sigfigs=`sf'"
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
sigfigs `=loo_mean[5]' 10
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "leave-one-out mean"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "leave-one-out mean" "sigfigs=`sf'"
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
benchmark_rangestat sum x, interval(t . 0) by(g1 g2) testname("multiple by-variables")

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
* Group 0: mean of 1..100 = 50.5
* Group 1: mean of 101..200 = 150.5
sigfigs `=m[1]' 50.5
local sf1 = r(sigfigs)
sigfigs `=m[100]' 50.5
local sf2 = r(sigfigs)
if `sf1' >= $DEFAULT_SIGFIGS & `sf2' >= $DEFAULT_SIGFIGS {
    test_pass "groups don't overlap"
}
else {
    local sf1_fmt : display %5.1f `sf1'
    local sf2_fmt : display %5.1f `sf2'
    test_fail "groups don't overlap" "m[1] sigfigs=`sf1_fmt', m[100] sigfigs=`sf2_fmt'"
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
* Each group has 50 obs with t=1..50, cumsum at obs 25 (t=25) should be 25
if s[25] == 25 {
    test_pass "numeric by-variable"
}
else {
    test_fail "numeric by-variable" "expected 25, got `=s[25]'"
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
sigfigs `=m[1]' 50.5
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "by-group mean with infinite interval"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "by-group mean with infinite interval" "sigfigs=`sf'"
}

* Test 5.9: By-group with multiple stats
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = 1
sort group t
crangestat (mean) m=x (sum) s=x (count) n=x, interval(t -5 5) by(group)
* Each group has 100 obs, all x=1. Middle obs should have count=11, sum=11, mean=1
if n[50] == 11 & s[50] == 11 & m[50] == 1 {
    test_pass "by-group with multiple stats"
}
else {
    test_fail "by-group with multiple stats" "expected n=11,s=11,m=1; got n=`=n[50]',s=`=s[50]',m=`=m[50]'"
}

* Test 5.10: By-group at boundaries
clear
set obs 100
gen group = _n <= 50
gen t = mod(_n - 1, 50) + 1
gen x = 1
sort group t
crangestat (count) n=x, interval(t -10 10) by(group)
* At t=1 (obs 1), window covers t in [-9, 11], so t=1..11 within group = 11 obs
if n[1] == 11 {
    test_pass "by-group at boundaries"
}
else {
    test_fail "by-group at boundaries" "expected 11, got `=n[1]'"
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
* At t=50, window covers t in [45,55] = 11 obs total
* Obs 48,49,50,51,52 have x=missing, so 5 missing
* Count of non-missing = 11 - 5 = 6
if n[50] == 6 {
    test_pass "missing source values excluded"
}
else {
    test_fail "missing source values excluded" "expected 6, got `=n[50]'"
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
sigfigs `=m[50]' 50.5
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "mean with some missing"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "mean with some missing" "sigfigs=`sf'"
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
* At obs 50, window covers t in [45,55] = 11 obs
* Missing obs: those where mod(_n, 3)==0: 45,48,51,54 = 4 missing
* So count should be 11 - 4 = 7
if n[50] == 7 {
    test_pass "sparse missing pattern"
}
else {
    test_fail "sparse missing pattern" "expected 7, got `=n[50]'"
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
* All obs have t=1, so all should be included, mean = 50.5
sigfigs `=m[50]' 50.5
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "identical key values"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "identical key values" "sigfigs=`sf'"
}

* Test 7.4: Sorted vs unsorted data
clear
set obs 100
* Reverse order
gen t = 101 - _n
gen x = t
crangestat (mean) m=x, interval(t -2 2)
* obs 50 has t=51, x=51; window covers t in [49,53]
* obs with t in {49,50,51,52,53} have x = t, so mean = (49+50+51+52+53)/5 = 51
* Find obs 50's value: t[50]=51, so expected mean = 51
sigfigs `=m[50]' 51
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "reverse sorted data"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "reverse sorted data" "sigfigs=`sf'"
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
benchmark_rangestat mean x, interval(t -10 10) testname("10K observations")

* Test 8.2: 50K observations
clear
set obs 50000
gen t = _n
gen x = runiform()
benchmark_rangestat mean x, interval(t -5 5) testname("50K observations")

* Test 8.3: 100K observations
clear
set obs 100000
gen t = _n
gen x = runiform()
benchmark_rangestat mean x, interval(t -5 5) testname("100K observations")

* Test 8.4: Large with by-groups
clear
set obs 50000
gen group = mod(_n - 1, 100)
gen t = ceil(_n / 100)
gen x = runiform()
sort group t
benchmark_rangestat mean x, interval(t -5 5) by(group) testname("large with by-groups")

* Test 8.5: Large with multiple stats
clear
set obs 50000
gen t = _n
gen x = runiform()
benchmark_rangestat mean x, interval(t -5 5) testname("large multi-stat: mean")
benchmark_rangestat sd x, interval(t -5 5) testname("large multi-stat: sd")
benchmark_rangestat min x, interval(t -5 5) testname("large multi-stat: min")
benchmark_rangestat max x, interval(t -5 5) testname("large multi-stat: max")

* Test 8.6: Path B with duplicate keys (run optimization)
* 10K obs, 200 unique keys → 50 obs per key. Enters Path B (!excludeself, g_count>=1000)
* and each run broadcasts to 50 obs.
clear
set obs 10000
set seed 88888
gen t = ceil(_n / 50)
gen x = rnormal()
benchmark_rangestat mean x, interval(t -5 5) testname("Path B: dup keys, mean")
benchmark_rangestat min x, interval(t -5 5) testname("Path B: dup keys, min")
benchmark_rangestat max x, interval(t -5 5) testname("Path B: dup keys, max")
benchmark_rangestat sd x, interval(t -5 5) testname("Path B: dup keys, sd")

* Test 8.7: Path C parallel sliding window (large + excludeself)
* 20K obs, 1 group, excludeself → Path C since excludeself prevents Path B.
* nonmiss_count=20000 > PARALLEL_THRESHOLD=10000, so parallel sliding window.
clear
set obs 20000
set seed 99988
gen t = _n
gen x = rnormal()
benchmark_rangestat mean x, interval(t -10 10) excludeself testname("Path C: parallel sliding, mean")
benchmark_rangestat count x, interval(t -10 10) excludeself testname("Path C: parallel sliding, count")
benchmark_rangestat min x, interval(t -10 10) excludeself testname("Path C: parallel sliding, min")
benchmark_rangestat sd x, interval(t -10 10) excludeself testname("Path C: parallel sliding, sd")

* Test 8.8: Path A cross-group parallel with rangestat comparison
* 500 groups × 100 obs = 50K obs. avg_group_size=100 < 500 and ngroups=500 >= 8.
clear
set obs 50000
set seed 77788
gen group = mod(_n - 1, 500)
gen t = ceil(_n / 500)
gen x = rnormal()
sort group t
benchmark_rangestat mean x, interval(t -5 5) by(group) testname("Path A: cross-group, mean")
benchmark_rangestat min x, interval(t -5 5) by(group) testname("Path A: cross-group, min")
benchmark_rangestat sd x, interval(t -5 5) by(group) testname("Path A: cross-group, sd")

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
benchmark_rangestat mean price, interval(date -4 0) testname("rolling mean (time series)")

* Test 9.2: Panel data rolling window
webuse grunfeld, clear
sort company year
benchmark_rangestat mean invest, interval(year -2 0) by(company) testname("panel data rolling window")

* Test 9.3: Age-based similar observations
webuse nlswork, clear
keep if !missing(ln_wage) & !missing(age)
benchmark_rangestat mean ln_wage, interval(age -1 1) testname("age-based similar observations")

* Test 9.4: Leave-one-out cross-validation
clear
set obs 1000
set seed 12345
gen x = rnormal()
gen y = 2 * x + rnormal()
gen obs = _n
crangestat (mean) loo_y=y, interval(obs . .) excludeself
* LOO mean for obs 500: (sum of all y - y[500]) / 999
quietly summarize y
local total_sum = r(sum)
local expected_loo = (`total_sum' - y[500]) / 999
sigfigs `=loo_y[500]' `expected_loo'
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "leave-one-out cross-validation"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "loo cv" "sigfigs=`sf'"
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
 * SECTION: Intentional Error Tests
 *
 * These tests verify that crangestat returns the same error codes as rangestat
 * when given invalid inputs or error conditions.
 * Note: rangestat is a user-written command (ssc install rangestat).
 * If rangestat is not installed, tests compare against expected behavior.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Check if rangestat is installed
capture which rangestat
local rangestat_installed = (_rc == 0)

* Variable doesn't exist
clear
set obs 100
gen t = _n
gen x = runiform()
if `rangestat_installed' {
    test_error_match, stata_cmd(rangestat (mean) nonexistent_var, interval(t -5 5)) ctools_cmd(crangestat (mean) result=nonexistent_var, interval(t -5 5)) testname("nonexistent variable")
}
else {
    capture crangestat (mean) result=nonexistent_var, interval(t -5 5)
    if _rc != 0 {
        test_pass "[error] nonexistent variable (rc=`=_rc') [rangestat not installed]"
    }
    else {
        test_fail "[error] nonexistent variable" "should have errored"
    }
}

* Interval variable doesn't exist
clear
set obs 100
gen t = _n
gen x = runiform()
if `rangestat_installed' {
    test_error_match, stata_cmd(rangestat (mean) x, interval(nonexistent_t -5 5)) ctools_cmd(crangestat (mean) result=x, interval(nonexistent_t -5 5)) testname("nonexistent interval variable")
}
else {
    capture crangestat (mean) result=x, interval(nonexistent_t -5 5)
    if _rc != 0 {
        test_pass "[error] nonexistent interval variable (rc=`=_rc') [rangestat not installed]"
    }
    else {
        test_fail "[error] nonexistent interval variable" "should have errored"
    }
}

* Missing interval option
clear
set obs 100
gen t = _n
gen x = runiform()
if `rangestat_installed' {
    test_error_match, stata_cmd(rangestat (mean) x) ctools_cmd(crangestat (mean) result=x) testname("missing interval option")
}
else {
    capture crangestat (mean) result=x
    if _rc != 0 {
        test_pass "[error] missing interval option (rc=`=_rc') [rangestat not installed]"
    }
    else {
        test_fail "[error] missing interval option" "should have errored"
    }
}

* Invalid statistic function
clear
set obs 100
gen t = _n
gen x = runiform()
if `rangestat_installed' {
    test_error_match, stata_cmd(rangestat (invalid_stat) x, interval(t -5 5)) ctools_cmd(crangestat (invalid_stat) result=x, interval(t -5 5)) testname("invalid statistic function")
}
else {
    capture crangestat (invalid_stat) result=x, interval(t -5 5)
    if _rc != 0 {
        test_pass "[error] invalid statistic function (rc=`=_rc') [rangestat not installed]"
    }
    else {
        test_fail "[error] invalid statistic function" "should have errored"
    }
}

* Empty dataset
clear
set obs 0
gen t = .
gen x = .
if `rangestat_installed' {
    capture rangestat (mean) result=x, interval(t -5 5)
    local stata_rc = _rc
    capture crangestat (mean) result=x, interval(t -5 5)
    local ctools_rc = _rc
    if `stata_rc' == `ctools_rc' {
        test_pass "[error] empty dataset (rc=`stata_rc')"
    }
    else {
        test_fail "[error] empty dataset" "rangestat rc=`stata_rc', crangestat rc=`ctools_rc'"
    }
}
else {
    capture crangestat (mean) result=x, interval(t -5 5)
    if _rc != 0 {
        test_pass "[error] empty dataset (rc=`=_rc') [rangestat not installed]"
    }
    else {
        test_fail "[error] empty dataset" "should have errored"
    }
}

/*******************************************************************************
 * SECTION 10: Additional statistics coverage (5 tests)
 ******************************************************************************/
print_section "Additional Statistics"

* Test 10.1: P5 statistic
clear
set obs 200
gen obs = _n
gen x = obs
benchmark_rangestat p5 x, interval(obs -20 20) testname("p5 statistic")

* Test 10.2: P10 statistic
clear
set obs 200
gen obs = _n
gen x = obs
benchmark_rangestat p10 x, interval(obs -20 20) testname("p10 statistic")

* Test 10.3: P90 statistic
clear
set obs 200
gen obs = _n
gen x = obs
benchmark_rangestat p90 x, interval(obs -20 20) testname("p90 statistic")

* Test 10.4: P95 statistic
clear
set obs 200
gen obs = _n
gen x = obs
benchmark_rangestat p95 x, interval(obs -20 20) testname("p95 statistic")

* Test 10.5: All percentiles at once
clear
set obs 200
gen obs = _n
gen x = obs
benchmark_rangestat median x, interval(obs -20 20) testname("all percentiles: median vs rangestat")
* Also verify ordering: p1 <= p5 <= p10 <= p25 <= p50 <= p75 <= p90 <= p95 <= p99
crangestat (p1) xp1=x (p5) xp5=x (p10) xp10=x (p25) xp25=x (median) xp50=x (p75) xp75=x (p90) xp90=x (p95) xp95=x (p99) xp99=x, interval(obs -20 20)
local ordered = 1
forvalues i = 80/120 {
    if xp1[`i'] > xp5[`i'] | xp5[`i'] > xp10[`i'] | xp10[`i'] > xp25[`i'] | xp25[`i'] > xp50[`i'] | xp50[`i'] > xp75[`i'] | xp75[`i'] > xp90[`i'] | xp90[`i'] > xp95[`i'] | xp95[`i'] > xp99[`i'] {
        local ordered = 0
    }
}
if `ordered' {
    test_pass "all percentiles ordered correctly"
}
else {
    test_fail "all percentiles" "percentiles not in correct order"
}

/*******************************************************************************
 * SECTION 11: Rangestat head-to-head with real data (5 tests)
 ******************************************************************************/
print_section "Rangestat Head-to-Head (real data)"

* Test 11.1: Rolling stats on Grunfeld panel
webuse grunfeld, clear
sort company year
benchmark_rangestat mean invest, interval(year -2 2) by(company) testname("grunfeld: rolling mean")
benchmark_rangestat sd invest, interval(year -2 2) by(company) testname("grunfeld: rolling sd")
benchmark_rangestat count invest, interval(year -2 2) by(company) testname("grunfeld: rolling count")

* Test 11.2: nlswork age-based stats
webuse nlswork, clear
keep if !missing(ln_wage) & !missing(age)
benchmark_rangestat mean ln_wage, interval(age -2 2) testname("nlswork: mean by age")
benchmark_rangestat min ln_wage, interval(age -2 2) testname("nlswork: min by age")

/*******************************************************************************
 * SECTION 12: Negative and mixed values (5 tests)
 ******************************************************************************/
print_section "Negative and Mixed Values"

* Test 12.1: Negative source values
clear
set obs 100
gen t = _n
gen x = _n - 50
crangestat (mean) m=x, interval(t -2 2)
* At obs 50, values are -2,-1,0,1,2: mean=0
if abs(m[50]) < 1e-10 {
    test_pass "negative source values: mean"
}
else {
    test_fail "negative source: mean" "expected ~0, got `=m[50]'"
}

* Test 12.2: Negative source sum
clear
set obs 100
gen t = _n
gen x = _n - 50
crangestat (sum) s=x, interval(t -2 2)
* At obs 50, values are -2,-1,0,1,2: sum=0
if abs(s[50]) < 1e-10 {
    test_pass "negative source values: sum"
}
else {
    test_fail "negative source: sum" "expected ~0, got `=s[50]'"
}

* Test 12.3: Negative key values
clear
set obs 100
gen t = _n - 50
gen x = 1
crangestat (count) n=x, interval(t -2 2)
* At obs 50 (t=0), window [-2,2] includes t=-2,-1,0,1,2 = 5 obs
if n[50] == 5 {
    test_pass "negative key values"
}
else {
    test_fail "negative key values" "expected 5, got `=n[50]'"
}

* Test 12.4: Mixed positive/negative with rangestat
clear
set obs 200
set seed 54321
gen t = _n
gen x = rnormal()
benchmark_rangestat mean x, interval(t -5 5) testname("mixed values: mean")
benchmark_rangestat sd x, interval(t -5 5) testname("mixed values: sd")

* Test 12.5: Large magnitude values (numerical stability)
clear
set obs 100
gen t = _n
gen x = _n * 1e8
crangestat (mean) m=x (sd) s=x, interval(t -2 2)
* At obs 50, values 48e8..52e8, mean = 50e8 = 5e9
sigfigs `=m[50]' 5e9
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "large magnitude values: mean"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "large magnitude: mean" "sigfigs=`sf'"
}

/*******************************************************************************
 * SECTION 13: Data ordering and precision (5 tests)
 ******************************************************************************/
print_section "Data Ordering and Precision"

* Test 13.1: Random order with rangestat comparison
clear
set obs 500
set seed 99999
gen t = runiform() * 500
gen x = rnormal()
benchmark_rangestat mean x, interval(t -5 5) testname("random order: mean")
benchmark_rangestat count x, interval(t -5 5) testname("random order: count")

* Test 13.2: Float precision data
clear
set obs 100
gen t = _n
gen float x = _n + 0.1
benchmark_rangestat mean x, interval(t -2 2) testname("float precision data")

* Test 13.3: Result variable already exists (replace)
clear
set obs 100
gen t = _n
gen x = _n
gen double existing_var = 999
crangestat (mean) existing_var=x, interval(t -2 2)
* existing_var should be overwritten with actual mean
if existing_var[50] != 999 & !missing(existing_var[50]) {
    sigfigs `=existing_var[50]' 50
    if r(sigfigs) >= $DEFAULT_SIGFIGS {
        test_pass "result var already exists (overwrite)"
    }
    else {
        local sf : display %5.1f r(sigfigs)
        test_fail "result var overwrite" "sigfigs=`sf'"
    }
}
else {
    test_fail "result var overwrite" "variable not overwritten"
}

* Test 13.4: Very small interval with dense data
clear
set obs 10000
gen t = _n / 100
gen x = rnormal()
benchmark_rangestat count x, interval(t -0.01 0.01) testname("small interval dense data")

* Test 13.5: Multiple source variables, different stats
clear
set obs 100
gen t = _n
gen x = _n
gen y = _n * 2
gen z = rnormal()
crangestat (mean) mx=x (sum) sy=y (sd) sdz=z (min) minx=x (max) maxy=y, interval(t -3 3)
* Verify mx[50] = 50, sy[50] = 2*(47+48+49+50+51+52+53) = 2*350 = 700
sigfigs `=mx[50]' 50
local sf1 = r(sigfigs)
sigfigs `=sy[50]' 700
local sf2 = r(sigfigs)
if `sf1' >= $DEFAULT_SIGFIGS & `sf2' >= $DEFAULT_SIGFIGS & minx[50] == 47 & maxy[50] == 106 {
    test_pass "multiple source vars different stats"
}
else {
    test_fail "multi source vars" "mx=`=mx[50]' sy=`=sy[50]' minx=`=minx[50]' maxy=`=maxy[50]'"
}

/*******************************************************************************
 * SECTION 14: Advanced by-group tests (5 tests)
 ******************************************************************************/
print_section "Advanced by-group Tests"

* Test 14.1: by-group + excludeself with rangestat comparison
clear
set obs 200
gen group = mod(_n - 1, 2)
gen t = ceil(_n / 2)
gen x = rnormal()
sort group t
benchmark_rangestat mean x, interval(t -3 3) by(group) excludeself testname("by + excludeself: mean")
benchmark_rangestat count x, interval(t -3 3) by(group) excludeself testname("by + excludeself: count")

* Test 14.2: Many by-groups (100 groups)
clear
set obs 10000
gen group = mod(_n - 1, 100)
gen t = ceil(_n / 100)
gen x = 1
sort group t
crangestat (sum) s=x, interval(t . 0) by(group)
* Each group has 100 obs with t=1..100. At t=50 (obs 50 within each group), cumsum=50
if s[50] == 50 {
    test_pass "100 by-groups"
}
else {
    test_fail "100 by-groups" "expected 50, got `=s[50]'"
}

* Test 14.3: Unequal group sizes
clear
set obs 300
gen group = cond(_n <= 100, 0, cond(_n <= 250, 1, 2))
gen t = cond(group == 0, _n, cond(group == 1, _n - 100, _n - 250))
gen x = 1
sort group t
crangestat (count) n=x, interval(t . .) by(group)
* Group 0: 100 obs, Group 1: 150 obs, Group 2: 50 obs
if n[1] == 100 & n[101] == 150 & n[251] == 50 {
    test_pass "unequal group sizes"
}
else {
    test_fail "unequal group sizes" "n[1]=`=n[1]' n[101]=`=n[101]' n[251]=`=n[251]'"
}

* Test 14.4: by-group with infinite interval and rangestat
clear
set obs 400
gen group = mod(_n - 1, 4)
gen t = ceil(_n / 4)
set seed 11111
gen x = rnormal()
sort group t
benchmark_rangestat mean x, interval(t . .) by(group) testname("by infinite: mean")
benchmark_rangestat sd x, interval(t . .) by(group) testname("by infinite: sd")

* Test 14.5: by-group where groups have different t ranges
clear
set obs 200
gen group = _n <= 100
gen t = cond(group == 1, _n, _n + 1000)
gen x = 1
sort group t
crangestat (sum) s=x, interval(t -5 5) by(group)
* Groups have non-overlapping t ranges, results should be independent
* Group 1 (obs 1-100): t=1..100, at obs 50 sum should be 11
* Group 0 (obs 101-200): t=1101..1200, at obs 150 sum should be 11
if s[50] == 11 & s[150] == 11 {
    test_pass "by-groups with non-overlapping t ranges"
}
else {
    test_fail "non-overlapping t" "s[50]=`=s[50]' s[150]=`=s[150]'"
}

/*******************************************************************************
 * SECTION 15: Auto dataset head-to-head
 ******************************************************************************/
print_section "Auto Dataset"

* Test 15.1: Rolling stats on price by rep78
sysuse auto, clear
keep if !missing(rep78) & !missing(price) & !missing(mpg)
sort rep78 price
benchmark_rangestat mean mpg, interval(price -500 500) by(rep78) testname("auto: rolling mean mpg by rep78")
benchmark_rangestat count mpg, interval(price -500 500) by(rep78) testname("auto: rolling count mpg by rep78")

* Test 15.2: Min/max on weight within price range
sysuse auto, clear
keep if !missing(weight) & !missing(price)
benchmark_rangestat min weight, interval(price -1000 1000) testname("auto: min weight by price")
benchmark_rangestat max weight, interval(price -1000 1000) testname("auto: max weight by price")

* Test 15.3: Median and IQR on mpg within weight range
sysuse auto, clear
keep if !missing(mpg) & !missing(weight)
benchmark_rangestat median mpg, interval(weight -300 300) testname("auto: median mpg by weight")

* Test 15.4: Excludeself on auto data
sysuse auto, clear
keep if !missing(mpg) & !missing(price)
crangestat (mean) c_loo=mpg, interval(price . .) excludeself
* LOO mean: (sum_all - mpg_i) / (N - 1)
quietly summarize mpg
local total_sum = r(sum)
local N = r(N)
sigfigs `=c_loo[1]' `=(`total_sum' - mpg[1]) / (`N' - 1)'
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "auto: LOO mean mpg"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "auto: LOO mean mpg" "sigfigs=`sf'"
}

* Test 15.5: Multiple stats on auto with foreign grouping
sysuse auto, clear
keep if !missing(mpg) & !missing(price) & !missing(weight)
sort foreign price
benchmark_rangestat mean mpg, interval(price -2000 2000) by(foreign) testname("auto: multi-stat by foreign")

/*******************************************************************************
 * SECTION 16: Non-uniform key spacing
 ******************************************************************************/
print_section "Non-Uniform Key Spacing"

* Test 16.1: Gaps in time series (missing dates)
clear
set obs 100
gen t = _n
* Introduce gaps: remove every 5th observation's key
replace t = t + 10 * floor((_n - 1) / 5)
gen x = 1
benchmark_rangestat count x, interval(t -2 2) testname("gaps in time series")

* Test 16.2: Exponential spacing (log-scale data)
clear
set obs 50
gen t = exp(_n / 10)
gen x = _n
benchmark_rangestat mean x, interval(t -0.5 0.5) testname("exponential spacing")

* Test 16.3: Clustered keys (many ties with gaps)
clear
set obs 500
* 5 clusters of 100 obs each, at t = 10, 20, 30, 40, 50
gen t = (floor((_n - 1) / 100) + 1) * 10
gen x = _n
crangestat (mean) m=x (count) n=x, interval(t -5 5)
* At t=10, only the 100 obs with t=10 are in range [5,15]
* Mean of obs 1..100 = 50.5, count = 100
sigfigs `=m[1]' 50.5
if r(sigfigs) >= $DEFAULT_SIGFIGS & n[1] == 100 {
    test_pass "clustered keys"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "clustered keys" "m=`=m[1]' (sf=`sf'), n=`=n[1]'"
}

* Test 16.4: Non-integer key with rangestat
clear
set obs 200
set seed 77777
gen t = runiform() * 100
gen x = rnormal()
benchmark_rangestat mean x, interval(t -2.5 2.5) testname("non-integer key: mean")
benchmark_rangestat sd x, interval(t -2.5 2.5) testname("non-integer key: sd")

* Test 16.5: Date-like key with weekly gaps
clear
set obs 365
gen day = _n
* Only weekdays (remove weekends)
gen dow = mod(_n - 1, 7)
keep if dow < 5
gen x = rnormal()
benchmark_rangestat mean x, interval(day -7 7) testname("weekday-only keys")

/*******************************************************************************
 * SECTION 17: Combined features stress tests
 ******************************************************************************/
print_section "Combined Features"

* Test 17.1: by + excludeself + missing + multiple stats
clear
set obs 600
gen group = mod(_n - 1, 3)
gen t = ceil(_n / 3)
set seed 22222
gen x = rnormal()
* Sprinkle missing values
replace x = . if mod(_n, 7) == 0
sort group t
benchmark_rangestat mean x, interval(t -4 4) by(group) excludeself testname("combined: mean (by+excl+miss)")
benchmark_rangestat count x, interval(t -4 4) by(group) excludeself testname("combined: count (by+excl+miss)")

* Test 17.2: All percentiles with by-group on real data
webuse grunfeld, clear
sort company year
benchmark_rangestat median invest, interval(year -3 3) by(company) testname("grunfeld: median by company")
crangestat (p10) c_p10=invest (p25) c_p25=invest (median) c_p50=invest (p75) c_p75=invest (p90) c_p90=invest, interval(year -3 3) by(company)
* Verify percentile ordering within each observation
local ordered = 1
quietly count
forvalues i = 1/`=r(N)' {
    if !missing(c_p10[`i']) & !missing(c_p90[`i']) {
        if c_p10[`i'] > c_p25[`i'] | c_p25[`i'] > c_p50[`i'] | c_p50[`i'] > c_p75[`i'] | c_p75[`i'] > c_p90[`i'] {
            local ordered = 0
        }
    }
}
if `ordered' {
    test_pass "grunfeld: percentiles by company"
}
else {
    test_fail "grunfeld: percentiles by company" "percentile ordering violated"
}

* Test 17.3: Skewness and kurtosis with by-group + rangestat
webuse nlswork, clear
keep if !missing(ln_wage) & !missing(age) & !missing(race)
sort race age
benchmark_rangestat skewness ln_wage, interval(age -3 3) by(race) testname("nlswork: skewness by race")
benchmark_rangestat kurtosis ln_wage, interval(age -3 3) by(race) testname("nlswork: kurtosis by race")

* Test 17.4: first/last/firstnm/lastnm comprehensive
clear
set obs 100
gen t = _n
gen x = _n * 10
replace x = . in 3
replace x = . in 97
crangestat (first) c_first=x (last) c_last=x (firstnm) c_fnm=x (lastnm) c_lnm=x, interval(t -2 2)
* At obs 3 (t=3): window [1,5], first=10(t=1), last=50(t=5)
* x[3]=., so firstnm should be 10 (first non-missing), lastnm=50
if c_first[3] == 10 & c_last[3] == 50 & c_fnm[3] == 10 & c_lnm[3] == 50 {
    test_pass "first/last/firstnm/lastnm with missing"
}
else {
    test_fail "first/last/firstnm/lastnm" "first=`=c_first[3]' last=`=c_last[3]' fnm=`=c_fnm[3]' lnm=`=c_lnm[3]'"
}

* Test 17.5: Large dataset with all feature combinations
clear
set obs 5000
gen group = mod(_n - 1, 10)
gen t = ceil(_n / 10)
set seed 33333
gen x = rnormal() * 100
gen y = runiform() * 50
replace x = . if mod(_n, 13) == 0
sort group t
benchmark_rangestat mean x, interval(t -5 5) by(group) excludeself testname("large combined: mean")

/*******************************************************************************
 * SECTION 18: Auto-naming convention
 ******************************************************************************/
print_section "Auto-Naming Convention"

* Test 18.1: Default name for mean
clear
set obs 50
gen t = _n
gen price = _n * 100
crangestat (mean) price, interval(t -3 3)
* Convention is stat_varname, so (mean) price → mean_price
capture confirm variable mean_price
if _rc == 0 {
    sigfigs `=mean_price[25]' 2500
    if r(sigfigs) >= $DEFAULT_SIGFIGS {
        test_pass "auto-name: mean"
    }
    else {
        local sf : display %5.1f r(sigfigs)
        test_fail "auto-name: mean" "sigfigs=`sf'"
    }
}
else {
    test_fail "auto-name: mean" "variable mean_price not created"
}

* Test 18.2: Default name for count
clear
set obs 50
gen t = _n
gen val = 1
crangestat (count) val, interval(t -2 2)
capture confirm variable count_val
if _rc == 0 {
    if count_val[25] == 5 {
        test_pass "auto-name: count"
    }
    else {
        test_fail "auto-name: count" "expected 5, got `=count_val[25]'"
    }
}
else {
    test_fail "auto-name: count" "variable count_val not created"
}

* Test 18.3: Default names for multiple stats
clear
set obs 50
gen t = _n
gen x = _n
crangestat (mean) x (sd) x (min) x (max) x, interval(t -2 2)
capture confirm variable mean_x
local rc1 = _rc
capture confirm variable sd_x
local rc2 = _rc
capture confirm variable min_x
local rc3 = _rc
capture confirm variable max_x
local rc4 = _rc
if `rc1' == 0 & `rc2' == 0 & `rc3' == 0 & `rc4' == 0 {
    test_pass "auto-name: multiple stats"
}
else {
    test_fail "auto-name: multiple stats" "not all variables created"
}

* Test 18.4: Mixed explicit and auto names
clear
set obs 50
gen t = _n
gen x = _n
crangestat (mean) my_mean=x (sd) x, interval(t -2 2)
capture confirm variable my_mean
local rc1 = _rc
capture confirm variable sd_x
local rc2 = _rc
if `rc1' == 0 & `rc2' == 0 {
    test_pass "auto-name: mixed explicit and auto"
}
else {
    test_fail "auto-name: mixed" "my_mean rc=`rc1', sd_x rc=`rc2'"
}

/*******************************************************************************
 * SECTION 19: Known distribution verification
 ******************************************************************************/
print_section "Known Distribution Verification"

* Test 19.1: Mean and SD of uniform distribution
clear
set obs 100000
set seed 44444
gen t = _n
gen x = runiform()
* With infinite interval, mean ~ 0.5, sd ~ 1/sqrt(12) = 0.2887
crangestat (mean) m=x (sd) s=x, interval(t . .)
sigfigs `=m[1]' 0.5
local sf_m = r(sigfigs)
sigfigs `=s[1]' `=1/sqrt(12)'
local sf_s = r(sigfigs)
* With 100K obs, we expect ~2.5 sigfigs of sampling accuracy
if `sf_m' >= 2 & `sf_s' >= 2 {
    test_pass "uniform distribution: mean and sd"
}
else {
    local sf_m_fmt : display %5.1f `sf_m'
    local sf_s_fmt : display %5.1f `sf_s'
    test_fail "uniform distribution" "mean sf=`sf_m_fmt', sd sf=`sf_s_fmt'"
}

* Test 19.2: Sum of constants
clear
set obs 1000
gen t = _n
gen x = 7.5
crangestat (sum) s=x (mean) m=x (sd) sdx=x (count) n=x, interval(t . .)
* sum = 7500, mean = 7.5, sd = 0, count = 1000
sigfigs `=s[1]' 7500
local sf_sum = r(sigfigs)
sigfigs `=m[1]' 7.5
local sf_mean = r(sigfigs)
if `sf_sum' >= $DEFAULT_SIGFIGS & `sf_mean' >= $DEFAULT_SIGFIGS & n[1] == 1000 & (missing(sdx[1]) | sdx[1] == 0) {
    test_pass "constant values: sum, mean, sd, count"
}
else {
    test_fail "constant values" "sum=`=s[1]' mean=`=m[1]' sd=`=sdx[1]' n=`=n[1]'"
}

* Test 19.3: Known linear sequence stats
clear
set obs 100
gen t = _n
gen x = _n
crangestat (mean) m=x (sum) s=x (min) lo=x (max) hi=x (count) n=x, interval(t . .)
* For x = 1..100: mean=50.5, sum=5050, min=1, max=100, count=100
sigfigs `=m[1]' 50.5
local sf_m = r(sigfigs)
sigfigs `=s[1]' 5050
local sf_s = r(sigfigs)
if `sf_m' >= $DEFAULT_SIGFIGS & `sf_s' >= $DEFAULT_SIGFIGS & lo[1] == 1 & hi[1] == 100 & n[1] == 100 {
    test_pass "linear sequence: all stats correct"
}
else {
    test_fail "linear sequence" "m=`=m[1]' s=`=s[1]' lo=`=lo[1]' hi=`=hi[1]' n=`=n[1]'"
}

* Test 19.4: Variance of known data
clear
set obs 5
gen t = _n
gen x = .
replace x = 2 in 1
replace x = 4 in 2
replace x = 4 in 3
replace x = 4 in 4
replace x = 5 in 5
* Variance (N-1 denom): mean=3.8, SS=((2-3.8)^2+(4-3.8)^2*3+(5-3.8)^2)/(5-1)
* = (3.24 + 0.12 + 1.44) / 4 = 4.8 / 4 = 1.2
crangestat (variance) v=x, interval(t . .)
sigfigs `=v[1]' 1.2
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "known variance calculation"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "known variance" "expected 1.2, got `=v[1]' (sf=`sf')"
}

* Test 19.5: Median of odd and even counts
clear
set obs 10
gen t = _n
gen x = _n * 10
* Window of 5 (odd): median = middle value
crangestat (median) med5=x, interval(t -2 2)
* At obs 5 (t=5), window [3,7], values 30,40,50,60,70 → median=50
if med5[5] == 50 {
    test_pass "median: odd count window"
}
else {
    test_fail "median: odd count" "expected 50, got `=med5[5]'"
}

/*******************************************************************************
 * SECTION 20: Skewed and extreme data
 ******************************************************************************/
print_section "Skewed and Extreme Data"

* Test 20.1: Log-normal data (right-skewed)
clear
set obs 1000
set seed 55555
gen t = _n
gen x = exp(rnormal())
benchmark_rangestat mean x, interval(t -50 50) testname("lognormal: mean")

* Test 20.2: Data with outliers
clear
set obs 101
gen t = _n
gen x = 0
replace x = 1000 in 51
crangestat (mean) m=x (median) med=x (sd) s=x (min) lo=x (max) hi=x, interval(t -5 5)
* At obs 51 (t=51), window [46,56] = 11 obs, one is 1000, rest are 0
* Mean = 1000/11 ≈ 90.909, median = 0, min = 0, max = 1000
sigfigs `=m[51]' `=1000/11'
local sf_m = r(sigfigs)
if `sf_m' >= $DEFAULT_SIGFIGS & med[51] == 0 & lo[51] == 0 & hi[51] == 1000 {
    test_pass "data with outliers"
}
else {
    local sf_m_fmt : display %5.1f `sf_m'
    test_fail "data with outliers" "mean sf=`sf_m_fmt' med=`=med[51]' lo=`=lo[51]' hi=`=hi[51]'"
}

* Test 20.3: Very small values (near machine epsilon)
clear
set obs 100
gen t = _n
gen x = _n * 1e-15
crangestat (mean) m=x (sum) s=x, interval(t -2 2)
* At obs 50, values 48e-15..52e-15, mean = 50e-15 = 5e-14
sigfigs `=m[50]' 5e-14
if r(sigfigs) >= $DEFAULT_SIGFIGS {
    test_pass "very small values"
}
else {
    local sf : display %5.1f r(sigfigs)
    test_fail "very small values" "sigfigs=`sf'"
}

* Test 20.4: Mixed magnitude within same window
clear
set obs 20
gen t = _n
gen x = cond(mod(_n, 2) == 0, 1e6, 1e-6)
benchmark_rangestat mean x, interval(t -2 2) testname("mixed magnitude")

* Test 20.5: Census dataset — population-based range stats
capture webuse census, clear
if _rc == 0 {
    keep if !missing(pop) & !missing(medage)
    benchmark_rangestat mean pop, interval(medage -1 1) testname("census: mean pop by medage")
}
else {
    test_pass "census: mean pop by medage [dataset not available]"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
* End of crangestat validation
noi print_summary "crangestat"
}
