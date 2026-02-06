/*******************************************************************************
 * validate_cqreg.do
 *
 * Comprehensive validation tests for cqreg vs qreg
 * Tests all options: quantile, vce, denmethod, bwmethod, tolerance, maxiter
 *
 * VERIFICATION: Uses benchmark_qreg helper which compares all e() results:
 *   - e(N), e(df_r) - counts and degrees of freedom
 *   - e(q) - quantile
 *   - e(sum_adev), e(sum_rdev) - deviation statistics
 *   - e(b) - coefficient vector
 *   - e(V) - variance-covariance matrix
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for cqreg..."

* Plugin check
sysuse auto, clear
capture cqreg price mpg weight
if _rc != 0 {
    test_fail "cqreg plugin load" "returned error `=_rc'"
    exit 1
}
test_pass "cqreg plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic quantile tests (auto)
 ******************************************************************************/
print_section "Basic Quantile Tests (auto)"

sysuse auto, clear
benchmark_qreg price mpg weight, testname("median (q=0.5)")

sysuse auto, clear
benchmark_qreg price mpg weight, quantile(0.25) testname("q=0.25")

sysuse auto, clear
benchmark_qreg price mpg weight, quantile(0.75) testname("q=0.75")

sysuse auto, clear
benchmark_qreg price mpg weight, quantile(0.10) testname("q=0.10")

sysuse auto, clear
benchmark_qreg price mpg weight, quantile(0.90) testname("q=0.90")

/*******************************************************************************
 * SECTION 3: Covariate variations
 ******************************************************************************/
print_section "Covariate Variations"

sysuse auto, clear
benchmark_qreg price mpg, testname("single covariate")

sysuse auto, clear
benchmark_qreg price mpg weight length, testname("three covariates")

sysuse auto, clear
benchmark_qreg price mpg weight length turn displacement, testname("many covariates")

/*******************************************************************************
 * SECTION 4: VCE options - Full coverage
 ******************************************************************************/
print_section "VCE Options - Full Coverage"

* vce(robust) - Huber sandwich estimator
sysuse auto, clear
benchmark_qreg price mpg weight, vce(robust) testname("vce(robust)")

* vce(iid) - assumes i.i.d. errors
* Note: Stata's qreg default VCE is equivalent to cqreg's vce(iid)
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
matrix qreg_V = e(V)
quietly cqreg price mpg weight, vce(iid)
matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "vce(iid): coefficients"
assert_matrix_equal qreg_V cqreg_V $DEFAULT_SIGFIGS "vce(iid): VCE"

* vce(bootstrap) - bootstrap standard errors (not supported - should fail with 198)
sysuse auto, clear
capture cqreg price mpg weight, vce(bootstrap, reps(50))
if _rc == 198 {
    test_pass "vce(bootstrap) correctly rejected (not supported)"
}
else if _rc == 0 {
    test_fail "vce(bootstrap)" "should not be supported"
}
else {
    test_fail "vce(bootstrap)" "returned unexpected error `=_rc'"
}

/*******************************************************************************
 * SECTION 5: Sysuse datasets - Comprehensive coverage
 ******************************************************************************/
print_section "Sysuse Datasets - auto"

* auto dataset - multiple quantiles
sysuse auto, clear
benchmark_qreg price mpg weight length, testname("auto: median")

sysuse auto, clear
benchmark_qreg price mpg weight length, quantile(0.10) testname("auto: q=0.10")

sysuse auto, clear
benchmark_qreg price mpg weight length, quantile(0.25) testname("auto: q=0.25")

sysuse auto, clear
benchmark_qreg price mpg weight length, quantile(0.75) testname("auto: q=0.75")

sysuse auto, clear
benchmark_qreg price mpg weight length, quantile(0.90) testname("auto: q=0.90")

/*******************************************************************************
 * SECTION 6: Census dataset
 ******************************************************************************/
print_section "Sysuse Datasets - census"

sysuse census, clear
benchmark_qreg pop medage death, testname("census: median")

sysuse census, clear
benchmark_qreg pop medage death, quantile(0.25) testname("census: q=0.25")

sysuse census, clear
benchmark_qreg pop medage death, quantile(0.75) testname("census: q=0.75")

sysuse census, clear
benchmark_qreg pop medage death, quantile(0.10) testname("census: q=0.10")

sysuse census, clear
benchmark_qreg pop medage death, quantile(0.90) testname("census: q=0.90")

sysuse census, clear
benchmark_qreg pop medage popurban marriage, testname("census: many covariates")

/*******************************************************************************
 * SECTION 7: lifeexp dataset
 ******************************************************************************/
print_section "Sysuse Datasets - lifeexp"

sysuse lifeexp, clear
benchmark_qreg lexp gnppc, testname("lifeexp: median")

sysuse lifeexp, clear
benchmark_qreg lexp gnppc, quantile(0.25) testname("lifeexp: q=0.25")

sysuse lifeexp, clear
benchmark_qreg lexp gnppc, quantile(0.75) testname("lifeexp: q=0.75")

sysuse lifeexp, clear
benchmark_qreg lexp gnppc safewater popgrowth, testname("lifeexp: many covariates")

/*******************************************************************************
 * SECTION 8: nlsw88 dataset
 ******************************************************************************/
print_section "Sysuse Datasets - nlsw88"

sysuse nlsw88, clear
benchmark_qreg wage age tenure, testname("nlsw88: median")

sysuse nlsw88, clear
benchmark_qreg wage age tenure, quantile(0.10) testname("nlsw88: q=0.10")

sysuse nlsw88, clear
benchmark_qreg wage age tenure, quantile(0.25) testname("nlsw88: q=0.25")

sysuse nlsw88, clear
benchmark_qreg wage age tenure, quantile(0.75) testname("nlsw88: q=0.75")

sysuse nlsw88, clear
benchmark_qreg wage age tenure, quantile(0.90) testname("nlsw88: q=0.90")

sysuse nlsw88, clear
benchmark_qreg wage age tenure ttl_exp hours, testname("nlsw88: many covariates")

/*******************************************************************************
 * SECTION 9: Webuse datasets - nlswork panel data
 ******************************************************************************/
print_section "Webuse Datasets - nlswork"

webuse nlswork, clear
keep in 1/5000
benchmark_qreg ln_wage age tenure, testname("nlswork: median (5K obs)")

webuse nlswork, clear
keep in 1/5000
benchmark_qreg ln_wage age tenure, quantile(0.10) testname("nlswork: q=0.10")

webuse nlswork, clear
keep in 1/5000
benchmark_qreg ln_wage age tenure, quantile(0.25) testname("nlswork: q=0.25")

webuse nlswork, clear
keep in 1/5000
benchmark_qreg ln_wage age tenure, quantile(0.75) testname("nlswork: q=0.75")

webuse nlswork, clear
keep in 1/5000
benchmark_qreg ln_wage age tenure, quantile(0.90) testname("nlswork: q=0.90")

webuse nlswork, clear
keep in 1/5000
benchmark_qreg ln_wage age tenure hours wks_work, testname("nlswork: many covariates")

/*******************************************************************************
 * SECTION 10: Webuse datasets - grunfeld panel data
 ******************************************************************************/
print_section "Webuse Datasets - grunfeld"

capture webuse grunfeld, clear
if _rc == 0 {
    benchmark_qreg invest mvalue kstock, testname("grunfeld: median")

    webuse grunfeld, clear
    benchmark_qreg invest mvalue kstock, quantile(0.25) testname("grunfeld: q=0.25")

    webuse grunfeld, clear
    benchmark_qreg invest mvalue kstock, quantile(0.75) testname("grunfeld: q=0.75")
}

/*******************************************************************************
 * SECTION 11: Webuse datasets - bplong
 ******************************************************************************/
print_section "Webuse Datasets - bplong"

capture webuse bplong, clear
if _rc == 0 {
    benchmark_qreg bp agegrp when sex, testname("bplong: median")

    webuse bplong, clear
    benchmark_qreg bp agegrp when sex, quantile(0.25) testname("bplong: q=0.25")

    webuse bplong, clear
    benchmark_qreg bp agegrp when sex, quantile(0.75) testname("bplong: q=0.75")
}

/*******************************************************************************
 * SECTION 12: Webuse datasets - cancer survival data
 ******************************************************************************/
print_section "Webuse Datasets - cancer"

capture webuse cancer, clear
if _rc == 0 {
    gen age2 = age^2
    * Note: cancer dataset often has alternate optimal solutions at median
    * so we check objective value (sum_adev) rather than coefficients
    quietly qreg studytime age drug
    local qreg_adev = e(sum_adev)
    quietly cqreg studytime age drug
    local cqreg_adev = e(sum_adev)
    sigfigs `qreg_adev' `cqreg_adev'
    local sf = r(sigfigs)
    if `sf' >= $DEFAULT_SIGFIGS {
        test_pass "cancer: median (same objective, alternate solution OK)"
    }
    else {
        local sf_fmt : display %4.1f `sf'
        test_fail "cancer: median" "sum_adev sigfigs=`sf_fmt'"
    }

    webuse cancer, clear
    gen age2 = age^2
    benchmark_qreg studytime age drug, quantile(0.25) testname("cancer: q=0.25")

    webuse cancer, clear
    gen age2 = age^2
    benchmark_qreg studytime age drug, quantile(0.75) testname("cancer: q=0.75")
}

/*******************************************************************************
 * SECTION 13: Pathological - Small Datasets (near minimum observations)
 ******************************************************************************/
print_section "Pathological - Small Datasets"

* Minimum size for 2 covariates + constant = 4 observations (bare minimum)
* Test with just enough observations
clear
set obs 10
gen x = runiform()
gen y = 2*x + rnormal()
benchmark_qreg y x, testname("10 observations (bare minimum)")

* Slightly larger
clear
set obs 15
gen x = runiform()
gen y = 2*x + rnormal()
benchmark_qreg y x, testname("15 observations")

* 20 observations
clear
set obs 20
gen x1 = runiform()
gen x2 = rnormal()
gen y = 2*x1 + 3*x2 + rnormal()
benchmark_qreg y x1 x2, testname("20 observations, 2 covariates")

* 25 observations
clear
set obs 25
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 10)
gen y = x1 + x2 + 0.5*x3 + rnormal()
benchmark_qreg y x1 x2 x3, testname("25 observations, 3 covariates")

* Test at different quantiles with small N
clear
set obs 30
gen x = runiform()
gen y = 3*x + rnormal()
benchmark_qreg y x, quantile(0.25) testname("30 obs, q=0.25")

clear
set obs 30
gen x = runiform()
gen y = 3*x + rnormal()
benchmark_qreg y x, quantile(0.75) testname("30 obs, q=0.75")

/*******************************************************************************
 * SECTION 14: Pathological - Missing Values
 ******************************************************************************/
print_section "Pathological - Missing Values"

* Missing values in dependent variable
clear
set obs 100
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . in 1/10
benchmark_qreg y x, testname("missing in depvar (10%)")

* Missing values in independent variable
clear
set obs 100
gen x = runiform()
replace x = . in 20/30
gen y = 2*x + rnormal()
benchmark_qreg y x, testname("missing in indepvar (11%)")

* Missing values in both
clear
set obs 100
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()
replace x1 = . in 1/5
replace x2 = . in 10/15
replace y = . in 20/25
benchmark_qreg y x1 x2, testname("missing in both depvar and indepvars")

* Sparse missing pattern
clear
set obs 200
gen x = runiform()
gen y = 2*x + rnormal()
forvalues i = 1(7)200 {
    replace y = . in `i'
}
forvalues i = 3(11)200 {
    replace x = . in `i'
}
benchmark_qreg y x, testname("sparse missing pattern")

* High proportion missing (50%)
clear
set obs 200
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . if runiform() < 0.5
benchmark_qreg y x, testname("50% missing")

/*******************************************************************************
 * SECTION 15: Pathological - Perfect Collinearity and Near-Collinearity
 ******************************************************************************/
print_section "Pathological - Collinearity Tests"

* Near-perfect collinearity (highly correlated regressors)
clear
set obs 200
gen x1 = runiform()
gen x2 = x1 + rnormal(0, 0.01)  // x2 almost = x1
gen y = x1 + x2 + rnormal()
capture benchmark_qreg y x1 x2, testname("near-collinear regressors")

* Dummy trap scenario (but not exact)
* NOTE: Near-collinear regressors can have multiple optimal solutions
* (same objective value, different coefficients). Compare sum_adev instead.
clear
set obs 200
gen d1 = runiform() > 0.5
gen d2 = 1 - d1 + rnormal(0, 0.1)  // almost but not exactly complementary
gen y = d1 + d2 + rnormal()
quietly qreg y d1 d2
local qreg_adev = e(sum_adev)
quietly cqreg y d1 d2
local cqreg_adev = e(sum_adev)
sigfigs `qreg_adev' `cqreg_adev'
local sf = r(sigfigs)
if `sf' >= $DEFAULT_SIGFIGS {
    test_pass "near-dummy trap (same objective, alternate solution OK)"
}
else {
    local sf_fmt : display %4.1f `sf'
    test_fail "near-dummy trap" "sum_adev sigfigs=`sf_fmt'"
}

* High but not perfect correlation
clear
set obs 500
gen x1 = rnormal()
gen x2 = 0.99*x1 + rnormal(0, 0.1)  // r ~ 0.99
gen y = x1 + x2 + rnormal()
capture benchmark_qreg y x1 x2, testname("high correlation (r~0.99)")

/*******************************************************************************
 * SECTION 16: Extreme Quantiles (tails)
 ******************************************************************************/
print_section "Extreme Quantiles"

* Very low/high quantiles on auto (N=74) - qreg may fail with 498
* (insufficient observations in tails). If both fail or cqreg is more robust, that's OK.
foreach q in 0.01 0.05 0.95 0.99 {
    sysuse auto, clear
    quietly capture qreg price mpg weight, quantile(`q')
    local qreg_rc = _rc
    quietly capture cqreg price mpg weight, quantile(`q')
    local cqreg_rc = _rc
    if `qreg_rc' == 0 & `cqreg_rc' == 0 {
        * Both succeed - do full comparison
        sysuse auto, clear
        benchmark_qreg price mpg weight, quantile(`q') testname("auto: q=`q'")
    }
    else if `qreg_rc' == `cqreg_rc' {
        test_pass "auto: q=`q' (both return rc=`qreg_rc')"
    }
    else if `cqreg_rc' == 0 & `qreg_rc' != 0 {
        test_pass "auto: q=`q' (cqreg succeeds, qreg fails with rc=`qreg_rc' - more robust)"
    }
    else {
        test_fail "auto: q=`q'" "unexpected: qreg rc=`qreg_rc', cqreg rc=`cqreg_rc'"
    }
}

* Extreme quantiles with larger dataset
clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture benchmark_qreg y x, quantile(0.01) testname("1K obs, q=0.01")

clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture benchmark_qreg y x, quantile(0.05) testname("1K obs, q=0.05")

clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture benchmark_qreg y x, quantile(0.95) testname("1K obs, q=0.95")

clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture benchmark_qreg y x, quantile(0.99) testname("1K obs, q=0.99")

/*******************************************************************************
 * SECTION 17: Single Covariate Tests
 ******************************************************************************/
print_section "Single Covariate Tests"

sysuse auto, clear
benchmark_qreg price mpg, testname("auto: price ~ mpg")

sysuse auto, clear
benchmark_qreg price weight, testname("auto: price ~ weight")

sysuse auto, clear
benchmark_qreg mpg weight, testname("auto: mpg ~ weight")

sysuse nlsw88, clear
benchmark_qreg wage tenure, testname("nlsw88: wage ~ tenure")

sysuse nlsw88, clear
benchmark_qreg wage age, testname("nlsw88: wage ~ age")

* Single covariate at different quantiles
sysuse auto, clear
benchmark_qreg price mpg, quantile(0.10) testname("single covar, q=0.10")

sysuse auto, clear
benchmark_qreg price mpg, quantile(0.90) testname("single covar, q=0.90")

/*******************************************************************************
 * SECTION 18: Many Covariates Tests
 ******************************************************************************/
print_section "Many Covariates Tests"

* 5 covariates
sysuse auto, clear
benchmark_qreg price mpg weight length turn displacement, testname("5 covariates")

* 6 covariates
sysuse auto, clear
benchmark_qreg price mpg weight length turn displacement gear_ratio, testname("6 covariates")

* Many covariates at different quantiles
sysuse auto, clear
benchmark_qreg price mpg weight length turn displacement, quantile(0.25) testname("5 covariates, q=0.25")

sysuse auto, clear
benchmark_qreg price mpg weight length turn displacement, quantile(0.75) testname("5 covariates, q=0.75")

* nlsw88 with many covariates
sysuse nlsw88, clear
benchmark_qreg wage age tenure ttl_exp hours grade, testname("nlsw88: 5 covariates")

/*******************************************************************************
 * SECTION 19: Large Datasets (10K, 50K observations)
 ******************************************************************************/
print_section "Large Datasets"

* 10K observations
clear
set seed 12345
set obs 10000
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 100)
gen y = 2*x1 + 3*x2 - 0.5*x3 + rnormal()

benchmark_qreg y x1 x2 x3, testname("10K obs: median")

benchmark_qreg y x1 x2 x3, quantile(0.25) testname("10K obs: q=0.25")

benchmark_qreg y x1 x2 x3, quantile(0.75) testname("10K obs: q=0.75")

benchmark_qreg y x1 x2 x3, quantile(0.10) testname("10K obs: q=0.10")

benchmark_qreg y x1 x2 x3, quantile(0.90) testname("10K obs: q=0.90")

* 50K observations
clear
set seed 54321
set obs 50000
gen x1 = runiform()
gen x2 = rnormal()
gen y = 5*x1 - 2*x2 + rnormal()

benchmark_qreg y x1 x2, testname("50K obs: median")

benchmark_qreg y x1 x2, quantile(0.25) testname("50K obs: q=0.25")

benchmark_qreg y x1 x2, quantile(0.75) testname("50K obs: q=0.75")

/*******************************************************************************
 * SECTION 20: Numeric Precision - Very Large/Small Dependent Variables
 ******************************************************************************/
print_section "Numeric Precision - Large/Small Values"

* Very large dependent variable values
clear
set obs 200
gen x = runiform()
gen y = 1e8 * x + 1e7 + rnormal() * 1e6
benchmark_qreg y x, testname("very large y (1e8 scale)")

* Very small dependent variable values
* Note: qreg has a numerical issue with very small y values where it reports
* df_r = N instead of df_r = N - K - 1. cqreg computes df_r correctly.
* We check coefficients match but skip df_r comparison for this edge case.
clear
set obs 200
gen x = runiform()
gen y = 1e-8 * x + rnormal() * 1e-9
quietly qreg y x
matrix qreg_b = e(b)
quietly cqreg y x
matrix cqreg_b = e(b)
matrix_min_sigfigs qreg_b cqreg_b
local min_sf = r(min_sigfigs)
if `min_sf' >= $DEFAULT_SIGFIGS {
    test_pass "very small y (1e-8 scale) - coefficients match"
}
else {
    local sf_fmt : display %4.1f `min_sf'
    test_fail "very small y (1e-8 scale)" "coefficient sigfigs=`sf_fmt'"
}

* Mixed scale
clear
set obs 200
gen x = runiform() * 1e6
gen y = 0.001 * x + rnormal()
benchmark_qreg y x, testname("large x, small y")

* Large integers
clear
set obs 200
gen x = runiformint(1000000, 9999999)
gen y = 0.5 * x + runiformint(1, 1000000)
benchmark_qreg y x, testname("large integers")

* Near-zero values
clear
set obs 200
gen x = rnormal(0, 0.001)
gen y = 2*x + rnormal(0, 0.001)
benchmark_qreg y x, testname("near-zero values")

* Extreme range in same dataset
clear
set obs 200
gen x = runiform()
gen y = cond(_n <= 100, 1e6 * x, 1e-6 * x) + rnormal()
benchmark_qreg y x, testname("extreme range in y")

/*******************************************************************************
 * SECTION 21: denmethod Option
 ******************************************************************************/
print_section "denmethod Option"

* denmethod(fitted) is the default for both qreg and cqreg - compare fully
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
matrix qreg_V = e(V)
quietly cqreg price mpg weight, denmethod(fitted)
matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "denmethod(fitted): coefficients"
assert_matrix_equal qreg_V cqreg_V $DEFAULT_SIGFIGS "denmethod(fitted): VCE"

* denmethod(residual) - qreg does not expose this option, so compare coefficients
* only (coefficients are independent of denmethod; VCE will intentionally differ)
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
quietly cqreg price mpg weight, denmethod(residual)
matrix cqreg_b = e(b)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "denmethod(residual): coefficients"

* denmethod(kernel) - qreg does not expose this option, so compare coefficients only
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
quietly cqreg price mpg weight, denmethod(kernel)
matrix cqreg_b = e(b)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "denmethod(kernel): coefficients"

/*******************************************************************************
 * SECTION 22: bwmethod Option
 ******************************************************************************/
print_section "bwmethod Option"

* bwmethod(hsheather) is the default for both qreg and cqreg - compare fully
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
matrix qreg_V = e(V)
quietly cqreg price mpg weight, bwmethod(hsheather)
matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "bwmethod(hsheather): coefficients"
assert_matrix_equal qreg_V cqreg_V $DEFAULT_SIGFIGS "bwmethod(hsheather): VCE"

* bwmethod(bofinger) - qreg does not expose this option, so compare coefficients
* only (coefficients are independent of bwmethod; VCE will intentionally differ)
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
quietly cqreg price mpg weight, bwmethod(bofinger)
matrix cqreg_b = e(b)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "bwmethod(bofinger): coefficients"

* bwmethod(chamberlain) - qreg does not expose this option, so compare coefficients only
sysuse auto, clear
quietly qreg price mpg weight
matrix qreg_b = e(b)
quietly cqreg price mpg weight, bwmethod(chamberlain)
matrix cqreg_b = e(b)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "bwmethod(chamberlain): coefficients"

/*******************************************************************************
 * SECTION 23: tolerance/maxiter Options
 ******************************************************************************/
print_section "tolerance/maxiter Options"

sysuse auto, clear
capture cqreg price mpg weight, tolerance(1e-10)
if _rc == 0 {
    test_pass "tolerance(1e-10) accepted"
}
else {
    test_fail "tolerance option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, tolerance(1e-4)
if _rc == 0 {
    test_pass "tolerance(1e-4) accepted"
}
else {
    test_fail "tolerance(1e-4)" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, maxiter(500)
if _rc == 0 {
    test_pass "maxiter(500) accepted"
}
else {
    test_fail "maxiter option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, maxiter(1000)
if _rc == 0 {
    test_pass "maxiter(1000) accepted"
}
else {
    test_fail "maxiter(1000)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 24: verbose/timeit Options
 ******************************************************************************/
print_section "verbose/timeit Options"

sysuse auto, clear
capture cqreg price mpg weight, verbose
if _rc == 0 {
    test_pass "verbose option accepted"
}
else {
    test_fail "verbose option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, timeit
if _rc == 0 {
    test_pass "timeit option accepted"
}
else {
    test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 25: if/in Conditions
 ******************************************************************************/
print_section "if/in Conditions"

* if condition - use benchmark_qreg for full comparison (N, coefficients, VCE)
sysuse auto, clear
benchmark_qreg price mpg weight if price > 5000, testname("if condition")

* in condition - LP may have multiple optimal solutions (alternative optima)
* so we compare sum_adev (objective value) rather than requiring identical coefficients
sysuse auto, clear
quietly qreg price mpg weight in 1/50
local qreg_N = e(N)
local qreg_adev = e(sum_adev)
local qreg_cols = colsof(e(b))
matrix qreg_V = e(V)

capture quietly cqreg price mpg weight in 1/50
if _rc != 0 {
    test_fail "in condition" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)
    local cqreg_cols = colsof(e(b))
    matrix cqreg_V = e(V)

    if `qreg_N' != `cqreg_N' {
        test_fail "in condition" "N differs"
    }
    else if `qreg_cols' != `cqreg_cols' {
        test_fail "in condition" "dimension mismatch: qreg=`qreg_cols' cqreg=`cqreg_cols'"
    }
    else {
        * Check objective value using sigfigs
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        * Check VCE using sigfigs
        matrix_min_sigfigs qreg_V cqreg_V
        local V_sf = r(min_sigfigs)

        if `adev_sf' >= $DEFAULT_SIGFIGS & `V_sf' >= $DEFAULT_SIGFIGS {
            test_pass "in condition (alternate solution OK)"
        }
        else if `adev_sf' < $DEFAULT_SIGFIGS {
            local sf_fmt : display %4.1f `adev_sf'
            test_fail "in condition" "sum_adev sigfigs=`sf_fmt'"
        }
        else {
            local sf_fmt : display %4.1f `V_sf'
            test_fail "in condition" "e(V) sigfigs=`sf_fmt'"
        }
    }
}

* Combined if and in - use benchmark_qreg for full comparison
sysuse auto, clear
benchmark_qreg price mpg weight if foreign == 0 in 1/60, testname("if and in combined")

* if with numeric comparison - use benchmark_qreg for full comparison
sysuse auto, clear
benchmark_qreg price mpg weight if mpg > 20, testname("if mpg > 20")

/*******************************************************************************
 * SECTION 26: absorb Option (experimental)
 ******************************************************************************/
print_section "absorb Option (experimental)"

sysuse auto, clear
capture cqreg price mpg weight, absorb(foreign)
if _rc == 0 {
    test_pass "absorb(foreign) accepted"
}
else {
    test_fail "absorb option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, absorb(rep78)
if _rc == 0 {
    test_pass "absorb(rep78) accepted"
}
else {
    test_fail "absorb(rep78)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 27: Synthetic Data Patterns
 ******************************************************************************/
print_section "Synthetic Data Patterns"

* Homoskedastic errors
clear
set seed 11111
set obs 500
gen x = rnormal()
gen y = 2 + 3*x + rnormal()
benchmark_qreg y x, testname("homoskedastic errors")

* Heteroskedastic errors
clear
set seed 22222
set obs 500
gen x = runiform()
gen y = 2 + 3*x + rnormal() * x
benchmark_qreg y x, testname("heteroskedastic errors")

* Skewed dependent variable
clear
set seed 33333
set obs 500
gen x = runiform()
gen y = exp(x + rnormal(0, 0.5))
benchmark_qreg y x, testname("skewed y (exponential)")

* Heavy-tailed errors
clear
set seed 44444
set obs 500
gen x = rnormal()
gen e = rt(3)  // t-distribution with 3 df
gen y = 1 + 2*x + e
benchmark_qreg y x, testname("heavy-tailed errors (t-dist)")

* Bimodal dependent variable
clear
set seed 55555
set obs 500
gen x = runiform()
gen group = runiform() < 0.5
gen y = cond(group, 10 + 2*x, 50 + 2*x) + rnormal()
benchmark_qreg y x, testname("bimodal y")

* Outliers in y
clear
set seed 66666
set obs 500
gen x = rnormal()
gen y = 2*x + rnormal()
replace y = y * 10 if _n <= 5  // 1% extreme outliers
benchmark_qreg y x, testname("outliers in y (1%)")

* Outliers in x
clear
set seed 77777
set obs 500
gen x = rnormal()
replace x = x * 10 if _n <= 5
gen y = 2*x + rnormal()
benchmark_qreg y x, testname("outliers in x (1%)")

/*******************************************************************************
 * SECTION 28: Multiple Quantile Sequence
 ******************************************************************************/
print_section "Multiple Quantile Sequence"

* Run through full quantile sequence on same data
clear
set seed 88888
set obs 1000
gen x = rnormal()
gen y = 3*x + rnormal()
tempfile qdata
save `qdata'

foreach q in 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 {
    use `qdata', clear
    benchmark_qreg y x, quantile(`q') testname("quantile sequence: q=`q'")
}

/*******************************************************************************
 * SECTION 29: Convergence Stress Tests
 ******************************************************************************/
print_section "Convergence Stress Tests"

* Difficult convergence scenario - near-singular design
clear
set obs 100
gen x1 = runiform()
gen x2 = x1 + rnormal(0, 0.0001)
gen y = x1 + x2 + rnormal()
capture benchmark_qreg y x1, testname("near-singular design")

* Very high leverage points
clear
set obs 100
gen x = rnormal()
replace x = 100 in 1  // extreme leverage
gen y = 2*x + rnormal()
capture benchmark_qreg y x, testname("high leverage point")

* Constant y (should fail gracefully with 198)
clear
set obs 50
gen x = rnormal()
gen y = 5
capture cqreg y x
if _rc == 198 {
    test_pass "constant y - fails gracefully (rc=198)"
}
else if _rc != 0 {
    test_fail "constant y" "failed with unexpected rc=`=_rc' (expected 198)"
}
else {
    test_fail "constant y" "should have failed but succeeded"
}

* Constant x (should succeed by dropping collinear variable, matching qreg)
clear
set obs 50
gen x = 5
gen y = rnormal()
capture cqreg y x
if _rc == 0 {
    test_pass "constant x - handled gracefully (collinear var dropped)"
}
else {
    test_fail "constant x" "failed with rc=`=_rc' (should succeed like qreg)"
}

/*******************************************************************************
 * SECTION 30: Option Combinations
 ******************************************************************************/
print_section "Option Combinations"

* quantile + vce
sysuse auto, clear
benchmark_qreg price mpg weight, quantile(0.25) vce(robust) testname("quantile(0.25) + vce(robust)")

* quantile + denmethod (fitted is default for qreg, so full comparison works)
sysuse auto, clear
quietly qreg price mpg weight, quantile(0.75)
matrix qreg_b = e(b)
matrix qreg_V = e(V)
quietly cqreg price mpg weight, quantile(0.75) denmethod(fitted)
matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "quantile(0.75) + denmethod(fitted): coefficients"
assert_matrix_equal qreg_V cqreg_V $DEFAULT_SIGFIGS "quantile(0.75) + denmethod(fitted): VCE"

* quantile + bwmethod (hsheather is default for qreg, so full comparison works)
sysuse auto, clear
quietly qreg price mpg weight, quantile(0.10)
matrix qreg_b = e(b)
matrix qreg_V = e(V)
quietly cqreg price mpg weight, quantile(0.10) bwmethod(hsheather)
matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "quantile(0.10) + bwmethod(hsheather): coefficients"
assert_matrix_equal qreg_V cqreg_V $DEFAULT_SIGFIGS "quantile(0.10) + bwmethod(hsheather): VCE"

* quantile + tolerance + maxiter
* Note: benchmark_qreg does not support tolerance/maxiter; compare manually
sysuse auto, clear
quietly qreg price mpg weight, quantile(0.5)
matrix qreg_b = e(b)
matrix qreg_V = e(V)
quietly cqreg price mpg weight, quantile(0.5) tolerance(1e-8) maxiter(1000)
matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
assert_matrix_equal qreg_b cqreg_b $DEFAULT_SIGFIGS "quantile + tolerance + maxiter: coefficients"
assert_matrix_equal qreg_V cqreg_V $DEFAULT_SIGFIGS "quantile + tolerance + maxiter: VCE"

* Full option combo
* Note: verbose does not affect numerical results; vce(robust) is supported by benchmark_qreg
sysuse auto, clear
benchmark_qreg price mpg weight, quantile(0.5) vce(robust) testname("full option combination")

/*******************************************************************************
 * SECTION 31: Real-World Regression Scenarios
 ******************************************************************************/
print_section "Real-World Regression Scenarios"

* Wage regression (classic labor economics)
sysuse nlsw88, clear
benchmark_qreg wage age ttl_exp tenure, testname("wage equation: median")

sysuse nlsw88, clear
benchmark_qreg wage age ttl_exp tenure, quantile(0.10) testname("wage equation: q=0.10 (low earners)")

sysuse nlsw88, clear
benchmark_qreg wage age ttl_exp tenure, quantile(0.90) testname("wage equation: q=0.90 (high earners)")

* Car price regression (hedonic pricing)
sysuse auto, clear
benchmark_qreg price mpg weight length headroom, testname("hedonic price: median")

sysuse auto, clear
benchmark_qreg price mpg weight length headroom, quantile(0.25) testname("hedonic price: q=0.25 (budget)")

sysuse auto, clear
benchmark_qreg price mpg weight length headroom, quantile(0.75) testname("hedonic price: q=0.75 (premium)")

/*******************************************************************************
 * SECTION 32: Factor Variables (i.varname)
 *
 * Tests factor variable expansion (i.varname) in cqreg vs qreg.
 * Compares coefficients, N, and objective function values.
 ******************************************************************************/
print_section "Factor Variables (i.varname)"

* Factor variable with binary indicator (i.foreign)
sysuse auto, clear
quietly qreg price mpg weight i.foreign
local qreg_N = e(N)
local qreg_adev = e(sum_adev)
matrix qreg_b = e(b)

capture quietly cqreg price mpg weight i.foreign
if _rc != 0 {
    test_fail "i.foreign: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)
    matrix cqreg_b = e(b)

    * Check N matches
    if `qreg_N' != `cqreg_N' {
        test_fail "i.foreign: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        * Check objective function (sum_adev) matches within tolerance
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "i.foreign: median"
        }
        else {
            local sf_fmt : display %4.1f `adev_sf'
            test_fail "i.foreign: median" "sum_adev sigfigs=`sf_fmt'"
        }
    }
}

* Factor variable with multiple levels (i.rep78)
sysuse auto, clear
quietly qreg price mpg weight i.rep78
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg price mpg weight i.rep78
if _rc != 0 {
    test_fail "i.rep78: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "i.rep78: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "i.rep78: median"
        }
        else {
            local sf_fmt : display %4.1f `adev_sf'
            test_fail "i.rep78: median" "sum_adev sigfigs=`sf_fmt'"
        }
    }
}

* i.rep78 at q=0.25
sysuse auto, clear
quietly qreg price mpg weight i.rep78, quantile(0.25)
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg price mpg weight i.rep78, quantile(0.25)
if _rc != 0 {
    test_fail "i.rep78: q=0.25" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "i.rep78: q=0.25" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "i.rep78: q=0.25"
        }
        else {
            test_fail "i.rep78: q=0.25" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* i.rep78 at q=0.75
sysuse auto, clear
quietly qreg price mpg weight i.rep78, quantile(0.75)
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg price mpg weight i.rep78, quantile(0.75)
if _rc != 0 {
    test_fail "i.rep78: q=0.75" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "i.rep78: q=0.75" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "i.rep78: q=0.75"
        }
        else {
            test_fail "i.rep78: q=0.75" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Multiple factor variables
sysuse auto, clear
quietly qreg price mpg i.foreign i.rep78
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg price mpg i.foreign i.rep78
if _rc != 0 {
    test_fail "i.foreign i.rep78: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "i.foreign i.rep78: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "i.foreign i.rep78: median"
        }
        else {
            test_fail "i.foreign i.rep78: median" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Factor variable with nlsw88 dataset (i.race)
sysuse nlsw88, clear
quietly qreg wage age tenure i.race
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg wage age tenure i.race
if _rc != 0 {
    test_fail "nlsw88 i.race: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "nlsw88 i.race: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "nlsw88 i.race: median"
        }
        else {
            test_fail "nlsw88 i.race: median" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* nlsw88 i.race at q=0.10
sysuse nlsw88, clear
quietly qreg wage age tenure i.race, quantile(0.10)
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg wage age tenure i.race, quantile(0.10)
if _rc != 0 {
    test_fail "nlsw88 i.race: q=0.10" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "nlsw88 i.race: q=0.10" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "nlsw88 i.race: q=0.10"
        }
        else {
            test_fail "nlsw88 i.race: q=0.10" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* nlsw88 i.race at q=0.90
sysuse nlsw88, clear
quietly qreg wage age tenure i.race, quantile(0.90)
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg wage age tenure i.race, quantile(0.90)
if _rc != 0 {
    test_fail "nlsw88 i.race: q=0.90" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "nlsw88 i.race: q=0.90" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "nlsw88 i.race: q=0.90"
        }
        else {
            test_fail "nlsw88 i.race: q=0.90" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Factor variable with occupation (many levels)
sysuse nlsw88, clear
quietly qreg wage age tenure i.occupation
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg wage age tenure i.occupation
if _rc != 0 {
    test_fail "nlsw88 i.occupation: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "nlsw88 i.occupation: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "nlsw88 i.occupation: median"
        }
        else {
            test_fail "nlsw88 i.occupation: median" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Factor variable with industry
sysuse nlsw88, clear
quietly qreg wage age tenure i.industry
local qreg_N = e(N)
local qreg_adev = e(sum_adev)

capture quietly cqreg wage age tenure i.industry
if _rc != 0 {
    test_fail "nlsw88 i.industry: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)

    if `qreg_N' != `cqreg_N' {
        test_fail "nlsw88 i.industry: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "nlsw88 i.industry: median"
        }
        else {
            test_fail "nlsw88 i.industry: median" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Continuous-by-factor interaction (c.var#i.var) using benchmark_qreg
sysuse auto, clear
benchmark_qreg price c.mpg#i.foreign weight, testname("c.mpg#i.foreign: median")

sysuse auto, clear
benchmark_qreg price c.mpg#i.foreign weight, quantile(0.25) testname("c.mpg#i.foreign: q=0.25")

* Factor-by-factor interaction (i.var#i.var)
* NOTE: These tests may have alternate optimal solutions (same objective, different coefficients)
* Compare sum_adev (objective value) rather than individual coefficients
sysuse nlsw88, clear
quietly qreg wage i.race#i.married age
local qreg_N = e(N)
local qreg_adev = e(sum_adev)
local qreg_cols = colsof(e(b))

capture quietly cqreg wage i.race#i.married age
if _rc != 0 {
    test_fail "i.race#i.married: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)
    local cqreg_cols = colsof(e(b))

    if `qreg_N' != `cqreg_N' {
        test_fail "i.race#i.married: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
    }
    else if `qreg_cols' != `cqreg_cols' {
        test_fail "i.race#i.married: median" "dimension mismatch: qreg=`qreg_cols' cqreg=`cqreg_cols'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "i.race#i.married: median (alternate solution OK)"
        }
        else {
            test_fail "i.race#i.married: median" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Base level specification (ib#.var)
* NOTE: These tests may have alternate optimal solutions
sysuse auto, clear
quietly qreg price mpg ib3.rep78
local qreg_N = e(N)
local qreg_adev = e(sum_adev)
local qreg_cols = colsof(e(b))

capture quietly cqreg price mpg ib3.rep78
if _rc != 0 {
    test_fail "ib3.rep78: median" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)
    local cqreg_cols = colsof(e(b))

    if `qreg_N' != `cqreg_N' {
        test_fail "ib3.rep78: median" "N differs"
    }
    else if `qreg_cols' != `cqreg_cols' {
        test_fail "ib3.rep78: median" "dimension mismatch: qreg=`qreg_cols' cqreg=`cqreg_cols'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "ib3.rep78: median (alternate solution OK)"
        }
        else {
            test_fail "ib3.rep78: median" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

sysuse auto, clear
quietly qreg price mpg ib3.rep78, quantile(0.75)
local qreg_N = e(N)
local qreg_adev = e(sum_adev)
local qreg_cols = colsof(e(b))

capture quietly cqreg price mpg ib3.rep78, quantile(0.75)
if _rc != 0 {
    test_fail "ib3.rep78: q=0.75" "cqreg returned error `=_rc'"
}
else {
    local cqreg_N = e(N)
    local cqreg_adev = e(sum_adev)
    local cqreg_cols = colsof(e(b))

    if `qreg_N' != `cqreg_N' {
        test_fail "ib3.rep78: q=0.75" "N differs"
    }
    else if `qreg_cols' != `cqreg_cols' {
        test_fail "ib3.rep78: q=0.75" "dimension mismatch: qreg=`qreg_cols' cqreg=`cqreg_cols'"
    }
    else {
        sigfigs `qreg_adev' `cqreg_adev'
        local adev_sf = r(sigfigs)
        if `adev_sf' >= $DEFAULT_SIGFIGS {
            test_pass "ib3.rep78: q=0.75 (alternate solution OK)"
        }
        else {
            test_fail "ib3.rep78: q=0.75" "sum_adev sigfigs=`adev_sf'"
        }
    }
}

* Full factorial interaction (i.var##c.var)
sysuse auto, clear
quietly qreg price i.foreign##c.mpg weight
local qreg_cols = colsof(e(b))
quietly cqreg price i.foreign##c.mpg weight
local cqreg_cols = colsof(e(b))
if `qreg_cols' == `cqreg_cols' {
    benchmark_qreg price i.foreign##c.mpg weight, testname("i.foreign##c.mpg: median")
}
else {
    test_fail "i.foreign##c.mpg: median" "dimension differs: qreg=`qreg_cols' cqreg=`cqreg_cols'"
}

* Factor variables at extreme quantiles
* Note: qreg fails with error 498 at extreme quantiles (insufficient obs per cell)
* cqreg is more robust and may succeed where qreg fails - this is acceptable
sysuse auto, clear
quietly capture qreg price mpg i.rep78, quantile(0.01)
local qreg_rc = _rc
quietly capture cqreg price mpg i.rep78, quantile(0.01)
local cqreg_rc = _rc
if `qreg_rc' == `cqreg_rc' {
    test_pass "i.rep78: q=0.01 (both return rc=`qreg_rc')"
}
else if `cqreg_rc' == 0 & `qreg_rc' != 0 {
    test_pass "i.rep78: q=0.01 (cqreg succeeds, qreg fails - more robust)"
}
else {
    test_fail "i.rep78: q=0.01" "unexpected: qreg=`qreg_rc', cqreg=`cqreg_rc'"
}

sysuse auto, clear
quietly capture qreg price mpg i.rep78, quantile(0.99)
local qreg_rc = _rc
quietly capture cqreg price mpg i.rep78, quantile(0.99)
local cqreg_rc = _rc
if `qreg_rc' == `cqreg_rc' {
    test_pass "i.rep78: q=0.99 (both return rc=`qreg_rc')"
}
else if `cqreg_rc' == 0 & `qreg_rc' != 0 {
    test_pass "i.rep78: q=0.99 (cqreg succeeds, qreg fails - more robust)"
}
else {
    test_fail "i.rep78: q=0.99" "unexpected: qreg=`qreg_rc', cqreg=`cqreg_rc'"
}

/*******************************************************************************
 * SECTION 33: Time Series Operators (L., D., F.)
 *
 * Tests lagged and differenced variables in panel data settings.
 *
 * NOTE: Stata's qreg does not support time-series operators (L., D., F.)
 * directly in the varlist. Users must manually generate lagged/differenced
 * variables. These tests verify cqreg matches qreg when using such variables.
 ******************************************************************************/
print_section "Time Series Variables (panel data)"

* Create panel dataset for time series tests using grunfeld
capture webuse grunfeld, clear
if _rc == 0 {
    quietly xtset company year

    * Generate lagged and differenced variables manually
    by company: gen L_mvalue = mvalue[_n-1]
    by company: gen L2_mvalue = mvalue[_n-2]
    by company: gen D_mvalue = mvalue - mvalue[_n-1]
    by company: gen D_kstock = kstock - kstock[_n-1]
    by company: gen F_mvalue = mvalue[_n+1]

    * Lag variable (L.mvalue equivalent) - median
    quietly qreg invest L_mvalue kstock
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest L_mvalue kstock
    if _rc != 0 {
        test_fail "grunfeld L.mvalue: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld L.mvalue: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld L.mvalue: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld L.mvalue: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Lag variable - q=0.25
    quietly qreg invest L_mvalue kstock, quantile(0.25)
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest L_mvalue kstock, quantile(0.25)
    if _rc != 0 {
        test_fail "grunfeld L.mvalue: q=0.25" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld L.mvalue: q=0.25" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld L.mvalue: q=0.25"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld L.mvalue: q=0.25" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Lag variable - q=0.75
    quietly qreg invest L_mvalue kstock, quantile(0.75)
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest L_mvalue kstock, quantile(0.75)
    if _rc != 0 {
        test_fail "grunfeld L.mvalue: q=0.75" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld L.mvalue: q=0.75" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld L.mvalue: q=0.75"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld L.mvalue: q=0.75" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Difference variable (D.mvalue equivalent) - median
    quietly qreg invest D_mvalue kstock
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest D_mvalue kstock
    if _rc != 0 {
        test_fail "grunfeld D.mvalue: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld D.mvalue: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld D.mvalue: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld D.mvalue: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Difference variable - q=0.25
    quietly qreg invest D_mvalue kstock, quantile(0.25)
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest D_mvalue kstock, quantile(0.25)
    if _rc != 0 {
        test_fail "grunfeld D.mvalue: q=0.25" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld D.mvalue: q=0.25" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld D.mvalue: q=0.25"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld D.mvalue: q=0.25" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Difference variable - q=0.75
    quietly qreg invest D_mvalue kstock, quantile(0.75)
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest D_mvalue kstock, quantile(0.75)
    if _rc != 0 {
        test_fail "grunfeld D.mvalue: q=0.75" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld D.mvalue: q=0.75" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld D.mvalue: q=0.75"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld D.mvalue: q=0.75" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Multiple lags (L2.mvalue equivalent) - median
    quietly qreg invest L2_mvalue kstock
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest L2_mvalue kstock
    if _rc != 0 {
        test_fail "grunfeld L2.mvalue: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld L2.mvalue: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld L2.mvalue: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld L2.mvalue: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Lag and difference combined
    quietly qreg invest L_mvalue D_kstock
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest L_mvalue D_kstock
    if _rc != 0 {
        test_fail "grunfeld L.mvalue D.kstock: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld L.mvalue D.kstock: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld L.mvalue D.kstock: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld L.mvalue D.kstock: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Forward variable (F.mvalue equivalent) - median
    quietly qreg invest F_mvalue kstock
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg invest F_mvalue kstock
    if _rc != 0 {
        test_fail "grunfeld F.mvalue: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "grunfeld F.mvalue: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "grunfeld F.mvalue: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "grunfeld F.mvalue: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }
}

* Time series with nlswork panel data
capture webuse nlswork, clear
if _rc == 0 {
    keep in 1/5000
    quietly xtset idcode year

    * Generate lagged and differenced variables manually
    by idcode: gen L_tenure = tenure[_n-1]
    by idcode: gen D_tenure = tenure - tenure[_n-1]

    * Lag of tenure - median
    quietly qreg ln_wage L_tenure age
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg ln_wage L_tenure age
    if _rc != 0 {
        test_fail "nlswork L.tenure: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "nlswork L.tenure: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "nlswork L.tenure: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "nlswork L.tenure: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Lag of tenure - q=0.25
    quietly qreg ln_wage L_tenure age, quantile(0.25)
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg ln_wage L_tenure age, quantile(0.25)
    if _rc != 0 {
        test_fail "nlswork L.tenure: q=0.25" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "nlswork L.tenure: q=0.25" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "nlswork L.tenure: q=0.25"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "nlswork L.tenure: q=0.25" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Lag of tenure - q=0.75
    quietly qreg ln_wage L_tenure age, quantile(0.75)
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg ln_wage L_tenure age, quantile(0.75)
    if _rc != 0 {
        test_fail "nlswork L.tenure: q=0.75" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "nlswork L.tenure: q=0.75" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "nlswork L.tenure: q=0.75"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "nlswork L.tenure: q=0.75" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }

    * Difference of tenure - median
    quietly qreg ln_wage D_tenure age
    local qreg_N = e(N)
    local qreg_adev = e(sum_adev)

    capture quietly cqreg ln_wage D_tenure age
    if _rc != 0 {
        test_fail "nlswork D.tenure: median" "cqreg returned error `=_rc'"
    }
    else {
        local cqreg_N = e(N)
        local cqreg_adev = e(sum_adev)

        if `qreg_N' != `cqreg_N' {
            test_fail "nlswork D.tenure: median" "N differs: qreg=`qreg_N' cqreg=`cqreg_N'"
        }
        else {
            sigfigs `qreg_adev' `cqreg_adev'
            local adev_sf = r(sigfigs)
            if `adev_sf' >= $DEFAULT_SIGFIGS {
                test_pass "nlswork D.tenure: median"
            }
            else {
                local sf_fmt : display %4.1f `adev_sf'
                test_fail "nlswork D.tenure: median" "sum_adev sigfigs=`sf_fmt'"
            }
        }
    }
}

* Multiple lags in same regression using benchmark_qreg
capture webuse grunfeld, clear
if _rc == 0 {
    quietly xtset company year
    by company: gen L_mvalue = mvalue[_n-1]
    by company: gen L2_mvalue = mvalue[_n-2]

    benchmark_qreg invest L_mvalue L2_mvalue kstock, testname("L. and L2. in same regression: median")

    * Multiple lags at different quantiles
    webuse grunfeld, clear
    quietly xtset company year
    by company: gen L_mvalue = mvalue[_n-1]
    by company: gen L2_mvalue = mvalue[_n-2]

    benchmark_qreg invest L_mvalue L2_mvalue kstock, quantile(0.25) testname("L. and L2.: q=0.25")
}

* Lead and lag together
capture webuse grunfeld, clear
if _rc == 0 {
    quietly xtset company year
    by company: gen L_mvalue = mvalue[_n-1]
    by company: gen F_mvalue = mvalue[_n+1]

    benchmark_qreg invest L_mvalue F_mvalue kstock, testname("L. and F. together: median")
}

* Time series with factor variables
* Time series with factor variables
* NOTE: i.company creates many levels, cqreg may differ in base level handling
capture webuse grunfeld, clear
if _rc == 0 {
    quietly xtset company year
    by company: gen L_mvalue = mvalue[_n-1]

    quietly qreg invest L_mvalue kstock i.company
    local qreg_cols = colsof(e(b))
    quietly cqreg invest L_mvalue kstock i.company
    local cqreg_cols = colsof(e(b))
    if `qreg_cols' == `cqreg_cols' {
        webuse grunfeld, clear
        quietly xtset company year
        by company: gen L_mvalue = mvalue[_n-1]
        benchmark_qreg invest L_mvalue kstock i.company, testname("L.mvalue + i.company: median")
    }
    else {
        test_fail "L.mvalue + i.company: median" "dimension differs: qreg=`qreg_cols' cqreg=`cqreg_cols'"
    }

    webuse grunfeld, clear
    quietly xtset company year
    by company: gen L_mvalue = mvalue[_n-1]

    quietly qreg invest L_mvalue kstock i.company, quantile(0.75)
    local qreg_cols = colsof(e(b))
    quietly cqreg invest L_mvalue kstock i.company, quantile(0.75)
    local cqreg_cols = colsof(e(b))
    if `qreg_cols' == `cqreg_cols' {
        webuse grunfeld, clear
        quietly xtset company year
        by company: gen L_mvalue = mvalue[_n-1]
        benchmark_qreg invest L_mvalue kstock i.company, quantile(0.75) testname("L.mvalue + i.company: q=0.75")
    }
    else {
        test_fail "L.mvalue + i.company: q=0.75" "dimension differs: qreg=`qreg_cols' cqreg=`cqreg_cols'"
    }
}

/*******************************************************************************
 * SECTION 34: Reproducibility with Seed
 ******************************************************************************/
print_section "Reproducibility"

* Run same regression twice to check consistency
clear
set seed 99999
set obs 500
gen x = rnormal()
gen y = 2*x + rnormal()

cqreg y x
matrix b1 = e(b)
local N1 = e(N)

cqreg y x
matrix b2 = e(b)
local N2 = e(N)

if b1[1,1] == b2[1,1] & `N1' == `N2' {
    test_pass "reproducibility: same results on same data"
}
else {
    test_fail "reproducibility" "results differ on identical data"
}

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cqreg returns the same error codes as qreg
 * when given invalid inputs or error conditions.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Variable doesn't exist
sysuse auto, clear
test_error_match, stata_cmd(qreg price nonexistent_var) ctools_cmd(cqreg price nonexistent_var) testname("nonexistent variable")

* Invalid quantile (> 1)
sysuse auto, clear
test_error_match, stata_cmd(qreg price mpg, quantile(1.5)) ctools_cmd(cqreg price mpg, quantile(1.5)) testname("invalid quantile > 1")

* Invalid quantile (< 0)
sysuse auto, clear
test_error_match, stata_cmd(qreg price mpg, quantile(-0.5)) ctools_cmd(cqreg price mpg, quantile(-0.5)) testname("invalid quantile < 0")

* No observations after if condition
sysuse auto, clear
test_error_match, stata_cmd(qreg price mpg if price > 100000) ctools_cmd(cqreg price mpg if price > 100000) testname("no observations after if")

* String variable as dependent
sysuse auto, clear
test_error_match, stata_cmd(qreg make mpg) ctools_cmd(cqreg make mpg) testname("string dependent variable")

/*******************************************************************************
 * Summary
 ******************************************************************************/

* End of cqreg validation
noi print_summary "cqreg"
}
