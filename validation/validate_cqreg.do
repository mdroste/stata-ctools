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

di as text ""
di as text "======================================================================"
di as text "              CQREG VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture cqreg price mpg weight
if _rc != 0 {
    noi test_fail "cqreg plugin load" "returned error `=_rc'"
    noi print_summary "cqreg"
    exit 1
}
noi test_pass "cqreg plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic quantile tests (auto)
 ******************************************************************************/
noi print_section "Basic Quantile Tests (auto)"

sysuse auto, clear
noi benchmark_qreg price mpg weight, testname("median (q=0.5)")

sysuse auto, clear
noi benchmark_qreg price mpg weight, quantile(0.25) testname("q=0.25")

sysuse auto, clear
noi benchmark_qreg price mpg weight, quantile(0.75) testname("q=0.75")

sysuse auto, clear
noi benchmark_qreg price mpg weight, quantile(0.10) testname("q=0.10")

sysuse auto, clear
noi benchmark_qreg price mpg weight, quantile(0.90) testname("q=0.90")

/*******************************************************************************
 * SECTION 3: Covariate variations
 ******************************************************************************/
noi print_section "Covariate Variations"

sysuse auto, clear
noi benchmark_qreg price mpg, testname("single covariate")

sysuse auto, clear
noi benchmark_qreg price mpg weight length, testname("three covariates")

sysuse auto, clear
noi benchmark_qreg price mpg weight length turn displacement, testname("many covariates")

/*******************************************************************************
 * SECTION 4: VCE options - Full coverage
 ******************************************************************************/
noi print_section "VCE Options - Full Coverage"

* vce(robust) - Huber sandwich estimator
sysuse auto, clear
noi benchmark_qreg price mpg weight, vce(robust) testname("vce(robust)")

* vce(iid) - assumes i.i.d. errors
sysuse auto, clear
capture cqreg price mpg weight, vce(iid)
if _rc == 0 {
    noi test_pass "vce(iid) accepted"
}
else {
    noi test_fail "vce(iid)" "returned error `=_rc'"
}

* vce(bootstrap) - bootstrap standard errors
sysuse auto, clear
capture cqreg price mpg weight, vce(bootstrap, reps(50))
if _rc == 0 {
    noi test_pass "vce(bootstrap) accepted"
}
else {
    noi test_fail "vce(bootstrap)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 5: Sysuse datasets - Comprehensive coverage
 ******************************************************************************/
noi print_section "Sysuse Datasets - auto"

* auto dataset - multiple quantiles
sysuse auto, clear
noi benchmark_qreg price mpg weight length, testname("auto: median")

sysuse auto, clear
noi benchmark_qreg price mpg weight length, quantile(0.10) testname("auto: q=0.10")

sysuse auto, clear
noi benchmark_qreg price mpg weight length, quantile(0.25) testname("auto: q=0.25")

sysuse auto, clear
noi benchmark_qreg price mpg weight length, quantile(0.75) testname("auto: q=0.75")

sysuse auto, clear
noi benchmark_qreg price mpg weight length, quantile(0.90) testname("auto: q=0.90")

/*******************************************************************************
 * SECTION 6: Census dataset
 ******************************************************************************/
noi print_section "Sysuse Datasets - census"

sysuse census, clear
noi benchmark_qreg pop medage death, testname("census: median")

sysuse census, clear
noi benchmark_qreg pop medage death, quantile(0.25) testname("census: q=0.25")

sysuse census, clear
noi benchmark_qreg pop medage death, quantile(0.75) testname("census: q=0.75")

sysuse census, clear
noi benchmark_qreg pop medage death, quantile(0.10) testname("census: q=0.10")

sysuse census, clear
noi benchmark_qreg pop medage death, quantile(0.90) testname("census: q=0.90")

sysuse census, clear
noi benchmark_qreg pop medage popurban marriage, testname("census: many covariates")

/*******************************************************************************
 * SECTION 7: lifeexp dataset
 ******************************************************************************/
noi print_section "Sysuse Datasets - lifeexp"

sysuse lifeexp, clear
noi benchmark_qreg lexp gnppc, testname("lifeexp: median")

sysuse lifeexp, clear
noi benchmark_qreg lexp gnppc, quantile(0.25) testname("lifeexp: q=0.25")

sysuse lifeexp, clear
noi benchmark_qreg lexp gnppc, quantile(0.75) testname("lifeexp: q=0.75")

sysuse lifeexp, clear
noi benchmark_qreg lexp gnppc safewater popgrowth, testname("lifeexp: many covariates")

/*******************************************************************************
 * SECTION 8: nlsw88 dataset
 ******************************************************************************/
noi print_section "Sysuse Datasets - nlsw88"

sysuse nlsw88, clear
noi benchmark_qreg wage age tenure, testname("nlsw88: median")

sysuse nlsw88, clear
noi benchmark_qreg wage age tenure, quantile(0.10) testname("nlsw88: q=0.10")

sysuse nlsw88, clear
noi benchmark_qreg wage age tenure, quantile(0.25) testname("nlsw88: q=0.25")

sysuse nlsw88, clear
noi benchmark_qreg wage age tenure, quantile(0.75) testname("nlsw88: q=0.75")

sysuse nlsw88, clear
noi benchmark_qreg wage age tenure, quantile(0.90) testname("nlsw88: q=0.90")

sysuse nlsw88, clear
noi benchmark_qreg wage age tenure ttl_exp hours, testname("nlsw88: many covariates")

/*******************************************************************************
 * SECTION 9: Webuse datasets - nlswork panel data
 ******************************************************************************/
noi print_section "Webuse Datasets - nlswork"

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, testname("nlswork: median (5K obs)")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, quantile(0.10) testname("nlswork: q=0.10")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, quantile(0.25) testname("nlswork: q=0.25")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, quantile(0.75) testname("nlswork: q=0.75")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, quantile(0.90) testname("nlswork: q=0.90")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure hours wks_work, testname("nlswork: many covariates")

/*******************************************************************************
 * SECTION 10: Webuse datasets - grunfeld panel data
 ******************************************************************************/
noi print_section "Webuse Datasets - grunfeld"

capture webuse grunfeld, clear
if _rc == 0 {
    noi benchmark_qreg invest mvalue kstock, testname("grunfeld: median")

    webuse grunfeld, clear
    noi benchmark_qreg invest mvalue kstock, quantile(0.25) testname("grunfeld: q=0.25")

    webuse grunfeld, clear
    noi benchmark_qreg invest mvalue kstock, quantile(0.75) testname("grunfeld: q=0.75")
}

/*******************************************************************************
 * SECTION 11: Webuse datasets - bplong
 ******************************************************************************/
noi print_section "Webuse Datasets - bplong"

capture webuse bplong, clear
if _rc == 0 {
    noi benchmark_qreg bp agegrp when sex, testname("bplong: median")

    webuse bplong, clear
    noi benchmark_qreg bp agegrp when sex, quantile(0.25) testname("bplong: q=0.25")

    webuse bplong, clear
    noi benchmark_qreg bp agegrp when sex, quantile(0.75) testname("bplong: q=0.75")
}

/*******************************************************************************
 * SECTION 12: Webuse datasets - cancer survival data
 ******************************************************************************/
noi print_section "Webuse Datasets - cancer"

capture webuse cancer, clear
if _rc == 0 {
    gen age2 = age^2
    noi benchmark_qreg studytime age drug, testname("cancer: median")

    webuse cancer, clear
    gen age2 = age^2
    noi benchmark_qreg studytime age drug, quantile(0.25) testname("cancer: q=0.25")

    webuse cancer, clear
    gen age2 = age^2
    noi benchmark_qreg studytime age drug, quantile(0.75) testname("cancer: q=0.75")
}

/*******************************************************************************
 * SECTION 13: Pathological - Small Datasets (near minimum observations)
 ******************************************************************************/
noi print_section "Pathological - Small Datasets"

* Minimum size for 2 covariates + constant = 4 observations (bare minimum)
* Test with just enough observations
clear
set obs 10
gen x = runiform()
gen y = 2*x + rnormal()
noi benchmark_qreg y x, testname("10 observations (bare minimum)")

* Slightly larger
clear
set obs 15
gen x = runiform()
gen y = 2*x + rnormal()
noi benchmark_qreg y x, testname("15 observations")

* 20 observations
clear
set obs 20
gen x1 = runiform()
gen x2 = rnormal()
gen y = 2*x1 + 3*x2 + rnormal()
noi benchmark_qreg y x1 x2, testname("20 observations, 2 covariates")

* 25 observations
clear
set obs 25
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 10)
gen y = x1 + x2 + 0.5*x3 + rnormal()
noi benchmark_qreg y x1 x2 x3, testname("25 observations, 3 covariates")

* Test at different quantiles with small N
clear
set obs 30
gen x = runiform()
gen y = 3*x + rnormal()
noi benchmark_qreg y x, quantile(0.25) testname("30 obs, q=0.25")

clear
set obs 30
gen x = runiform()
gen y = 3*x + rnormal()
noi benchmark_qreg y x, quantile(0.75) testname("30 obs, q=0.75")

/*******************************************************************************
 * SECTION 14: Pathological - Missing Values
 ******************************************************************************/
noi print_section "Pathological - Missing Values"

* Missing values in dependent variable
clear
set obs 100
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . in 1/10
noi benchmark_qreg y x, testname("missing in depvar (10%)")

* Missing values in independent variable
clear
set obs 100
gen x = runiform()
replace x = . in 20/30
gen y = 2*x + rnormal()
noi benchmark_qreg y x, testname("missing in indepvar (11%)")

* Missing values in both
clear
set obs 100
gen x1 = runiform()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal()
replace x1 = . in 1/5
replace x2 = . in 10/15
replace y = . in 20/25
noi benchmark_qreg y x1 x2, testname("missing in both depvar and indepvars")

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
noi benchmark_qreg y x, testname("sparse missing pattern")

* High proportion missing (50%)
clear
set obs 200
gen x = runiform()
gen y = 2*x + rnormal()
replace y = . if runiform() < 0.5
noi benchmark_qreg y x, testname("50% missing")

/*******************************************************************************
 * SECTION 15: Pathological - Perfect Collinearity and Near-Collinearity
 ******************************************************************************/
noi print_section "Pathological - Collinearity Tests"

* Near-perfect collinearity (highly correlated regressors)
clear
set obs 200
gen x1 = runiform()
gen x2 = x1 + rnormal(0, 0.01)  // x2 almost = x1
gen y = x1 + x2 + rnormal()
capture noi benchmark_qreg y x1 x2, testname("near-collinear regressors")

* Dummy trap scenario (but not exact)
clear
set obs 200
gen d1 = runiform() > 0.5
gen d2 = 1 - d1 + rnormal(0, 0.1)  // almost but not exactly complementary
gen y = d1 + d2 + rnormal()
capture noi benchmark_qreg y d1 d2, testname("near-dummy trap")

* High but not perfect correlation
clear
set obs 500
gen x1 = rnormal()
gen x2 = 0.99*x1 + rnormal(0, 0.1)  // r ~ 0.99
gen y = x1 + x2 + rnormal()
capture noi benchmark_qreg y x1 x2, testname("high correlation (r~0.99)")

/*******************************************************************************
 * SECTION 16: Extreme Quantiles (tails)
 ******************************************************************************/
noi print_section "Extreme Quantiles"

* Very low quantiles
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.01)
if _rc == 0 {
    noi test_pass "q=0.01 accepted"
}
else {
    noi test_fail "q=0.01" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, quantile(0.05)
if _rc == 0 {
    noi test_pass "q=0.05 accepted"
}
else {
    noi test_fail "q=0.05" "returned error `=_rc'"
}

* Very high quantiles
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.95)
if _rc == 0 {
    noi test_pass "q=0.95 accepted"
}
else {
    noi test_fail "q=0.95" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, quantile(0.99)
if _rc == 0 {
    noi test_pass "q=0.99 accepted"
}
else {
    noi test_fail "q=0.99" "returned error `=_rc'"
}

* Extreme quantiles with larger dataset
clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture noi benchmark_qreg y x, quantile(0.01) testname("1K obs, q=0.01")

clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture noi benchmark_qreg y x, quantile(0.05) testname("1K obs, q=0.05")

clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture noi benchmark_qreg y x, quantile(0.95) testname("1K obs, q=0.95")

clear
set obs 1000
gen x = runiform()
gen y = 2*x + rnormal()
capture noi benchmark_qreg y x, quantile(0.99) testname("1K obs, q=0.99")

/*******************************************************************************
 * SECTION 17: Single Covariate Tests
 ******************************************************************************/
noi print_section "Single Covariate Tests"

sysuse auto, clear
noi benchmark_qreg price mpg, testname("auto: price ~ mpg")

sysuse auto, clear
noi benchmark_qreg price weight, testname("auto: price ~ weight")

sysuse auto, clear
noi benchmark_qreg mpg weight, testname("auto: mpg ~ weight")

sysuse nlsw88, clear
noi benchmark_qreg wage tenure, testname("nlsw88: wage ~ tenure")

sysuse nlsw88, clear
noi benchmark_qreg wage age, testname("nlsw88: wage ~ age")

* Single covariate at different quantiles
sysuse auto, clear
noi benchmark_qreg price mpg, quantile(0.10) testname("single covar, q=0.10")

sysuse auto, clear
noi benchmark_qreg price mpg, quantile(0.90) testname("single covar, q=0.90")

/*******************************************************************************
 * SECTION 18: Many Covariates Tests
 ******************************************************************************/
noi print_section "Many Covariates Tests"

* 5 covariates
sysuse auto, clear
noi benchmark_qreg price mpg weight length turn displacement, testname("5 covariates")

* 6 covariates
sysuse auto, clear
noi benchmark_qreg price mpg weight length turn displacement gear_ratio, testname("6 covariates")

* Many covariates at different quantiles
sysuse auto, clear
noi benchmark_qreg price mpg weight length turn displacement, quantile(0.25) testname("5 covariates, q=0.25")

sysuse auto, clear
noi benchmark_qreg price mpg weight length turn displacement, quantile(0.75) testname("5 covariates, q=0.75")

* nlsw88 with many covariates
sysuse nlsw88, clear
noi benchmark_qreg wage age tenure ttl_exp hours grade, testname("nlsw88: 5 covariates")

/*******************************************************************************
 * SECTION 19: Large Datasets (10K, 50K observations)
 ******************************************************************************/
noi print_section "Large Datasets"

* 10K observations
clear
set seed 12345
set obs 10000
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 100)
gen y = 2*x1 + 3*x2 - 0.5*x3 + rnormal()

noi benchmark_qreg y x1 x2 x3, testname("10K obs: median")

noi benchmark_qreg y x1 x2 x3, quantile(0.25) testname("10K obs: q=0.25")

noi benchmark_qreg y x1 x2 x3, quantile(0.75) testname("10K obs: q=0.75")

noi benchmark_qreg y x1 x2 x3, quantile(0.10) testname("10K obs: q=0.10")

noi benchmark_qreg y x1 x2 x3, quantile(0.90) testname("10K obs: q=0.90")

* 50K observations
clear
set seed 54321
set obs 50000
gen x1 = runiform()
gen x2 = rnormal()
gen y = 5*x1 - 2*x2 + rnormal()

noi benchmark_qreg y x1 x2, testname("50K obs: median")

noi benchmark_qreg y x1 x2, quantile(0.25) testname("50K obs: q=0.25")

noi benchmark_qreg y x1 x2, quantile(0.75) testname("50K obs: q=0.75")

/*******************************************************************************
 * SECTION 20: Numeric Precision - Very Large/Small Dependent Variables
 ******************************************************************************/
noi print_section "Numeric Precision - Large/Small Values"

* Very large dependent variable values
clear
set obs 200
gen x = runiform()
gen y = 1e8 * x + 1e7 + rnormal() * 1e6
noi benchmark_qreg y x, testname("very large y (1e8 scale)")

* Very small dependent variable values
clear
set obs 200
gen x = runiform()
gen y = 1e-8 * x + rnormal() * 1e-9
noi benchmark_qreg y x, testname("very small y (1e-8 scale)")

* Mixed scale
clear
set obs 200
gen x = runiform() * 1e6
gen y = 0.001 * x + rnormal()
noi benchmark_qreg y x, testname("large x, small y")

* Large integers
clear
set obs 200
gen x = runiformint(1000000, 9999999)
gen y = 0.5 * x + runiformint(1, 1000000)
noi benchmark_qreg y x, testname("large integers")

* Near-zero values
clear
set obs 200
gen x = rnormal(0, 0.001)
gen y = 2*x + rnormal(0, 0.001)
noi benchmark_qreg y x, testname("near-zero values")

* Extreme range in same dataset
clear
set obs 200
gen x = runiform()
gen y = cond(_n <= 100, 1e6 * x, 1e-6 * x) + rnormal()
noi benchmark_qreg y x, testname("extreme range in y")

/*******************************************************************************
 * SECTION 21: denmethod Option
 ******************************************************************************/
noi print_section "denmethod Option"

sysuse auto, clear
capture cqreg price mpg weight, denmethod(fitted)
if _rc == 0 {
    noi test_pass "denmethod(fitted) accepted"
}
else {
    noi test_fail "denmethod(fitted)" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, denmethod(residual)
if _rc == 0 {
    noi test_pass "denmethod(residual) accepted"
}
else {
    noi test_fail "denmethod(residual)" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, denmethod(kernel)
if _rc == 0 {
    noi test_pass "denmethod(kernel) accepted"
}
else {
    noi test_fail "denmethod(kernel)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 22: bwmethod Option
 ******************************************************************************/
noi print_section "bwmethod Option"

sysuse auto, clear
capture cqreg price mpg weight, bwmethod(hsheather)
if _rc == 0 {
    noi test_pass "bwmethod(hsheather) accepted"
}
else {
    noi test_fail "bwmethod(hsheather)" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, bwmethod(bofinger)
if _rc == 0 {
    noi test_pass "bwmethod(bofinger) accepted"
}
else {
    noi test_fail "bwmethod(bofinger)" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, bwmethod(chamberlain)
if _rc == 0 {
    noi test_pass "bwmethod(chamberlain) accepted"
}
else {
    noi test_fail "bwmethod(chamberlain)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 23: tolerance/maxiter Options
 ******************************************************************************/
noi print_section "tolerance/maxiter Options"

sysuse auto, clear
capture cqreg price mpg weight, tolerance(1e-10)
if _rc == 0 {
    noi test_pass "tolerance(1e-10) accepted"
}
else {
    noi test_fail "tolerance option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, tolerance(1e-4)
if _rc == 0 {
    noi test_pass "tolerance(1e-4) accepted"
}
else {
    noi test_fail "tolerance(1e-4)" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, maxiter(500)
if _rc == 0 {
    noi test_pass "maxiter(500) accepted"
}
else {
    noi test_fail "maxiter option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, maxiter(1000)
if _rc == 0 {
    noi test_pass "maxiter(1000) accepted"
}
else {
    noi test_fail "maxiter(1000)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 24: verbose/timeit Options
 ******************************************************************************/
noi print_section "verbose/timeit Options"

sysuse auto, clear
capture cqreg price mpg weight, verbose
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, timeit
if _rc == 0 {
    noi test_pass "timeit option accepted"
}
else {
    noi test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 25: if/in Conditions
 ******************************************************************************/
noi print_section "if/in Conditions"

* if condition
sysuse auto, clear
qreg price mpg weight if price > 5000
local qreg_N = e(N)

cqreg price mpg weight if price > 5000
local cqreg_N = e(N)

if `qreg_N' == `cqreg_N' {
    noi test_pass "if condition: N matches"
}
else {
    noi test_fail "if condition" "N differs"
}

* in condition
sysuse auto, clear
qreg price mpg weight in 1/50
local qreg_N = e(N)

cqreg price mpg weight in 1/50
local cqreg_N = e(N)

if `qreg_N' == `cqreg_N' {
    noi test_pass "in condition: N matches"
}
else {
    noi test_fail "in condition" "N differs"
}

* Combined if and in
sysuse auto, clear
qreg price mpg weight if foreign == 0 in 1/60
local qreg_N = e(N)

cqreg price mpg weight if foreign == 0 in 1/60
local cqreg_N = e(N)

if `qreg_N' == `cqreg_N' {
    noi test_pass "if and in combined: N matches"
}
else {
    noi test_fail "if and in combined" "N differs"
}

* if with numeric comparison
sysuse auto, clear
qreg price mpg weight if mpg > 20
local qreg_N = e(N)

cqreg price mpg weight if mpg > 20
local cqreg_N = e(N)

if `qreg_N' == `cqreg_N' {
    noi test_pass "if mpg > 20: N matches"
}
else {
    noi test_fail "if mpg > 20" "N differs"
}

/*******************************************************************************
 * SECTION 26: absorb Option (experimental)
 ******************************************************************************/
noi print_section "absorb Option (experimental)"

sysuse auto, clear
capture cqreg price mpg weight, absorb(foreign)
if _rc == 0 {
    noi test_pass "absorb(foreign) accepted"
}
else {
    noi test_fail "absorb option" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, absorb(rep78)
if _rc == 0 {
    noi test_pass "absorb(rep78) accepted"
}
else {
    noi test_fail "absorb(rep78)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 27: Synthetic Data Patterns
 ******************************************************************************/
noi print_section "Synthetic Data Patterns"

* Homoskedastic errors
clear
set seed 11111
set obs 500
gen x = rnormal()
gen y = 2 + 3*x + rnormal()
noi benchmark_qreg y x, testname("homoskedastic errors")

* Heteroskedastic errors
clear
set seed 22222
set obs 500
gen x = runiform()
gen y = 2 + 3*x + rnormal() * x
noi benchmark_qreg y x, testname("heteroskedastic errors")

* Skewed dependent variable
clear
set seed 33333
set obs 500
gen x = runiform()
gen y = exp(x + rnormal(0, 0.5))
noi benchmark_qreg y x, testname("skewed y (exponential)")

* Heavy-tailed errors
clear
set seed 44444
set obs 500
gen x = rnormal()
gen e = rt(3)  // t-distribution with 3 df
gen y = 1 + 2*x + e
noi benchmark_qreg y x, testname("heavy-tailed errors (t-dist)")

* Bimodal dependent variable
clear
set seed 55555
set obs 500
gen x = runiform()
gen group = runiform() < 0.5
gen y = cond(group, 10 + 2*x, 50 + 2*x) + rnormal()
noi benchmark_qreg y x, testname("bimodal y")

* Outliers in y
clear
set seed 66666
set obs 500
gen x = rnormal()
gen y = 2*x + rnormal()
replace y = y * 10 if _n <= 5  // 1% extreme outliers
noi benchmark_qreg y x, testname("outliers in y (1%)")

* Outliers in x
clear
set seed 77777
set obs 500
gen x = rnormal()
replace x = x * 10 if _n <= 5
gen y = 2*x + rnormal()
noi benchmark_qreg y x, testname("outliers in x (1%)")

/*******************************************************************************
 * SECTION 28: Multiple Quantile Sequence
 ******************************************************************************/
noi print_section "Multiple Quantile Sequence"

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
    noi benchmark_qreg y x, quantile(`q') testname("quantile sequence: q=`q'")
}

/*******************************************************************************
 * SECTION 29: Convergence Stress Tests
 ******************************************************************************/
noi print_section "Convergence Stress Tests"

* Difficult convergence scenario - near-singular design
clear
set obs 100
gen x1 = runiform()
gen x2 = x1 + rnormal(0, 0.0001)
gen y = x1 + x2 + rnormal()
capture noi benchmark_qreg y x1, testname("near-singular design")

* Very high leverage points
clear
set obs 100
gen x = rnormal()
replace x = 100 in 1  // extreme leverage
gen y = 2*x + rnormal()
capture noi benchmark_qreg y x, testname("high leverage point")

* Constant y (should fail gracefully)
clear
set obs 50
gen x = rnormal()
gen y = 5
capture cqreg y x
if _rc != 0 {
    noi test_pass "constant y - fails gracefully (rc=`=_rc')"
}
else {
    noi test_fail "constant y" "should have failed"
}

* Constant x (should fail gracefully)
clear
set obs 50
gen x = 5
gen y = rnormal()
capture cqreg y x
if _rc != 0 {
    noi test_pass "constant x - fails gracefully (rc=`=_rc')"
}
else {
    noi test_fail "constant x" "should have failed"
}

/*******************************************************************************
 * SECTION 30: Option Combinations
 ******************************************************************************/
noi print_section "Option Combinations"

* quantile + vce
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.25) vce(robust)
if _rc == 0 {
    noi test_pass "quantile(0.25) + vce(robust)"
}
else {
    noi test_fail "quantile + vce combo" "returned error `=_rc'"
}

* quantile + denmethod
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.75) denmethod(fitted)
if _rc == 0 {
    noi test_pass "quantile(0.75) + denmethod(fitted)"
}
else {
    noi test_fail "quantile + denmethod combo" "returned error `=_rc'"
}

* quantile + bwmethod
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.10) bwmethod(hsheather)
if _rc == 0 {
    noi test_pass "quantile(0.10) + bwmethod(hsheather)"
}
else {
    noi test_fail "quantile + bwmethod combo" "returned error `=_rc'"
}

* quantile + tolerance + maxiter
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.5) tolerance(1e-8) maxiter(1000)
if _rc == 0 {
    noi test_pass "quantile + tolerance + maxiter"
}
else {
    noi test_fail "multiple options combo" "returned error `=_rc'"
}

* Full option combo
sysuse auto, clear
capture cqreg price mpg weight, quantile(0.5) vce(robust) tolerance(1e-8) verbose
if _rc == 0 {
    noi test_pass "full option combination"
}
else {
    noi test_fail "full option combo" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 31: Real-World Regression Scenarios
 ******************************************************************************/
noi print_section "Real-World Regression Scenarios"

* Wage regression (classic labor economics)
sysuse nlsw88, clear
noi benchmark_qreg wage age ttl_exp tenure, testname("wage equation: median")

sysuse nlsw88, clear
noi benchmark_qreg wage age ttl_exp tenure, quantile(0.10) testname("wage equation: q=0.10 (low earners)")

sysuse nlsw88, clear
noi benchmark_qreg wage age ttl_exp tenure, quantile(0.90) testname("wage equation: q=0.90 (high earners)")

* Car price regression (hedonic pricing)
sysuse auto, clear
noi benchmark_qreg price mpg weight length headroom, testname("hedonic price: median")

sysuse auto, clear
noi benchmark_qreg price mpg weight length headroom, quantile(0.25) testname("hedonic price: q=0.25 (budget)")

sysuse auto, clear
noi benchmark_qreg price mpg weight length headroom, quantile(0.75) testname("hedonic price: q=0.75 (premium)")

/*******************************************************************************
 * SECTION 32: Reproducibility with Seed
 ******************************************************************************/
noi print_section "Reproducibility"

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
    noi test_pass "reproducibility: same results on same data"
}
else {
    noi test_fail "reproducibility" "results differ on identical data"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "cqreg"

if $TESTS_FAILED > 0 {
    exit 1
}

}
