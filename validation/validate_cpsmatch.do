/*******************************************************************************
 * validate_cpsmatch.do
 *
 * Comprehensive validation tests for cpsmatch vs psmatch2
 * Tests all matching methods, options, and edge cases
 *
 * VERIFICATION: Uses benchmark_psmatch2 helper which compares:
 *   - r(n_treated), r(n_controls), r(n_matched) - sample sizes
 *   - r(att), r(att_se) - treatment effect estimates
 *   - _pscore variable - propensity scores
 *
 * All comparisons use significant figures (default 7 sigfigs).
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for cpsmatch..."

/*******************************************************************************
 * SECTION 1: Plugin functionality check
 ******************************************************************************/
print_section "Plugin Check"

* Check psmatch2 is installed
capture which psmatch2
if _rc != 0 {
    noi di as error "psmatch2 not installed - skipping comparison tests"
    noi di as error "Install with: ssc install psmatch2"
    test_fail "psmatch2 availability" "not installed"
    print_summary "cpsmatch"
    exit 1
}
test_pass "psmatch2 available"

* Create simple test data
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.5 + 0.3*x1 + 0.2*x2 + rnormal() > 0.5)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

capture noisily cpsmatch treat x1 x2
if _rc != 0 {
    local errc = _rc
    test_fail "cpsmatch plugin load" "plugin returned error `errc'"
    print_summary "cpsmatch"
    exit 1
}
test_pass "cpsmatch plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic functionality (creates expected variables)
 ******************************************************************************/
print_section "Basic Functionality"

* Test: Creates expected variables
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.5 + 0.3*x1 + 0.2*x2 + rnormal() > 0.5)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

cpsmatch treat x1 x2

* Store r() values immediately before they get overwritten
local cps_n_treated = r(n_treated)
local cps_n_controls = r(n_controls)
local cps_n_matched = r(n_matched)

capture confirm numeric variable _pscore
if _rc == 0 {
    test_pass "creates _pscore variable"
}
else {
    test_fail "creates _pscore" "variable not created"
}

capture confirm numeric variable _weight
if _rc == 0 {
    test_pass "creates _weight variable"
}
else {
    test_fail "creates _weight" "variable not created"
}

capture confirm numeric variable _support
if _rc == 0 {
    test_pass "creates _support variable"
}
else {
    test_fail "creates _support" "variable not created"
}

* Test: Propensity scores are between 0 and 1
quietly sum _pscore
local pmin = r(min)
local pmax = r(max)
if `pmin' >= 0 & `pmax' <= 1 {
    test_pass "pscores in [0,1] range"
}
else {
    test_fail "pscores range" "min=`pmin', max=`pmax'"
}

* Test: Returns expected scalars
if `cps_n_treated' > 0 & `cps_n_controls' > 0 & `cps_n_matched' > 0 {
    test_pass "returns sample size scalars"
}
else {
    test_fail "sample scalars" "n_treated/n_controls/n_matched missing"
}

/*******************************************************************************
 * SECTION 3: Comparison with psmatch2 - Nearest Neighbor Matching
 ******************************************************************************/
print_section "Nearest Neighbor Matching (vs psmatch2)"

* Test: Basic NN matching
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("NN: basic (1 neighbor)")

* Test: 3 neighbors
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(3) testname("NN: 3 neighbors")

* Test: 5 neighbors
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(5) testname("NN: 5 neighbors")

* Test: NN with caliper
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) caliper(0.1) testname("NN: with caliper(0.1)")

* Test: NN with common support
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) common testname("NN: with common support")

* Test: NN without replacement
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) noreplacement testname("NN: without replacement")

/*******************************************************************************
 * SECTION 4: Comparison with psmatch2 - Kernel Matching
 ******************************************************************************/
print_section "Kernel Matching (vs psmatch2)"

* Test: Kernel matching (default epan)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel testname("kernel: epanechnikov (default)")

* Test: Kernel matching (normal/gaussian)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(normal) testname("kernel: normal")

* Test: Kernel matching (biweight)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(biweight) testname("kernel: biweight")

* Test: Kernel matching (uniform)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(uniform) testname("kernel: uniform")

* Test: Kernel matching (tricube)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(tricube) testname("kernel: tricube")

* Test: Kernel with custom bandwidth
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel bwidth(0.1) testname("kernel: custom bwidth(0.1)")

/*******************************************************************************
 * SECTION 5: Comparison with psmatch2 - Radius Matching
 ******************************************************************************/
print_section "Radius Matching (vs psmatch2)"

* Test: Radius matching
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) radius caliper(0.1) testname("radius: caliper(0.1)")

* Test: Radius with wider caliper
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) radius caliper(0.25) testname("radius: caliper(0.25)")

/*******************************************************************************
 * SECTION 6: Comparison with psmatch2 - Estimation Methods
 ******************************************************************************/
print_section "Estimation Methods (vs psmatch2)"

* Test: Logit estimation
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) logit testname("logit estimation")

* Test: Probit estimation (explicit)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) probit testname("probit estimation")

/*******************************************************************************
 * SECTION 7: Comparison with psmatch2 - Large Datasets
 ******************************************************************************/
print_section "Large Datasets (vs psmatch2)"

* Test: 10K observations
clear
set obs 10000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("10K obs: basic NN")

* Test: 10K with kernel
clear
set obs 10000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel testname("10K obs: kernel")

/*******************************************************************************
 * SECTION 8: Comparison with psmatch2 - Multiple Covariates
 ******************************************************************************/
print_section "Multiple Covariates (vs psmatch2)"

* Test: Single covariate
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen treat = (0.5*x1 + rnormal() > 0)
gen y = 2*treat + x1 + rnormal()

benchmark_psmatch2 treat x1, outcome(y) testname("1 covariate")

* Test: 5 covariates
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen treat = (0.1*x1 + 0.1*x2 + 0.1*x3 + rnormal() > 0)
gen y = 2*treat + x1 + x2 + rnormal()

benchmark_psmatch2 treat x1 x2 x3 x4 x5, outcome(y) testname("5 covariates")

/*******************************************************************************
 * SECTION 9: Comparison with psmatch2 - Unbalanced Treatment
 ******************************************************************************/
print_section "Unbalanced Treatment (vs psmatch2)"

* Test: 20% treated
clear
set obs 1000
set seed 12345
gen x1 = rnormal()
gen treat = runiform() < 0.2
gen y = 2*treat + x1 + rnormal()

benchmark_psmatch2 treat x1, outcome(y) testname("20% treated")

* Test: 80% treated
clear
set obs 1000
set seed 12345
gen x1 = rnormal()
gen treat = runiform() < 0.8
gen y = 2*treat + x1 + rnormal()

benchmark_psmatch2 treat x1, outcome(y) testname("80% treated")

/*******************************************************************************
 * SECTION 10: Comparison with psmatch2 - Combined Options
 ******************************************************************************/
print_section "Combined Options (vs psmatch2)"

* Test: NN + caliper + common support
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(3) caliper(0.1) common testname("NN + caliper + common")

* Test: Kernel + common support
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel common testname("kernel + common")

/*******************************************************************************
 * SECTION 11: if/in Conditions
 ******************************************************************************/
print_section "if/in Conditions"

* Test: if condition
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen treat = (0.3*x1 + rnormal() > 0)
gen y = 2*treat + x1 + rnormal()
gen group = mod(_n, 2)

capture cpsmatch treat x1 if group == 1
if _rc == 0 {
    test_pass "if condition"
}
else {
    local errc = _rc
    test_fail "if condition" "rc=`errc'"
}

* Test: in range
capture cpsmatch treat x1 in 1/250
if _rc == 0 {
    test_pass "in range"
}
else {
    local errc = _rc
    test_fail "in range" "rc=`errc'"
}

/*******************************************************************************
 * SECTION 12: Edge Cases
 ******************************************************************************/
print_section "Edge Cases"

* Test: Small sample (50 obs)
clear
set obs 50
set seed 12345
gen x1 = rnormal()
gen treat = (0.5*x1 + rnormal() > 0)
gen y = 2*treat + x1 + rnormal()

capture cpsmatch treat x1
if _rc == 0 {
    test_pass "small sample (50 obs)"
}
else {
    local errc = _rc
    test_fail "small sample" "rc=`errc'"
}

* Test: Very tight caliper (may drop many matches)
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen treat = (0.3*x1 + rnormal() > 0)
gen y = 2*treat + x1 + rnormal()

capture cpsmatch treat x1, caliper(0.001)
if _rc == 0 {
    test_pass "very tight caliper"
}
else {
    local errc = _rc
    test_fail "tight caliper" "rc=`errc'"
}

* Test: Poor overlap with common support
clear
set obs 200
set seed 12345
gen x1 = cond(_n <= 100, rnormal() - 3, rnormal() + 3)
gen treat = _n <= 100
gen y = 2*treat + x1 + rnormal()

cpsmatch treat x1, common
local off = r(n_off_support)
if `off' > 0 {
    test_pass "poor overlap detected"
}
else {
    test_fail "poor overlap" "expected off-support observations"
}

/*******************************************************************************
 * SECTION 13: cpsmatch-only Options (no psmatch2 comparison)
 ******************************************************************************/
print_section "cpsmatch-only Options"

* Test: verbose option
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen treat = (0.3*x1 + rnormal() > 0)

capture cpsmatch treat x1, verbose
if _rc == 0 {
    test_pass "verbose option"
}
else {
    local errc = _rc
    test_fail "verbose" "rc=`errc'"
}

* Test: threads option
capture cpsmatch treat x1, threads(2)
if _rc == 0 {
    test_pass "threads(2) option"
}
else {
    local errc = _rc
    test_fail "threads" "rc=`errc'"
}

/*******************************************************************************
 * SECTION 14: ATT Verification
 ******************************************************************************/
print_section "ATT Verification"

* Test: ATT sign and magnitude (true effect = 2)
clear
set obs 2000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

cpsmatch treat x1 x2, outcome(y)
local att = r(att)
local att_se = r(att_se)
local att_t = r(att_t)

* With large sample, ATT should be close to 2 (within 1.5 of true value)
if abs(`att' - 2) < 1.5 {
    test_pass "ATT magnitude reasonable (near true=2)"
}
else {
    test_fail "ATT magnitude" "expected ~2, got `att'"
}

* Test: ATT returns SE and t-stat
if `att_se' != . & `att_se' > 0 {
    test_pass "ATT SE computed and positive"
}
else {
    test_fail "ATT SE" "missing or non-positive"
}

if `att_t' != . {
    test_pass "ATT t-stat computed"
}
else {
    test_fail "ATT t-stat" "missing"
}

/*******************************************************************************
 * SECTION 15: Manual Verification
 ******************************************************************************/
print_section "Manual Verification"

* Test: Pscore matches probit prediction
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

cpsmatch treat x1 x2
probit treat x1 x2
predict pscore_manual, pr
quietly gen diff = abs(_pscore - pscore_manual)
quietly sum diff
local maxdiff = r(max)
if `maxdiff' < 0.0001 {
    test_pass "pscore matches probit prediction"
}
else {
    test_fail "pscore match" "max diff = `maxdiff'"
}
drop pscore_manual diff

* Test: _support is binary
quietly count if _support != 0 & _support != 1 & !missing(_support)
local nonbinary = r(N)
if `nonbinary' == 0 {
    test_pass "_support is binary"
}
else {
    test_fail "_support binary" "`nonbinary' non-binary values"
}

* Test: Weights are non-negative
quietly sum _weight
local minw = r(min)
if `minw' >= 0 {
    test_pass "weights are non-negative"
}
else {
    test_fail "weights non-negative" "min=`minw'"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cpsmatch"
}
