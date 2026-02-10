/*******************************************************************************
 * validate_cpsmatch.do
 *
 * Validation tests for cpsmatch propensity score matching
 * Compares cpsmatch results with psmatch2 for all matching methods
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
    noi print_summary "cpsmatch"
    exit 1
}
test_pass "psmatch2 available"

* Plugin test
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
    noi print_summary "cpsmatch"
    exit 1
}
test_pass "cpsmatch plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic functionality and return values
 ******************************************************************************/
print_section "Basic Functionality & Return Values"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.5 + 0.3*x1 + 0.2*x2 + rnormal() > 0.5)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Run psmatch2 first to get ground-truth values for comparison
psmatch2 treat x1 x2, outcome(y)
local psm2_att = r(att)
local psm2_att_se = r(seatt)
quietly count if _treated == 1 & _support == 1
local psm2_n_treated = r(N)
quietly count if _treated == 0 & _support == 1
local psm2_n_controls = r(N)
quietly count if _weight != . & _weight > 0 & _treated == 1
local psm2_n_matched = r(N)

* Clean up psmatch2 variables before running cpsmatch
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

cpsmatch treat x1 x2, outcome(y)

* Store ALL r() values before they get overwritten
local cps_N = r(N)
local cps_n_treated = r(n_treated)
local cps_n_controls = r(n_controls)
local cps_n_matched = r(n_matched)
local cps_n_off_support = r(n_off_support)
local cps_common_min = r(common_min)
local cps_common_max = r(common_max)
local cps_att = r(att)
local cps_att_se = r(att_se)
local cps_att_t = r(att_t)
local cps_method = "`r(method)'"
local cps_treatvar = "`r(treatvar)'"
local cps_outcome = "`r(outcome)'"

* Test variable creation
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

sum _pscore, meanonly
local pmin = r(min)
local pmax = r(max)
if `pmin' >= 0 & `pmax' <= 1 {
    test_pass "pscores in [0,1] range"
}
else {
    test_fail "pscores range" "min=`pmin', max=`pmax'"
}

* Test ALL scalar return values
local all_scalars_ok = 1

* r(N) - should equal number of observations
if `cps_N' != 500 {
    local all_scalars_ok = 0
}

* r(n_treated) + r(n_controls) should equal N (approximately, excluding missings)
if `cps_n_treated' <= 0 | `cps_n_controls' <= 0 {
    local all_scalars_ok = 0
}

* r(n_matched) should be positive
if `cps_n_matched' <= 0 {
    local all_scalars_ok = 0
}

* r(n_off_support) should be non-negative
if `cps_n_off_support' < 0 | `cps_n_off_support' == . {
    local all_scalars_ok = 0
}

* r(common_min) and r(common_max) should be valid pscores
if `cps_common_min' < 0 | `cps_common_min' > 1 | `cps_common_max' < 0 | `cps_common_max' > 1 {
    local all_scalars_ok = 0
}

* r(att) should be defined (we specified outcome)
if `cps_att' == . {
    local all_scalars_ok = 0
}

* r(att_se) should be positive
if `cps_att_se' <= 0 | `cps_att_se' == . {
    local all_scalars_ok = 0
}

* r(att_t) should be defined
if `cps_att_t' == . {
    local all_scalars_ok = 0
}

if `all_scalars_ok' {
    test_pass "returns all scalar values (N, n_treated, n_controls, n_matched, n_off_support, common_min, common_max, att, att_se, att_t)"
}
else {
    test_fail "scalar returns" "N=`cps_N' n_treated=`cps_n_treated' n_controls=`cps_n_controls' n_matched=`cps_n_matched' n_off=`cps_n_off_support' att=`cps_att'"
}

* Compare scalar return values against psmatch2 ground truth
assert_scalar_equal `cps_n_treated' `psm2_n_treated' $DEFAULT_SIGFIGS "r(n_treated) matches psmatch2"
assert_scalar_equal `cps_n_controls' `psm2_n_controls' $DEFAULT_SIGFIGS "r(n_controls) matches psmatch2"
assert_scalar_equal `cps_n_matched' `psm2_n_matched' $DEFAULT_SIGFIGS "r(n_matched) matches psmatch2"
assert_scalar_equal `cps_att' `psm2_att' $DEFAULT_SIGFIGS "r(att) matches psmatch2"

* Test local return values
local all_locals_ok = 1

* r(method) should be "nearest" (default)
if "`cps_method'" != "nearest" {
    local all_locals_ok = 0
}

* r(treatvar) should be "treat"
if "`cps_treatvar'" != "treat" {
    local all_locals_ok = 0
}

* r(outcome) should be "y"
if "`cps_outcome'" != "y" {
    local all_locals_ok = 0
}

if `all_locals_ok' {
    test_pass "returns all local values (method, treatvar, outcome)"
}
else {
    test_fail "local returns" "method=`cps_method' treatvar=`cps_treatvar' outcome=`cps_outcome'"
}

/*******************************************************************************
 * SECTION 3: psmatch2 comparison - Nearest Neighbor
 ******************************************************************************/
print_section "Nearest Neighbor (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("NN: vs psmatch2")

/*******************************************************************************
 * SECTION 4: psmatch2 comparison - Multiple Neighbors
 ******************************************************************************/
print_section "Multiple Neighbors (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(3) testname("NN(3): vs psmatch2")

/*******************************************************************************
 * SECTION 5: psmatch2 comparison - Kernel Matching
 ******************************************************************************/
print_section "Kernel Matching (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel testname("Kernel: vs psmatch2")

/*******************************************************************************
 * SECTION 6: psmatch2 comparison - Radius Matching
 ******************************************************************************/
print_section "Radius Matching (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) radius caliper(0.1) testname("Radius: vs psmatch2")

/*******************************************************************************
 * SECTION 7: psmatch2 comparison - Common Support
 ******************************************************************************/
print_section "Common Support (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) common testname("Common support: vs psmatch2")

/*******************************************************************************
 * SECTION 8: psmatch2 comparison - Without Replacement
 ******************************************************************************/
print_section "Without Replacement (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) noreplacement testname("Noreplacement: vs psmatch2")

/*******************************************************************************
 * SECTION 9: psmatch2 comparison - Logit
 ******************************************************************************/
print_section "Logit Estimation (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) logit testname("Logit: vs psmatch2")

/*******************************************************************************
 * SECTION 10: Large datasets
 ******************************************************************************/
print_section "Large Datasets"

clear
set obs 10000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("10K obs: vs psmatch2")

/*******************************************************************************
 * SECTION 11: if/in Conditions
 ******************************************************************************/
print_section "if/in Conditions"

clear
set obs 1000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()
gen group = mod(_n, 2)

benchmark_psmatch2 treat x1 x2 if group == 1, outcome(y) testname("if condition: vs psmatch2")
benchmark_psmatch2 treat x1 x2 in 1/400, outcome(y) testname("in condition: vs psmatch2")
benchmark_psmatch2 treat x1 x2 if group == 1 in 1/600, outcome(y) testname("if/in combined: vs psmatch2")

/*******************************************************************************
 * SECTION 13: Pre-estimated Propensity Score
 ******************************************************************************/
print_section "Pre-estimated Propensity Score"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Estimate propensity score externally
probit treat x1 x2
predict double pscore_ext, pr

benchmark_psmatch2 treat, pscore(pscore_ext) outcome(y) testname("pscore(): vs psmatch2")

/*******************************************************************************
 * SECTION 14: Kernel Types
 ******************************************************************************/
print_section "Kernel Types"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(normal) testname("kernel(normal): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(biweight) testname("kernel(biweight): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(uniform) testname("kernel(uniform): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(tricube) testname("kernel(tricube): vs psmatch2")

/*******************************************************************************
 * SECTION 15: Bandwidth
 ******************************************************************************/
print_section "Bandwidth"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel bwidth(0.08) testname("bwidth(0.08): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) kernel bwidth(0.03) testname("bwidth(0.03): vs psmatch2")

/*******************************************************************************
 * SECTION 16: Ties
 ******************************************************************************/
print_section "Ties"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) ties testname("ties: vs psmatch2")

/*******************************************************************************
 * SECTION 17: Descending
 ******************************************************************************/
print_section "Descending"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) noreplacement descending testname("descending: vs psmatch2")

/*******************************************************************************
 * SECTION 18: Combined Options
 ******************************************************************************/
print_section "Combined Options"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(3) common logit testname("neighbor(3) common logit: vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) kernel common testname("kernel common: vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) radius caliper(0.05) common testname("radius caliper(0.05) common: vs psmatch2")

/*******************************************************************************
 * SECTION 19: cpsmatch-only Features
 ******************************************************************************/
print_section "cpsmatch-only Features"

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
    test_fail "verbose option" "rc=`=_rc'"
}

capture drop _pscore _weight _treated _support _nn _id

capture cpsmatch treat x1, threads(2)
if _rc == 0 {
    test_pass "threads(2) option"
}
else {
    test_fail "threads option" "rc=`=_rc'"
}

/*******************************************************************************
 * SECTION 20: Default Epanechnikov Kernel (vs psmatch2)
 ******************************************************************************/
print_section "Epanechnikov Kernel (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) kernel kerneltype(epan) testname("kernel(epan): vs psmatch2")

/*******************************************************************************
 * SECTION 21: NN with Caliper (vs psmatch2)
 ******************************************************************************/
print_section "NN with Caliper (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) caliper(0.1) testname("NN caliper(0.1): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) caliper(0.01) testname("NN caliper(0.01): vs psmatch2")

/*******************************************************************************
 * SECTION 22: Multiple Covariates (vs psmatch2)
 ******************************************************************************/
print_section "Multiple Covariates (vs psmatch2)"

clear
set obs 1000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = runiform()
gen x5 = rnormal() * 2
gen treat = (0.2*x1 + 0.15*x2 + 0.1*x3 - 0.1*x4 + 0.05*x5 + rnormal() > 0)
gen y = 3*treat + x1 + 0.5*x2 - 0.3*x3 + x4 + rnormal()

benchmark_psmatch2 treat x1 x2 x3 x4 x5, outcome(y) testname("5 covariates: vs psmatch2")
benchmark_psmatch2 treat x1 x2 x3 x4 x5, outcome(y) kernel testname("5 covariates kernel: vs psmatch2")

/*******************************************************************************
 * SECTION 23: Many Neighbors (vs psmatch2)
 ******************************************************************************/
print_section "Many Neighbors (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(5) testname("NN(5): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(10) testname("NN(10): vs psmatch2")

/*******************************************************************************
 * SECTION 24: Unbalanced Treatment Groups (vs psmatch2)
 ******************************************************************************/
print_section "Unbalanced Treatment (vs psmatch2)"

* Many treated, few controls
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > -0.5)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("unbalanced (many treated): vs psmatch2")

* Few treated, many controls
clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 1.0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("unbalanced (few treated): vs psmatch2")

/*******************************************************************************
 * SECTION 25: Weight Variable Comparison (vs psmatch2)
 ******************************************************************************/
print_section "Weight Variable Comparison (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Comprehensive comparison via benchmark helper (checks _weight for controls too)
benchmark_psmatch2 treat x1 x2, outcome(y) testname("NN: weight comparison via benchmark")

/*******************************************************************************
 * SECTION 26: pscore() with Different Methods (vs psmatch2)
 ******************************************************************************/
print_section "pscore() with Methods (vs psmatch2)"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

* Estimate propensity score externally
logit treat x1 x2
predict double pscore_ext, pr

benchmark_psmatch2 treat, pscore(pscore_ext) outcome(y) kernel testname("pscore() kernel: vs psmatch2")
benchmark_psmatch2 treat, pscore(pscore_ext) outcome(y) radius caliper(0.1) testname("pscore() radius: vs psmatch2")

/*******************************************************************************
 * SECTION 27: Larger Dataset (vs psmatch2)
 ******************************************************************************/
print_section "Larger Dataset (vs psmatch2)"

clear
set obs 25000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + 0.1*x3 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 - 0.3*x3 + rnormal()

benchmark_psmatch2 treat x1 x2 x3, outcome(y) testname("25K obs: vs psmatch2")
benchmark_psmatch2 treat x1 x2 x3, outcome(y) kernel testname("25K kernel: vs psmatch2")

/*******************************************************************************
 * SECTION 28: benchmark_psmatch2 Helper Tests
 ******************************************************************************/
print_section "benchmark_psmatch2 Helper"

* Use the benchmark_psmatch2 helper for comprehensive comparison

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) testname("benchmark: NN default")

preserve
benchmark_psmatch2 treat x1 x2, outcome(y) neighbor(3) testname("benchmark: NN(3)")
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

preserve
benchmark_psmatch2 treat x1 x2, outcome(y) kernel testname("benchmark: kernel default")
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

preserve
benchmark_psmatch2 treat x1 x2, outcome(y) radius caliper(0.1) testname("benchmark: radius 0.1")
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

preserve
benchmark_psmatch2 treat x1 x2, outcome(y) logit testname("benchmark: logit")
restore
capture drop _pscore _weight _treated _support _nn _id _n1 _pdif

/*******************************************************************************
 * SECTION 29: Small Sample Edge Cases
 ******************************************************************************/
print_section "Small Sample Edge Cases"

* Minimal sample
clear
set obs 30
set seed 12345
gen x1 = rnormal()
gen treat = (_n > 15)
gen y = 2*treat + x1 + rnormal()

benchmark_psmatch2 treat x1, outcome(y) testname("small sample (N=30): vs psmatch2")

* Single covariate
clear
set obs 200
set seed 12345
gen x1 = rnormal()
gen treat = (0.5*x1 + rnormal() > 0)
gen y = 2*treat + x1 + rnormal()

benchmark_psmatch2 treat x1, outcome(y) testname("single covariate: vs psmatch2")

/*******************************************************************************
 * SECTION 30: Combined Options with benchmark_psmatch2
 ******************************************************************************/
print_section "Additional Combined Options"

clear
set obs 500
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen treat = (0.3*x1 + 0.2*x2 + rnormal() > 0)
gen y = 2*treat + x1 + 0.5*x2 + rnormal()

benchmark_psmatch2 treat x1 x2, outcome(y) noreplacement ties testname("noreplacement ties: vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) kernel bwidth(0.2) testname("kernel bwidth(0.2): vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) radius caliper(0.1) common logit testname("radius common logit: vs psmatch2")
benchmark_psmatch2 treat x1 x2, outcome(y) noreplacement descending common ties testname("descending common ties: vs psmatch2")

/*******************************************************************************
 * SECTION: Intentional Error Tests
 *
 * These tests verify that cpsmatch returns the same error codes as psmatch2
 * when given invalid inputs or error conditions.
 * Note: psmatch2 is a user-written command (ssc install psmatch2).
 * If psmatch2 is not installed, tests compare against expected behavior.
 ******************************************************************************/
print_section "Intentional Error Tests"

* Check if psmatch2 is installed
capture which psmatch2
local psmatch2_installed = (_rc == 0)

* Variable doesn't exist (expected rc=111: variable not found)
clear
set obs 100
gen treat = runiform() > 0.5
gen x = rnormal()
if `psmatch2_installed' {
    test_error_match, stata_cmd(psmatch2 treat nonexistent_var) ctools_cmd(cpsmatch treat nonexistent_var) testname("nonexistent variable")
}
else {
    capture cpsmatch treat nonexistent_var
    local ctools_rc = _rc
    if `ctools_rc' == 111 {
        test_pass "[error] nonexistent variable (rc=`ctools_rc') [psmatch2 not installed]"
    }
    else if `ctools_rc' != 0 {
        test_fail "[error] nonexistent variable" "expected rc=111, got rc=`ctools_rc' [psmatch2 not installed]"
    }
    else {
        test_fail "[error] nonexistent variable" "should have errored [psmatch2 not installed]"
    }
}

* Treatment variable doesn't exist (expected rc=111: variable not found)
clear
set obs 100
gen x = rnormal()
if `psmatch2_installed' {
    test_error_match, stata_cmd(psmatch2 nonexistent_treat x) ctools_cmd(cpsmatch nonexistent_treat x) testname("nonexistent treatment variable")
}
else {
    capture cpsmatch nonexistent_treat x
    local ctools_rc = _rc
    if `ctools_rc' == 111 {
        test_pass "[error] nonexistent treatment variable (rc=`ctools_rc') [psmatch2 not installed]"
    }
    else if `ctools_rc' != 0 {
        test_fail "[error] nonexistent treatment variable" "expected rc=111, got rc=`ctools_rc' [psmatch2 not installed]"
    }
    else {
        test_fail "[error] nonexistent treatment variable" "should have errored [psmatch2 not installed]"
    }
}

* String treatment variable (expected rc=109: type mismatch, or rc=198: syntax)
clear
set obs 100
gen str10 treat = cond(runiform() > 0.5, "yes", "no")
gen x = rnormal()
if `psmatch2_installed' {
    test_error_match, stata_cmd(psmatch2 treat x) ctools_cmd(cpsmatch treat x) testname("string treatment variable")
}
else {
    capture cpsmatch treat x
    local ctools_rc = _rc
    if `ctools_rc' == 109 | `ctools_rc' == 198 {
        test_pass "[error] string treatment variable (rc=`ctools_rc') [psmatch2 not installed]"
    }
    else if `ctools_rc' != 0 {
        test_fail "[error] string treatment variable" "expected rc=109 or 198, got rc=`ctools_rc' [psmatch2 not installed]"
    }
    else {
        test_fail "[error] string treatment variable" "should have errored [psmatch2 not installed]"
    }
}

* No covariates specified (expected rc=198: syntax error, or rc=100)
clear
set obs 100
gen treat = runiform() > 0.5
if `psmatch2_installed' {
    test_error_match, stata_cmd(psmatch2 treat) ctools_cmd(cpsmatch treat) testname("no covariates specified")
}
else {
    capture cpsmatch treat
    local ctools_rc = _rc
    if `ctools_rc' == 198 | `ctools_rc' == 100 {
        test_pass "[error] no covariates specified (rc=`ctools_rc') [psmatch2 not installed]"
    }
    else if `ctools_rc' != 0 {
        test_fail "[error] no covariates specified" "expected rc=198 or 100, got rc=`ctools_rc' [psmatch2 not installed]"
    }
    else {
        test_fail "[error] no covariates specified" "should have errored [psmatch2 not installed]"
    }
}

* All treated or all control (expected rc=459: no observations, or rc=198, or rc=2000+)
clear
set obs 100
gen treat = 1
gen x = rnormal()
if `psmatch2_installed' {
    test_error_match, stata_cmd(psmatch2 treat x) ctools_cmd(cpsmatch treat x) testname("all treated observations")
}
else {
    capture cpsmatch treat x
    local ctools_rc = _rc
    if `ctools_rc' == 459 | `ctools_rc' == 198 | `ctools_rc' == 2000 | `ctools_rc' == 2001 | `ctools_rc' == 480 {
        test_pass "[error] all treated observations (rc=`ctools_rc') [psmatch2 not installed]"
    }
    else if `ctools_rc' != 0 {
        test_fail "[error] all treated observations" "expected rc=459/198/2000/2001/480, got rc=`ctools_rc' [psmatch2 not installed]"
    }
    else {
        test_fail "[error] all treated observations" "should have errored [psmatch2 not installed]"
    }
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
noi print_summary "cpsmatch"
}
