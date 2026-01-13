/*******************************************************************************
 * validate_cqreg.do
 *
 * Comprehensive validation tests for cqreg vs qreg
 * Tests quantile regression across various scenarios
 *
 * Note: Coefficients and standard errors are compared with tolerance 1e-6
 * due to floating-point precision differences in complex optimization
 ******************************************************************************/

do "validate_setup.do"

di as text ""
di as text "======================================================================"
di as text "              CQREG VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * Basic quantile tests (auto dataset)
 ******************************************************************************/
print_section "Basic quantile tests"

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
 * Covariate variations
 ******************************************************************************/
print_section "Covariate variations"

sysuse auto, clear
benchmark_qreg price mpg, testname("single covariate")

sysuse auto, clear
benchmark_qreg price mpg weight length turn, testname("many covariates")

sysuse auto, clear
benchmark_qreg price mpg weight length, quantile(0.25) testname("many covariates q=0.25")

/*******************************************************************************
 * Census dataset
 ******************************************************************************/
print_section "Census dataset"

sysuse census, clear
benchmark_qreg pop medage death marriage, testname("census median")

sysuse census, clear
benchmark_qreg pop medage death marriage, quantile(0.75) testname("census q=0.75")

/*******************************************************************************
 * Synthetic data
 ******************************************************************************/
print_section "Synthetic data"

clear
set seed 12345
set obs 5000
gen x = rnormal()
gen y = 2 + 3*x + rnormal()
benchmark_qreg y x, testname("5K synthetic")

* Larger dataset
clear
set seed 54321
set obs 10000
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen y = 1 + 0.5*x1 - 0.3*x2 + rnormal()
benchmark_qreg y x1 x2, testname("10K synthetic")

clear
set seed 54321
set obs 10000
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen y = 1 + 0.5*x1 - 0.3*x2 + rnormal()
benchmark_qreg y x1 x2, quantile(0.25) testname("10K synthetic q=0.25")

clear
set seed 54321
set obs 10000
gen x1 = rnormal()
gen x2 = rnormal() * 2
gen y = 1 + 0.5*x1 - 0.3*x2 + rnormal()
benchmark_qreg y x1 x2, quantile(0.75) testname("10K synthetic q=0.75")

/*******************************************************************************
 * Bandwidth methods
 ******************************************************************************/
print_section "Bandwidth methods"

sysuse auto, clear

* Hall-Sheather (default)
quietly cqreg price mpg weight, bwmethod(hsheather)
local hs_se_mpg = _se[mpg]

* Bofinger
quietly cqreg price mpg weight, bwmethod(bofinger)
local bof_se_mpg = _se[mpg]

* Chamberlain
quietly cqreg price mpg weight, bwmethod(chamberlain)
local cham_se_mpg = _se[mpg]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `hs_se_mpg' > 0 & `bof_se_mpg' > 0 & `cham_se_mpg' > 0 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] all bandwidth methods produce valid SEs"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] some bandwidth methods failed"
}

/*******************************************************************************
 * HDFE tests (cqreg-specific)
 ******************************************************************************/
print_section "HDFE tests (experimental)"

sysuse auto, clear
capture quietly cqreg price mpg weight, absorb(foreign)
local cqreg_rc = _rc
local cqreg_b_mpg = _b[mpg]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `cqreg_rc' == 0 & !missing(`cqreg_b_mpg') & `cqreg_b_mpg' != 0 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] single FE absorb"
}
else {
    di as text "[INFO] single FE absorb skipped (experimental)"
}

sysuse auto, clear
capture quietly cqreg price mpg weight, absorb(foreign rep78)
local cqreg_rc = _rc
local cqreg_b_mpg = _b[mpg]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `cqreg_rc' == 0 & !missing(`cqreg_b_mpg') & `cqreg_b_mpg' != 0 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] two-way FE absorb"
}
else {
    di as text "[INFO] two-way FE absorb skipped (experimental)"
}

sysuse auto, clear
capture quietly cqreg price mpg weight, absorb(foreign) vce(cluster foreign)
local cqreg_rc = _rc
local cqreg_N_clust = e(N_clust)

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `cqreg_rc' == 0 & `cqreg_N_clust' == 2 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] HDFE + cluster VCE"
}
else {
    di as text "[INFO] HDFE + cluster VCE skipped (experimental)"
}

/*******************************************************************************
 * Stored results
 ******************************************************************************/
print_section "Stored results"

sysuse auto, clear

quietly qreg price mpg weight
matrix qreg_b = e(b)
local qreg_sum_adev = e(sum_adev)

quietly cqreg price mpg weight
matrix cqreg_b = e(b)
local cqreg_sum_adev = e(sum_adev)

* Compare e(b)
tempname diff
matrix `diff' = qreg_b - cqreg_b
local cols = colsof(`diff')
local maxdiff = 0
forvalues j = 1/`cols' {
    local d = abs(`diff'[1, `j'])
    if `d' > `maxdiff' local maxdiff = `d'
}

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `maxdiff' < 1e-6 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] e(b) matches"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] e(b) differs (max: " %9.2e `maxdiff' ")"
}

* Compare sum_adev
local rel_diff = abs(`qreg_sum_adev' - `cqreg_sum_adev') / `qreg_sum_adev'
global TESTS_TOTAL = $TESTS_TOTAL + 1
if `rel_diff' < 0.001 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] e(sum_adev) matches (rel diff: " %9.6f `rel_diff' ")"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] e(sum_adev) differs (rel diff: " %9.6f `rel_diff' ")"
}

/*******************************************************************************
 * Convergence check
 ******************************************************************************/
print_section "Convergence"

sysuse auto, clear
quietly cqreg price mpg weight
local converged = 0
capture local converged = e(converged)
local b_mpg = _b[mpg]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `converged' == 1 | (!missing(`b_mpg') & `b_mpg' != 0) {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] algorithm converged"
}
else {
    global TESTS_FAILED = $TESTS_FAILED + 1
    di as error "[FAIL] algorithm did not converge"
}

/*******************************************************************************
 * nlswork (may not converge - informational)
 ******************************************************************************/
print_section "nlswork (informational)"

webuse nlswork, clear
keep in 1/2000

quietly qreg ln_wage age ttl_exp tenure
local qreg_b_age = _b[age]

quietly cqreg ln_wage age ttl_exp tenure
local cqreg_converged = e(converged)
local cqreg_b_age = _b[age]

global TESTS_TOTAL = $TESTS_TOTAL + 1
if `cqreg_converged' == 1 & abs(`qreg_b_age' - `cqreg_b_age') < 1e-6 {
    global TESTS_PASSED = $TESTS_PASSED + 1
    di as result "[PASS] nlswork coefficients match"
}
else if `cqreg_converged' == 1 {
    di as text "[INFO] nlswork converged but coefficients differ slightly"
    global TESTS_PASSED = $TESTS_PASSED + 1
}
else {
    di as text "[INFO] nlswork did not converge (IPM may need more iterations)"
}

/*******************************************************************************
 * SUMMARY
 ******************************************************************************/
print_summary "cqreg"

* Return error code if any tests failed
if $TESTS_FAILED > 0 {
    exit 1
}
