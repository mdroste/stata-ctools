/*******************************************************************************
 * validate_cqreg.do
 *
 * Comprehensive validation tests for cqreg vs qreg
 * Tests all options: quantile, vce, denmethod, bwmethod, tolerance, maxiter
 ******************************************************************************/

do "validate_setup.do"

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
 * SECTION 4: VCE options
 ******************************************************************************/
noi print_section "VCE Options"

sysuse auto, clear
noi benchmark_qreg price mpg weight, vce(robust) testname("vce(robust)")

* Note: cqreg may have different VCE options than Stata's qreg
sysuse auto, clear
capture cqreg price mpg weight, vce(iid)
if _rc == 0 {
    noi test_pass "vce(iid) accepted"
}
else {
    noi test_fail "vce(iid)" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 5: Census dataset
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
noi benchmark_qreg pop medage death, testname("census median")

sysuse census, clear
noi benchmark_qreg pop medage death, quantile(0.25) testname("census q=0.25")

sysuse census, clear
noi benchmark_qreg pop medage death, quantile(0.75) testname("census q=0.75")

/*******************************************************************************
 * SECTION 6: nlswork panel data
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, testname("nlswork median")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, quantile(0.25) testname("nlswork q=0.25")

webuse nlswork, clear
keep in 1/5000
noi benchmark_qreg ln_wage age tenure, quantile(0.75) testname("nlswork q=0.75")

/*******************************************************************************
 * SECTION 7: Large dataset
 ******************************************************************************/
noi print_section "Large Dataset"

clear
set seed 12345
set obs 10000
gen x1 = runiform()
gen x2 = rnormal()
gen x3 = runiformint(1, 100)
gen y = 2*x1 + 3*x2 - 0.5*x3 + rnormal()

noi benchmark_qreg y x1 x2 x3, testname("10K median")

noi benchmark_qreg y x1 x2 x3, quantile(0.25) testname("10K q=0.25")

noi benchmark_qreg y x1 x2 x3, quantile(0.75) testname("10K q=0.75")

/*******************************************************************************
 * SECTION 8: denmethod option
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
 * SECTION 9: bwmethod option
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
 * SECTION 10: tolerance/maxiter options
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
capture cqreg price mpg weight, maxiter(500)
if _rc == 0 {
    noi test_pass "maxiter(500) accepted"
}
else {
    noi test_fail "maxiter option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: verbose/timeit options
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
 * SECTION 12: if/in conditions
 ******************************************************************************/
noi print_section "if/in Conditions"

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

/*******************************************************************************
 * SECTION 13: absorb option (experimental)
 ******************************************************************************/
noi print_section "absorb Option (experimental)"

sysuse auto, clear
capture cqreg price mpg weight, absorb(foreign)
if _rc == 0 {
    noi test_pass "absorb option accepted"
}
else {
    noi test_fail "absorb option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 14: Extreme quantiles
 ******************************************************************************/
noi print_section "Extreme Quantiles"

sysuse auto, clear
capture cqreg price mpg weight, quantile(0.05)
if _rc == 0 {
    noi test_pass "q=0.05 accepted"
}
else {
    noi test_fail "q=0.05" "returned error `=_rc'"
}

sysuse auto, clear
capture cqreg price mpg weight, quantile(0.95)
if _rc == 0 {
    noi test_pass "q=0.95 accepted"
}
else {
    noi test_fail "q=0.95" "returned error `=_rc'"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "cqreg"

if $TESTS_FAILED > 0 {
    exit 1
}

}
