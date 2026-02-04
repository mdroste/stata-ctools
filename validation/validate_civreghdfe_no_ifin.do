/*******************************************************************************
 * validate_civreghdfe_no_ifin.do
 *
 * Same as validate_civreghdfe.do but WITHOUT if/in tests (Section 13)
 * Used to isolate whether the crash is caused by if/in filtering code path.
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running validation tests for civreghdfe (NO if/in tests)..."

* Check if ivreghdfe is installed
capture which ivreghdfe
if _rc != 0 {
    test_fail "ivreghdfe installation" "ivreghdfe not found"
    exit 0
}
test_pass "ivreghdfe is installed"

* Plugin check
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)
if _rc != 0 {
    test_fail "civreghdfe plugin load" "returned error `=_rc'"
    exit 1
}
test_pass "civreghdfe plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic IV tests (auto)
 ******************************************************************************/
print_section "Basic IV Tests (auto)"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("single endog, single instr")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("single endog, two instr")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) testname("with exog var")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length) turn displacement, absorb(foreign) testname("with multiple exog")

/*******************************************************************************
 * SECTION 3: Multiple endogenous variables
 ******************************************************************************/
print_section "Multiple Endogenous Variables"

sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement), absorb(foreign) testname("two endog")

sysuse auto, clear
benchmark_ivreghdfe price (mpg weight = length turn displacement headroom), absorb(foreign) testname("two endog, more instr")

/*******************************************************************************
 * SECTION 4: VCE options
 ******************************************************************************/
print_section "VCE Options"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("vce(robust)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("vce(cluster)")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("multi-instr + robust")

/*******************************************************************************
 * SECTION 5: Two-way fixed effects
 ******************************************************************************/
print_section "Two-Way Fixed Effects"

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) vce(robust) testname("two-way + robust")

/*******************************************************************************
 * SECTION 6: Census dataset
 ******************************************************************************/
print_section "Census Dataset"

sysuse census, clear
benchmark_ivreghdfe pop (medage = death), absorb(region) testname("census basic")

sysuse census, clear
benchmark_ivreghdfe pop (medage = death marriage), absorb(region) testname("census two instr")

sysuse census, clear
benchmark_ivreghdfe pop (medage = death) divorce, absorb(region) testname("census with exog")

/*******************************************************************************
 * SECTION 7: Panel data (nlswork)
 ******************************************************************************/
print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) testname("nlswork basic")

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age ttl_exp), absorb(idcode) testname("nlswork two instr")

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(robust) testname("nlswork robust")

webuse nlswork, clear
keep in 1/5000
benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(cluster idcode) testname("nlswork cluster")

/*******************************************************************************
 * SECTION 8: Large dataset
 ******************************************************************************/
print_section "Large Dataset"

clear
set seed 12345
set obs 20000
gen id = runiformint(1, 200)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* Use benchmark_ivreghdfe for full comparison (N, coefficients, VCE)
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) testname("20K dataset")

/*******************************************************************************
 * SECTION 9: first option
 ******************************************************************************/
print_section "first Option"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) first
if _rc == 0 {
    test_pass "first option accepted"
}
else {
    test_fail "first option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 10: small option
 ******************************************************************************/
print_section "small Option"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) small
if _rc == 0 {
    test_pass "small option accepted"
}
else {
    test_fail "small option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: tolerance/maxiter options
 ******************************************************************************/
print_section "tolerance/maxiter Options"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10)
if _rc == 0 {
    test_pass "tolerance(1e-10) accepted"
}
else {
    test_fail "tolerance option" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) maxiter(1000)
if _rc == 0 {
    test_pass "maxiter(1000) accepted"
}
else {
    test_fail "maxiter option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 12: verbose/timeit options
 ******************************************************************************/
print_section "verbose/timeit Options"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
if _rc == 0 {
    test_pass "verbose option accepted"
}
else {
    test_fail "verbose option" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) timeit
if _rc == 0 {
    test_pass "timeit option accepted"
}
else {
    test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 13: if/in conditions - SKIPPED
 ******************************************************************************/
print_section "if/in Conditions - SKIPPED"
noi di as text "  [SKIPPED] if/in tests excluded from this validation file"

/*******************************************************************************
 * SECTION 14: Coefficient comparison
 ******************************************************************************/
print_section "Coefficient Comparison"

sysuse auto, clear

ivreghdfe price (mpg = weight length), absorb(foreign)
matrix ivreghdfe_b = e(b)

civreghdfe price (mpg = weight length), absorb(foreign)
matrix civreghdfe_b = e(b)

* Compare coefficients using sigfigs
matrix_min_sigfigs ivreghdfe_b civreghdfe_b
local min_sf = r(min_sigfigs)

if `min_sf' >= $DEFAULT_SIGFIGS {
    local sf_fmt : display %4.1f `min_sf'
    test_pass "coefficients match (sigfigs=`sf_fmt')"
}
else {
    local sf_fmt : display %4.1f `min_sf'
    test_fail "coefficients" "sigfigs=`sf_fmt' below threshold"
}

/*******************************************************************************
 * SECTION 15: Weights
 ******************************************************************************/
print_section "Weights"

sysuse auto, clear
capture civreghdfe price (mpg = weight) [aw=weight], absorb(foreign)
if _rc == 0 {
    test_pass "aweight accepted"
}
else {
    test_fail "aweight" "returned error `=_rc'"
}

sysuse auto, clear
gen int fw = ceil(mpg/5)
capture civreghdfe price (mpg = weight) [fw=fw], absorb(foreign)
if _rc == 0 {
    test_pass "fweight accepted"
}
else {
    test_fail "fweight" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight) [pw=weight], absorb(foreign)
if _rc == 0 {
    test_pass "pweight accepted"
}
else {
    test_fail "pweight" "returned error `=_rc'"
}

/*******************************************************************************
 * Summary - Keep it short for debugging
 ******************************************************************************/

* Cleanup to prevent state corruption before next validation script
capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

* End of civreghdfe validation (no if/in)
noi print_summary "civreghdfe (no if/in)"
}
