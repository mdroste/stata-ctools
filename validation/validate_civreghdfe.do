/*******************************************************************************
 * validate_civreghdfe.do
 *
 * Comprehensive validation tests for civreghdfe vs ivreghdfe
 * Tests all options: absorb, vce, weights, first, small, tolerance, maxiter
 *
 * Requires ivreghdfe to be installed (ssc install ivreghdfe)
 ******************************************************************************/

* Load setup (works from project root or validation dir)
capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

di as text ""
di as text "======================================================================"
di as text "              CIVREGHDFE VALIDATION TEST SUITE"
di as text "======================================================================"

/*******************************************************************************
 * Check if ivreghdfe is installed
 ******************************************************************************/
noi print_section "Prerequisite Check"

capture which ivreghdfe
if _rc != 0 {
    noi di as error "ivreghdfe not installed. Please run: ssc install ivreghdfe"
    noi di as error "Skipping civreghdfe validation tests"
    noi test_fail "ivreghdfe installation" "ivreghdfe not found"
    noi print_summary "civreghdfe"
    exit 0
}
noi test_pass "ivreghdfe is installed"

/*******************************************************************************
 * SECTION 1: Plugin check
 ******************************************************************************/
noi print_section "Plugin Check"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)
if _rc != 0 {
    noi test_fail "civreghdfe plugin load" "returned error `=_rc'"
    noi print_summary "civreghdfe"
    exit 1
}
noi test_pass "civreghdfe plugin loads and runs"

/*******************************************************************************
 * SECTION 2: Basic IV tests (auto)
 ******************************************************************************/
noi print_section "Basic IV Tests (auto)"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("single endog, single instr")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) testname("single endog, two instr")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight) turn, absorb(foreign) testname("with exog var")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight length) turn displacement, absorb(foreign) testname("with multiple exog")

/*******************************************************************************
 * SECTION 3: Multiple endogenous variables
 ******************************************************************************/
noi print_section "Multiple Endogenous Variables"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg weight = length turn displacement), absorb(foreign) testname("two endog")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg weight = length turn displacement headroom), absorb(foreign) testname("two endog, more instr")

/*******************************************************************************
 * SECTION 4: VCE options
 ******************************************************************************/
noi print_section "VCE Options"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("vce(robust)")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("vce(cluster)")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) testname("multi-instr + robust")

/*******************************************************************************
 * SECTION 5: Two-way fixed effects
 ******************************************************************************/
noi print_section "Two-Way Fixed Effects"

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) testname("two-way FE")

sysuse auto, clear
noi benchmark_ivreghdfe price (mpg = weight), absorb(foreign rep78) vce(robust) testname("two-way + robust")

/*******************************************************************************
 * SECTION 6: Census dataset
 ******************************************************************************/
noi print_section "Census Dataset"

sysuse census, clear
noi benchmark_ivreghdfe pop (medage = death), absorb(region) testname("census basic")

sysuse census, clear
noi benchmark_ivreghdfe pop (medage = death marriage), absorb(region) testname("census two instr")

sysuse census, clear
noi benchmark_ivreghdfe pop (medage = death) divorce, absorb(region) testname("census with exog")

/*******************************************************************************
 * SECTION 7: Panel data (nlswork)
 ******************************************************************************/
noi print_section "Panel Data (nlswork)"

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) testname("nlswork basic")

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age ttl_exp), absorb(idcode) testname("nlswork two instr")

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(robust) testname("nlswork robust")

webuse nlswork, clear
keep in 1/5000
noi benchmark_ivreghdfe ln_wage (tenure = age), absorb(idcode) vce(cluster idcode) testname("nlswork cluster")

/*******************************************************************************
 * SECTION 8: Large dataset
 ******************************************************************************/
noi print_section "Large Dataset"

clear
set seed 12345
set obs 20000
gen id = runiformint(1, 200)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id)
matrix ivreghdfe_b = e(b)
local ivreghdfe_N = e(N)

civreghdfe y (x_endog = z1 z2) x_exog, absorb(id)
matrix civreghdfe_b = e(b)
local civreghdfe_N = e(N)

if `ivreghdfe_N' == `civreghdfe_N' {
    noi test_pass "20K dataset: N matches"
}
else {
    noi test_fail "20K dataset" "N differs"
}

/*******************************************************************************
 * SECTION 9: first option
 ******************************************************************************/
noi print_section "first Option"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) first
if _rc == 0 {
    noi test_pass "first option accepted"
}
else {
    noi test_fail "first option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 10: small option
 ******************************************************************************/
noi print_section "small Option"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) small
if _rc == 0 {
    noi test_pass "small option accepted"
}
else {
    noi test_fail "small option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 11: tolerance/maxiter options
 ******************************************************************************/
noi print_section "tolerance/maxiter Options"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10)
if _rc == 0 {
    noi test_pass "tolerance(1e-10) accepted"
}
else {
    noi test_fail "tolerance option" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) maxiter(1000)
if _rc == 0 {
    noi test_pass "maxiter(1000) accepted"
}
else {
    noi test_fail "maxiter option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 12: verbose/timeit options
 ******************************************************************************/
noi print_section "verbose/timeit Options"

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
if _rc == 0 {
    noi test_pass "verbose option accepted"
}
else {
    noi test_fail "verbose option" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) timeit
if _rc == 0 {
    noi test_pass "timeit option accepted"
}
else {
    noi test_fail "timeit option" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 13: if/in conditions
 ******************************************************************************/
noi print_section "if/in Conditions"

sysuse auto, clear
ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign)
local ivreghdfe_N = e(N)

civreghdfe price (mpg = weight) if price > 5000, absorb(foreign)
local civreghdfe_N = e(N)

if `ivreghdfe_N' == `civreghdfe_N' {
    noi test_pass "if condition: N matches"
}
else {
    noi test_fail "if condition" "N differs"
}

sysuse auto, clear
ivreghdfe price (mpg = weight) in 1/50, absorb(foreign)
local ivreghdfe_N = e(N)

civreghdfe price (mpg = weight) in 1/50, absorb(foreign)
local civreghdfe_N = e(N)

if `ivreghdfe_N' == `civreghdfe_N' {
    noi test_pass "in condition: N matches"
}
else {
    noi test_fail "in condition" "N differs"
}

/*******************************************************************************
 * SECTION 14: Coefficient comparison
 ******************************************************************************/
noi print_section "Coefficient Comparison"

sysuse auto, clear

ivreghdfe price (mpg = weight length), absorb(foreign)
matrix ivreghdfe_b = e(b)

civreghdfe price (mpg = weight length), absorb(foreign)
matrix civreghdfe_b = e(b)

tempname diff
matrix `diff' = ivreghdfe_b - civreghdfe_b
local cols = colsof(`diff')
local maxdiff = 0
forvalues j = 1/`cols' {
    local d = abs(`diff'[1, `j'])
    if `d' > `maxdiff' local maxdiff = `d'
}

if `maxdiff' < 1e-6 {
    noi test_pass "coefficients match (maxdiff=`maxdiff')"
}
else {
    noi test_fail "coefficients" "maxdiff=`maxdiff' exceeds tolerance"
}

/*******************************************************************************
 * SECTION 15: Weights
 ******************************************************************************/
noi print_section "Weights"

sysuse auto, clear
capture civreghdfe price (mpg = weight) [aw=weight], absorb(foreign)
if _rc == 0 {
    noi test_pass "aweight accepted"
}
else {
    noi test_fail "aweight" "returned error `=_rc'"
}

sysuse auto, clear
gen int fw = ceil(mpg/5)
capture civreghdfe price (mpg = weight) [fw=fw], absorb(foreign)
if _rc == 0 {
    noi test_pass "fweight accepted"
}
else {
    noi test_fail "fweight" "returned error `=_rc'"
}

sysuse auto, clear
capture civreghdfe price (mpg = weight) [pw=weight], absorb(foreign)
if _rc == 0 {
    noi test_pass "pweight accepted"
}
else {
    noi test_fail "pweight" "returned error `=_rc'"
}

/*******************************************************************************
 * SECTION 16: Display options
 ******************************************************************************/
noi print_section "Display Options"

* Test level()
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) level(90)
if e(level) == 90 {
    noi test_pass "level(90) sets correct confidence level"
}
else {
    noi test_fail "level(90)" "e(level) = `e(level)' instead of 90"
}

* Test noheader - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader
if _rc == 0 {
    noi test_pass "noheader option accepted"
}
else {
    noi test_fail "noheader" "returned error `=_rc'"
}

* Test nofooter - just verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nofooter
if _rc == 0 {
    noi test_pass "nofooter option accepted"
}
else {
    noi test_fail "nofooter" "returned error `=_rc'"
}

* Test nooutput - verify it runs and stores results
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) nooutput
if _rc == 0 & e(N) > 0 {
    noi test_pass "nooutput option accepted and stores e(N)"
}
else {
    noi test_fail "nooutput" "returned error or missing e(N)"
}

* Test title() - verify it runs without error
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) title("Custom Title Test")
if _rc == 0 {
    noi test_pass "title() option accepted"
}
else {
    noi test_fail "title()" "returned error `=_rc'"
}

* Test depname() - verify it sets e(depvar)
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) depname("outcome_var")
if "`e(depvar)'" == "outcome_var" {
    noi test_pass "depname() sets e(depvar) correctly"
}
else {
    noi test_fail "depname()" "e(depvar) = `e(depvar)' instead of outcome_var"
}

* Test noid - verify underid test is suppressed but e(idstat) still stored
sysuse auto, clear
civreghdfe price (mpg = weight length), absorb(foreign) noid
if e(idstat) != . {
    noi test_pass "noid suppresses display but e(idstat) still computed"
}
else {
    noi test_fail "noid" "e(idstat) not stored"
}

* Test combined options
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) noheader nofooter level(99)
if _rc == 0 & e(level) == 99 {
    noi test_pass "combined display options work together"
}
else {
    noi test_fail "combined options" "returned error or wrong level"
}

/*******************************************************************************
 * Summary
 ******************************************************************************/

noi print_summary "civreghdfe"

if $TESTS_FAILED > 0 {
    exit 1
}

}
