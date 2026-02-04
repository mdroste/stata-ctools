/*******************************************************************************
 * test_bisect_crash.do - Bisect to find which tests cause the crash
 *
 * Run sections incrementally, then clear all to see when crash occurs.
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running sections 1-13 then clear all..."

* Check if ivreghdfe is installed
capture which ivreghdfe
if _rc != 0 {
    test_fail "ivreghdfe installation" "ivreghdfe not found"
    exit 0
}

* Plugin check
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign)
if _rc != 0 {
    exit 1
}

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
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) testname("20K dataset")

/*******************************************************************************
 * SECTION 9-12: Options tests
 ******************************************************************************/
print_section "Options Tests"
sysuse auto, clear
capture civreghdfe price (mpg = weight), absorb(foreign) first
capture civreghdfe price (mpg = weight), absorb(foreign) small
capture civreghdfe price (mpg = weight), absorb(foreign) tolerance(1e-10)
capture civreghdfe price (mpg = weight), absorb(foreign) maxiter(1000)
capture civreghdfe price (mpg = weight), absorb(foreign) verbose
capture civreghdfe price (mpg = weight), absorb(foreign) timeit

/*******************************************************************************
 * SECTION 13: if/in conditions
 ******************************************************************************/
print_section "if/in Conditions"
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) if price > 5000, absorb(foreign) testname("if condition")
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight) in 1/50, absorb(foreign) testname("in condition")

noi di as text _n "Sections 1-13 complete. Now trying clear all..."

* Cleanup
capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

}

* Now clear all - this should trigger the crash if the bug exists
di as text _n "Now running: clear all"
clear all

di as text "SUCCESS: No crash occurred after sections 1-13!"
