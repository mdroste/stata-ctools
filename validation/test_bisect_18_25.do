/*******************************************************************************
 * test_bisect_18_25.do - Test sections 18-25 (later sections only)
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

quietly {

noi di as text "Running sections 18-25 then clear all..."

/*******************************************************************************
 * SECTION 18: Diagnostic test options
 ******************************************************************************/
print_section "Diagnostic Test Options"
scalar drop _all
matrix drop _all
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3)
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) orthog(z3)
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) noid

drop x_endog y
gen x_endog1 = 0.5*z1 + 0.3*z2 + rnormal()
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
gen y = 2*x_endog1 + 1.5*x_endog2 + 1.0*x_exog + rnormal()

capture civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2)
civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) endogtest(x_endog2)
capture civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3)
civreghdfe y (x_endog1 x_endog2 = z1 z2 z3) x_exog, absorb(firm) redundant(z3)

/*******************************************************************************
 * SECTION 19: partial() option
 ******************************************************************************/
print_section "partial() option"
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog1 = runiform()
gen x_exog2 = rnormal()
gen y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + rnormal()

capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm)
local b_full = _b[x_endog]
civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2, absorb(firm) partial(x_exog2)
local b_partial = _b[x_endog]
gen x_exog3 = rnormal()
replace y = 2*x_endog + 1.0*x_exog1 + 0.5*x_exog2 + 0.3*x_exog3 + rnormal()
capture civreghdfe y (x_endog = z1 z2) x_exog1 x_exog2 x_exog3, absorb(firm) partial(x_exog2 x_exog3)

/*******************************************************************************
 * SECTION 20: ffirst option
 ******************************************************************************/
print_section "ffirst option"
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.0*x_exog + rnormal()

capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst
civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) ffirst
gen x_endog2 = 0.4*z2 + 0.3*z3 + rnormal()
replace y = 2*x_endog + 1.5*x_endog2 + 1.0*x_exog + rnormal()
civreghdfe y (x_endog x_endog2 = z1 z2 z3) x_exog, absorb(firm) ffirst

/*******************************************************************************
 * SECTION 21: Save options
 ******************************************************************************/
print_section "Save options"
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.0*x_exog + rnormal()

capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst
capture estimates restore _civreghdfe_main
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) savefirst savefprefix(test_)
capture estimates restore test_main
capture civreghdfe y (x_endog = z1 z2) x_exog, absorb(firm) saverf

/*******************************************************************************
 * SECTION 22: Factor variables
 ******************************************************************************/
print_section "Factor Variables"
sysuse auto, clear
capture civreghdfe price (mpg = weight) i.rep78, absorb(foreign)
capture ivreghdfe price (mpg = weight) i.rep78, absorb(foreign)
if _rc == 0 {
    matrix ivreghdfe_b = e(b)
    local ivreghdfe_N = e(N)
    civreghdfe price (mpg = weight) i.rep78, absorb(foreign)
}
sysuse auto, clear
capture civreghdfe price (mpg = i.rep78), absorb(foreign)
sysuse auto, clear
capture civreghdfe price (mpg = weight length) i.rep78 i.foreign, absorb(headroom)
sysuse auto, clear
capture civreghdfe price (mpg = weight) i.rep78##c.turn, absorb(foreign)
sysuse auto, clear
capture civreghdfe price (mpg = weight length) ib3.rep78, absorb(foreign)
sysuse auto, clear
capture civreghdfe price (mpg = weight) c.turn#i.foreign, absorb(rep78)
webuse nlswork, clear
keep in 1/5000
capture civreghdfe ln_wage (tenure = age) i.race#i.union, absorb(idcode)
sysuse auto, clear
capture ivreghdfe price (mpg = weight length) c.turn#i.foreign, absorb(rep78)
if _rc == 0 {
    local ivreghdfe_N = e(N)
    capture civreghdfe price (mpg = weight length) c.turn#i.foreign, absorb(rep78)
}

/*******************************************************************************
 * SECTION 23: Time series operators
 ******************************************************************************/
print_section "Time Series Operators"
clear
set seed 54321
set obs 500
gen id = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
tsset id time
gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

capture civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = z1 z2) D.x_exog, absorb(id)
capture civreghdfe y (x_endog = F.z1 z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = L2.z1 z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = L(1/2).z1) x_exog, absorb(id)
capture ivreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id)
if _rc == 0 {
    local ivreghdfe_N = e(N)
    civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id)
}
capture civreghdfe y (x_endog = L.z1 D.z2) x_exog, absorb(id)
capture civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(robust)
capture civreghdfe y (x_endog = L.z1 z2) x_exog, absorb(id) vce(cluster id)

/*******************************************************************************
 * SECTION 24: IV Estimator Options
 ******************************************************************************/
print_section "IV Estimator Options"
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) liml testname("liml")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(1) testname("fuller(1)")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) fuller(4) testname("fuller(4)")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(0.5) testname("kclass(0.5)")
ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) kclass(1)
local kclass_b = _b[x_endog]
ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm)
local tsls_b = _b[x_endog]
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s testname("gmm2s")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(robust) testname("gmm2s + robust")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue testname("cue")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) coviv testname("coviv")
matrix b0 = (2, 1.5)
capture civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) b0(b0)

noi di as text _n "Sections 18-24 complete. Now trying clear all..."

capture estimates drop _all
capture scalar drop _all
capture matrix drop _all
capture clear

}

di as text _n "Now running: clear all"
clear all

di as text "SUCCESS: No crash after sections 18-24!"
