/*******************************************************************************
 * test_failure.do
 * Test benchmark_ivreghdfe calls (runs both ivreghdfe + civreghdfe)
 ******************************************************************************/

clear all
set more off
adopath ++ "build"

capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

di as text "Testing benchmark_ivreghdfe calls..."

quietly {

* Basic tests
sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) testname("basic")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(robust) testname("robust")

sysuse auto, clear
benchmark_ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) testname("cluster")

* HAC tests
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

benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(bartlett) bw(2) testname("kernel")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2) testname("dkraay")
benchmark_ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer testname("kiefer")

* IV estimator tests
clear
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) liml testname("liml")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s testname("gmm2s")
benchmark_ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) cue testname("cue")

}

di as text "benchmark_ivreghdfe tests done. Running clear all..."
clear all
di as text "NO CRASH"
