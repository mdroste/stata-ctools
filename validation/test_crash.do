* Test crash: find which VCE option causes heap corruption
clear all
set more off
adopath ++ "build"
adopath ++ "../build"

* Generate test data once
set seed 12345
clear
set obs 200
gen id = mod(_n - 1, 20) + 1
gen time = ceil(_n / 20)
gen x_exog = rnormal()
gen z1 = rnormal()
gen z2 = rnormal()
gen x_endog = 0.5 * z1 + 0.3 * z2 + rnormal()
gen y = 1 + 0.5 * x_endog + 0.3 * x_exog + rnormal()
xtset id time
tempfile testdata
save `testdata'

* Run civreghdfe with various VCE options that might cause issues
di "Test 1: basic..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id)

di "Test 2: robust..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) vce(robust)

di "Test 3: cluster..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) vce(cluster id)

di "Test 4: two-way cluster..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) vce(cluster id time)

di "Test 5: kernel bartlett..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(bartlett) bw(2)

di "Test 6: kernel parzen..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)

di "Test 7: dkraay..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2)

di "Test 8: kiefer..."
quietly civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer

di "All VCE tests complete"

* Now run cdecode
di "Running cdecode..."
sysuse auto, clear
cdecode foreign, gen(foreign_str)

di "All completed successfully!"
