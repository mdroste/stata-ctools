* Test noreturn option for civreghdfe
* This tests the hypothesis that SF_scal_save/SF_mat_store are causing heap corruption

clear all
set more off
adopath ++ "build"

di "Testing civreghdfe with noreturn option..."

* Generate test data
set seed 12345
set obs 500
gen id = ceil(_n / 10)
gen time = mod(_n - 1, 10) + 1
gen x1 = rnormal()
gen x2 = rnormal()
gen z1 = rnormal()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.5*rnormal()
gen y = 1 + 0.5*x_endog + 0.3*x1 + 0.2*x2 + rnormal()

* Run civreghdfe with noreturn option
di "Running civreghdfe with noreturn..."
civreghdfe y (x_endog = z1 z2) x1 x2, absorb(id) vce(cluster id) noreturn

di "civreghdfe with noreturn completed successfully"

* Now run cdecode to test if heap is corrupted
di "Running cdecode to test heap integrity..."
gen id_label = "group_" + string(id)
encode id_label, gen(id_encoded)
cdecode id_encoded, gen(id_decoded)

di "cdecode completed successfully - no heap corruption!"

* Clean up
drop id_label id_encoded id_decoded

di ""
di "=============================================="
di "TEST PASSED: noreturn option works correctly"
di "=============================================="
