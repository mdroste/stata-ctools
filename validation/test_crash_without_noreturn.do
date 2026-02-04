* Test WITHOUT noreturn - this should crash if the hypothesis is correct
* Run this to verify the crash still happens without noreturn

clear all
set more off
adopath ++ "build"

di "=============================================="
di "Testing civreghdfe WITHOUT noreturn option"
di "Running 50 iterations (expect crash)..."
di "=============================================="

* Generate test data
set seed 12345
set obs 1000
gen id = ceil(_n / 20)
gen time = mod(_n - 1, 20) + 1
gen x1 = rnormal()
gen x2 = rnormal()
gen z1 = rnormal()
gen z2 = rnormal()
gen z3 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + 0.5*rnormal()
gen x_endog2 = 0.3*z1 + 0.4*z2 + 0.3*z3 + 0.5*rnormal()
gen y = 1 + 0.5*x_endog + 0.3*x_endog2 + 0.2*x1 + 0.1*x2 + rnormal()
gen weight = 0.5 + uniform()

* Run 50 iterations WITHOUT noreturn
forval i = 1/50 {
    if mod(`i', 10) == 0 {
        di "Iteration `i'..."
    }

    * Vary specifications to stress test - same as noreturn test but without noreturn
    if mod(`i', 5) == 0 {
        quietly civreghdfe y (x_endog = z1 z2) x1, absorb(id)
    }
    else if mod(`i', 5) == 1 {
        quietly civreghdfe y (x_endog = z1 z2 z3) x1 x2, absorb(id) vce(robust)
    }
    else if mod(`i', 5) == 2 {
        quietly civreghdfe y (x_endog = z1 z2) x1, absorb(id) vce(cluster id)
    }
    else if mod(`i', 5) == 3 {
        quietly civreghdfe y (x_endog x_endog2 = z1 z2 z3) x1, absorb(id)
    }
    else {
        quietly civreghdfe y (x_endog = z1 z2) x1 x2, absorb(id time)
    }
}

di "All 50 civreghdfe iterations completed WITHOUT noreturn"

* Now run cdecode to test if heap is corrupted
di ""
di "Testing heap integrity with cdecode..."
gen id_label = "group_" + string(id)
encode id_label, gen(id_encoded)
cdecode id_encoded, gen(id_decoded)

di ""
di "If you see this, the crash did NOT happen (unexpected)"
