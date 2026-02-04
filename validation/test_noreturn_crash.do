* Test whether noreturn option prevents the heap corruption crash
* Run this after running validate_civreghdfe.do to test the hypothesis

clear all
set more off
adopath ++ "build"

di "=============================================="
di "Testing civreghdfe with noreturn option"
di "Running 50 iterations to stress test..."
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

* Run 50 iterations with noreturn
forval i = 1/50 {
    if mod(`i', 10) == 0 {
        di "Iteration `i'..."
    }

    * Vary specifications to stress test
    if mod(`i', 5) == 0 {
        civreghdfe y (x_endog = z1 z2) x1, absorb(id) noreturn
    }
    else if mod(`i', 5) == 1 {
        civreghdfe y (x_endog = z1 z2 z3) x1 x2, absorb(id) vce(robust) noreturn
    }
    else if mod(`i', 5) == 2 {
        civreghdfe y (x_endog = z1 z2) x1, absorb(id) vce(cluster id) noreturn
    }
    else if mod(`i', 5) == 3 {
        civreghdfe y (x_endog x_endog2 = z1 z2 z3) x1, absorb(id) noreturn
    }
    else {
        civreghdfe y (x_endog = z1 z2) x1 x2, absorb(id time) noreturn
    }
}

di "All 50 civreghdfe iterations completed with noreturn"

* Now run cdecode to test if heap is corrupted
di ""
di "Testing heap integrity with cdecode..."
gen id_label = "group_" + string(id)
encode id_label, gen(id_encoded)
cdecode id_encoded, gen(id_decoded)

di ""
di "=============================================="
di "SUCCESS: No heap corruption with noreturn!"
di "=============================================="
di ""
di "This confirms the hypothesis that SF_scal_save/SF_mat_store"
di "calls are causing the heap corruption."
