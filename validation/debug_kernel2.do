* Debug: verify which HAC code path is used
clear all
adopath ++ "build"
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

* Check the scalar values
di _n "Checking civreghdfe scalars before call..."

* Run civreghdfe with verbose to see debug output
civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3) verbose

* Check what scalars were set
di _n "Scalars set:"
di "  __civreghdfe_hac_panel = " scalar(__civreghdfe_hac_panel)
di "  __civreghdfe_kernel = " scalar(__civreghdfe_kernel)
di "  __civreghdfe_bw = " scalar(__civreghdfe_bw)
