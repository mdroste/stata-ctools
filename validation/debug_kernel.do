* Diagnostic script for kernel HAC VCE issues
clear all

* Load ctools from build directory
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

di _n "=========================================="
di "TEST 1: Parzen kernel bw(3)"
di "=========================================="

* Run ivreghdfe
qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
matrix V_ivreg = e(V)
scalar widstat_ivreg = e(widstat)
scalar idstat_ivreg = e(idstat)
scalar idp_ivreg = e(idp)

di "ivreghdfe results:"
di "  V[1,1] = " V_ivreg[1,1]
di "  V[2,2] = " V_ivreg[2,2]
di "  widstat = " widstat_ivreg
di "  idstat = " idstat_ivreg
di "  idp = " idp_ivreg

* Run civreghdfe
qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
matrix V_civreg = e(V)
scalar widstat_civreg = e(widstat)
scalar idstat_civreg = e(idstat)
scalar idp_civreg = e(idp)

di _n "civreghdfe results:"
di "  V[1,1] = " V_civreg[1,1]
di "  V[2,2] = " V_civreg[2,2]
di "  widstat = " widstat_civreg
di "  idstat = " idstat_civreg
di "  idp = " idp_civreg

di _n "Comparison (civreghdfe - ivreghdfe):"
di "  V[1,1] diff = " (V_civreg[1,1] - V_ivreg[1,1])
di "  V[2,2] diff = " (V_civreg[2,2] - V_ivreg[2,2])
di "  widstat diff = " (widstat_civreg - widstat_ivreg)
di "  idstat diff = " (idstat_civreg - idstat_ivreg)
di "  V[1,1] ratio = " (V_civreg[1,1] / V_ivreg[1,1])
di "  V[2,2] ratio = " (V_civreg[2,2] / V_ivreg[2,2])

di _n "=========================================="
di "TEST 2: Quadratic Spectral kernel bw(2)"
di "=========================================="

* Run ivreghdfe
qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(qs) bw(2)
matrix V_ivreg = e(V)
scalar widstat_ivreg = e(widstat)
scalar idstat_ivreg = e(idstat)
scalar idp_ivreg = e(idp)
scalar F_ivreg = e(F)

di "ivreghdfe results:"
di "  V[1,1] = " V_ivreg[1,1]
di "  V[2,2] = " V_ivreg[2,2]
di "  widstat = " widstat_ivreg
di "  idstat = " idstat_ivreg
di "  idp = " idp_ivreg
di "  F = " F_ivreg

* Run civreghdfe
qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(qs) bw(2)
matrix V_civreg = e(V)
scalar widstat_civreg = e(widstat)
scalar idstat_civreg = e(idstat)
scalar idp_civreg = e(idp)
scalar F_civreg = e(F)

di _n "civreghdfe results:"
di "  V[1,1] = " V_civreg[1,1]
di "  V[2,2] = " V_civreg[2,2]
di "  widstat = " widstat_civreg
di "  idstat = " idstat_civreg
di "  idp = " idp_civreg
di "  F = " F_civreg

di _n "Comparison (civreghdfe - ivreghdfe):"
di "  V[1,1] diff = " (V_civreg[1,1] - V_ivreg[1,1])
di "  V[2,2] diff = " (V_civreg[2,2] - V_ivreg[2,2])
di "  widstat diff = " (widstat_civreg - widstat_ivreg)
di "  idstat diff = " (idstat_civreg - idstat_ivreg)
di "  F diff = " (F_civreg - F_ivreg)
di "  V[1,1] ratio = " (V_civreg[1,1] / V_ivreg[1,1])
di "  V[2,2] ratio = " (V_civreg[2,2] / V_ivreg[2,2])

di _n "=========================================="
di "TEST 3: Robust only (no kernel) - baseline comparison"
di "=========================================="

* Run ivreghdfe
qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) vce(robust)
matrix V_ivreg = e(V)
scalar widstat_ivreg = e(widstat)

di "ivreghdfe results:"
di "  V[1,1] = " V_ivreg[1,1]
di "  V[2,2] = " V_ivreg[2,2]
di "  widstat = " widstat_ivreg

* Run civreghdfe
qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) vce(robust)
matrix V_civreg = e(V)
scalar widstat_civreg = e(widstat)

di _n "civreghdfe results:"
di "  V[1,1] = " V_civreg[1,1]
di "  V[2,2] = " V_civreg[2,2]
di "  widstat = " widstat_civreg

di _n "Comparison (civreghdfe - ivreghdfe):"
di "  V[1,1] ratio = " (V_civreg[1,1] / V_ivreg[1,1])
di "  V[2,2] ratio = " (V_civreg[2,2] / V_ivreg[2,2])

di _n "=========================================="
di "TEST 4: Bartlett bw(2) - for comparison"
di "=========================================="

* Run ivreghdfe
qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(bartlett) bw(2)
matrix V_ivreg = e(V)
scalar widstat_ivreg = e(widstat)

di "ivreghdfe results:"
di "  V[1,1] = " V_ivreg[1,1]
di "  V[2,2] = " V_ivreg[2,2]
di "  widstat = " widstat_ivreg

* Run civreghdfe
qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(bartlett) bw(2)
matrix V_civreg = e(V)
scalar widstat_civreg = e(widstat)

di _n "civreghdfe results:"
di "  V[1,1] = " V_civreg[1,1]
di "  V[2,2] = " V_civreg[2,2]
di "  widstat = " widstat_civreg

di _n "Comparison (civreghdfe - ivreghdfe):"
di "  V[1,1] ratio = " (V_civreg[1,1] / V_ivreg[1,1])
di "  V[2,2] ratio = " (V_civreg[2,2] / V_ivreg[2,2])
