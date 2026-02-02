* Test HAC: compare panel vs no-panel structure
clear all
adopath ++ "build"
set seed 54321

* Create panel data: 2 panels of 100 obs each
set obs 200
gen id = ceil(_n / 100)
gen time = mod(_n - 1, 100) + 1
tsset id time

gen z1 = runiform()
gen z2 = rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

di _n "=========================================="
di "TEST A: Panel data (2 panels x 100 time periods)"
di "=========================================="

* Run ivreghdfe
qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
matrix V_ivreg = e(V)
scalar widstat_ivreg = e(widstat)
di "ivreghdfe: V[1,1]=" V_ivreg[1,1] " V[2,2]=" V_ivreg[2,2] " widstat=" widstat_ivreg

* Run civreghdfe
qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
matrix V_civreg = e(V)
scalar widstat_civreg = e(widstat)
di "civreghdfe: V[1,1]=" V_civreg[1,1] " V[2,2]=" V_civreg[2,2] " widstat=" widstat_civreg
di "Ratio: V[1,1]=" (V_civreg[1,1] / V_ivreg[1,1]) " V[2,2]=" (V_civreg[2,2] / V_ivreg[2,2])

di _n "=========================================="
di "TEST B: Same data, 1 panel x 200 time periods (restructured)"
di "=========================================="

* Restructure: single panel
replace id = 1
replace time = _n
tsset id time

* Run ivreghdfe
qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
matrix V_ivreg = e(V)
scalar widstat_ivreg = e(widstat)
di "ivreghdfe: V[1,1]=" V_ivreg[1,1] " V[2,2]=" V_ivreg[2,2] " widstat=" widstat_ivreg

* Run civreghdfe
qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
matrix V_civreg = e(V)
scalar widstat_civreg = e(widstat)
di "civreghdfe: V[1,1]=" V_civreg[1,1] " V[2,2]=" V_civreg[2,2] " widstat=" widstat_civreg
di "Ratio: V[1,1]=" (V_civreg[1,1] / V_ivreg[1,1]) " V[2,2]=" (V_civreg[2,2] / V_ivreg[2,2])
