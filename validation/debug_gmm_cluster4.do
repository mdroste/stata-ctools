* Debug GMM2S + cluster - fresh comparison after fix
adopath ++ "build/"
clear all
set seed 12345
set obs 500
gen firm = ceil(_n / 10)
gen z1 = runiform()
gen z2 = rnormal()
gen z3 = runiform() + 0.1*rnormal()
gen x_endog = 0.5*z1 + 0.3*z2 + 0.2*z3 + rnormal()
gen x_exog = runiform()
gen y = 2*x_endog + 1.5*x_exog + rnormal()

* ivreghdfe GMM2S + cluster
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local iv_V11 = e(V)[1,1]

* civreghdfe GMM2S + cluster
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_V11 = e(V)[1,1]

di "========== After Fix =========="
di "ivreghdfe V[1,1] = " %20.15f `iv_V11'
di "civreghdfe V[1,1] = " %20.15f `c_V11'
di "Diff = " %20.15e (`c_V11' - `iv_V11')
di "Ratio = " %20.15f (`c_V11' / `iv_V11')
