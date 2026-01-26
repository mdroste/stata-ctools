* Debug GMM2S + cluster - check coefficients
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
local iv_b1 = _b[x_endog]
local iv_b2 = _b[x_exog]
local iv_V11 = e(V)[1,1]
local iv_V22 = e(V)[2,2]

* civreghdfe GMM2S + cluster
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_b1 = _b[x_endog]
local c_b2 = _b[x_exog]
local c_V11 = e(V)[1,1]
local c_V22 = e(V)[2,2]

di "========== Coefficient Comparison =========="
di "b[x_endog]: ivreghdfe = " %20.15f `iv_b1'
di "            civreghdfe = " %20.15f `c_b1'
di "            diff = " %20.15e (`c_b1' - `iv_b1')
di ""
di "b[x_exog]:  ivreghdfe = " %20.15f `iv_b2'
di "            civreghdfe = " %20.15f `c_b2'
di "            diff = " %20.15e (`c_b2' - `iv_b2')

di _n "========== VCE Comparison =========="
di "V[1,1]: ivreghdfe = " %20.15f `iv_V11'
di "        civreghdfe = " %20.15f `c_V11'
di "        diff = " %20.15e (`c_V11' - `iv_V11')
di ""
di "V[2,2]: ivreghdfe = " %20.15f `iv_V22'
di "        civreghdfe = " %20.15f `c_V22'
di "        diff = " %20.15e (`c_V22' - `iv_V22')

di _n "========== Ratio Analysis =========="
di "V[1,1] ratio = " %20.15f (`c_V11' / `iv_V11')
di "V[2,2] ratio = " %20.15f (`c_V22' / `iv_V22')
di "Both ratios should be close if it's a simple scaling issue"
