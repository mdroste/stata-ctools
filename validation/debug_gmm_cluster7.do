* Debug GMM2S + cluster - check nested adjustment
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

* ivreghdfe with nested=0 (force non-nested)
di "========== Testing nested adjustment =========="

* ivreghdfe default (nested detection)
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local iv_V11 = e(V)[1,1]
local iv_nested = e(df_a_nested)
di "ivreghdfe (default):"
di "  V[1,1] = " %20.15f `iv_V11'
di "  df_a_nested = `iv_nested'"

* civreghdfe
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_V11 = e(V)[1,1]
di "civreghdfe:"
di "  V[1,1] = " %20.15f `c_V11'

di _n "Diff = " %20.15e (`c_V11' - `iv_V11')

* Check what multiplier we need
local needed = `iv_V11' / `c_V11'
di "Needed multiplier = " %20.15f `needed'

* Check if multiplying by (N-K)/(N-K-1) helps
local N = 500
local K = 2
local factor1 = (`N'-`K') / (`N'-`K'-1)
di _n "Potential factors:"
di "(N-K)/(N-K-1) = " %15.12f `factor1'
di "Multiplied VCE would be: " %20.15f (`c_V11' * `factor1')

* Check factor that would exactly match
di _n "To match exactly, factor = " %15.12f `needed'
di "This corresponds to 1 + " %15.12e (`needed' - 1)
