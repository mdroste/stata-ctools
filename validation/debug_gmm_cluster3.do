* Debug GMM2S + cluster DOF adjustment
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

* Get baseline values
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local iv_V11 = e(V)[1,1]

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_V11 = e(V)[1,1]

local N = 500
local K = 2
local df_a = 50
local G = 50
local df_r = `N' - `K' - `df_a'

di "========== DOF Adjustment Analysis =========="
di "N = `N', K = `K', df_a = `df_a', G = `G', df_r = `df_r'"
di ""

* Current formula: (N-1)/df_r * G/(G-1)
local current_adj = ((`N'-1)/`df_r') * (`G'/(`G'-1))
di "Current: (N-1)/df_r * G/(G-1) = " %12.10f `current_adj'

* Alternative: N/df_r * G/(G-1)
local alt_adj = (`N'/`df_r') * (`G'/(`G'-1))
di "Alt:     N/df_r * G/(G-1)     = " %12.10f `alt_adj'

di ""
di "Ratio of adjustments: " %12.10f (`alt_adj' / `current_adj')
di "Needed ratio (iv/c): " %12.10f (`iv_V11' / `c_V11')
di ""
di "If we use N/df_r * G/(G-1):"
di "  New c_V11 = c_V11 * (alt/current) = " %15.12f (`c_V11' * `alt_adj' / `current_adj')
di "  ivreghdfe V11 =                     " %15.12f `iv_V11'
di "  Diff = " %15.12e ((`c_V11' * `alt_adj' / `current_adj') - `iv_V11')
