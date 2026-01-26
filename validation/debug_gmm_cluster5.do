* Debug GMM2S + cluster - check potential small factors
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

* Get values
qui ivreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local iv_V11 = e(V)[1,1]
local iv_df_r = e(df_r)

qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_V11 = e(V)[1,1]
local c_df_r = e(df_r)

local N = 500
local K = 2
local df_a = 50
local G = 50

di "========== Current Comparison =========="
di "ivreghdfe V[1,1] = " %20.15f `iv_V11'
di "civreghdfe V[1,1] = " %20.15f `c_V11'
di "Diff = " %20.15e (`c_V11' - `iv_V11')
di "Ratio c/iv = " %20.15f (`c_V11' / `iv_V11')
di "Needed multiplier = " %20.15f (`iv_V11' / `c_V11')

di _n "========== e(df_r) values =========="
di "ivreghdfe e(df_r) = `iv_df_r'"
di "civreghdfe e(df_r) = `c_df_r'"

di _n "========== Potential small corrections =========="
local needed = `iv_V11' / `c_V11'
di "Needed multiplier = " %15.12f `needed'
di "N/(N-1) = " %15.12f (500/499)
di "(N+1)/N = " %15.12f (501/500)
di "(G+1)/G = " %15.12f (51/50)
di "df_r/(df_r-1) = " %15.12f (448/447)

* Check if nested adjustment matters
di _n "With df_r = G - 1 = 49 (ivreghdfe reports):"
local iv_adj = (`N' / 49.0) * (`G' / (`G' - 1.0))
di "DOF adj using df_r=49: " %15.12f `iv_adj'

local our_df_r = `N' - `K' - `df_a'
local our_adj = (`N' / `our_df_r') * (`G' / (`G' - 1.0))
di "DOF adj using df_r=448: " %15.12f `our_adj'

di _n "Ratio of adjustments: " %15.12f (`iv_adj' / `our_adj')
