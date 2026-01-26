* Debug GMM2S + cluster VCE difference
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
local iv_V22 = e(V)[2,2]
local iv_N = e(N)
local iv_df_r = e(df_r)
local iv_N_clust = e(N_clust)

* civreghdfe GMM2S + cluster
qui civreghdfe y (x_endog = z1 z2 z3) x_exog, absorb(firm) gmm2s vce(cluster firm)
local c_V11 = e(V)[1,1]
local c_V22 = e(V)[2,2]
local c_N = e(N)
local c_df_r = e(df_r)
local c_N_clust = e(N_clust)

di "========== VCE Comparison =========="
di "V[1,1]: ivreghdfe = " %20.15f `iv_V11' ", civreghdfe = " %20.15f `c_V11'
di "        diff = " %20.15e (`c_V11' - `iv_V11')
di "        ratio = " %20.15f (`c_V11' / `iv_V11')
di ""
di "V[2,2]: ivreghdfe = " %20.15f `iv_V22' ", civreghdfe = " %20.15f `c_V22'
di "        diff = " %20.15e (`c_V22' - `iv_V22')
di "        ratio = " %20.15f (`c_V22' / `iv_V22')

di _n "========== Metadata =========="
di "N:       ivreghdfe = `iv_N', civreghdfe = `c_N'"
di "df_r:    ivreghdfe = `iv_df_r', civreghdfe = `c_df_r'"
di "N_clust: ivreghdfe = `iv_N_clust', civreghdfe = `c_N_clust'"

di _n "========== DOF Adjustment Analysis =========="
local G = `iv_N_clust'
local N = `iv_N'
local K = 2
local df_a = 50

* Current formula: (N-1)/df_r * G/(G-1) where df_r = N - K - df_a
local df_r = `N' - `K' - `df_a'
local current_adj = ((`N'-1)/`df_r') * (`G'/(`G'-1))
di "Current DOF adj: (N-1)/df_r * G/(G-1) = " %20.15f `current_adj'

* What if df_a is not included in df_r?
local df_r_no_fe = `N' - `K'
local alt_adj1 = ((`N'-1)/`df_r_no_fe') * (`G'/(`G'-1))
di "Alt adj (no FE in df_r): (N-1)/(N-K) * G/(G-1) = " %20.15f `alt_adj1'

* What adjustment would make them equal?
local needed_ratio = `iv_V11' / `c_V11'
di "Needed ratio (iv/c) = " %20.15f `needed_ratio'

* Check if ivreghdfe uses nested adjustment (FE nested in cluster)
di _n "If FE nested in cluster, df_a might be treated differently"
local df_r_nested = `N' - `K'  // FE absorbed but not counted
local nested_adj = ((`N'-1)/`df_r_nested') * (`G'/(`G'-1))
di "Nested adj: (N-1)/(N-K) * G/(G-1) = " %20.15f `nested_adj'
