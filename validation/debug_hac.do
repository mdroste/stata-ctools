* Debug HAC kernel differences between civreghdfe and ivreghdfe
* Test different T values to check if ratio scales with T/(T-1)
clear all
adopath + "build"

foreach T in 5 10 20 {
    di _n _n "{hline 60}"
    di "TESTING WITH T = `T' TIME PERIODS"
    di "{hline 60}"

    clear
    set seed 54321
    local N_panels = 50
    local N = `N_panels' * `T'
    set obs `N'
    gen id = ceil(_n / `T')
    gen time = mod(_n - 1, `T') + 1
    tsset id time

    gen z1 = runiform()
    gen z2 = rnormal()
    gen x_endog = 0.5*z1 + 0.3*z2 + rnormal()
    gen x_exog = runiform()
    gen y = 2*x_endog + 1.5*x_exog + rnormal()

    * Run dkraay(2) with both
    qui ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2)
    local iv_widstat = e(widstat)
    local iv_N_clust = e(N_clust)
    local iv_idstat = e(idstat)

    qui civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2)
    local c_widstat = e(widstat)
    local c_N_clust = e(N_clust)
    local c_idstat = e(idstat)

    local ratio_w = `iv_widstat' / `c_widstat'
    local ratio_id = `iv_idstat' / `c_idstat'
    local expected = (`T' - 1) / `T'

    di "  T = `T', N_clust = `iv_N_clust' (iv) / `c_N_clust' (civ)"
    di "  widstat: iv=" %10.4f `iv_widstat' " civ=" %10.4f `c_widstat' " ratio=" %8.6f `ratio_w'
    di "  idstat:  iv=" %10.4f `iv_idstat' " civ=" %10.4f `c_idstat' " ratio=" %8.6f `ratio_id'
    di "  (T-1)/T = " %8.6f `expected'
    di "  widstat matches (T-1)/T? " cond(abs(`ratio_w' - `expected') < 0.01, "YES", "NO")
    di "  idstat matches?  " cond(abs(`ratio_id' - 1.0) < 0.01, "YES", "NO")
}

* Test dkraay(2)
di _n "{hline 60}"
di "TEST: dkraay(2)"
di "{hline 60}"

di _n "=== ivreghdfe ==="
ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2)
local iv_b = _b[x_endog]
local iv_idstat = e(idstat)
local iv_idp = e(idp)
local iv_widstat = e(widstat)
local iv_V11 = e(V)[1,1]
local iv_N = e(N)
local iv_N_clust = e(N_clust)

di _n "Results:"
di "  N = `iv_N', N_clust = `iv_N_clust'"
di "  b[x_endog] = `iv_b'"
di "  idstat = `iv_idstat'"
di "  idp = `iv_idp'"
di "  widstat = `iv_widstat'"
di "  V[1,1] = `iv_V11'"

di _n "=== civreghdfe ==="
civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2)
local c_b = _b[x_endog]
local c_idstat = e(idstat)
local c_idp = e(idp)
local c_widstat = e(widstat)
local c_V11 = e(V)[1,1]
local c_N = e(N)
local c_N_clust = e(N_clust)

di _n "Results:"
di "  N = `c_N', N_clust = `c_N_clust'"
di "  b[x_endog] = `c_b'"
di "  idstat = `c_idstat'"
di "  idp = `c_idp'"
di "  widstat = `c_widstat'"
di "  V[1,1] = `c_V11'"

di _n "=== Comparison ==="
di "  b ratio (iv/c) = " %8.6f `iv_b'/`c_b'
di "  idstat ratio = " %8.6f `iv_idstat'/`c_idstat'
di "  widstat ratio = " %8.6f `iv_widstat'/`c_widstat'
di "  V[1,1] ratio = " %8.6f `iv_V11'/`c_V11'

* Test parzen kernel
di _n _n "{hline 60}"
di "TEST: kernel(parzen) bw(3)"
di "{hline 60}"

di _n "=== ivreghdfe ==="
ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
local iv_b = _b[x_endog]
local iv_idstat = e(idstat)
local iv_widstat = e(widstat)
local iv_V11 = e(V)[1,1]

di _n "Results:"
di "  b[x_endog] = `iv_b'"
di "  idstat = `iv_idstat'"
di "  widstat = `iv_widstat'"
di "  V[1,1] = `iv_V11'"

di _n "=== civreghdfe ==="
civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kernel(parzen) bw(3)
local c_b = _b[x_endog]
local c_idstat = e(idstat)
local c_widstat = e(widstat)
local c_V11 = e(V)[1,1]

di _n "Results:"
di "  b[x_endog] = `c_b'"
di "  idstat = `c_idstat'"
di "  widstat = `c_widstat'"
di "  V[1,1] = `c_V11'"

di _n "=== Comparison ==="
di "  b ratio (iv/c) = " %8.6f `iv_b'/`c_b'
di "  idstat ratio = " %8.6f `iv_idstat'/`c_idstat'
di "  widstat ratio = " %8.6f `iv_widstat'/`c_widstat'
di "  V[1,1] ratio = " %8.6f `iv_V11'/`c_V11'

* Test kiefer
di _n _n "{hline 60}"
di "TEST: kiefer"
di "{hline 60}"

di _n "=== ivreghdfe ==="
ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer
local iv_b = _b[x_endog]
local iv_idstat = e(idstat)
local iv_widstat = e(widstat)
local iv_sargan = e(sargan)
local iv_V11 = e(V)[1,1]

di _n "Results:"
di "  b[x_endog] = `iv_b'"
di "  idstat = `iv_idstat'"
di "  widstat = `iv_widstat'"
di "  sargan = `iv_sargan'"
di "  V[1,1] = `iv_V11'"

di _n "=== civreghdfe ==="
civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) kiefer
local c_b = _b[x_endog]
local c_idstat = e(idstat)
local c_widstat = e(widstat)
local c_sargan = e(sargan)
local c_V11 = e(V)[1,1]

di _n "Results:"
di "  b[x_endog] = `c_b'"
di "  idstat = `c_idstat'"
di "  widstat = `c_widstat'"
di "  sargan = `c_sargan'"
di "  V[1,1] = `c_V11'"

di _n "=== Comparison ==="
di "  b ratio (iv/c) = " %8.6f `iv_b'/`c_b'
di "  idstat ratio = " %8.6f `iv_idstat'/`c_idstat'
di "  widstat ratio = " %8.6f `iv_widstat'/`c_widstat'
di "  sargan ratio = " %8.6f `iv_sargan'/`c_sargan'
di "  V[1,1] ratio = " %8.6f `iv_V11'/`c_V11'
