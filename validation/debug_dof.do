* Debug DOF differences for dkraay
clear all
adopath + "build"

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

* Run with verbose output
di _n "=== ivreghdfe dkraay(2) ==="
ivreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2)

* Show key e() scalars
di _n "Key scalars from ivreghdfe:"
di "  e(N) = " e(N)
di "  e(N_clust) = " e(N_clust)
di "  e(df_r) = " e(df_r)
di "  e(df_m) = " e(df_m)
di "  e(widstat) = " e(widstat)
di "  e(idstat) = " e(idstat)
capture di "  e(sdofminus) = " e(sdofminus)
capture di "  e(dofminus) = " e(dofminus)
capture di "  e(partial_ct) = " e(partial_ct)
capture di "  e(partialcons) = " e(partialcons)

* The widstat formula is: rkf = chi2/(N-1) * dof * (G-1)/G / L
* Solve for chi2: chi2 = widstat * (N-1) * L * G / (dof * (G-1))
* If dof = N - K - sdofminus, we can back out chi2

* Try different dof values to see which gives consistent chi2
local widstat = e(widstat)
local N = e(N)
local G = e(N_clust)
local L = 2

foreach assumed_dof in 496 446 447 448 449 450 {
    local chi2_implied = `widstat' * (`N'-1) * `L' * `G' / (`assumed_dof' * (`G'-1))
    di "If dof = `assumed_dof': implied chi2 = " %10.4f `chi2_implied'
}

di _n "=== civreghdfe dkraay(2) ==="
civreghdfe y (x_endog = z1 z2) x_exog, absorb(id) dkraay(2) verbose
di _n "civreghdfe widstat = " e(widstat)
di "civreghdfe idstat = " e(idstat)
