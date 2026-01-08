clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid

local N = 69
local q = 0.5

* First, let's see what Stata's kdensity gives
kdensity resid, at(0) generate(den_at_0) n(1) nograph
summ den_at_0
local kde_f0 = r(mean)
display "Stata kdensity f(0): " %12.10f `kde_f0'
display "kdensity sparsity: " %12.4f (1/`kde_f0')
drop den_at_0

* Now let's try different bandwidth methods for kdensity
* Silverman rule of thumb
summ resid
local sd = r(sd)
local iqr_raw = .
_pctile resid, percentile(25 75)
local iqr_raw = r(r2) - r(r1)
local h_silverman = 0.9 * min(`sd', `iqr_raw'/1.349) * `N'^(-0.2)
display "Silverman bandwidth: " %12.4f `h_silverman'

* Now estimate f(0) with Epanechnikov kernel manually
local sum = 0
forvalues i = 1/`N' {
    local r = resid[`i']
    local u = `r' / `h_silverman'
    if abs(`u') < 1 {
        local k = 0.75 * (1 - `u'^2)
        local sum = `sum' + `k'
    }
}
local f0_manual = `sum' / (`N' * `h_silverman')
display "Manual Epan f(0) with Silverman h: " %12.10f `f0_manual'
display "Manual Epan sparsity: " %12.4f (1/`f0_manual')

* Try with Hall-Sheather bandwidth in residual units
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h_hs = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

* Convert h_hs (probability) to residual scale using IQR
* h_resid ≈ h_prob * IQR / 0.5
_pctile resid, percentile(`=(`q'-`h_hs')*100' `=(`q'+`h_hs')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local h_resid = (`r_hi' - `r_lo')  // This is the full width

display _n "HS bandwidth (prob): " %12.8f `h_hs'
display "r_lo, r_hi: " %10.2f `r_lo' " " %10.2f `r_hi'
display "h_resid (width): " %12.4f `h_resid'

* Estimate f(0) with this bandwidth
local sum2 = 0
forvalues i = 1/`N' {
    local r = resid[`i']
    local u = `r' / `h_resid'
    if abs(`u') < 1 {
        local k = 0.75 * (1 - `u'^2)
        local sum2 = `sum2' + `k'
    }
}
local f0_hs = `sum2' / (`N' * `h_resid')
display "KDE f(0) with h_resid: " %12.10f `f0_hs'
display "KDE sparsity: " %12.4f (1/`f0_hs')

* Get qreg implied sparsity
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)
local v11 = e(V)[1,1]
local implied_s = sqrt(`v11' / (0.25 * XtX_inv[1,1]))
display _n "qreg implied sparsity: " %12.4f `implied_s'
display "implied f(0): " %12.10f (1/`implied_s')

* What bandwidth would give the qreg sparsity?
* f(0) = (1/nh) * sum K(r_i/h)
* For Epanechnikov: K(u) = 0.75*(1-u^2) for |u|<1
* At f(0), only residuals within [-h, h] contribute
* Roughly, if distribution is ~uniform near 0, f(0) ≈ 0.75 * 2h / (n*h) = 1.5/n
* Not quite right...

* Actually let's just try different bandwidths
display _n "=== TESTING DIFFERENT BANDWIDTHS ==="
foreach mult in 0.5 0.7 0.8 0.9 1.0 1.1 1.2 1.5 2.0 {
    local h_test = `h_resid' * `mult'
    local sum3 = 0
    forvalues i = 1/`N' {
        local r = resid[`i']
        local u = `r' / `h_test'
        if abs(`u') < 1 {
            local k = 0.75 * (1 - `u'^2)
            local sum3 = `sum3' + `k'
        }
    }
    local f0_test = `sum3' / (`N' * `h_test')
    local sparsity_test = 1/`f0_test'
    display "h_mult=" %4.2f `mult' " h=" %8.1f `h_test' " f(0)=" %10.8f `f0_test' " sparsity=" %10.2f `sparsity_test'
}
