clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

* Try different VCE options in qreg to understand the formula
display "=== QREG WITH DIFFERENT VCE OPTIONS ==="

qreg price mpg weight, quantile(0.5)
matrix V_iid = e(V)
display "IID SE[mpg]: " %12.6f sqrt(V_iid[1,1])

* Try with vce(robust)
qreg price mpg weight, quantile(0.5) vce(robust)
matrix V_robust = e(V)
display "Robust SE[mpg]: " %12.6f sqrt(V_robust[1,1])

* Try with bootstrap
qreg price mpg weight, quantile(0.5) vce(boot, reps(100))
matrix V_boot = e(V)
display "Bootstrap SE[mpg]: " %12.6f sqrt(V_boot[1,1])

* Compute IQR-based sparsity estimate (Stata's traditional method)
* Stata docs say: "sparsity is estimated using the method of Koenker (2005)"
* This typically involves fitted density with Gaussian kernel

* Let's check fitted values and residuals
qreg price mpg weight, quantile(0.5)
predict double xb, xb
predict double resid, resid

* Sort by fitted values
sort xb

* Koenker's method: estimate density of Y conditional on X
* Using local polynomial or kernel methods

* For now, let's verify the formula: V = s² * τ(1-τ) * (X'X)^{-1}
* And check if there's a scale factor

gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)

local N = 69
local K = 3
local q = 0.5
local tau_factor = `q' * (1 - `q')

display _n "=== FORMULA VERIFICATION ==="

* Get qreg's VCE
local v11 = V_iid[1,1]
local xtxinv11 = XtX_inv[1,1]

* If V = s² * τ(1-τ) * (X'X)^{-1}, then s² = V / (τ(1-τ) * (X'X)^{-1})
local implied_s2 = `v11' / (`tau_factor' * `xtxinv11')
local implied_s = sqrt(`implied_s2')

display "V[1,1] = " %12.4f `v11'
display "(X'X)^{-1}[1,1] = " %12.8f `xtxinv11'
display "Implied s = " %12.4f `implied_s'

* Now let's compute sparsity using fitted Gaussian kernel density
* Stata may fit a local density at the median of residuals

* Get Hall-Sheather bandwidth
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h_prob = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

* Stata's fitted density uses Gaussian kernel
* f(0) = (1/nh) * sum_i K((r_i - 0)/h)
* where K is Gaussian: (1/sqrt(2π)) * exp(-u²/2)

summ resid
local sd = r(sd)
_pctile resid, percentile(25 75)
local iqr = r(r2) - r(r1)

* Standard bandwidth for Gaussian KDE
local h_gauss = 1.06 * min(`sd', `iqr'/1.349) * `N'^(-0.2)
display _n "Gaussian KDE bandwidth: " %12.4f `h_gauss'

* Compute Gaussian KDE at 0
local sum_gauss = 0
forvalues i = 1/`N' {
    local r = resid[`i']
    local u = `r' / `h_gauss'
    local k = exp(-0.5 * `u'^2) / sqrt(2 * _pi)
    local sum_gauss = `sum_gauss' + `k'
}
local f0_gauss = `sum_gauss' / (`N' * `h_gauss')
local sparsity_gauss = 1/`f0_gauss'
display "Gaussian f(0) = " %12.10f `f0_gauss'
display "Gaussian sparsity = " %12.4f `sparsity_gauss'
display "Ratio to qreg: " %8.4f (`sparsity_gauss' / `implied_s')

* Try with different bandwidth scales
display _n "=== GAUSSIAN KDE WITH SCALED BANDWIDTH ==="
foreach mult in 0.5 0.75 1.0 1.25 1.5 2.0 2.5 3.0 {
    local h = `h_gauss' * `mult'
    local sum = 0
    forvalues i = 1/`N' {
        local r = resid[`i']
        local u = `r' / `h'
        local k = exp(-0.5 * `u'^2) / sqrt(2 * _pi)
        local sum = `sum' + `k'
    }
    local f0 = `sum' / (`N' * `h')
    local sparsity = 1/`f0'
    display "mult=" %4.2f `mult' " h=" %8.1f `h' " sparsity=" %10.2f `sparsity' " ratio=" %6.4f (`sparsity'/`implied_s')
}
