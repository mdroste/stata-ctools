clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid

local N = 69
local q = 0.5

* Get qreg implied sparsity
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)
local v11 = e(V)[1,1]
local implied_s = sqrt(`v11' / (0.25 * XtX_inv[1,1]))

display "qreg implied sparsity: " %12.4f `implied_s'

* Compute Silverman bandwidth for Gaussian kernel
summ resid
local sd = r(sd)
_pctile resid, percentile(25 75)
local iqr = r(r2) - r(r1)

local h_gauss = 1.06 * min(`sd', `iqr'/1.349) * `N'^(-0.2)
display "Gaussian KDE bandwidth: " %12.4f `h_gauss'

* Compute Gaussian KDE at 0
local sum_gauss = 0
forvalues i = 1/`N' {
    local r = resid[`i']
    local u = `r' / `h_gauss'
    local k = exp(-0.5 * `u'^2) / sqrt(2 * c(pi))
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
        local k = exp(-0.5 * `u'^2) / sqrt(2 * c(pi))
        local sum = `sum' + `k'
    }
    local f0 = `sum' / (`N' * `h')
    local sparsity = 1/`f0'
    display "mult=" %4.2f `mult' " h=" %8.1f `h' " sparsity=" %10.2f `sparsity' " ratio=" %6.4f (`sparsity'/`implied_s')
}
