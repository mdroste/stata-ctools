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
local implied_f0 = 1/`implied_s'

display "qreg implied sparsity: " %12.4f `implied_s'
display "qreg implied f(0): " %12.10f `implied_f0'

* Hall-Sheather bandwidth
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h_prob = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

* Get residual quantiles
_pctile resid, percentile(`=(`q'-`h_prob')*100' `=(`q'+`h_prob')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local h_resid = (`r_hi' - `r_lo')

display _n "h_prob = " %12.8f `h_prob'
display "h_resid = " %12.4f `h_resid'

* Compute difference quotient sparsity
local sparsity_dq = (`r_hi' - `r_lo') / (2 * `h_prob')
display "DQ sparsity = " %12.4f `sparsity_dq'

* Standard deviation of residuals
summ resid
local sd = r(sd)
display _n "SD of residuals: " %12.4f `sd'

* Silverman optimal bandwidth
_pctile resid, percentile(25 75)
local iqr = r(r2) - r(r1)
local h_silverman = 0.9 * min(`sd', `iqr'/1.349) * `N'^(-0.2)
display "Silverman bandwidth: " %12.4f `h_silverman'

* Use Silverman bandwidth for KDE
* Epanechnikov kernel: K(u) = 0.75*(1-u^2) for |u|<1
local sum_kde = 0
forvalues i = 1/`N' {
    local r = resid[`i']
    local u = `r' / `h_silverman'
    if abs(`u') < 1 {
        local k = 0.75 * (1 - `u'^2)
        local sum_kde = `sum_kde' + `k'
    }
}
local f0_silverman = `sum_kde' / (`N' * `h_silverman')
display _n "=== KDE with Silverman bandwidth ==="
display "f(0) = " %12.10f `f0_silverman'
display "sparsity = " %12.4f (1/`f0_silverman')

* Now try with different bandwidths to find what matches qreg
display _n "=== SEARCHING FOR MATCHING BANDWIDTH ==="
display "Target f(0) = " %12.10f `implied_f0'

foreach h in 500 1000 1500 2000 2500 3000 3500 4000 5000 {
    local sum_kde = 0
    forvalues i = 1/`N' {
        local r = resid[`i']
        local u = `r' / `h'
        if abs(`u') < 1 {
            local k = 0.75 * (1 - `u'^2)
            local sum_kde = `sum_kde' + `k'
        }
    }
    local f0 = `sum_kde' / (`N' * `h')
    local sparsity = 1/`f0'
    display "h=" %5.0f `h' " f(0)=" %12.10f `f0' " sparsity=" %10.2f `sparsity' " ratio=" %8.4f (`sparsity'/`implied_s')
}

* Try the h_resid value
display _n "=== WITH h_resid ==="
local h = `h_resid'
local sum_kde = 0
forvalues i = 1/`N' {
    local r = resid[`i']
    local u = `r' / `h'
    if abs(`u') < 1 {
        local k = 0.75 * (1 - `u'^2)
        local sum_kde = `sum_kde' + `k'
    }
}
local f0 = `sum_kde' / (`N' * `h')
display "h_resid=" %8.2f `h_resid' " f(0)=" %12.10f `f0' " sparsity=" %10.2f (1/`f0')
