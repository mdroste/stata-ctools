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

* Check normal-scale estimator: sparsity = σ * sqrt(2π)
summ resid
local sd = r(sd)
local sparsity_normal = `sd' * sqrt(2 * c(pi))
display _n "=== NORMAL-SCALE ESTIMATOR ==="
display "SD of residuals: " %12.4f `sd'
display "Normal-scale sparsity: " %12.4f `sparsity_normal'
display "Ratio to qreg: " %8.4f (`sparsity_normal' / `implied_s')

* Try IQR-based scale (more robust)
_pctile resid, percentile(25 75)
local iqr = r(r2) - r(r1)
local sigma_iqr = `iqr' / 1.349  // 1.349 = 2*Φ^{-1}(0.75) for normal
local sparsity_iqr = `sigma_iqr' * sqrt(2 * c(pi))
display _n "=== IQR-BASED SCALE ESTIMATOR ==="
display "IQR: " %12.4f `iqr'
display "sigma_iqr: " %12.4f `sigma_iqr'
display "IQR sparsity: " %12.4f `sparsity_iqr'
display "Ratio to qreg: " %8.4f (`sparsity_iqr' / `implied_s')

* Try MAD-based scale
_pctile resid, percentile(50)
local med = r(r1)
gen abs_dev = abs(resid - `med')
summ abs_dev
_pctile abs_dev, percentile(50)
local mad = r(r1)
local sigma_mad = `mad' / 0.6745  // 0.6745 = Φ^{-1}(0.75) for normal
local sparsity_mad = `sigma_mad' * sqrt(2 * c(pi))
display _n "=== MAD-BASED SCALE ESTIMATOR ==="
display "MAD: " %12.4f `mad'
display "sigma_mad: " %12.4f `sigma_mad'
display "MAD sparsity: " %12.4f `sparsity_mad'
display "Ratio to qreg: " %8.4f (`sparsity_mad' / `implied_s')

* Try DQ method but with different bandwidth formulas
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')

* Hall-Sheather bandwidth
local h_hs = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

* Bofinger bandwidth  
local h_bof = ((9/2) * `phi_zq'^4 * (2*`z_q'^2 + 1)^2 / `N')^(1/5)

display _n "=== DIFFERENCE QUOTIENT SPARSITY ==="
display "h_hs = " %12.8f `h_hs'
display "h_bof = " %12.8f `h_bof'

foreach bw in hs bof {
    local h = `h_`bw''
    local q_lo = `q' - `h'
    local q_hi = `q' + `h'
    _pctile resid, percentile(`=`q_lo'*100' `=`q_hi'*100')
    local r_lo = r(r1)
    local r_hi = r(r2)
    local sparsity = (`r_hi' - `r_lo') / (2 * `h')
    display "`bw' DQ sparsity: " %12.4f `sparsity' " ratio=" %8.4f (`sparsity'/`implied_s')
}

* Key insight: what if Stata uses a fitted quantile function rather than residuals?
* The conditional density f_Y|X(ξ(τ)|X) at the τ-quantile
* For IID VCE, this is constant across X, estimated from residuals

* Check what bandwidth would make DQ match qreg
* sparsity = (r_hi - r_lo) / (2h)
* So h = (r_hi - r_lo) / (2 * sparsity)
_pctile resid, percentile(`=(`q'-`h_hs')*100' `=(`q'+`h_hs')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local h_implied = (`r_hi' - `r_lo') / (2 * `implied_s')
display _n "=== FINDING MATCHING BANDWIDTH ==="
display "Width (r_hi - r_lo): " %12.4f (`r_hi' - `r_lo')
display "h that would match qreg: " %12.8f `h_implied'
display "h_hs: " %12.8f `h_hs'
display "Ratio h_implied/h_hs: " %8.4f (`h_implied' / `h_hs')
