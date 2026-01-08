clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid
predict double yhat, xb

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

* Hall-Sheather bandwidth
local z_alpha = invnormal(0.975)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')
local h = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

display _n "Hall-Sheather h: " %12.8f `h'

* Try DQ on fitted values instead of residuals
display _n "=== DQ ON FITTED VALUES ==="
_pctile yhat, percentile(`=(`q'-`h')*100' `=(`q'+`h')*100')
local y_lo = r(r1)
local y_hi = r(r2)
local sparsity_yhat = (`y_hi' - `y_lo') / (2 * `h')
display "yhat DQ sparsity: " %12.4f `sparsity_yhat' " ratio=" %8.4f (`sparsity_yhat'/`implied_s')

* Try DQ on Y (depvar)
display _n "=== DQ ON Y (price) ==="
_pctile price, percentile(`=(`q'-`h')*100' `=(`q'+`h')*100')
local p_lo = r(r1)
local p_hi = r(r2)
local sparsity_y = (`p_hi' - `p_lo') / (2 * `h')
display "Y DQ sparsity: " %12.4f `sparsity_y' " ratio=" %8.4f (`sparsity_y'/`implied_s')

* Try normal-scale with df correction
display _n "=== NORMAL-SCALE WITH DF CORRECTION ==="
summ resid
local sd = r(sd)  // This uses n-1 divisor

* If variance is sum(r²)/(n-1), then for qreg we might need n-k divisor
* sd_nk = sd * sqrt((n-1)/(n-k))
local sd_nk = `sd' * sqrt(68/66)
local sparsity_nk = `sd_nk' * sqrt(2 * c(pi))
display "SD (n-1 divisor): " %12.4f `sd'
display "SD (n-k divisor): " %12.4f `sd_nk'
display "Sparsity (n-k): " %12.4f `sparsity_nk' " ratio=" %8.4f (`sparsity_nk'/`implied_s')

* Or maybe the original sum squares with just n divisor?
* Then convert: sd_n = sd * sqrt((n-1)/n)
local sd_n = `sd' * sqrt(68/69)
local sparsity_n = `sd_n' * sqrt(2 * c(pi))
display "SD (n divisor): " %12.4f `sd_n'
display "Sparsity (n): " %12.4f `sparsity_n' " ratio=" %8.4f (`sparsity_n'/`implied_s')

* Check if Stata is using bootstrapped SE or some other adjustment
* The ratio qreg/normal_scale = 6947/6464 = 1.0747
* What mathematical factor gives this?
display _n "=== INVESTIGATION ==="
display "qreg/normal = " %8.4f (`implied_s' / (`sd' * sqrt(2*c(pi))))
display "This factor of ~1.075 might come from:"
display "  sqrt(N/(N-K)) = " %8.4f sqrt(69/66)
display "  N/(N-K) = " %8.4f (69/66)
display "  (N-1)/(N-K) = " %8.4f (68/66)

* Actually, let me try: sparsity = sd * sqrt(2π) * some_factor
* Where some_factor = 1.0747
* factor² = 1.155
* Maybe related to kurtosis or other moment?

* Let me check if Stata uses Bofinger bandwidth with normal-scale
local h_bof = ((9/2) * `phi_zq'^4 * (2*`z_q'^2 + 1)^2 / `N')^(1/5)
local sparsity_bof_norm = `sd' * sqrt(2 * c(pi)) * (`h_bof' / `h')
display _n "Bofinger/HS ratio: " %8.4f (`h_bof' / `h')
display "Normal-scale * Bof/HS: " %12.4f `sparsity_bof_norm' " ratio=" %8.4f (`sparsity_bof_norm'/`implied_s')
