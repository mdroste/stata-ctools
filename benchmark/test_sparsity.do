clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

* Run qreg and get density estimate
qreg price mpg weight, quantile(0.5)
local qreg_f = e(f)
local qreg_sum_adev = e(sum_adev)
local qreg_sum_rdev = e(sum_rdev)

display "qreg f(0) = " %12.8f `qreg_f'
display "qreg sparsity (1/f) = " %12.8f 1/`qreg_f'

* Get residuals
predict double r, resid
summ r

* Manual sparsity calculation using Stata's method
local N = 69
local q = 0.5
local alpha = 0.05

* Hall-Sheather bandwidth
local z_alpha = invnormal(1 - `alpha'/2)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')

local h = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

display "Hall-Sheather h = " %12.8f `h'
display "q-h = " %12.8f (`q' - `h') " q+h = " %12.8f (`q' + `h')

* Get residual quantiles at q-h and q+h
_pctile r, percentile(`=(`q'-`h')*100' `=(`q'+`h')*100')
local r_lo = r(r1)
local r_hi = r(r2)

display "r_lo = " %12.6f `r_lo'
display "r_hi = " %12.6f `r_hi'
display "2h = " %12.8f (2*`h')

* Difference quotient sparsity
local sparsity_dq = (`r_hi' - `r_lo') / (2 * `h')
display "Sparsity (diff quotient) = " %12.4f `sparsity_dq'

* Compare to qreg's sparsity
display "qreg's sparsity = " %12.4f (1/`qreg_f')
display "Ratio (DQ/qreg) = " %12.8f (`sparsity_dq' / (1/`qreg_f'))
