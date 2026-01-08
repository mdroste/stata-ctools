clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

* Run qreg and examine what it stores
qreg price mpg weight, quantile(0.5)
ereturn list

* Get the VCE matrix
matrix list e(V)

* Compute X'X inverse manually
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)
matrix list XtX_inv

* IID VCE formula: V = s^2 * tau*(1-tau) * (X'X)^{-1}
* So: s^2 = V[1,1] / (tau*(1-tau) * (X'X)^{-1}[1,1])
local v11 = e(V)[1,1]
local xtxinv11 = XtX_inv[1,1]
local tau = 0.5
local scale = `tau' * (1 - `tau')
local implied_s2 = `v11' / (`scale' * `xtxinv11')
local implied_s = sqrt(`implied_s2')

display _n "====================================="
display "BACK-CALCULATED SPARSITY FROM QREG"
display "====================================="
display "V[1,1] = " %12.6f `v11'
display "(X'X)^{-1}[1,1] = " %12.8f `xtxinv11'
display "tau*(1-tau) = " %6.4f `scale'
display "implied s^2 = " %12.2f `implied_s2'
display "implied s = " %12.4f `implied_s'
display _n "This is qreg's effective sparsity."

* Now compute difference quotient sparsity
predict double resid, resid
local N = 69
local alpha = 0.05
local z_alpha = invnormal(1 - `alpha'/2)
local z_q = invnormal(`tau')
local phi_zq = normalden(`z_q')

local h = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

_pctile resid, percentile(`=(`tau'-`h')*100' `=(`tau'+`h')*100')
local r_lo = r(r1)
local r_hi = r(r2)

local sparsity_dq = (`r_hi' - `r_lo') / (2 * `h')

display _n "====================================="
display "DIFFERENCE QUOTIENT SPARSITY"
display "====================================="
display "h = " %12.8f `h'
display "r_lo = " %12.4f `r_lo'
display "r_hi = " %12.4f `r_hi'
display "sparsity_dq = " %12.4f `sparsity_dq'

display _n "====================================="
display "COMPARISON"
display "====================================="
display "qreg implied sparsity: " %12.4f `implied_s'
display "DQ sparsity:           " %12.4f `sparsity_dq'
display "Ratio (qreg/DQ):       " %12.6f (`implied_s' / `sparsity_dq')

* What multiplier would make DQ match qreg?
local multiplier = `implied_s' / `sparsity_dq'
display _n "To match qreg, DQ sparsity needs multiplier: " %8.6f `multiplier'
