clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

* Run qreg with default (fitted) method
qreg price mpg weight, quantile(0.5)
local se_mpg_fitted = _se[mpg]
local se_wt_fitted = _se[weight]
display "=== FITTED METHOD (default) ==="
display "SE[mpg]: " %12.6f `se_mpg_fitted'
display "SE[weight]: " %12.6f `se_wt_fitted'

* Run qreg with residual method
qreg price mpg weight, quantile(0.5) vce(iid, residual)
local se_mpg_resid = _se[mpg]
local se_wt_resid = _se[weight]
display _n "=== RESIDUAL METHOD ==="
display "SE[mpg]: " %12.6f `se_mpg_resid'
display "SE[weight]: " %12.6f `se_wt_resid'
display "sparsity: " %12.4f e(sparsity)
display "bandwidth: " %12.8f e(bwidth)

* Run cqreg
cqreg price mpg weight, quantile(0.5)
local se_mpg_cqreg = _se[mpg]
local se_wt_cqreg = _se[weight]
display _n "=== CQREG ==="
display "SE[mpg]: " %12.6f `se_mpg_cqreg'
display "SE[weight]: " %12.6f `se_wt_cqreg'

display _n "=== RATIOS (cqreg / qreg) ==="
display "vs fitted method: " %8.4f (`se_mpg_cqreg' / `se_mpg_fitted')
display "vs residual method: " %8.4f (`se_mpg_cqreg' / `se_mpg_resid')
