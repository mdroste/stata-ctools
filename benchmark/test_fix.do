clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display "============================================="
display "TESTING CQREG VCE FIX"
display "============================================="

* Run qreg with residual method (what we're matching)
qreg price mpg weight, quantile(0.5) vce(iid, residual)
matrix V_qreg_resid = e(V)
local se_mpg_qreg_resid = _se[mpg]
local se_wt_qreg_resid = _se[weight]
local sparsity_qreg = e(sparsity)
local bw_qreg = e(bwidth)

display _n "=== QREG VCE(IID, RESIDUAL) ==="
display "SE[mpg]: " %12.6f `se_mpg_qreg_resid'
display "SE[weight]: " %12.6f `se_wt_qreg_resid'
display "Sparsity: " %12.4f `sparsity_qreg'
display "Bandwidth: " %12.8f `bw_qreg'

* Run qreg with default (fitted) method
quietly qreg price mpg weight, quantile(0.5)
local se_mpg_qreg_fitted = _se[mpg]
local se_wt_qreg_fitted = _se[weight]

display _n "=== QREG VCE(IID, FITTED) [default] ==="
display "SE[mpg]: " %12.6f `se_mpg_qreg_fitted'
display "SE[weight]: " %12.6f `se_wt_qreg_fitted'

* Run cqreg
cqreg price mpg weight, quantile(0.5) verbose
local se_mpg_cqreg = _se[mpg]
local se_wt_cqreg = _se[weight]

display _n "=== CQREG ==="
display "SE[mpg]: " %12.6f `se_mpg_cqreg'
display "SE[weight]: " %12.6f `se_wt_cqreg'

display _n "============================================="
display "COMPARISON (cqreg vs qreg residual method)"
display "============================================="
display "SE[mpg] ratio: " %8.6f (`se_mpg_cqreg' / `se_mpg_qreg_resid')
display "SE[weight] ratio: " %8.6f (`se_wt_cqreg' / `se_wt_qreg_resid')

local diff_mpg = abs(`se_mpg_cqreg' - `se_mpg_qreg_resid')
local diff_wt = abs(`se_wt_cqreg' - `se_wt_qreg_resid')

display _n "Absolute differences:"
display "  SE[mpg]: " %12.8f `diff_mpg'
display "  SE[weight]: " %12.10f `diff_wt'

if `diff_mpg' < 0.001 & `diff_wt' < 0.00001 {
    display _n "{result}SUCCESS: cqreg matches qreg vce(iid, residual)!{txt}"
} 
else {
    display _n "{error}MISMATCH: cqreg differs from qreg{txt}"
}

display _n "============================================="
display "NOTE: cqreg matches qreg vce(iid, residual)."
display "For qreg's default (fitted) method, ratio is:"
display "  cqreg/fitted = " %8.4f (`se_mpg_cqreg' / `se_mpg_qreg_fitted')
display "============================================="
