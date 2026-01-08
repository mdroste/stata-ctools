clear all
set more off
adopath + "../build"
sysuse auto, clear

* qreg results
quietly qreg price mpg, vce(robust)
local qreg_robust = _se[mpg]

quietly qreg price mpg, vce(iid, residual)
local qreg_iid_resid = _se[mpg]
local qreg_sparsity_resid = e(sparsity)

quietly qreg price mpg
local qreg_iid_fitted = _se[mpg]
local qreg_sparsity_fitted = e(sparsity)

* cqreg results
quietly cqreg price mpg, vce(robust) denmethod(residual)
local cqreg_robust_resid = _se[mpg]

quietly cqreg price mpg, vce(iid) denmethod(residual)
local cqreg_iid_resid = _se[mpg]
local cqreg_sparsity_resid = e(sparsity)

display "=== Summary ==="
display "qreg robust:        SE[mpg] = " %8.4f `qreg_robust'
display "qreg IID fitted:    SE[mpg] = " %8.4f `qreg_iid_fitted' "  (sparsity=" %6.1f `qreg_sparsity_fitted' ")"
display "qreg IID residual:  SE[mpg] = " %8.4f `qreg_iid_resid' "  (sparsity=" %6.1f `qreg_sparsity_resid' ")"
display ""
display "cqreg robust (residual): SE[mpg] = " %8.4f `cqreg_robust_resid'
display "cqreg IID (residual):    SE[mpg] = " %8.4f `cqreg_iid_resid' "  (sparsity=" %6.1f `cqreg_sparsity_resid' ")"

display _n "Ratios:"
display "qreg robust / qreg IID residual = " %6.4f `qreg_robust'/`qreg_iid_resid'
display "cqreg robust / cqreg IID = " %6.4f `cqreg_robust_resid'/`cqreg_iid_resid'
