clear all
set more off
adopath + "../build"
sysuse auto, clear

display "============================================================"
display "Comparing qreg and cqreg IID VCE methods"
display "============================================================"

* qreg IID fitted (default)
quietly qreg price mpg
local qreg_iid_fitted_se = _se[mpg]
local qreg_iid_fitted_sp = e(sparsity)
local qreg_iid_fitted_bw = e(bwidth)

* qreg IID residual
quietly qreg price mpg, vce(iid, residual)
local qreg_iid_resid_se = _se[mpg]
local qreg_iid_resid_sp = e(sparsity)
local qreg_iid_resid_bw = e(bwidth)

* cqreg IID fitted (default)
quietly cqreg price mpg, vce(iid)
local cqreg_iid_fitted_se = _se[mpg]
local cqreg_iid_fitted_sp = e(sparsity)
local cqreg_iid_fitted_bw = e(bwidth)

* cqreg IID residual
quietly cqreg price mpg, vce(iid) denmethod(residual)
local cqreg_iid_resid_se = _se[mpg]
local cqreg_iid_resid_sp = e(sparsity)
local cqreg_iid_resid_bw = e(bwidth)

display ""
display "=== qreg results ==="
display "IID fitted (default):  SE[mpg] = " %8.4f `qreg_iid_fitted_se' "  sparsity = " %8.2f `qreg_iid_fitted_sp' "  bwidth = " %8.6f `qreg_iid_fitted_bw'
display "IID residual:          SE[mpg] = " %8.4f `qreg_iid_resid_se' "  sparsity = " %8.2f `qreg_iid_resid_sp' "  bwidth = " %8.6f `qreg_iid_resid_bw'

display ""
display "=== cqreg results ==="
display "IID fitted (default):  SE[mpg] = " %8.4f `cqreg_iid_fitted_se' "  sparsity = " %8.2f `cqreg_iid_fitted_sp' "  bwidth = " %8.6f `cqreg_iid_fitted_bw'
display "IID residual:          SE[mpg] = " %8.4f `cqreg_iid_resid_se' "  sparsity = " %8.2f `cqreg_iid_resid_sp' "  bwidth = " %8.6f `cqreg_iid_resid_bw'

display ""
display "=== Comparison ==="
display "IID fitted:  cqreg/qreg SE ratio = " %6.4f `cqreg_iid_fitted_se'/`qreg_iid_fitted_se'
display "IID fitted:  cqreg/qreg sparsity ratio = " %6.4f `cqreg_iid_fitted_sp'/`qreg_iid_fitted_sp'
display "IID residual: cqreg/qreg SE ratio = " %6.4f `cqreg_iid_resid_se'/`qreg_iid_resid_se'
display "IID residual: cqreg/qreg sparsity ratio = " %6.4f `cqreg_iid_resid_sp'/`qreg_iid_resid_sp'

* Target: ratios should be 1.0000 for exact match
