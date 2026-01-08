clear all
set more off
adopath + "../build"
sysuse auto, clear

display "============================================================"
display "Comprehensive comparison: qreg vs cqreg"
display "============================================================"

* Test 1: Default (should be IID fitted)
display ""
display "=== Test 1: Default settings ==="
quietly qreg price mpg
local qreg_default_se = _se[mpg]
local qreg_default_sp = e(sparsity)
local qreg_default_den = e(denmethod)

quietly cqreg price mpg
local cqreg_default_se = _se[mpg]
local cqreg_default_sp = e(sparsity)
local cqreg_default_den = e(denmethod)

display "qreg default:  SE=" %8.4f `qreg_default_se' " sparsity=" %8.2f `qreg_default_sp' " denmethod=" "`qreg_default_den'"
display "cqreg default: SE=" %8.4f `cqreg_default_se' " sparsity=" %8.2f `cqreg_default_sp' " denmethod=" "`cqreg_default_den'"
display "Ratio: " %6.4f `cqreg_default_se'/`qreg_default_se'

* Test 2: IID with residual method
display ""
display "=== Test 2: IID with residual method ==="
quietly qreg price mpg, vce(iid, residual)
local qreg_resid_se = _se[mpg]
local qreg_resid_sp = e(sparsity)

quietly cqreg price mpg, vce(iid) denmethod(residual)
local cqreg_resid_se = _se[mpg]
local cqreg_resid_sp = e(sparsity)

display "qreg iid resid:  SE=" %8.4f `qreg_resid_se' " sparsity=" %8.2f `qreg_resid_sp'
display "cqreg iid resid: SE=" %8.4f `cqreg_resid_se' " sparsity=" %8.2f `cqreg_resid_sp'
display "Ratio: " %6.4f `cqreg_resid_se'/`qreg_resid_se'

* Test 3: Robust VCE
display ""
display "=== Test 3: Robust VCE ==="
quietly qreg price mpg, vce(robust)
local qreg_robust_se = _se[mpg]
local qreg_robust_sp = e(sparsity)

quietly cqreg price mpg, vce(robust)
local cqreg_robust_se = _se[mpg]
local cqreg_robust_sp = e(sparsity)

display "qreg robust:  SE=" %8.4f `qreg_robust_se' " sparsity=" %8.2f `qreg_robust_sp'
display "cqreg robust: SE=" %8.4f `cqreg_robust_se' " sparsity=" %8.2f `cqreg_robust_sp'
display "Ratio: " %6.4f `cqreg_robust_se'/`qreg_robust_se'

* Test 4: Multiple regressors
display ""
display "=== Test 4: Multiple regressors ==="
quietly qreg price mpg weight
local qreg_mpg = _se[mpg]
local qreg_wgt = _se[weight]

quietly cqreg price mpg weight
local cqreg_mpg = _se[mpg]
local cqreg_wgt = _se[weight]

display "qreg:  SE[mpg]=" %8.4f `qreg_mpg' " SE[weight]=" %8.4f `qreg_wgt'
display "cqreg: SE[mpg]=" %8.4f `cqreg_mpg' " SE[weight]=" %8.4f `cqreg_wgt'
display "Ratios: mpg=" %6.4f `cqreg_mpg'/`qreg_mpg' " weight=" %6.4f `cqreg_wgt'/`qreg_wgt'

* Test 5: Different quantile
display ""
display "=== Test 5: 75th percentile ==="
quietly qreg price mpg, quantile(0.75)
local qreg_q75_se = _se[mpg]
local qreg_q75_sp = e(sparsity)

quietly cqreg price mpg, quantile(0.75)
local cqreg_q75_se = _se[mpg]
local cqreg_q75_sp = e(sparsity)

display "qreg Q75:  SE=" %8.4f `qreg_q75_se' " sparsity=" %8.2f `qreg_q75_sp'
display "cqreg Q75: SE=" %8.4f `cqreg_q75_se' " sparsity=" %8.2f `cqreg_q75_sp'
display "Ratio: " %6.4f `cqreg_q75_se'/`qreg_q75_se'

display ""
display "============================================================"
display "Summary: All ratios should be 1.0000 for exact match"
display "============================================================"
