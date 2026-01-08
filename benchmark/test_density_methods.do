* Test density estimation methods for cqreg VCE
* Compare cqreg with different denmethod options to qreg

clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display _n "========================================"
display "=== DENSITY METHOD COMPARISON TEST ==="
display "========================================" _n

display _n "Dataset: auto (N=" _N "), quantile=0.5" _n

* ========================================
* IID VCE Comparison
* ========================================
display _n "=== IID VCE COMPARISON ===" _n

* qreg default (fitted method)
quietly qreg price mpg weight, quantile(0.5)
display "qreg vce(iid) [fitted method - default]:"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_qreg_fitted = _se[mpg]
local se_wt_qreg_fitted = _se[weight]

* qreg residual method
quietly qreg price mpg weight, quantile(0.5) vce(iid, residual)
display _n "qreg vce(iid, residual):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_qreg_resid = _se[mpg]
local se_wt_qreg_resid = _se[weight]

* cqreg fitted method (default)
quietly cqreg price mpg weight, quantile(0.5) vce(iid) denmethod(fitted)
display _n "cqreg vce(iid) denmethod(fitted):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_cqreg_fitted = _se[mpg]
local se_wt_cqreg_fitted = _se[weight]

* cqreg residual method
quietly cqreg price mpg weight, quantile(0.5) vce(iid) denmethod(residual)
display _n "cqreg vce(iid) denmethod(residual):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_cqreg_resid = _se[mpg]
local se_wt_cqreg_resid = _se[weight]

* Check matches
display _n "=== IID VCE Match Check ===" _n
local ratio_fitted = `se_mpg_cqreg_fitted' / `se_mpg_qreg_fitted'
local ratio_resid = `se_mpg_cqreg_resid' / `se_mpg_qreg_resid'
display "Fitted method ratio (cqreg/qreg): " %6.4f `ratio_fitted'
display "Residual method ratio (cqreg/qreg): " %6.4f `ratio_resid'

if abs(`ratio_resid' - 1) < 0.001 {
    display _n "RESIDUAL METHOD MATCHES qreg vce(iid, residual)"
}
else {
    display _n "WARNING: Residual method mismatch"
}

* ========================================
* Robust VCE Comparison
* ========================================
display _n _n "=== ROBUST VCE COMPARISON ===" _n

* qreg robust
quietly qreg price mpg weight, quantile(0.5) vce(robust)
display "qreg vce(robust):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_qreg_robust = _se[mpg]
local se_wt_qreg_robust = _se[weight]

* cqreg robust with fitted method (should match qreg robust)
quietly cqreg price mpg weight, quantile(0.5) vce(robust) denmethod(fitted)
display _n "cqreg vce(robust) denmethod(fitted):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_cqreg_robust_fitted = _se[mpg]
local se_wt_cqreg_robust_fitted = _se[weight]

* cqreg robust with residual method
quietly cqreg price mpg weight, quantile(0.5) vce(robust) denmethod(residual)
display _n "cqreg vce(robust) denmethod(residual):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
local se_mpg_cqreg_robust_resid = _se[mpg]
local se_wt_cqreg_robust_resid = _se[weight]

* Check matches
display _n "=== Robust VCE Match Check ===" _n
local ratio_robust = `se_mpg_cqreg_robust_fitted' / `se_mpg_qreg_robust'
display "Robust fitted method ratio (cqreg/qreg): " %6.4f `ratio_robust'

if abs(`ratio_robust' - 1) < 0.10 {
    display _n "ROBUST VCE MATCHES qreg vce(robust) (within 10%)"
}
else {
    display _n "WARNING: Robust VCE mismatch - ratio: " %6.4f `ratio_robust'
}

* ========================================
* Cluster VCE (unique to cqreg)
* ========================================
display _n _n "=== CLUSTER VCE (cqreg only) ===" _n

quietly cqreg price mpg weight, quantile(0.5) vce(cluster rep78)
display "cqreg vce(cluster rep78):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
display "  N_clust    = " e(N_clust)

display _n "Note: Stata's qreg does not support cluster VCE" _n

* ========================================
* Summary
* ========================================
display _n "========================================"
display "=== SUMMARY ==="
display "========================================" _n

display "IID methods:"
display "  qreg fitted:     SE[mpg]=" %8.2f `se_mpg_qreg_fitted'
display "  qreg residual:   SE[mpg]=" %8.2f `se_mpg_qreg_resid'
display "  cqreg fitted:    SE[mpg]=" %8.2f `se_mpg_cqreg_fitted'
display "  cqreg residual:  SE[mpg]=" %8.2f `se_mpg_cqreg_resid'

display _n "Robust methods:"
display "  qreg robust:     SE[mpg]=" %8.2f `se_mpg_qreg_robust'
display "  cqreg fitted:    SE[mpg]=" %8.2f `se_mpg_cqreg_robust_fitted'
display "  cqreg residual:  SE[mpg]=" %8.2f `se_mpg_cqreg_robust_resid'
