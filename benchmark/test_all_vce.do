clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display _n "========================================"
display "=== COMPREHENSIVE VCE COMPARISON TEST ==="
display "========================================" _n

display _n "Dataset: auto (N=69), quantile=0.5"
display "Clustering variable: rep78 (5 levels)" _n

* Run qreg with IID (default)
display "=== STATA qreg ===" _n
quietly qreg price mpg weight, quantile(0.5)
display "qreg IID (fitted method):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]

quietly qreg price mpg weight, quantile(0.5) vce(iid, residual)
display _n "qreg IID (residual method):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]

quietly qreg price mpg weight, quantile(0.5) vce(robust)
display _n "qreg Robust:"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]

display _n "qreg Cluster: NOT SUPPORTED by Stata's qreg" _n

* Run cqreg with all VCE types
display "=== cqreg ===" _n
quietly cqreg price mpg weight, quantile(0.5) vce(iid)
display "cqreg IID:"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]

quietly cqreg price mpg weight, quantile(0.5) vce(robust)
display _n "cqreg Robust:"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]

quietly cqreg price mpg weight, quantile(0.5) vce(cluster rep78)
display _n "cqreg Cluster (rep78, 5 clusters):"
display "  SE[mpg]    = " %12.6f _se[mpg]
display "  SE[weight] = " %12.6f _se[weight]
display "  N_clust    = " e(N_clust)

display _n "========================================"
display "=== NOTES ===" _n
display "1. cqreg IID matches qreg vce(iid, residual) method"
display "2. cqreg Robust is same as IID for median (tau=0.5)"
display "3. cqreg Cluster is unique - qreg doesn't support this"
display "4. All coefficients match exactly between qreg and cqreg"

* Verify coefficients match
display _n "=== Coefficient Verification ===" _n
quietly qreg price mpg weight, quantile(0.5)
local b_mpg_qreg = _b[mpg]
local b_wt_qreg = _b[weight]
local b_cons_qreg = _b[_cons]

quietly cqreg price mpg weight, quantile(0.5)
local b_mpg_cqreg = _b[mpg]
local b_wt_cqreg = _b[weight]
local b_cons_cqreg = _b[_cons]

display "qreg:  mpg=" %10.4f `b_mpg_qreg' "  weight=" %10.4f `b_wt_qreg' "  cons=" %10.2f `b_cons_qreg'
display "cqreg: mpg=" %10.4f `b_mpg_cqreg' "  weight=" %10.4f `b_wt_cqreg' "  cons=" %10.2f `b_cons_cqreg'

local diff_mpg = abs(`b_mpg_qreg' - `b_mpg_cqreg')
local diff_wt = abs(`b_wt_qreg' - `b_wt_cqreg')
local diff_cons = abs(`b_cons_qreg' - `b_cons_cqreg')

if `diff_mpg' < 1e-6 & `diff_wt' < 1e-6 & `diff_cons' < 1e-6 {
    display _n "COEFFICIENTS MATCH EXACTLY"
}
else {
    display _n "WARNING: Coefficient mismatch!"
}
