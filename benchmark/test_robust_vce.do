clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display _n "=== TESTING ROBUST VCE ===" _n

* First run qreg with default VCE (iid) to get baseline
quietly qreg price mpg weight, quantile(0.5)
local se_mpg_iid = _se[mpg]
local se_wt_iid = _se[weight]
display "qreg IID VCE:"
display "  SE[mpg]    = " %12.6f `se_mpg_iid'
display "  SE[weight] = " %12.6f `se_wt_iid'

* Now run with robust VCE
quietly qreg price mpg weight, quantile(0.5) vce(robust)
local se_mpg_robust = _se[mpg]
local se_wt_robust = _se[weight]
display _n "qreg Robust VCE:"
display "  SE[mpg]    = " %12.6f `se_mpg_robust'
display "  SE[weight] = " %12.6f `se_wt_robust'

* Check ratio robust/iid
display _n "Ratio (robust/iid):"
display "  mpg:    " %8.4f (`se_mpg_robust'/`se_mpg_iid')
display "  weight: " %8.4f (`se_wt_robust'/`se_wt_iid')

* Now run cqreg with IID VCE
quietly cqreg price mpg weight, quantile(0.5) vce(iid)
local se_mpg_cqreg_iid = _se[mpg]
local se_wt_cqreg_iid = _se[weight]
display _n "cqreg IID VCE:"
display "  SE[mpg]    = " %12.6f `se_mpg_cqreg_iid'
display "  SE[weight] = " %12.6f `se_wt_cqreg_iid'

* Now run cqreg with robust VCE
quietly cqreg price mpg weight, quantile(0.5) vce(robust)
local se_mpg_cqreg_robust = _se[mpg]
local se_wt_cqreg_robust = _se[weight]
display _n "cqreg Robust VCE:"
display "  SE[mpg]    = " %12.6f `se_mpg_cqreg_robust'
display "  SE[weight] = " %12.6f `se_wt_cqreg_robust'

* Compare cqreg robust to qreg robust
display _n "Ratio cqreg/qreg (robust VCE):"
display "  mpg:    " %12.6f (`se_mpg_cqreg_robust'/`se_mpg_robust')
display "  weight: " %12.6f (`se_wt_cqreg_robust'/`se_wt_robust')

* Also check what the variance matrix looks like
display _n "=== Full ereturn list from qreg robust ==="
quietly qreg price mpg weight, quantile(0.5) vce(robust)
ereturn list
matrix list e(V), format(%12.6f)
