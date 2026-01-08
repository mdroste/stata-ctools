clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display "=== ROBUST VCE ===" 
quietly qreg price mpg weight, quantile(0.5) vce(robust)
local se_mpg_qreg = _se[mpg]
local se_wt_qreg = _se[weight]

display "qreg vce(robust):"
display "  SE[mpg] = " %12.6f `se_mpg_qreg'
display "  SE[weight] = " %12.6f `se_wt_qreg'

quietly cqreg price mpg weight, quantile(0.5) vce(robust)
local se_mpg_cqreg = _se[mpg]
local se_wt_cqreg = _se[weight]

display _n "cqreg vce(robust):"
display "  SE[mpg] = " %12.6f `se_mpg_cqreg'
display "  SE[weight] = " %12.6f `se_wt_cqreg'

display _n "Ratios (cqreg/qreg):"
display "  mpg: " %8.4f (`se_mpg_cqreg'/`se_mpg_qreg')
display "  weight: " %8.4f (`se_wt_cqreg'/`se_wt_qreg')
