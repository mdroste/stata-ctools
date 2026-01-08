clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display "=== CLUSTER VCE ===" 
quietly qreg price mpg weight, quantile(0.5) vce(cluster rep78)
local se_mpg_qreg = _se[mpg]

display "qreg vce(cluster rep78):"
display "  SE[mpg] = " %12.6f `se_mpg_qreg'

quietly cqreg price mpg weight, quantile(0.5) vce(cluster rep78)
local se_mpg_cqreg = _se[mpg]

display _n "cqreg vce(cluster rep78):"
display "  SE[mpg] = " %12.6f `se_mpg_cqreg'

display _n "Ratio: " %8.4f (`se_mpg_cqreg'/`se_mpg_qreg')
