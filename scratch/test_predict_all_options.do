* Test all predict options match between reghdfe and creghdfe
* Options: xb, stdp, residuals, d, xbd, dresiduals

clear all
set more off

adopath + "../build"

sysuse auto, clear

di _n "=== Testing All Predict Options ===" _n

* Run reghdfe with resid option (required for most predict options)
reghdfe price mpg weight, absorb(foreign trunk) resid
predict yhat_r_xb, xb
predict yhat_r_stdp, stdp
predict yhat_r_resid, residuals
predict yhat_r_d, d
predict yhat_r_xbd, xbd
predict yhat_r_dresid, dresiduals

di "reghdfe predictions created"

* Run creghdfe with resid option
creghdfe price mpg weight, absorb(foreign trunk) resid
predict yhat_c_xb, xb
predict yhat_c_stdp, stdp
predict yhat_c_resid, residuals
predict yhat_c_d, d
predict yhat_c_xbd, xbd
predict yhat_c_dresid, dresiduals

di "creghdfe predictions created"

* Compare each option
di _n "{hline 60}"
di "Prediction Option Comparison (reghdfe vs creghdfe)"
di "{hline 60}"

* XB
gen diff_xb = yhat_r_xb - yhat_c_xb
qui sum diff_xb
local max_xb = max(abs(r(min)), abs(r(max)))
di "XB:          Max abs diff = " %12.6e `max_xb' _col(50) cond(`max_xb' < 1e-8, "PASS", "FAIL")

* STDP
gen diff_stdp = yhat_r_stdp - yhat_c_stdp
qui sum diff_stdp
local max_stdp = max(abs(r(min)), abs(r(max)))
di "STDP:        Max abs diff = " %12.6e `max_stdp' _col(50) cond(`max_stdp' < 1e-8, "PASS", "FAIL")

* Residuals
gen diff_resid = yhat_r_resid - yhat_c_resid
qui sum diff_resid
local max_resid = max(abs(r(min)), abs(r(max)))
di "Residuals:   Max abs diff = " %12.6e `max_resid' _col(50) cond(`max_resid' < 1e-8, "PASS", "FAIL")

* D (fixed effects)
gen diff_d = yhat_r_d - yhat_c_d
qui sum diff_d
local max_d = max(abs(r(min)), abs(r(max)))
di "D:           Max abs diff = " %12.6e `max_d' _col(50) cond(`max_d' < 1e-8, "PASS", "FAIL")

* XBD (xb + d)
gen diff_xbd = yhat_r_xbd - yhat_c_xbd
qui sum diff_xbd
local max_xbd = max(abs(r(min)), abs(r(max)))
di "XBD:         Max abs diff = " %12.6e `max_xbd' _col(50) cond(`max_xbd' < 1e-8, "PASS", "FAIL")

* DResiduals (d + residuals)
gen diff_dresid = yhat_r_dresid - yhat_c_dresid
qui sum diff_dresid
local max_dresid = max(abs(r(min)), abs(r(max)))
di "DResiduals:  Max abs diff = " %12.6e `max_dresid' _col(50) cond(`max_dresid' < 1e-8, "PASS", "FAIL")

di "{hline 60}"

* Show summary stats for each prediction type
di _n "Summary of predictions:"
summarize yhat_r_* yhat_c_*

di _n "Test complete."
