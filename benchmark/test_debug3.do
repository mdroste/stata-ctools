clear all
set more off
sysuse auto, clear

quietly qreg price mpg
predict resid, residuals
gen fitted = price - resid

display "=== Data Ranges ==="
sum price, detail
local y_iqr = r(p75) - r(p25)
display "Y (price) IQR = " `y_iqr'

sum fitted, detail
local fit_iqr = r(p75) - r(p25)
display "Fitted IQR = " `fit_iqr'

sum resid, detail
local resid_iqr = r(p75) - r(p25)
display "Residual IQR = " `resid_iqr'

local N = 74
local q = 0.5
local h_prob = 0.23141524

display _n "=== Bandwidth in different scales ==="
display "h_prob = " `h_prob'

* Using residual quantiles
_pctile resid, p(`=(0.5-`h_prob')*100' `=(0.5+`h_prob')*100')
local h_resid = (r(r2) - r(r1)) / 2
display "h from residual quantiles = " `h_resid'

* Using fitted value quantiles
_pctile fitted, p(`=(0.5-`h_prob')*100' `=(0.5+`h_prob')*100')
local h_fitted = (r(r2) - r(r1)) / 2
display "h from fitted quantiles = " `h_fitted'

* Using Y quantiles
_pctile price, p(`=(0.5-`h_prob')*100' `=(0.5+`h_prob')*100')
local h_y = (r(r2) - r(r1)) / 2
display "h from Y quantiles = " `h_y'

* What sparsity would difference quotient give with residual h?
local s_dq = (r(r2) - r(r1)) / (2 * `h_prob')
* Wait, that's not right. Let me use residuals.
_pctile resid, p(`=(0.5-`h_prob')*100' `=(0.5+`h_prob')*100')
local r_lo = r(r1)
local r_hi = r(r2)
local s_dq = (`r_hi' - `r_lo') / (2 * `h_prob')
display _n "Difference quotient sparsity = " `s_dq'

* What does Stata report?
display _n "=== Stata reported values ==="
quietly qreg price mpg, vce(iid, residual)
display "qreg IID residual sparsity = " e(sparsity)

quietly qreg price mpg
display "qreg IID fitted sparsity = " e(sparsity)
