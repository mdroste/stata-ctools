* Test predict functionality
* Compare reghdfe vs creghdfe predictions

clear all
set more off

adopath + "../build"

sysuse auto, clear

* Test 1: Simple case with single FE
di _n "=== Test 1: Single FE ===" _n

* Run reghdfe first
reghdfe price mpg, absorb(foreign) resid

* Store reghdfe coefficients and prediction
local b_mpg_reghdfe = _b[mpg]
local b_cons_reghdfe = _b[_cons]
predict yhat_reghdfe1, xb

di "reghdfe: b[mpg] = " %12.8f `b_mpg_reghdfe'
di "reghdfe: b[_cons] = " %12.8f `b_cons_reghdfe'

* Run creghdfe
creghdfe price mpg, absorb(foreign)

* Store creghdfe coefficients and prediction
local b_mpg_creghdfe = _b[mpg]
local b_cons_creghdfe = _b[_cons]
predict yhat_creghdfe1, xb

di "creghdfe: b[mpg] = " %12.8f `b_mpg_creghdfe'
di "creghdfe: b[_cons] = " %12.8f `b_cons_creghdfe'

* Compare coefficients
di _n "Coefficient differences:"
di "diff b[mpg] = " %12.8e (`b_mpg_reghdfe' - `b_mpg_creghdfe')
di "diff b[_cons] = " %12.8e (`b_cons_reghdfe' - `b_cons_creghdfe')

* Compare xb predictions
gen diff_yhat1 = yhat_reghdfe1 - yhat_creghdfe1

di _n "Summary of XB prediction differences (Test 1):"
summarize diff_yhat1

qui summarize diff_yhat1
local max_diff = max(abs(r(min)), abs(r(max)))
di "Max abs diff in yhat: " %12.6e `max_diff'


* Test 2: Two FEs
di _n "=== Test 2: Two FEs ===" _n

drop yhat_* diff_*

reghdfe price mpg, absorb(foreign trunk) resid
local b_mpg_reghdfe2 = _b[mpg]
local b_cons_reghdfe2 = _b[_cons]
predict yhat_reghdfe2, xb

di "reghdfe: b[mpg] = " %12.8f `b_mpg_reghdfe2'
di "reghdfe: b[_cons] = " %12.8f `b_cons_reghdfe2'

creghdfe price mpg, absorb(foreign trunk)
local b_mpg_creghdfe2 = _b[mpg]
local b_cons_creghdfe2 = _b[_cons]
predict yhat_creghdfe2, xb

di "creghdfe: b[mpg] = " %12.8f `b_mpg_creghdfe2'
di "creghdfe: b[_cons] = " %12.8f `b_cons_creghdfe2'

* Compare coefficients
di _n "Coefficient differences:"
di "diff b[mpg] = " %12.8e (`b_mpg_reghdfe2' - `b_mpg_creghdfe2')
di "diff b[_cons] = " %12.8e (`b_cons_reghdfe2' - `b_cons_creghdfe2')

gen diff_yhat2 = yhat_reghdfe2 - yhat_creghdfe2

di _n "Summary of XB prediction differences (Test 2):"
summarize diff_yhat2

qui summarize diff_yhat2
local max_diff2 = max(abs(r(min)), abs(r(max)))
di "Max abs diff in yhat: " %12.6e `max_diff2'

* Show some example predictions
di _n "Sample predictions (first 10 obs):"
list make price mpg yhat_reghdfe2 yhat_creghdfe2 diff_yhat2 in 1/10

di _n "Test complete."
