* Test predict functionality - xbd option
* Compare reghdfe vs creghdfe predictions including fixed effects

clear all
set more off

adopath + "../build"

sysuse auto, clear

di _n "=== Test: XBD (fitted values including fixed effects) ===" _n

* Run reghdfe with resid option to enable xbd
reghdfe price mpg, absorb(foreign trunk) resid
predict yhat_reghdfe_xbd, xbd
predict yhat_reghdfe_xb, xb
predict yhat_reghdfe_d, d
predict yhat_reghdfe_resid, residuals

di "reghdfe predictions created:"
summarize yhat_reghdfe_*

* Now run creghdfe
creghdfe price mpg, absorb(foreign trunk)

* Try predict with xbd option - likely will fail
capture predict yhat_creghdfe_xbd, xbd
if _rc {
    di as error "xbd option not available for creghdfe (rc=" _rc ")"
    di as text "This is expected - creghdfe does not store fixed effects"
}
else {
    di "creghdfe xbd prediction created"
}

* xb should work
predict yhat_creghdfe_xb, xb

* Compare the xb predictions
gen diff_xb = yhat_reghdfe_xb - yhat_creghdfe_xb

di _n "XB prediction comparison:"
summarize diff_xb
qui sum diff_xb
local max_diff = max(abs(r(min)), abs(r(max)))
di "Max abs diff in xb: " %12.6e `max_diff'

* Show what reghdfe's xbd gives vs xb
di _n "reghdfe xbd vs xb comparison:"
gen diff_xbd_xb = yhat_reghdfe_xbd - yhat_reghdfe_xb
summarize diff_xbd_xb
di "This difference represents the absorbed fixed effects"

* Check: xbd should equal y - residuals
gen check = price - yhat_reghdfe_resid - yhat_reghdfe_xbd
di _n "Verification: y - resid - xbd should be 0:"
summarize check

di _n "Test complete."
