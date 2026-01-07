* Check e(b) matrix structure
clear all
set more off

adopath + "../build"

sysuse auto, clear

di _n "=== Check e(b) matrix structure ===" _n

* Run reghdfe
reghdfe price mpg, absorb(foreign) resid

di "reghdfe e(b) matrix:"
matrix list e(b), format(%20.15f)

di _n "reghdfe coefficients from _b[]:"
di "_b[mpg] = " %20.15f _b[mpg]
di "_b[_cons] = " %20.15f _b[_cons]

predict yhat_r, xb

* Now run creghdfe
creghdfe price mpg, absorb(foreign)

di _n "creghdfe e(b) matrix:"
matrix list e(b), format(%20.15f)

di _n "creghdfe coefficients from _b[]:"
di "_b[mpg] = " %20.15f _b[mpg]
di "_b[_cons] = " %20.15f _b[_cons]

predict yhat_c, xb

* Compare predictions
gen diff = yhat_r - yhat_c
qui sum diff
di _n "Prediction difference:"
di "  Mean: " %12.6e r(mean)
di "  Min:  " %12.6e r(min)
di "  Max:  " %12.6e r(max)

* Manual check
di _n "=== Manual calculation check ==="
* For first observation (mpg = 22)
local mpg1 = mpg[1]
di "First obs mpg = " `mpg1'

* reghdfe prediction for first obs
di "reghdfe yhat[1] = " %20.15f yhat_r[1]
* creghdfe prediction for first obs
di "creghdfe yhat[1] = " %20.15f yhat_c[1]

* Manual using creghdfe coefficients
local manual_c = _b[mpg] * `mpg1' + _b[_cons]
di "Manual (creghdfe coefs): " %20.15f `manual_c'

* Run reghdfe again to get its coefficients for manual calc
reghdfe price mpg, absorb(foreign) resid
local manual_r = _b[mpg] * `mpg1' + _b[_cons]
di "Manual (reghdfe coefs): " %20.15f `manual_r'

di _n "yhat_r[1] - manual_r = " %12.6e (yhat_r[1] - `manual_r')

di _n "Test complete."
