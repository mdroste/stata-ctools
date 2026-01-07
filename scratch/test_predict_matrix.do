* Check matrix storage issue
clear all
set more off

adopath + "../build"

sysuse auto, clear

di _n "=== Check e(b) matrix after creghdfe ===" _n

creghdfe price mpg, absorb(foreign)

di "Immediately after creghdfe:"
di "_b[mpg] = " %20.15f _b[mpg]
di "_b[_cons] = " %20.15f _b[_cons]

matrix list e(b), format(%20.15f)

* Check what _predict sees
predict yhat_c, xb

* Manual calculation
local mpg1 = mpg[1]
local manual = _b[mpg] * `mpg1' + _b[_cons]

di _n "For obs 1 (mpg=" `mpg1' "):"
di "  Manual calc: " %20.15f `manual'
di "  predict xb:  " %20.15f yhat_c[1]
di "  Difference:  " %12.6e (yhat_c[1] - `manual')

* Try manually computing for all obs
gen manual_pred = _b[mpg] * mpg + _b[_cons]
gen diff = yhat_c - manual_pred

summarize diff
di _n "Max abs difference: " %12.6e max(abs(diff), abs(diff))

* Check e(V) matrix
di _n "e(V) matrix:"
matrix list e(V), format(%12.6f)

* Check variable names in e(b)
local names : colfullnames e(b)
di _n "Column names in e(b): `names'"

* Check if there's an issue with how _predict works
* Create our own xb manually using matrix operations
matrix b = e(b)
di _n "Direct matrix access:"
di "  b[1,1] (mpg coef) = " %20.15f b[1,1]
di "  b[1,2] (_cons) = " %20.15f b[1,2]

gen manual2 = b[1,1] * mpg + b[1,2]
gen diff2 = yhat_c - manual2
summarize diff2

di _n "Test complete."
