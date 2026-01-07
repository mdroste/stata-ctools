* Check precision issue in predictions
clear all
set more off

adopath + "../build"

sysuse auto, clear

di _n "=== Precision investigation ===" _n

* First check the mpg variable storage type
describe mpg
di "mpg[1] = " %20.15f mpg[1]

* Run creghdfe
creghdfe price mpg, absorb(foreign)

* Check what type _predict creates
predict yhat_c, xb
describe yhat_c
di "yhat_c[1] = " %20.15f yhat_c[1]

* Now run reghdfe and check its predict type
reghdfe price mpg, absorb(foreign) resid
predict yhat_r, xb
describe yhat_r
di "yhat_r[1] = " %20.15f yhat_r[1]

* Generate manual predictions as double
gen double manual_double = _b[mpg] * mpg + _b[_cons]
describe manual_double

* Compare yhat_c to yhat_r
gen diff = yhat_r - yhat_c
summarize diff
di "Max diff: " %12.6e max(abs(diff))

* The issue is that mpg is stored as int
* When creghdfe posts e(b), what type is mpg treated as?
di _n "=== Direct comparison of predictions ==="
list make mpg yhat_r yhat_c diff in 1/5

di _n "Test complete."
