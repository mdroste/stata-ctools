* Detailed predict test to find source of differences
clear all
set more off

adopath + "../build"

sysuse auto, clear

di _n "=== Detailed Investigation of Predict Differences ===" _n

* Run reghdfe
reghdfe price mpg, absorb(foreign) resid

* Store everything
local b_mpg_r = _b[mpg]
local b_cons_r = _b[_cons]
di "reghdfe: b[mpg] = " %20.15f `b_mpg_r'
di "reghdfe: b[_cons] = " %20.15f `b_cons_r'

predict yhat_r, xb

* Manually compute xb using reghdfe coefficients
gen yhat_r_manual = `b_mpg_r' * mpg + `b_cons_r'

* Check if reghdfe's predict matches manual computation
gen diff_r = yhat_r - yhat_r_manual
summarize diff_r
di "reghdfe: Max diff between predict and manual: " %12.6e max(abs(diff_r))

* Now run creghdfe
creghdfe price mpg, absorb(foreign)

local b_mpg_c = _b[mpg]
local b_cons_c = _b[_cons]
di _n "creghdfe: b[mpg] = " %20.15f `b_mpg_c'
di "creghdfe: b[_cons] = " %20.15f `b_cons_c'

predict yhat_c, xb

* Manually compute xb using creghdfe coefficients
gen yhat_c_manual = `b_mpg_c' * mpg + `b_cons_c'

* Check if creghdfe's predict matches manual computation
gen diff_c = yhat_c - yhat_c_manual
summarize diff_c
di "creghdfe: Max diff between predict and manual: " %12.6e max(abs(diff_c))

* Compare the manual predictions using each set of coefficients
gen diff_manual = yhat_r_manual - yhat_c_manual
summarize diff_manual
di "Diff in manual predictions (coef difference): " %12.6e max(abs(diff_manual))

* Compare the actual predictions
gen diff_pred = yhat_r - yhat_c
summarize diff_pred
di "Diff in actual predictions: " %12.6e max(abs(diff_pred))

* Key question: are the coefficient differences responsible?
di _n "=== Coefficient Comparison ==="
di "b[mpg] diff = " %20.15e (`b_mpg_r' - `b_mpg_c')
di "b[_cons] diff = " %20.15e (`b_cons_r' - `b_cons_c')

* What is the actual mpg range?
summarize mpg

* Effect of mpg coef diff on predictions at different mpg values
di _n "=== Effect of coefficient differences ==="
di "At mpg=12: pred diff from mpg coef = " %12.6e (12 * (`b_mpg_r' - `b_mpg_c'))
di "At mpg=41: pred diff from mpg coef = " %12.6e (41 * (`b_mpg_r' - `b_mpg_c'))
di "Const diff contribution = " %12.6e (`b_cons_r' - `b_cons_c')

* Show predictions for specific observations
di _n "=== Sample Predictions ==="
list make mpg yhat_r yhat_c diff_pred in 1/5

di _n "Test complete."
