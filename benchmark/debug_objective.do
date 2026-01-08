* Debug: compute objective values at different solutions
clear all
set more off
adopath + "../build"

sysuse auto, clear

local q = 0.5  // median

* Compute objective function rho_q(r) = sum_i |r_i| * (q - I(r_i < 0))
* For median (q=0.5), this is 0.5 * sum |r_i|

* Step 1: OLS solution
quietly regress price mpg
matrix ols_b = e(b)
predict double r_ols, resid
gen double obj_ols_i = abs(r_ols) * (`q' - (r_ols < 0))
quietly summarize obj_ols_i, meanonly
local obj_ols = r(sum)
di "OLS coefficients: mpg=" %12.6f ols_b[1,1] " _cons=" %12.6f ols_b[1,2]
di "OLS objective: " %14.4f `obj_ols'
drop r_ols obj_ols_i

* Step 2: qreg solution
quietly qreg price mpg
matrix qreg_b = e(b)
local qreg_adev = e(sum_adev)
predict double r_qreg, resid
gen double obj_qreg_i = abs(r_qreg) * (`q' - (r_qreg < 0))
quietly summarize obj_qreg_i, meanonly
local obj_qreg = r(sum)
di _newline "QREG coefficients: mpg=" %12.6f qreg_b[1,1] " _cons=" %12.6f qreg_b[1,2]
di "QREG objective (computed): " %14.4f `obj_qreg'
di "QREG e(sum_adev):          " %14.4f `qreg_adev'
drop r_qreg obj_qreg_i

* Step 3: cqreg solution
quietly cqreg price mpg
matrix cqreg_b = e(b)
local cqreg_adev = e(sum_adev)
predict double r_cqreg, resid
gen double obj_cqreg_i = abs(r_cqreg) * (`q' - (r_cqreg < 0))
quietly summarize obj_cqreg_i, meanonly
local obj_cqreg = r(sum)
di _newline "CQREG coefficients: mpg=" %12.6f cqreg_b[1,1] " _cons=" %12.6f cqreg_b[1,2]
di "CQREG objective (computed): " %14.4f `obj_cqreg'
di "CQREG e(sum_adev):          " %14.4f `cqreg_adev'
drop r_cqreg obj_cqreg_i

di _newline "{hline 60}"
di "SUMMARY"
di "{hline 60}"
di "Objective at OLS:   " %14.4f `obj_ols'
di "Objective at QREG:  " %14.4f `obj_qreg' " (should be lowest)"
di "Objective at CQREG: " %14.4f `obj_cqreg' " (should match qreg)"
di "{hline 60}"
di "Improvement OLS -> QREG: " %10.4f `obj_ols' - `obj_qreg'
di "Improvement OLS -> CQREG: " %10.4f `obj_ols' - `obj_cqreg'
di "Gap CQREG - QREG (should be ~0): " %10.4f `obj_cqreg' - `obj_qreg'

* Compute objective manually at qreg coefficients
di _newline "{hline 60}"
di "VERIFY: Manual objective at qreg coefficients"
di "{hline 60}"
gen double yhat_qreg = qreg_b[1,1] * mpg + qreg_b[1,2]
gen double r_manual = price - yhat_qreg
gen double obj_manual_i = abs(r_manual) * (`q' - (r_manual < 0))
quietly summarize obj_manual_i, meanonly
di "Manual objective at qreg beta: " %14.4f r(sum)
drop yhat_qreg r_manual obj_manual_i

* Compute objective manually at cqreg coefficients
di _newline "VERIFY: Manual objective at cqreg coefficients"
gen double yhat_cqreg = cqreg_b[1,1] * mpg + cqreg_b[1,2]
gen double r_manual = price - yhat_cqreg
gen double obj_manual_i = abs(r_manual) * (`q' - (r_manual < 0))
quietly summarize obj_manual_i, meanonly
di "Manual objective at cqreg beta: " %14.4f r(sum)
drop yhat_cqreg r_manual obj_manual_i
