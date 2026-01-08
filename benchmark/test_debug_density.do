* Debug density estimation
clear all
set more off
sysuse auto, clear

quietly qreg price mpg, quantile(0.5) vce(robust)
display "qreg robust:"
display "  bwidth = " e(bwidth)
display "  N = " e(N)

* Get the residuals
predict resid, residuals
summarize resid, detail

* Check how many observations are near 0
count if abs(resid) < 500
count if abs(resid) < 1000
count if abs(resid) < 2000

* What are the quantiles of residuals?
_pctile resid, p(25 50 75)
display "r_25 = " r(r1)
display "r_50 = " r(r2)
display "r_75 = " r(r3)

* Bandwidth in probability units
local h = e(bwidth)
display "h (prob) = " `h'

* q-h and q+h
local q_lo = 0.5 - `h'
local q_hi = 0.5 + `h'
display "q_lo = " `q_lo'
display "q_hi = " `q_hi'

* Corresponding residual quantiles
_pctile resid, p(`=`q_lo'*100' `=`q_hi'*100')
display "r_lo = " r(r1)
display "r_hi = " r(r2)

local h_resid = (r(r2) - r(r1)) / 2
display "h_resid = " `h_resid'

* Density at 0 using this bandwidth
* f(0) ≈ P(|r| < h_resid) / (2 * h_resid)
count if abs(resid) < `h_resid'
local n_near = r(N)
local f0 = `n_near' / (2 * `h_resid' * 74)
display "Approx f(0) = " `f0'
display "Approx sparsity = " 1/`f0'
