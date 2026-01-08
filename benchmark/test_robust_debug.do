clear all
set more off
adopath + "../build"
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

display _n "=== DEBUGGING ROBUST VCE ===" _n

* First check what sparsity and bandwidth values we get from cqreg IID
quietly cqreg price mpg weight, quantile(0.5) vce(iid)
local sparsity = e(sparsity)
local bandwidth = e(bwidth)
di "cqreg IID - sparsity: " %12.4f `sparsity'
di "cqreg IID - bandwidth: " %12.6f `bandwidth'

* Compute expected h_residual
local h_residual = `bandwidth' * `sparsity'
di "h_residual = bandwidth * sparsity = " %12.4f `h_residual'

* Check residual distribution from qreg
quietly qreg price mpg weight, quantile(0.5)
predict double resid, residuals
summarize resid, detail
di _n "Observations with |resid| < h_residual: "
count if abs(resid) < `h_residual'

* Check basis observations (residual ≈ 0)
di _n "Basis observations (|resid| < 1e-8): "
count if abs(resid) < 1e-8

* Check residuals in different ranges
di _n "Residual distribution:"
di "  |resid| < 100:    " _continue
count if abs(resid) < 100
di "  |resid| < 500:    " _continue
count if abs(resid) < 500
di "  |resid| < 1000:   " _continue
count if abs(resid) < 1000
di "  |resid| < 2000:   " _continue
count if abs(resid) < 2000

* Now run cqreg with verbose to see more detail
di _n "=== Running cqreg vce(robust) with verbose ==="
cqreg price mpg weight, quantile(0.5) vce(robust) verbose
