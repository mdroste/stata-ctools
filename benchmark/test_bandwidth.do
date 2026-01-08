clear all
set more off
sysuse auto, clear
drop if missing(price, mpg, weight, rep78)

qreg price mpg weight, quantile(0.5)
predict double resid, resid

local N = 69
local q = 0.5
local alpha = 0.05
local z_alpha = invnormal(1 - `alpha'/2)
local z_q = invnormal(`q')
local phi_zq = normalden(`z_q')

* Hall-Sheather bandwidth
local h_hs = `N'^(-1/3) * `z_alpha'^(2/3) * ((1.5 * `phi_zq'^2) / (2*`z_q'^2 + 1))^(1/3)

* Bofinger bandwidth
local h_bof = ((9/2) * `phi_zq'^4 * (2*`z_q'^2 + 1)^2 / `N')^(1/5)

* Chamberlain bandwidth
local h_cham = `z_alpha' * sqrt(`q' * (1 - `q') / `N') / `phi_zq'

display "=== BANDWIDTH METHODS ==="
display "Hall-Sheather h: " %12.8f `h_hs'
display "Bofinger h:      " %12.8f `h_bof'
display "Chamberlain h:   " %12.8f `h_cham'

* Test with each bandwidth
foreach bw in hs bof cham {
    local h = `h_`bw''
    
    local q_lo = `q' - `h'
    local q_hi = `q' + `h'
    if `q_lo' < 0 local q_lo = 0
    if `q_hi' > 1 local q_hi = 1
    
    _pctile resid, percentile(`=`q_lo'*100' `=`q_hi'*100')
    local r_lo = r(r1)
    local r_hi = r(r2)
    
    local sparsity = (`r_hi' - `r_lo') / (`q_hi' - `q_lo')
    display "`bw' sparsity: " %12.4f `sparsity'
}

* The implied qreg sparsity
gen cons = 1
mkmat mpg weight cons, matrix(X)
matrix XtX = X' * X
matrix XtX_inv = syminv(XtX)
local v11 = e(V)[1,1]
local implied_s = sqrt(`v11' / (0.25 * XtX_inv[1,1]))
display _n "qreg implied sparsity: " %12.4f `implied_s'

* Check what residual quantile Stata is using 
* The quantile at median should be 0
_pctile resid, percentile(50)
display _n "Median residual: " %12.4f r(r1)

* Check IQR-based scaling
_pctile resid, percentile(25 50 75)
local iqr = r(r3) - r(r1)
display "IQR of residuals: " %12.4f `iqr'

* Stata's fitted density at median - let's use kdensity
kdensity resid, at(0) generate(resid_den) n(1)
summ resid_den
local kde_f0 = r(mean)
display _n "KDE density at 0: " %12.8f `kde_f0'
display "KDE sparsity (1/f): " %12.4f (1/`kde_f0')
