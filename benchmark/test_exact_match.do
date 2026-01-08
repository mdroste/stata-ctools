* Test exact match between cqreg and qreg
clear all
set more off
adopath + "../build"
sysuse auto, clear

display _n "=== EXACT MATCH TEST ===" _n

* Run qreg with robust VCE
quietly qreg price mpg, quantile(0.5) vce(robust)
display "qreg vce(robust):"
display "  beta[mpg]  = " %16.10f _b[mpg]
display "  beta[cons] = " %16.10f _b[_cons]
display "  SE[mpg]    = " %16.10f _se[mpg]
display "  SE[cons]   = " %16.10f _se[_cons]
matrix V_qreg = e(V)
display "  V[1,1]     = " %16.10e V_qreg[1,1]
display "  V[2,2]     = " %16.10e V_qreg[2,2]
local qreg_se_mpg = _se[mpg]
local qreg_se_cons = _se[_cons]

* Check what density method qreg used
display _n "  qreg denmethod = " e(denmethod)
display "  qreg bwidth    = " e(bwidth)
display "  qreg sparsity  = " e(sparsity)
display "  qreg f_r       = " e(f_r)

* Run cqreg with default settings (should be fitted)
display _n "cqreg (default):"
quietly cqreg price mpg, quantile(0.5) vce(robust)
display "  beta[mpg]  = " %16.10f _b[mpg]
display "  beta[cons] = " %16.10f _b[_cons]
display "  SE[mpg]    = " %16.10f _se[mpg]
display "  SE[cons]   = " %16.10f _se[_cons]
matrix V_cqreg = e(V)
display "  V[1,1]     = " %16.10e V_cqreg[1,1]
display "  V[2,2]     = " %16.10e V_cqreg[2,2]
display _n "  cqreg denmethod = " e(denmethod)
display "  cqreg bwidth    = " e(bwidth)
display "  cqreg sparsity  = " e(sparsity)

local cqreg_se_mpg = _se[mpg]

* Calculate ratios
display _n "=== COMPARISON ==="
display "SE[mpg] ratio (cqreg/qreg): " %12.10f `cqreg_se_mpg'/`qreg_se_mpg'
display "V[1,1] ratio:               " %12.10f V_cqreg[1,1]/V_qreg[1,1]

* Also check IID comparison
display _n _n "=== IID COMPARISON ==="
quietly qreg price mpg, quantile(0.5)
display "qreg vce(iid) [default fitted]:"
display "  SE[mpg]    = " %16.10f _se[mpg]
display "  sparsity   = " e(sparsity)
display "  bwidth     = " e(bwidth)

quietly qreg price mpg, quantile(0.5) vce(iid, residual)
display _n "qreg vce(iid, residual):"
display "  SE[mpg]    = " %16.10f _se[mpg]
display "  sparsity   = " e(sparsity)

quietly cqreg price mpg, quantile(0.5) vce(iid) denmethod(residual)
display _n "cqreg vce(iid) denmethod(residual):"
display "  SE[mpg]    = " %16.10f _se[mpg]
display "  sparsity   = " e(sparsity)
