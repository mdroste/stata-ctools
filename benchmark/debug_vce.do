* Debug VCE computation
clear all
set more off
adopath + "../build"

sysuse auto, clear

* Run qreg
qreg price mpg weight
matrix V_qreg = e(V)
di "qreg VCE:"
matrix list V_qreg

di ""
di "qreg sparsity = " e(sparsity)
di "qreg N = " e(N)

* Compute (X'X)^(-1)
matrix accum XtX = mpg weight
matrix XtX_inv = syminv(XtX)
di ""
di "(X'X)^(-1):"
matrix list XtX_inv

* What scale does qreg use?
* V = scale * (X'X)^(-1)
* scale = V[1,1] / XtX_inv[1,1]
local scale = V_qreg[1,1] / XtX_inv[1,1]
di ""
di "Implied scale = " `scale'

* What should scale be according to formula?
* scale = sparsity^2 * q*(1-q) / n
local expected_scale = e(sparsity)^2 * 0.25 / 74
di "Expected scale (sparsity^2 * 0.25 / 74) = " `expected_scale'

* Ratio
di "Ratio = " `scale' / `expected_scale'
