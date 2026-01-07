* Compare sparsity estimation between qreg and cqreg
clear all
set more off
adopath + "../build"

sysuse auto, clear

* Run qreg and check all e() returns
qreg price mpg weight
ereturn list

* What is the qreg sparsity?
di "qreg sparsity = " e(sparsity)
di "qreg bandwidth = " e(bwidth)

* Our cqreg
cqreg price mpg weight, verbose
di "cqreg sparsity = " e(sparsity)
di "cqreg bandwidth = " e(bwidth)

* Compare standard errors directly
di ""
di "SE comparison:"
di "qreg  mpg SE = " sqrt(V_qreg[1,1])
di "cqreg mpg SE = " sqrt(e(V)[1,1])
