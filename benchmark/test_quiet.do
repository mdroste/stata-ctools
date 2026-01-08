* Test quietly
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di "Testing quietly..."
quietly cqreg y x1 x2
di "Done - quietly worked"
matrix list e(b)
