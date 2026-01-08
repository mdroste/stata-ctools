* Test quietly prefix
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di "About to run quietly cqreg..."
quietly cqreg y x1 x2
di "Done - SUCCESS"
matrix list e(b)
