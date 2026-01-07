* Test at N=10000 with verbose output
clear all
set more off
adopath + "../build"

set seed 12345
set obs 10000
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 1 + 2*x1 - x2 + 0.5*x3 + rnormal()*2

di "Testing cqreg at N=10000..."
timer clear 1
timer on 1
cqreg y x1 x2 x3, verbose
timer off 1
quietly timer list 1
di "cqreg time: " r(t1) " seconds"

di "Done!"
