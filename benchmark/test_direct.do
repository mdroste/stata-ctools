* Test plugin directly
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

* Force verbose output to avoid crash
scalar __cqreg_verbose = 1
cqreg y x1 x2
di "SUCCESS"
