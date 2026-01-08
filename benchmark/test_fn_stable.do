* Test FN solver stability
clear all
set more off
adopath + "../build"

di "Testing FN solver stability..."

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

quietly qreg y x1 x2
matrix b_qreg = e(b)

cqreg y x1 x2, verbose
matrix b_cqreg = e(b)

di _newline "Results:"
di "qreg:  " %9.6f b_qreg[1,1] "  " %9.6f b_qreg[1,2] "  " %9.6f b_qreg[1,3]
di "cqreg: " %9.6f b_cqreg[1,1] "  " %9.6f b_cqreg[1,2] "  " %9.6f b_cqreg[1,3]

local diff = max(abs(b_qreg[1,1]-b_cqreg[1,1]), abs(b_qreg[1,2]-b_cqreg[1,2]), abs(b_qreg[1,3]-b_cqreg[1,3]))
di "Max diff: " `diff'
