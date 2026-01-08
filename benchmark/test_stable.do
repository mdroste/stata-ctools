* Test stable IPM solver
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di "=== Testing stable IPM ==="

quietly qreg y x1 x2
matrix b_qreg = e(b)
di "qreg:  " %9.6f b_qreg[1,1] "  " %9.6f b_qreg[1,2] "  " %9.6f b_qreg[1,3]

cqreg y x1 x2
matrix b_cqreg = e(b)
di "cqreg: " %9.6f b_cqreg[1,1] "  " %9.6f b_cqreg[1,2] "  " %9.6f b_cqreg[1,3]

local diff = max(abs(b_qreg[1,1]-b_cqreg[1,1]), abs(b_qreg[1,2]-b_cqreg[1,2]), abs(b_qreg[1,3]-b_cqreg[1,3]))
di "Max diff: " `diff'

if `diff' < 0.0001 {
    di "PASS - coefficients match perfectly"
}
else {
    di "FAIL"
}
