* Quick test
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
di "qreg: " r(t1) "s"
matrix b_qreg = e(b)

timer clear 2
timer on 2
quietly cqreg y x1 x2
timer off 2
quietly timer list 2
di "cqreg: " r(t2) "s"
matrix b_cqreg = e(b)

di "qreg coefs:  " b_qreg[1,1] " " b_qreg[1,2] " " b_qreg[1,3]
di "cqreg coefs: " b_cqreg[1,1] " " b_cqreg[1,2] " " b_cqreg[1,3]

local diff1 = abs(b_qreg[1,1] - b_cqreg[1,1])
local diff2 = abs(b_qreg[1,2] - b_cqreg[1,2])
local diff3 = abs(b_qreg[1,3] - b_cqreg[1,3])
di "Diffs: " `diff1' " " `diff2' " " `diff3'

if `diff1' < 0.01 & `diff2' < 0.01 & `diff3' < 0.01 {
    di "PASS - coefficients match"
}
else {
    di "FAIL - coefficients differ"
}
