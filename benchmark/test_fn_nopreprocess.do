* Test Frisch-Newton solver without preprocessing
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "Testing Frisch-Newton Solver (No Preprocessing)"
di "{hline 60}"

* Small test first
set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di _newline "=== N=1000 ==="

timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
di "qreg:  " r(t1) "s"
matrix b_qreg = e(b)

timer clear 2
timer on 2
cqreg y x1 x2, verbose nopreprocess
timer off 2
quietly timer list 2
di "cqreg: " r(t2) "s"
matrix b_cqreg = e(b)

di _newline "Coefficients:"
di "qreg:  x1=" %9.6f b_qreg[1,1] "  x2=" %9.6f b_qreg[1,2] "  cons=" %9.6f b_qreg[1,3]
di "cqreg: x1=" %9.6f b_cqreg[1,1] "  x2=" %9.6f b_cqreg[1,2] "  cons=" %9.6f b_cqreg[1,3]

local diff1 = abs(b_qreg[1,1] - b_cqreg[1,1])
local diff2 = abs(b_qreg[1,2] - b_cqreg[1,2])
local diff3 = abs(b_qreg[1,3] - b_cqreg[1,3])
di "Diffs: " `diff1' "  " `diff2' "  " `diff3'

if `diff1' < 0.1 & `diff2' < 0.1 & `diff3' < 0.1 {
    di "PASS - coefficients close"
}
else {
    di "FAIL - coefficients differ"
}
