* Test smoothed IPM only
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "Testing Smoothed IPM Solver"
di "{hline 60}"

* Generate test data
set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di _newline "=== N=1000 ==="

* Baseline: qreg
timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
di "qreg:         " r(t1) "s"
matrix b_qreg = e(b)
di "  beta: " %9.6f b_qreg[1,1] "  " %9.6f b_qreg[1,2] "  " %9.6f b_qreg[1,3]

* Test smoothed IPM (nopreprocess=2)
di _newline "Testing smoothed IPM (nopreprocess(2))..."
timer clear 2
timer on 2
cqreg y x1 x2, verbose nopreprocess(2)
timer off 2
quietly timer list 2
di "cqreg (IPM):  " r(t2) "s"
matrix b_ipm = e(b)
di "  beta: " %9.6f b_ipm[1,1] "  " %9.6f b_ipm[1,2] "  " %9.6f b_ipm[1,3]

local diff1 = abs(b_qreg[1,1] - b_ipm[1,1])
local diff2 = abs(b_qreg[1,2] - b_ipm[1,2])
local diff3 = abs(b_qreg[1,3] - b_ipm[1,3])
di "  Diffs: " `diff1' "  " `diff2' "  " `diff3'

if `diff1' < 0.01 & `diff2' < 0.01 & `diff3' < 0.01 {
    di "  PASS - coefficients match qreg"
}
else {
    di "  FAIL - coefficients differ from qreg"
}
