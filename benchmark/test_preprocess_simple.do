* Simple test of preprocessing algorithm
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "Simple preprocessing test"
di "{hline 60}"

* Start small
set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 1 + 2*x1 - x2 + 0.5*x3 + rnormal()*2

di _newline "=== N=1000 ==="

di "Running qreg..."
timer clear 1
timer on 1
qreg y x1 x2 x3, nolog
timer off 1
quietly timer list 1
di "qreg time: " r(t1) "s"
matrix b_qreg = e(b)

di _newline "Running cqreg with verbose..."
timer clear 2
timer on 2
cqreg y x1 x2 x3, verbose
timer off 2
quietly timer list 2
di "cqreg time: " r(t2) "s"
matrix b_cqreg = e(b)

di _newline "Comparing coefficients:"
di "qreg:  " b_qreg[1,1] " " b_qreg[1,2] " " b_qreg[1,3] " " b_qreg[1,4]
di "cqreg: " b_cqreg[1,1] " " b_cqreg[1,2] " " b_cqreg[1,3] " " b_cqreg[1,4]

di _newline "Done!"
