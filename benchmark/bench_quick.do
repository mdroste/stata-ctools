* Quick benchmark FN solver
clear all
set more off
adopath + "../build"

set seed 12345
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 1 + 2*x1 - x2 + rnormal()*2

di "=== N=1000 ==="

timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
di "qreg: " r(t1) "s"

timer clear 2
timer on 2
cqreg y x1 x2, nopreprocess(1)
timer off 2
quietly timer list 2
di "FN:   " r(t2) "s"
di "Speedup: " (r(t1)/r(t2)) "x"
