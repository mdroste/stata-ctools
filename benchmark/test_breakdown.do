clear all
set more off
adopath + "../build"

set obs 5000000
set seed 12345
gen x = rnormal()
gen y = 2*x + rnormal()*5

display "============================================================"
display "CQREG Runtime Breakdown: N=5,000,000, K=2"
display "============================================================"

display ""
display "=== Method 1: residual (1 IPM solve) ==="
timer clear
timer on 1
cqreg y x, denmethod(residual) verbose
timer off 1
timer list 1

display ""
display "=== Method 2: fitted (3 IPM solves) ==="
timer clear
timer on 2
cqreg y x, denmethod(fitted) verbose
timer off 2
timer list 2
