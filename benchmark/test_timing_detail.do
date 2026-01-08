clear all
set more off
adopath + "../build"

set obs 5000000
set seed 12345
gen x = rnormal()
gen y = 2*x + rnormal()*5

display "=== CQREG TIMING BREAKDOWN ==="
cqreg y x, denmethod(residual)
