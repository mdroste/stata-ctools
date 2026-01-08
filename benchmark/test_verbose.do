clear all
set more off
adopath + "../build"

set obs 5000000
set seed 12345
gen x = rnormal()
gen y = 2*x + rnormal()*5

display "=== RESIDUAL METHOD (1 solve) ==="
cqreg y x, denmethod(residual) verbose

display ""
display "=== FITTED METHOD (3 solves) ==="
cqreg y x, denmethod(fitted) verbose
