clear all
set more off
adopath + "../build"

set obs 5000000
set seed 12345
gen x = rnormal()
gen y = 2*x + rnormal()*5

display "=== Residual method (1 QR solve) ==="
timer clear
timer on 1
quietly cqreg y x, denmethod(residual)
timer off 1
quietly timer list 1
display "Total time: " %5.2f r(t1) " seconds"

display ""
display "=== Fitted method (3 QR solves for Siddiqui) ==="
timer clear
timer on 2
quietly cqreg y x, denmethod(fitted)
timer off 2
quietly timer list 2
display "Total time: " %5.2f r(t2) " seconds"

display ""
display "=== qreg default ==="
timer clear
timer on 3
quietly qreg y x
timer off 3
quietly timer list 3
display "Total time: " %5.2f r(t3) " seconds"
