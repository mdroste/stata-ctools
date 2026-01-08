clear all
set more off
adopath + "../build"

* Create large dataset
set obs 1000000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 2*x1 + 3*x2 - x3 + rnormal()*5

display "============================================================"
display "Performance test: N=1,000,000"
display "============================================================"

* Time qreg
timer clear
timer on 1
quietly qreg y x1 x2 x3
timer off 1
quietly timer list 1
local qreg_time = r(t1)
display "qreg time: " %6.3f `qreg_time' " seconds"

* Time cqreg with fitted (uses Siddiqui - 3 QR solves)
timer clear
timer on 2
quietly cqreg y x1 x2 x3
timer off 2
quietly timer list 2
local cqreg_fitted_time = r(t2)
display "cqreg (fitted) time: " %6.3f `cqreg_fitted_time' " seconds"

* Time cqreg with residual (1 QR solve)
timer clear
timer on 3
quietly cqreg y x1 x2 x3, denmethod(residual)
timer off 3
quietly timer list 3
local cqreg_resid_time = r(t3)
display "cqreg (residual) time: " %6.3f `cqreg_resid_time' " seconds"

display ""
display "Speedup vs qreg:"
display "  cqreg (fitted):   " %5.2f `qreg_time'/`cqreg_fitted_time' "x"
display "  cqreg (residual): " %5.2f `qreg_time'/`cqreg_resid_time' "x"
