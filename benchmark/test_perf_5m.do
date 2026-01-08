clear all
set more off
adopath + "../build"

* Create large dataset - 5 million obs, 2 variables (matching user's setup)
set obs 5000000
set seed 12345
gen x = rnormal()
gen y = 2*x + rnormal()*5

display "============================================================"
display "Performance test: N=5,000,000, K=2"
display "============================================================"

* Time qreg
timer clear
timer on 1
quietly qreg y x
timer off 1
quietly timer list 1
local qreg_time = r(t1)
display "qreg time: " %6.3f `qreg_time' " seconds"

* Time cqreg with residual (1 QR solve - fastest)
timer clear
timer on 2
quietly cqreg y x, denmethod(residual)
timer off 2
quietly timer list 2
local cqreg_resid_time = r(t2)
display "cqreg (residual) time: " %6.3f `cqreg_resid_time' " seconds"

* Time cqreg with fitted (3 QR solves for Siddiqui)
timer clear
timer on 3
quietly cqreg y x
timer off 3
quietly timer list 3
local cqreg_fitted_time = r(t3)
display "cqreg (fitted) time: " %6.3f `cqreg_fitted_time' " seconds"

display ""
display "Speedup vs qreg:"
display "  cqreg (residual): " %5.2f `qreg_time'/`cqreg_resid_time' "x"
display "  cqreg (fitted):   " %5.2f `qreg_time'/`cqreg_fitted_time' "x"
