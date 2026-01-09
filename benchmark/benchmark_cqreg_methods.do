* Benchmark: cqreg fitted vs residual density methods
* Purpose: Compare performance of optimized sparsity estimation methods

clear all
set more off
adopath + "build"

set seed 12345

di _newline "{hline 78}"
di "BENCHMARK: cqreg density estimation methods on 5 million observations"
di "{hline 78}"

* Generate dataset
di "Generating data..."
timer clear 1
timer on 1
quietly set obs 5000000
quietly gen double x = rnormal()
quietly gen double y = 2 + 3*x + rnormal()*2
timer off 1
quietly timer list 1
di "Data generation time: " %9.3f r(t1) " seconds"

di _newline "{hline 78}"
di "Dataset: N = " %12.0fc _N ", Variables: y, x"
di "{hline 78}"

* Benchmark qreg (reference)
di _newline "{hline 78}"
di "QREG (Stata built-in, default fitted method)"
di "{hline 78}"
timer clear 2
timer on 2
quietly qreg y x
timer off 2
quietly timer list 2
local qreg_time = r(t2)
di "qreg time: " %9.3f `qreg_time' " seconds"
matrix qreg_b = e(b)
local qreg_se_x = sqrt(e(V)[1,1])
local qreg_se_cons = sqrt(e(V)[2,2])

* Benchmark cqreg with FITTED method (default)
di _newline "{hline 78}"
di "CQREG with FITTED method (Siddiqui, default)"
di "{hline 78}"
timer clear 3
timer on 3
cqreg y x, verbose
timer off 3
quietly timer list 3
local cqreg_fitted_time = r(t3)
di "cqreg fitted time: " %9.3f `cqreg_fitted_time' " seconds"
matrix cqreg_fitted_b = e(b)
local cqreg_fitted_se_x = sqrt(e(V)[1,1])
local cqreg_fitted_se_cons = sqrt(e(V)[2,2])

* Benchmark cqreg with RESIDUAL method
di _newline "{hline 78}"
di "CQREG with RESIDUAL method (quickselect optimized)"
di "{hline 78}"
timer clear 4
timer on 4
cqreg y x, denmethod(residual) verbose
timer off 4
quietly timer list 4
local cqreg_residual_time = r(t4)
di "cqreg residual time: " %9.3f `cqreg_residual_time' " seconds"
matrix cqreg_residual_b = e(b)
local cqreg_residual_se_x = sqrt(e(V)[1,1])
local cqreg_residual_se_cons = sqrt(e(V)[2,2])

* Summary comparison
di _newline "{hline 78}"
di "PERFORMANCE SUMMARY"
di "{hline 78}"
di "                           qreg      cqreg(fitted)  cqreg(residual)"
di "{hline 70}"
di "Time (seconds):        " %9.3f `qreg_time' "       " %9.3f `cqreg_fitted_time' "       " %9.3f `cqreg_residual_time'
di "Speedup vs qreg:                     " %9.2f `qreg_time'/`cqreg_fitted_time' "x          " %9.2f `qreg_time'/`cqreg_residual_time' "x"
di ""
di "Coefficient x:         " %12.6f qreg_b[1,1] "  " %12.6f cqreg_fitted_b[1,1] "  " %12.6f cqreg_residual_b[1,1]
di "SE(x):                 " %12.6f `qreg_se_x' "  " %12.6f `cqreg_fitted_se_x' "  " %12.6f `cqreg_residual_se_x'
di ""
di "Coefficient _cons:     " %12.6f qreg_b[1,2] "  " %12.6f cqreg_fitted_b[1,2] "  " %12.6f cqreg_residual_b[1,2]
di "SE(_cons):             " %12.6f `qreg_se_cons' "  " %12.6f `cqreg_fitted_se_cons' "  " %12.6f `cqreg_residual_se_cons'
di ""
di "SE difference from qreg:"
di "  fitted vs qreg x:    " %12.6e `cqreg_fitted_se_x' - `qreg_se_x' " (" %6.3f 100*(`cqreg_fitted_se_x'-`qreg_se_x')/`qreg_se_x' "%)"
di "  residual vs qreg x:  " %12.6e `cqreg_residual_se_x' - `qreg_se_x' " (" %6.3f 100*(`cqreg_residual_se_x'-`qreg_se_x')/`qreg_se_x' "%)"

di _newline "{hline 78}"
di "END BENCHMARK"
di "{hline 78}"
