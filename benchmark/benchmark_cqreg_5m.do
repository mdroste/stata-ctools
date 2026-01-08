* Benchmark: cqreg vs qreg on 5 million observations with 2 variables
* Purpose: Identify performance bottlenecks in cqreg

clear all
set more off
adopath + "build"

* Set seed for reproducibility
set seed 12345

di _newline "{hline 78}"
di "BENCHMARK: cqreg vs qreg on 5 million observations"
di "{hline 78}"

* Generate dataset with 5 million observations
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

* Benchmark qreg
di _newline "{hline 78}"
di "QREG (Stata built-in)"
di "{hline 78}"
timer clear 2
timer on 2
quietly qreg y x
timer off 2
quietly timer list 2
local qreg_time = r(t2)
di "qreg time: " %9.3f `qreg_time' " seconds"

matrix qreg_b = e(b)
matrix qreg_V = e(V)
local qreg_N = e(N)
local qreg_sum_adev = e(sum_adev)
local qreg_sum_rdev = e(sum_rdev)
local qreg_r2_p = e(r2_p)

di "Coefficients:"
di "  x:     " %12.6f qreg_b[1,1]
di "  _cons: " %12.6f qreg_b[1,2]
di "Sum of absolute deviations: " %12.4f `qreg_sum_adev'
di "Pseudo R2: " %12.6f `qreg_r2_p'

* Benchmark cqreg with verbose timing
di _newline "{hline 78}"
di "CQREG (C-accelerated) with verbose timing"
di "{hline 78}"
timer clear 3
timer on 3
cqreg y x, verbose
timer off 3
quietly timer list 3
local cqreg_time = r(t3)
di "cqreg total time: " %9.3f `cqreg_time' " seconds"

matrix cqreg_b = e(b)
matrix cqreg_V = e(V)
local cqreg_N = e(N)
local cqreg_sum_adev = e(sum_adev)
local cqreg_sum_rdev = e(sum_rdev)
local cqreg_r2_p = e(r2_p)

* Run cqreg a second time for warm cache comparison
di _newline "{hline 78}"
di "CQREG (warm cache, second run)"
di "{hline 78}"
timer clear 4
timer on 4
quietly cqreg y x, verbose
timer off 4
quietly timer list 4
local cqreg_time2 = r(t4)
di "cqreg time (warm): " %9.3f `cqreg_time2' " seconds"

* Summary comparison
di _newline "{hline 78}"
di "PERFORMANCE SUMMARY"
di "{hline 78}"
di "              qreg         cqreg        speedup"
di "{hline 60}"
di "Time:     " %9.3f `qreg_time' "s    " %9.3f `cqreg_time' "s    " %9.2f `qreg_time'/`cqreg_time' "x"
di ""
di "Coefficient Comparison (should match closely):"
di "              qreg           cqreg          diff"
di "{hline 60}"
di "x:      " %12.6f qreg_b[1,1] "  " %12.6f cqreg_b[1,1] "  " %12.6e qreg_b[1,1]-cqreg_b[1,1]
di "_cons:  " %12.6f qreg_b[1,2] "  " %12.6f cqreg_b[1,2] "  " %12.6e qreg_b[1,2]-cqreg_b[1,2]
di ""
di "Standard Error Comparison:"
di "              qreg           cqreg          diff"
di "{hline 60}"
di "x:      " %12.6f sqrt(qreg_V[1,1]) "  " %12.6f sqrt(cqreg_V[1,1]) "  " %12.6e sqrt(qreg_V[1,1])-sqrt(cqreg_V[1,1])
di "_cons:  " %12.6f sqrt(qreg_V[2,2]) "  " %12.6f sqrt(cqreg_V[2,2]) "  " %12.6e sqrt(qreg_V[2,2])-sqrt(cqreg_V[2,2])
di ""
di "Model Statistics:"
di "              qreg           cqreg          diff"
di "{hline 60}"
di "N:      " %12.0f `qreg_N' "  " %12.0f `cqreg_N' "  " %12.0f `qreg_N'-`cqreg_N'
di "sum_adev:" %12.4f `qreg_sum_adev' "  " %12.4f `cqreg_sum_adev' "  " %12.4e `qreg_sum_adev'-`cqreg_sum_adev'
di "Pseudo R2:" %12.6f `qreg_r2_p' "  " %12.6f `cqreg_r2_p' "  " %12.6e `qreg_r2_p'-`cqreg_r2_p'

di _newline "{hline 78}"
di "END BENCHMARK"
di "{hline 78}"
