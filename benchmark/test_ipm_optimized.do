* Test optimized IPM implementation
clear all
set more off
adopath + "../build"

* =========================================================
* Test 1: Small dataset - verify correctness
* =========================================================

sysuse auto, clear

di _newline
di as text "{hline 70}"
di as text "Test 1: Correctness check (auto dataset, N=74)"
di as text "{hline 70}"

* Run qreg
di _newline
di as result "=== Stata's qreg ===" _newline
qreg price mpg weight, nolog
matrix b_qreg = e(b)
scalar qreg_obj = e(sum_adev)

* Run cqreg
di _newline
di as result "=== ctools cqreg (optimized IPM) ===" _newline
cqreg price mpg weight
matrix b_cqreg = e(b)

* Compare coefficients
di _newline
di as text "{hline 70}"
di as text "Coefficient Comparison:"
di as text "{hline 70}"
di as text "Variable     " _col(20) "qreg" _col(35) "cqreg" _col(50) "Diff"

local names : colnames b_qreg
local i = 1
local max_diff = 0
local pass = 1
foreach v of local names {
    local diff = b_cqreg[1,`i'] - b_qreg[1,`i']
    local abs_diff = abs(`diff')
    if `abs_diff' > `max_diff' {
        local max_diff = `abs_diff'
    }
    if `abs_diff' > 1e-4 {
        local pass = 0
    }
    di as text "`v'" _col(20) as result %12.6f b_qreg[1,`i'] _col(35) as result %12.6f b_cqreg[1,`i'] _col(50) as result %12.6f `diff'
    local i = `i' + 1
}

di _newline
di as text "Maximum absolute difference: " as result `max_diff'
if `pass' {
    di as result "PASS: All coefficients match within 1e-4"
}
else {
    di as error "FAIL: Some coefficients differ by more than 1e-4"
}

* =========================================================
* Test 2: Performance benchmark with larger dataset
* =========================================================

di _newline(2)
di as text "{hline 70}"
di as text "Test 2: Performance benchmark"
di as text "{hline 70}"

* Generate larger dataset
clear
set seed 12345
local N = 100000

di as text "Generating test data with N = " as result `N' as text " observations..."
quietly {
    set obs `N'
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen x3 = rnormal()
    gen x4 = rnormal()
    gen x5 = rnormal()
    gen y = 1 + 2*x1 - 1.5*x2 + 0.5*x3 + x4 - 0.25*x5 + rnormal()*2 + (runiform() - 0.5)*abs(x1)*3
}

* Benchmark qreg
di _newline
di as result "=== Stata's qreg (N=`N') ===" _newline
timer clear 1
timer on 1
quietly qreg y x1 x2 x3 x4 x5
timer off 1
quietly timer list 1
local qreg_time = r(t1)
di as text "Time: " as result %9.3f `qreg_time' as text " seconds"
matrix b_qreg = e(b)

* Benchmark cqreg
di _newline
di as result "=== ctools cqreg (optimized IPM, N=`N') ===" _newline
timer clear 2
timer on 2
cqreg y x1 x2 x3 x4 x5
timer off 2
quietly timer list 2
local cqreg_time = r(t2)
di as text "Time: " as result %9.3f `cqreg_time' as text " seconds"
matrix b_cqreg = e(b)

* Compare
di _newline
di as text "{hline 70}"
di as text "Comparison:"
di as text "{hline 70}"
di as text "qreg time:     " as result %9.3f `qreg_time' as text " seconds"
di as text "cqreg time:    " as result %9.3f `cqreg_time' as text " seconds"
di as text "Speedup:       " as result %9.2f (`qreg_time'/`cqreg_time') as text "x"

* Verify coefficients match
local names : colnames b_qreg
local i = 1
local max_diff = 0
foreach v of local names {
    local diff = abs(b_cqreg[1,`i'] - b_qreg[1,`i'])
    if `diff' > `max_diff' {
        local max_diff = `diff'
    }
    local i = `i' + 1
}
di as text "Max coef diff: " as result %9.6f `max_diff'

* =========================================================
* Test 3: Very large dataset
* =========================================================

di _newline(2)
di as text "{hline 70}"
di as text "Test 3: Large-scale benchmark"
di as text "{hline 70}"

* Generate larger dataset
clear
set seed 12345
local N = 500000

di as text "Generating test data with N = " as result `N' as text " observations..."
quietly {
    set obs `N'
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen x3 = rnormal()
    gen y = 1 + 2*x1 - 1.5*x2 + 0.5*x3 + rnormal()*2 + (runiform() - 0.5)*abs(x1)*3
}

* Benchmark cqreg only (qreg is too slow)
di _newline
di as result "=== ctools cqreg (N=`N') ===" _newline
timer clear 3
timer on 3
cqreg y x1 x2 x3, timeit
timer off 3
quietly timer list 3
local cqreg_time = r(t3)
di as text "Total cqreg time: " as result %9.3f `cqreg_time' as text " seconds"
di as text "Throughput: " as result %9.0fc (`N'/`cqreg_time') as text " obs/sec"

di _newline
di as text "{hline 70}"
di as text "All tests completed!"
di as text "{hline 70}"
