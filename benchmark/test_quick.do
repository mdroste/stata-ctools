* Quick performance test
clear all
set more off
adopath + "../build"

* Test correctness first
sysuse auto, clear
di "Testing correctness..."
qreg price mpg weight, nolog
matrix b_qreg = e(b)

cqreg price mpg weight
matrix b_cqreg = e(b)

* Check coefficients match
local max_diff = 0
forval i = 1/3 {
    local diff = abs(b_cqreg[1,`i'] - b_qreg[1,`i'])
    if `diff' > `max_diff' {
        local max_diff = `diff'
    }
}
di "Max coefficient diff: " `max_diff'
if `max_diff' < 1e-4 {
    di "PASS: Coefficients match"
}
else {
    di "FAIL: Coefficients differ"
}

* Small benchmark (10k obs)
clear
set seed 12345
set obs 10000
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 1 + 2*x1 - x2 + 0.5*x3 + rnormal()*2

di _newline
di "{hline 50}"
di "Benchmark N=10,000"
di "{hline 50}"

* qreg
timer clear 1
timer on 1
quietly qreg y x1 x2 x3
timer off 1
quietly timer list 1
di "qreg time:  " %6.3f r(t1) " seconds"

* cqreg
timer clear 2
timer on 2
quietly cqreg y x1 x2 x3
timer off 2
quietly timer list 2
di "cqreg time: " %6.3f r(t2) " seconds"

di "Speedup: " %6.2f (r(t1)/r(t2)) "x"

di _newline
di "Test completed!"
