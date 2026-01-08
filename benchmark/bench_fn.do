* Benchmark FN solver vs qreg
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "Benchmark: Frisch-Newton vs qreg"
di "{hline 60}"

foreach n in 1000 5000 10000 50000 {
    clear
    set seed 12345
    set obs `n'
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2

    di _newline "=== N=`n' ==="

    timer clear 1
    timer on 1
    quietly qreg y x1 x2
    timer off 1
    quietly timer list 1
    local t_qreg = r(t1)
    matrix b_qreg = e(b)

    timer clear 2
    timer on 2
    quietly cqreg y x1 x2, nopreprocess(1)
    timer off 2
    quietly timer list 2
    local t_fn = r(t2)
    matrix b_fn = e(b)

    local diff1 = abs(b_qreg[1,1] - b_fn[1,1])
    local diff2 = abs(b_qreg[1,2] - b_fn[1,2])
    local diff3 = abs(b_qreg[1,3] - b_fn[1,3])
    local max_diff = max(`diff1', `diff2', `diff3')

    di "qreg:  " %6.3f `t_qreg' "s"
    di "FN:    " %6.3f `t_fn' "s"
    di "Speedup: " %5.1f (`t_qreg'/`t_fn') "x"
    di "Max coef diff: " %8.5f `max_diff'
}
