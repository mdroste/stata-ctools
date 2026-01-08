* Benchmark FN solver - accuracy and performance
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "FN Solver Benchmark - Accuracy and Performance"
di "{hline 60}"

* Test accuracy and timing at various sizes
foreach n in 1000 5000 10000 50000 {
    clear
    set seed 12345
    set obs `n'
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2

    di _newline "N=`n'"

    * Time qreg
    timer clear 1
    timer on 1
    quietly qreg y x1 x2
    timer off 1
    matrix b1 = e(b)

    * Time cqreg
    timer clear 2
    timer on 2
    quietly cqreg y x1 x2
    timer off 2
    matrix b2 = e(b)

    * Compute max absolute difference
    local diff = max(abs(b1[1,1]-b2[1,1]), abs(b1[1,2]-b2[1,2]), abs(b1[1,3]-b2[1,3]))

    quietly timer list 1
    local t1 = r(t1)
    quietly timer list 2
    local t2 = r(t2)

    local speedup = `t1'/`t2'

    di "  qreg: " %7.4f `t1' "s  cqreg: " %7.4f `t2' "s  Speedup: " %5.2f `speedup' "x  MaxDiff: " %9.6f `diff'
}

di _newline "{hline 60}"
di "Done"
