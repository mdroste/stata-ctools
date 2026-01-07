* Test preprocessing algorithm for cqreg
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "Testing preprocessing algorithm (Chernozhukov et al. 2020)"
di "{hline 60}"

* Test at multiple sizes
foreach N in 1000 5000 10000 25000 50000 {
    clear
    set seed 12345
    set obs `N'
    quietly {
        gen x1 = rnormal()
        gen x2 = rnormal()
        gen x3 = rnormal()
        gen y = 1 + 2*x1 - x2 + 0.5*x3 + rnormal()*2
    }

    * qreg (baseline)
    timer clear 1
    timer on 1
    quietly qreg y x1 x2 x3
    timer off 1
    quietly timer list 1
    local qreg_time = r(t1)
    matrix b_qreg = e(b)

    * cqreg (with preprocessing)
    timer clear 2
    timer on 2
    quietly cqreg y x1 x2 x3
    timer off 2
    quietly timer list 2
    local cqreg_time = r(t2)
    matrix b_cqreg = e(b)

    * Check correctness
    local max_diff = 0
    forval i = 1/4 {
        local diff = abs(b_cqreg[1,`i'] - b_qreg[1,`i'])
        if `diff' > `max_diff' {
            local max_diff = `diff'
        }
    }

    local speedup = `qreg_time' / `cqreg_time'
    local correct = "FAIL"
    if `max_diff' < 0.001 {
        local correct = "PASS"
    }

    di "N=" %6.0fc `N' "  qreg=" %6.3f `qreg_time' "s  cqreg=" %6.3f `cqreg_time' "s  speedup=" %5.1f `speedup' "x  " "`correct'"
}

di "{hline 60}"
di "Test completed!"
