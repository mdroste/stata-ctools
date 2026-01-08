clear all
set more off
adopath + "../build"

* Test with different sizes to understand scaling
foreach n in 100000 500000 1000000 2000000 5000000 {
    clear
    set obs `n'
    set seed 12345
    gen x = rnormal()
    gen y = 2*x + rnormal()*5
    
    timer clear
    timer on 1
    quietly qreg y x
    timer off 1
    quietly timer list 1
    local qreg_time = r(t1)
    
    timer clear
    timer on 2
    cqreg y x, denmethod(residual)
    timer off 2
    quietly timer list 2
    local cqreg_time = r(t2)
    local cqreg_iter = e(iterations)
    
    display "N=`n': qreg=" %5.2f `qreg_time' "s  cqreg=" %5.2f `cqreg_time' "s  iterations=`cqreg_iter'  speedup=" %4.2f `qreg_time'/`cqreg_time' "x"
}
