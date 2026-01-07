* Scaling benchmark - sequential to avoid hangs
clear all
set more off
adopath + "../build"

di "{hline 60}"
di "cqreg vs qreg Scaling Benchmark"
di "{hline 60}"

* N = 1000
clear
set seed 12345
set obs 1000
quietly {
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2
}
timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
local qreg1 = r(t1)
timer clear 2
timer on 2
quietly cqreg y x1 x2
timer off 2
quietly timer list 2
local cqreg1 = r(t2)
di "N=1000:   qreg=" %6.3f `qreg1' "s  cqreg=" %6.3f `cqreg1' "s  speedup=" %5.1f (`qreg1'/`cqreg1') "x"

* N = 5000
clear
set seed 12345
set obs 5000
quietly {
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2
}
timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
local qreg2 = r(t1)
timer clear 2
timer on 2
quietly cqreg y x1 x2
timer off 2
quietly timer list 2
local cqreg2 = r(t2)
di "N=5000:   qreg=" %6.3f `qreg2' "s  cqreg=" %6.3f `cqreg2' "s  speedup=" %5.1f (`qreg2'/`cqreg2') "x"

* N = 10000
clear
set seed 12345
set obs 10000
quietly {
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2
}
timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
local qreg3 = r(t1)
timer clear 2
timer on 2
quietly cqreg y x1 x2
timer off 2
quietly timer list 2
local cqreg3 = r(t2)
di "N=10000:  qreg=" %6.3f `qreg3' "s  cqreg=" %6.3f `cqreg3' "s  speedup=" %5.1f (`qreg3'/`cqreg3') "x"

* N = 25000
clear
set seed 12345
set obs 25000
quietly {
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2
}
timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
local qreg4 = r(t1)
timer clear 2
timer on 2
quietly cqreg y x1 x2
timer off 2
quietly timer list 2
local cqreg4 = r(t2)
di "N=25000:  qreg=" %6.3f `qreg4' "s  cqreg=" %6.3f `cqreg4' "s  speedup=" %5.1f (`qreg4'/`cqreg4') "x"

* N = 50000
clear
set seed 12345
set obs 50000
quietly {
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 1 + 2*x1 - x2 + rnormal()*2
}
timer clear 1
timer on 1
quietly qreg y x1 x2
timer off 1
quietly timer list 1
local qreg5 = r(t1)
timer clear 2
timer on 2
quietly cqreg y x1 x2
timer off 2
quietly timer list 2
local cqreg5 = r(t2)
di "N=50000:  qreg=" %6.3f `qreg5' "s  cqreg=" %6.3f `cqreg5' "s  speedup=" %5.1f (`qreg5'/`cqreg5') "x"

di "{hline 60}"
di "Benchmark completed!"
