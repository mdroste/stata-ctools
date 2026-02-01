/****************************************************************************
 * benchmark_cqreg_perf.do - Performance comparison of cqreg vs qreg
 ****************************************************************************/

clear all
set more off

* Add build directory to adopath
capture adopath ++ "build"
capture adopath ++ "../build"

* Generate test data
local N = 50000
local K = 10

di as text "=============================================="
di as text "  cqreg Performance Benchmark"
di as text "=============================================="
di as text ""
di as text "Dataset: N = `N' observations, K = `K' regressors"
di as text ""

set seed 12345
set obs `N'

* Generate regressors
forvalues k = 1/`K' {
    gen x`k' = rnormal()
}

* Generate y with heterogeneous variance
gen y = 1 + x1 + 0.5*x2 + 0.3*x3 + (1 + 0.5*abs(x1))*rnormal()

* Warm up
quietly cqreg y x*, q(0.5)
quietly qreg y x*, q(0.5)

* Benchmark at different quantiles
local quantiles "0.25 0.50 0.75 0.90"

di as text "----------------------------------------------"
di as text "  Performance Results (seconds)"
di as text "----------------------------------------------"
di as text "  Quantile    qreg       cqreg      Speedup"
di as text "----------------------------------------------"

foreach q in `quantiles' {
    * Time qreg
    timer clear 1
    timer on 1
    forvalues i = 1/3 {
        quietly qreg y x*, q(`q')
    }
    timer off 1
    quietly timer list 1
    local t_qreg = r(t1) / 3

    * Time cqreg
    timer clear 2
    timer on 2
    forvalues i = 1/3 {
        quietly cqreg y x*, q(`q')
    }
    timer off 2
    quietly timer list 2
    local t_cqreg = r(t2) / 3

    local speedup = `t_qreg' / `t_cqreg'
    di as text "  `q'         " %6.3f `t_qreg' "     " %6.3f `t_cqreg' "     " %5.2f `speedup' "x"
}

di as text "----------------------------------------------"

* Larger dataset test
di as text ""
di as text "Testing with N = 100000..."
clear
local N = 100000
set obs `N'
forvalues k = 1/`K' {
    gen x`k' = rnormal()
}
gen y = 1 + x1 + 0.5*x2 + (1 + 0.5*abs(x1))*rnormal()

timer clear 1
timer on 1
quietly qreg y x*, q(0.5)
timer off 1
quietly timer list 1
local t_qreg = r(t1)

timer clear 2
timer on 2
quietly cqreg y x*, q(0.5)
timer off 2
quietly timer list 2
local t_cqreg = r(t2)

local speedup = `t_qreg' / `t_cqreg'
di as text "  N=100k      " %6.3f `t_qreg' "     " %6.3f `t_cqreg' "     " %5.2f `speedup' "x"

di as text ""
di as text "=============================================="
