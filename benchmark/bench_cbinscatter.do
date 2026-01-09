* Benchmark cbinscatter performance at scale
* Tests parallelization benefits on large datasets

clear all
set more off

cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath + "build"

di "======================================"
di "cbinscatter Large-Scale Benchmark"
di "======================================"

* Test sizes
local sizes "10000 50000 100000 500000 1000000"

foreach N of local sizes {
    di _n "======================================"
    di "N = `N' observations"
    di "======================================"

    * Generate test data
    clear
    set obs `N'
    gen x = rnormal()
    gen y = 2*x + rnormal()
    gen control1 = rnormal()
    gen control2 = rnormal()
    gen group = mod(_n, 5) + 1
    gen wt = runiform() + 0.5

    * Test 1: Basic binscatter
    di _n "Test 1: Basic (no controls)"
    timer clear 1
    timer on 1
    cbinscatter y x, nquantiles(20) nograph
    timer off 1
    timer list 1

    * Test 2: With controls
    di _n "Test 2: With 2 controls"
    timer clear 2
    timer on 2
    cbinscatter y x, controls(control1 control2) nquantiles(20) nograph
    timer off 2
    timer list 2

    * Test 3: With by()
    di _n "Test 3: With by(group)"
    timer clear 3
    timer on 3
    cbinscatter y x, by(group) nquantiles(20) nograph
    timer off 3
    timer list 3

    * Test 4: With weights
    di _n "Test 4: With aweights"
    timer clear 4
    timer on 4
    cbinscatter y x [aw=wt], nquantiles(20) nograph
    timer off 4
    timer list 4

    * Test 5: Line fitting
    di _n "Test 5: With linetype(linear)"
    timer clear 5
    timer on 5
    cbinscatter y x, nquantiles(20) linetype(linear) nograph
    timer off 5
    timer list 5

    * Test 6: Combined
    di _n "Test 6: Controls + by + weights + linetype"
    timer clear 6
    timer on 6
    cbinscatter y x [aw=wt], controls(control1 control2) by(group) nquantiles(20) linetype(linear) nograph
    timer off 6
    timer list 6
}

di _n "======================================"
di "Benchmark complete"
di "======================================"
