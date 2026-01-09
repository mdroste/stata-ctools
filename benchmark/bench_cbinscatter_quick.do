* Quick cbinscatter benchmark
* Tests bin computation performance at 5M observations

clear all
set more off

cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath + "build"

di "========================================"
di "cbinscatter Performance Benchmark"
di "========================================"

local N = 5000000
di "Generating `N' observations..."

clear
set obs `N'
gen x = rnormal()
gen y = 2*x + rnormal()
gen wt = runiform() + 0.5

di _n "Test 1: Basic binscatter (no weights)"
timer clear 1
timer on 1
cbinscatter y x, nquantiles(20) nograph
timer off 1
timer list 1

di _n "Test 2: With weights"
timer clear 2
timer on 2
cbinscatter y x [aw=wt], nquantiles(20) nograph
timer off 2
timer list 2

di _n "Test 3: With linear fit"
timer clear 3
timer on 3
cbinscatter y x, nquantiles(20) linetype(linear) nograph
timer off 3
timer list 3

di _n "Test 4: Verbose timing breakdown"
cbinscatter y x, nquantiles(20) verbose nograph

di _n "========================================"
di "Benchmark complete"
di "========================================"
