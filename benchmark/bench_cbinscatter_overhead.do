* Benchmark cbinscatter overhead at large scale
* Tests total time vs plugin time to identify .ado overhead

clear all
set more off

cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath + "build"

di "========================================"
di "cbinscatter Overhead Benchmark"
di "========================================"

* Test at 25M observations
local N = 25000000
di "Generating `N' observations..."

clear
set obs `N'
gen x = rnormal()
gen y = 2*x + rnormal()

di _n "Test 1: No graph"
timer clear 1
timer on 1
cbinscatter y x, nquantiles(20) nograph verbose
timer off 1
timer list 1

di _n "Test 2: With graph"
timer clear 2
timer on 2
cbinscatter y x, nquantiles(20) linetype(linear)
timer off 2
timer list 2

di _n "Overhead = Total Stata time - Plugin time"
di "If overhead is large, it's in the .ado file"

di _n "========================================"
di "Benchmark complete"
di "========================================"
