clear all
set more off

* Create large test data: N=10M, unique i=1M, unique t=10
local N = 10000000
local num_i = 1000000
local num_t = 10

di "Creating test data: N=`N', unique i=`num_i', unique t=`num_t'"
timer clear 1
timer on 1
quietly {
    set obs `N'
    gen double y = rnormal()
    gen double x = rnormal()
    gen long i = mod(_n-1, `num_i') + 1
    gen int t = mod(_n-1, `num_t') + 1
}
timer off 1
timer list 1

* Run creghdfe with verbose timing
di _n "Running creghdfe absorb(i t) vce(cluster i)..."
creghdfe y x, absorb(i t) vce(cluster i) verbose
