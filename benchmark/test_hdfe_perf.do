clear all
set more off

* Create test data: N=1M, unique i=100K, unique t=10
local N = 1000000
local num_i = 100000
local num_t = 10

di "Creating test data: N=`N', unique i=`num_i', unique t=`num_t'"
quietly {
    set obs `N'
    gen double y = rnormal()
    gen double x = rnormal()
    gen long i = mod(_n-1, `num_i') + 1
    gen int t = mod(_n-1, `num_t') + 1
}
di "Data created."

* Run creghdfe with verbose timing
di _n "Running creghdfe absorb(i t)..."
creghdfe y x, absorb(i t) verbose
