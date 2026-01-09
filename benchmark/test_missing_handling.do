* Test that cbinscatter correctly handles missing values
* This verifies the C plugin's missing value validation

clear all
set more off

cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath + "build"

di "========================================"
di "Missing Value Handling Test"
di "========================================"

* Create data with missing values
clear
set obs 1000
gen x = rnormal()
gen y = 2*x + rnormal()
gen control1 = rnormal()
gen wt = runiform() + 0.5

* Introduce missing values
replace y = . in 1/50
replace x = . in 51/100
replace control1 = . in 101/150
replace wt = . in 151/175

di _n "Data summary:"
di "Total obs: " _N
count if missing(y)
di "Missing y: " r(N)
count if missing(x)
di "Missing x: " r(N)
count if missing(control1)
di "Missing control1: " r(N)
count if missing(wt)
di "Missing wt: " r(N)

* Test 1: Basic - should drop missing y and x
di _n "Test 1: Basic (should use 900 obs)"
cbinscatter y x, nquantiles(10) nograph verbose
di "N used: " e(N)
assert e(N) == 900

* Test 2: With controls - should drop missing controls too
di _n "Test 2: With controls (should use 850 obs)"
cbinscatter y x, controls(control1) nquantiles(10) nograph verbose
di "N used: " e(N)
assert e(N) == 850

* Test 3: With weights - should drop missing weights too
di _n "Test 3: With weights (should use 875 obs)"
cbinscatter y x [aw=wt], nquantiles(10) nograph verbose
di "N used: " e(N)
assert e(N) == 875

di _n "========================================"
di "All missing value tests passed!"
di "========================================"
