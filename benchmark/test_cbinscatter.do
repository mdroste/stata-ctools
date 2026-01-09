* Test cbinscatter functionality
* Part of the ctools test suite

clear all
set more off

di "======================================"
di "cbinscatter Test Suite"
di "======================================"

* Change to project root
cd "/Users/Mike/Documents/GitHub/stata-ctools"

* Ensure build directory is in adopath
adopath + "build"

* Test 1: Basic usage with auto data
di _n "Test 1: Basic cbinscatter with auto data"
di "----------------------------------------"
sysuse auto, clear
describe

* Simple binscatter of price vs mpg
cbinscatter price mpg, nquantiles(10) verbose nograph

di "N = " e(N)
di "nquantiles = " e(nquantiles)
di "num_groups = " e(num_groups)
matrix list e(bindata)

* Test 2: With controls
di _n "Test 2: cbinscatter with controls"
di "----------------------------------------"
sysuse auto, clear
cbinscatter price mpg, controls(weight length) nquantiles(10) verbose nograph

di "N = " e(N)
matrix list e(bindata)

* Test 3: With absorb
di _n "Test 3: cbinscatter with absorb"
di "----------------------------------------"
sysuse auto, clear
* Create a non-missing rep78 for more observations
replace rep78 = 3 if rep78 == .
cbinscatter price mpg, absorb(foreign rep78) nquantiles(10) verbose nograph

di "N = " e(N)
matrix list e(bindata)

* Test 4: With by()
di _n "Test 4: cbinscatter with by()"
di "----------------------------------------"
sysuse auto, clear
cbinscatter price mpg, by(foreign) nquantiles(10) verbose nograph

di "N = " e(N)
di "num_groups = " e(num_groups)
matrix list e(bindata)

* Test 5: With weights
di _n "Test 5: cbinscatter with aweights"
di "----------------------------------------"
sysuse auto, clear
gen wt = _n
cbinscatter price mpg [aweight=wt], nquantiles(10) verbose nograph

di "N = " e(N)
matrix list e(bindata)

* Test 6: All options combined
di _n "Test 6: Combined options"
di "----------------------------------------"
sysuse auto, clear
replace rep78 = 3 if rep78 == .
cbinscatter price mpg, controls(weight) absorb(rep78) by(foreign) nquantiles(8) linetype(linear) verbose nograph

di "N = " e(N)
di "num_groups = " e(num_groups)
matrix list e(bindata)
matrix list e(coefs)

* Test 7: Discrete mode
di _n "Test 7: Discrete mode"
di "----------------------------------------"
sysuse auto, clear
cbinscatter price rep78 if rep78 != ., discrete verbose nograph

di "N = " e(N)
matrix list e(bindata)

* Test 8: Generate graph (uncomment to view)
di _n "Test 8: Generate graph"
di "----------------------------------------"
sysuse auto, clear
cbinscatter price mpg, nquantiles(20) linetype(linear) title("Price vs MPG") verbose
graph export "benchmark/cbinscatter_test.png", replace

di _n "======================================"
di "All tests completed successfully!"
di "======================================"
