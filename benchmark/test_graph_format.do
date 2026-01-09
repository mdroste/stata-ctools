* Test cbinscatter graph format
* Compare visual output to verify binscatter-style formatting

clear all
set more off

cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath + "build"

* Load test data
sysuse auto, clear

* Test 1: Basic plot (single group, with linear fit)
di "Test 1: Basic plot with linear fit"
cbinscatter price mpg, nquantiles(10) linetype(linear)
graph export "benchmark/graph_test1_basic.png", replace

* Test 2: No fit line
di "Test 2: No fit line"
cbinscatter price mpg, nquantiles(10) linetype(none)
graph export "benchmark/graph_test2_nofit.png", replace

* Test 3: Quadratic fit
di "Test 3: Quadratic fit"
cbinscatter price mpg, nquantiles(10) linetype(qfit)
graph export "benchmark/graph_test3_qfit.png", replace

* Test 4: By groups
di "Test 4: By groups"
cbinscatter price mpg, by(foreign) nquantiles(10) linetype(linear)
graph export "benchmark/graph_test4_by.png", replace

* Test 5: Custom titles
di "Test 5: Custom titles"
cbinscatter price mpg, nquantiles(10) linetype(linear) ///
    title("Price vs MPG") ytitle("Price (USD)") xtitle("Miles per Gallon")
graph export "benchmark/graph_test5_titles.png", replace

di "Graph format tests complete!"
