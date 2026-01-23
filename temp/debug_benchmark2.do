* Debug benchmark_import - why varnames(1) not applied

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* Load helpers
quietly do "validation/validate_setup.do"

* Test with all numeric headers file
file open fh using "temp/all_numeric_headers.csv", write replace
file write fh "1,2,3,4,5" _n
file write fh "10,20,30,40,50" _n
file write fh "11,21,31,41,51" _n
file close fh

di _n "=== Test: benchmark_import with all numeric headers ==="

* Run benchmark_import with debug output
di "Running benchmark_import..."
benchmark_import using "temp/all_numeric_headers.csv", testname("all numeric headers")

* Also test webuse datasets directly
di _n "=== Test: webuse airline ==="
di "Stata:"
webuse airline, clear
cexport delimited using "temp/airline.csv", replace
import delimited using "temp/airline.csv", clear varnames(1)
di "Stata N=" _N " K=" c(k)

di _n "cimport:"
cimport delimited using "temp/airline.csv", clear
di "cimport N=" _N " K=" c(k)
