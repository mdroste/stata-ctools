* Debug benchmark_import varnames(1) fix

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* Load helpers
quietly do "validation/validate_setup.do"

* Test with all numeric headers file
file open fh using "temp/numeric_headers.csv", write replace
file write fh "1,2,3,4,5" _n
file write fh "10,20,30,40,50" _n
file write fh "11,21,31,41,51" _n
file close fh

di _n "=== DEBUG: All numeric headers test ==="
di "File contents:"
type "temp/numeric_headers.csv"

di _n "Test 1: Stata with default options"
import delimited using "temp/numeric_headers.csv", clear
di "N = " _N ", K = " c(k)
desc, short

di _n "Test 2: Stata with varnames(1)"
import delimited using "temp/numeric_headers.csv", clear varnames(1)
di "N = " _N ", K = " c(k)
desc, short

di _n "Test 3: cimport with default options"
cimport delimited using "temp/numeric_headers.csv", clear
di "N = " _N ", K = " c(k)
desc, short

* Now test through benchmark_import
di _n "=== Running benchmark_import ==="
benchmark_import using "temp/numeric_headers.csv", testname("test numeric headers via benchmark")
