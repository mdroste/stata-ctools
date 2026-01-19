* Benchmark permutation performance with large dataset (10M obs, 10+ vars)
* This tests the use case the user cares about most

clear all
set more off

* Add build directory to adopath
adopath + "build"

* Test with 10M observations, 10 variables
local size = 10000000
local nvars = 10

di as text ""
di as text "========================================"
di as text "Dataset size: " `size' " observations"
di as text "Variables:    " `nvars' " numeric variables"
di as text "========================================"

clear
set seed 12345
quietly set obs `size'

* Create multiple variables
forval i = 1/`nvars' {
    quietly gen x`i' = runiform()
}

* Run csort with verbose to see timing breakdown
di as text ""
di as text "Sorting with IPS4O (default):"
csort x1, verbose

di as text ""
di as text "Sort time:    " %8.4f _csort_time_sort " sec"
di as text "Permute time: " %8.4f _csort_time_permute " sec"
di as text "Ratio (permute/sort): " %5.2f (_csort_time_permute / _csort_time_sort)

di as text ""
di as text "========================================"
di as text "Benchmark complete!"
di as text "========================================"
