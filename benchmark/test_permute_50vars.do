* Benchmark permutation with 10M obs, 50 variables

clear all
set more off
adopath + "build"

local size = 10000000
local nvars = 50

di as text ""
di as text "========================================"
di as text "Dataset size: " `size' " observations"
di as text "Variables:    " `nvars' " numeric variables"
di as text "========================================"

clear
set seed 12345
quietly set obs `size'

forval i = 1/`nvars' {
    quietly gen x`i' = runiform()
}

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
