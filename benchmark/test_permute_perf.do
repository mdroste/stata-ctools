* Benchmark permutation performance with large dataset
* Tests the optimized scatter-based permutation vs original gather

clear all
set more off

* Add build directory to adopath
adopath + "build"

* Test with different dataset sizes
foreach size in 100000 500000 1000000 {
    di as text ""
    di as text "========================================"
    di as text "Dataset size: " `size' " observations"
    di as text "========================================"

    clear
    set seed 12345
    quietly set obs `size'

    * Create multiple variables to stress the permutation
    forval i = 1/10 {
        quietly gen x`i' = runiform()
    }
    quietly gen str20 name = "name" + string(_n)

    * Run csort with verbose to see timing breakdown
    di as text ""
    di as text "Sorting with IPS4O (default):"
    csort x1, verbose

    di as text ""
    di as text "Sort time:    " %8.4f _csort_time_sort " sec"
    di as text "Permute time: " %8.4f _csort_time_permute " sec"
    di as text "Ratio (permute/sort): " %5.2f (_csort_time_permute / _csort_time_sort)

    * Also test with LSD radix
    di as text ""
    di as text "Sorting with LSD radix:"
    csort x2, verbose alg(lsd)

    di as text ""
    di as text "Sort time:    " %8.4f _csort_time_sort " sec"
    di as text "Permute time: " %8.4f _csort_time_permute " sec"
    di as text "Ratio (permute/sort): " %5.2f (_csort_time_permute / _csort_time_sort)
}

di as text ""
di as text "========================================"
di as text "Benchmark complete!"
di as text "========================================"
