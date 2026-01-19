* Test csort split timing (sort vs permutation)
* This tests the new timing breakdown showing sort computation vs permutation application

clear all
set more off

* Add build directory to adopath
adopath + "build"

* Generate test data
set seed 12345
set obs 100000

* Generate random data
gen id = _n
gen x = runiform()
gen y = rnormal()
gen str20 name = "name" + string(_n)

* Test 1: Default IPS4O algorithm with verbose timing
di as text ""
di as text "========================================"
di as text "Test 1: IPS4O (default) with verbose"
di as text "========================================"
csort x, verbose

* Check that timing scalars exist
di as text ""
di as text "Timing scalars:"
di as text "  _csort_time_load = " _csort_time_load
di as text "  _csort_time_sort = " _csort_time_sort
di as text "  _csort_time_permute = " _csort_time_permute
di as text "  _csort_time_store = " _csort_time_store
di as text "  _csort_time_cleanup = " _csort_time_cleanup
di as text "  _csort_time_total = " _csort_time_total

* Verify permutation time is > 0 for IPS4O
assert _csort_time_permute > 0

* Test 2: LSD radix sort with verbose timing
di as text ""
di as text "========================================"
di as text "Test 2: LSD radix sort with verbose"
di as text "========================================"
csort y, verbose alg(lsd)

* Check that timing scalars exist
di as text ""
di as text "Timing scalars:"
di as text "  _csort_time_load = " _csort_time_load
di as text "  _csort_time_sort = " _csort_time_sort
di as text "  _csort_time_permute = " _csort_time_permute
di as text "  _csort_time_store = " _csort_time_store
di as text "  _csort_time_cleanup = " _csort_time_cleanup
di as text "  _csort_time_total = " _csort_time_total

* Verify permutation time is > 0 for LSD
assert _csort_time_permute > 0

* Test 3: Merge sort (uses combined timing - permutation time should be 0)
di as text ""
di as text "========================================"
di as text "Test 3: Merge sort with verbose"
di as text "========================================"
csort id, verbose alg(merge)

* Check that timing scalars exist
di as text ""
di as text "Timing scalars:"
di as text "  _csort_time_load = " _csort_time_load
di as text "  _csort_time_sort = " _csort_time_sort
di as text "  _csort_time_permute = " _csort_time_permute
di as text "  _csort_time_store = " _csort_time_store
di as text "  _csort_time_cleanup = " _csort_time_cleanup
di as text "  _csort_time_total = " _csort_time_total

* For merge sort, permutation time should be 0 (combined with sort)
assert _csort_time_permute == 0

di as text ""
di as text "========================================"
di as text "All tests passed!"
di as text "========================================"
