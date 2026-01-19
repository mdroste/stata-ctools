* Test csort split timing (sort vs permutation) for ALL algorithms
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

* Test 1: IPS4O (default)
di as text ""
di as text "========================================"
di as text "Test 1: IPS4O (default)"
di as text "========================================"
csort x, verbose
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

* Test 2: LSD radix sort
di as text ""
di as text "========================================"
di as text "Test 2: LSD radix sort"
di as text "========================================"
csort y, verbose alg(lsd)
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

* Test 3: MSD radix sort
di as text ""
di as text "========================================"
di as text "Test 3: MSD radix sort"
di as text "========================================"
csort x, verbose alg(msd)
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

* Test 4: Timsort
di as text ""
di as text "========================================"
di as text "Test 4: Timsort"
di as text "========================================"
csort y, verbose alg(timsort)
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

* Test 5: Sample sort
di as text ""
di as text "========================================"
di as text "Test 5: Sample sort"
di as text "========================================"
csort x, verbose alg(sample)
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

* Test 6: Merge sort
di as text ""
di as text "========================================"
di as text "Test 6: Merge sort"
di as text "========================================"
csort id, verbose alg(merge)
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

* Test 7: Counting sort (integer data)
di as text ""
di as text "========================================"
di as text "Test 7: Counting sort"
di as text "========================================"
csort id, verbose alg(counting)
assert _csort_time_permute > 0
di as text "  PASS: permutation time = " %8.4f _csort_time_permute

di as text ""
di as text "========================================"
di as text "All 7 algorithm tests passed!"
di as text "========================================"
