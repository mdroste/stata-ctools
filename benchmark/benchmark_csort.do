* Benchmark: csort vs sort with 100 million observations, 5 double variables
* Single key sort on var1

clear all
set more off

* Set up timing
timer clear

display as text ""
display as text "==========================================="
display as text "Benchmark: 100M obs, 5 double vars, 1 key"
display as text "==========================================="
display as text ""

* Generate test data with 5 double variables from standard normal
display as text "Generating 100 million observations with 5 double variables..."
set seed 42
set obs 100000000

gen double var1 = rnormal()
gen double var2 = rnormal()
gen double var3 = rnormal()
gen double var4 = rnormal()
gen double var5 = rnormal()

display as text "Data generated. Starting benchmarks..."
display as text ""

* ---------------------------------------------
* Benchmark 1: Built-in Stata sort
* ---------------------------------------------
display as text "------------------------------------------"
display as text "Test 1: Built-in Stata sort (on var1)"
display as text "------------------------------------------"

* Randomize the data first (use fixed seed for reproducibility)
set seed 123
gen double sort_order = runiform()
sort sort_order
drop sort_order

* Time the sort
timer on 1
sort var1
timer off 1

* Save the Stata-sorted result for comparison
tempfile stata_sorted
quietly save `stata_sorted'

display as text "Stata sort completed."
display as text ""

* ---------------------------------------------
* Benchmark 2: csort (radix sort)
* ---------------------------------------------
display as text "------------------------------------------"
display as text "Test 2: csort (radix sort on var1)"
display as text "------------------------------------------"

* Regenerate data with same seeds
clear
set seed 42
set obs 100000000

gen double var1 = rnormal()
gen double var2 = rnormal()
gen double var3 = rnormal()
gen double var4 = rnormal()
gen double var5 = rnormal()

* Randomize with same seed as before
set seed 123
gen double sort_order = runiform()
sort sort_order
drop sort_order

* Add adopath for csort (need to be in project root for plugin path)
cd ..
adopath + "."

* Time csort with verbose output
timer on 2
csort var1, verbose
timer off 2

display as text ""
display as text "csort completed."
display as text ""

* ---------------------------------------------
* Verification: Compare datasets
* ---------------------------------------------
display as text "------------------------------------------"
display as text "Verification: Comparing sorted datasets"
display as text "------------------------------------------"

capture cf _all using `stata_sorted'
if _rc == 0 {
    display as text "PASS: Datasets are identical"
    local verification = "PASS"
}
else {
    display as error "FAIL: Datasets differ"
    local verification = "FAIL"
}

display as text ""

* ---------------------------------------------
* Summary
* ---------------------------------------------
display as text "==========================================="
display as text "SUMMARY"
display as text "==========================================="
display as text ""

timer list 1
timer list 2

display as text ""
display as text "Stata sort time:  " as result %9.3f r(t1) as text " seconds"
display as text "csort time:       " as result %9.3f r(t2) as text " seconds"
display as text ""

if r(t2) < r(t1) {
    local speedup = r(t1) / r(t2)
    display as text "csort is " as result %4.2f `speedup' as text "x faster than Stata sort"
}
else {
    local slowdown = r(t2) / r(t1)
    display as text "csort is " as result %4.2f `slowdown' as text "x slower than Stata sort"
}

display as text ""
display as text "Verification: `verification' - Datasets are identical"
display as text ""
