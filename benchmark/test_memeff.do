/*
    test_memeff.do
    Test and benchmark memory-efficient sorting mode

    Compares standard csort vs memory-efficient (memeff) mode.
    Memory-efficient mode is best for:
    - Wide datasets (many columns)
    - Limited memory situations
    - Large datasets that don't fit in RAM with standard approach
*/

clear all
set more off

* Add ctools build directory to adopath
* Use the directory where this script is located
local script_dir = c(pwd)
adopath + "`script_dir'/../build"
di "Added to adopath: `script_dir'/../build"

* Create test dataset with many columns (simulates wide dataset)
local nobs = 500000
local ncols = 100

di as text ""
di as text "=============================================="
di as text "Memory-Efficient Sort Test"
di as text "=============================================="
di as text "Observations: `nobs'"
di as text "Columns: `ncols' (1 sort key + `=`ncols'-1' data columns)"
di as text ""

* Generate sort key and many data columns
set seed 12345
set obs `nobs'

* Sort key (random integers)
gen sort_key = ceil(runiform() * 10000)

* Generate many data columns to simulate wide dataset
forvalues i = 1/`=`ncols'-1' {
    gen col`i' = runiform()
}

* Save original for verification
tempfile original
save "`original'"

* ============================================
* Test 1: Standard sort
* ============================================
di as text "----------------------------------------------"
di as text "Test 1: Standard csort"
di as text "----------------------------------------------"

use "`original'", clear

timer clear 1
timer on 1
csort sort_key, timeit
timer off 1

quietly timer list 1
local time_standard = r(t1)
di as text "Total time: " as result %8.4f `time_standard' " sec"

* Save sorted data for comparison
tempfile sorted_standard
save "`sorted_standard'"

* ============================================
* Test 2: Memory-efficient sort
* ============================================
di as text ""
di as text "----------------------------------------------"
di as text "Test 2: Memory-efficient csort (memeff)"
di as text "----------------------------------------------"

use "`original'", clear

timer clear 2
timer on 2
csort sort_key, timeit memeff
timer off 2

quietly timer list 2
local time_memeff = r(t2)
di as text "Total time: " as result %8.4f `time_memeff' " sec"

* Save sorted data for comparison
tempfile sorted_memeff
save "`sorted_memeff'"

* ============================================
* Verify correctness
* ============================================
di as text ""
di as text "----------------------------------------------"
di as text "Verification"
di as text "----------------------------------------------"

* Compare the two sorted datasets
use "`sorted_standard'", clear
tempfile std
save "`std'"

use "`sorted_memeff'", clear
cf _all using "`std'"

if _rc == 0 {
    di as result "PASS: Both methods produce identical results"
}
else {
    di as error "FAIL: Results differ between methods"
}

* ============================================
* Summary
* ============================================
di as text ""
di as text "=============================================="
di as text "Summary"
di as text "=============================================="
di as text "Standard sort:  " as result %8.4f `time_standard' " sec"
di as text "Memeff sort:    " as result %8.4f `time_memeff' " sec"
local speedup = `time_standard' / `time_memeff'
if `speedup' >= 1 {
    di as text "Memeff is " as result %5.2f `speedup' "x" as text " faster"
}
else {
    local slowdown = 1 / `speedup'
    di as text "Memeff is " as result %5.2f `slowdown' "x" as text " slower"
}
di as text ""
di as text "Note: Memory-efficient mode trades some speed for"
di as text "dramatically lower memory usage. Best for wide datasets"
di as text "or when memory is constrained."
di as text "=============================================="
