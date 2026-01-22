*! Test large m:1 merge - compare merge vs cmerge

clear all
set more off
adopath ++ "build"

timer clear

* === Create using dataset: 1M unique ids ===
di "Creating using dataset with 1M unique ids..."
clear
set obs 1000000
gen long id = _n
gen double using_val = runiform()
tempfile using_data
save `using_data'

* === Create master dataset: 25M observations ===
di "Creating master dataset with 25M observations..."
clear
set obs 25000000
gen long id = mod(_n - 1, 1000000) + 1  // ids 1 to 1M, each appears 25 times
gen double master_val = runiform()
tempfile master_data
save `master_data'

* === Run Stata's merge ===
di _n "Running Stata merge..."
use `master_data', clear
timer on 1
merge m:1 id using `using_data'
timer off 1

* Store results
local stata_N = _N
local stata_N1 = r(N_1)
local stata_N2 = r(N_2)
local stata_N3 = r(N_3)
local stata_N4 = r(N_4)
local stata_N5 = r(N_5)

* Save merged data for comparison
tempfile stata_merged
save `stata_merged'

di "Stata merge results:"
di "  N = `stata_N'"
di "  N_1 (master only) = `stata_N1'"
di "  N_2 (using only) = `stata_N2'"
di "  N_3 (matched) = `stata_N3'"
di "  N_4 (missing updated) = `stata_N4'"
di "  N_5 (nonmissing conflict) = `stata_N5'"
tab _merge

* === Run cmerge ===
di _n "Running cmerge..."
use `master_data', clear
timer on 2
cmerge m:1 id using `using_data'
timer off 2

* Store results
local cmerge_N = _N
local cmerge_N1 = r(N_1)
local cmerge_N2 = r(N_2)
local cmerge_N3 = r(N_3)
local cmerge_N4 = r(N_4)
local cmerge_N5 = r(N_5)

di "cmerge results:"
di "  N = `cmerge_N'"
di "  N_1 (master only) = `cmerge_N1'"
di "  N_2 (using only) = `cmerge_N2'"
di "  N_3 (matched) = `cmerge_N3'"
di "  N_4 (missing updated) = `cmerge_N4'"
di "  N_5 (nonmissing conflict) = `cmerge_N5'"
tab _merge

* Save cmerge data
tempfile cmerge_merged
save `cmerge_merged'

* === Compare results ===
di _n "=== COMPARISON ==="
di "N:   Stata=`stata_N' vs cmerge=`cmerge_N'"
di "N_1: Stata=`stata_N1' vs cmerge=`cmerge_N1'"
di "N_2: Stata=`stata_N2' vs cmerge=`cmerge_N2'"
di "N_3: Stata=`stata_N3' vs cmerge=`cmerge_N3'"

local match = (`stata_N' == `cmerge_N') & (`stata_N1' == `cmerge_N1') & (`stata_N2' == `cmerge_N2') & (`stata_N3' == `cmerge_N3')
if `match' {
    di _n "COUNTS MATCH"
}
else {
    di _n "COUNTS DO NOT MATCH!"
}

* === Verify match counts ===
di _n "=== VERIFICATION ==="
local test_pass = (`stata_N' == `cmerge_N') & (`cmerge_N1' == 0) & (`cmerge_N3' == 25000000)
if `test_pass' {
    di "TEST PASSED: cmerge matches Stata merge results"
    di "  Both: N = 25,000,000"
    di "  Both: Matched = 25,000,000"
    di "  Both: Master only = 0"
}
else {
    di "TEST FAILED: Results differ"
    di "  Stata:  N=`stata_N', matched=25000000"
    di "  cmerge: N=`cmerge_N', matched=`cmerge_N3', master_only=`cmerge_N1'"
}

* === Timing ===
di _n "=== TIMING ==="
timer list 1
timer list 2
