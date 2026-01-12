* Comprehensive test of ALL variable creation methods
clear all
set more off

local classpath "/Users/Mike/Documents/GitHub/stata-ctools"

* Create dataset WITH data - 100M obs, 20 existing vars
di "Creating test dataset..."
set obs 100000000
forval i = 1/20 {
    qui gen x`i' = rnormal()
}
gen g = ceil(runiform()*100000)

di "Dataset has " c(k) " variables, " _N " observations"
di ""

* ============================================================================
* Test 1: Stata gen (sequential)
* ============================================================================
di "Test 1: Stata gen sequential (24 vars)..."
timer clear 1
timer on 1
qui gen double v0 = .
qui gen double v1 = .
qui gen double v2 = .
qui gen double v3 = .
qui gen double v4 = .
qui gen double v5 = .
qui gen double v6 = .
qui gen double v7 = .
qui gen double v8 = .
qui gen double v9 = .
qui gen double v10 = .
qui gen double v11 = .
qui gen double v12 = .
qui gen double v13 = .
qui gen double v14 = .
qui gen double v15 = .
qui gen double v16 = .
qui gen double v17 = .
qui gen double v18 = .
qui gen double v19 = .
qui gen double v20 = .
qui gen double v21 = .
qui gen double v22 = .
qui gen double v23 = .
timer off 1
timer list 1
drop v0-v23

* ============================================================================
* Test 2: Mata st_addvar (with initialization)
* ============================================================================
di _n "Test 2: Mata st_addvar (24 vars, initializes to missing)..."
timer clear 2
timer on 2
mata: st_addvar(("double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double"), ("v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23"))
timer off 2
timer list 2
drop v0-v23

* ============================================================================
* Test 3: Mata _st_addvar (without initialization)
* ============================================================================
di _n "Test 3: Mata _st_addvar (24 vars, no initialization)..."
timer clear 3
timer on 3
mata: _st_addvar(("double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double"), ("v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23"))
timer off 3
timer list 3
drop v0-v23

* ============================================================================
* Test 4: Java AddVarBatch
* ============================================================================
di _n "Test 4: Java AddVarBatch (24 vars)..."
local types "double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double"
local names "v0+v1+v2+v3+v4+v5+v6+v7+v8+v9+v10+v11+v12+v13+v14+v15+v16+v17+v18+v19+v20+v21+v22+v23"
timer clear 4
timer on 4
javacall AddVarBatch create, classpath(`classpath') args(`types' `names')
timer off 4
timer list 4
drop v0-v23

* ============================================================================
* Test 5: Mata loop calling st_addvar one at a time
* ============================================================================
di _n "Test 5: Mata loop st_addvar (24 individual calls)..."
timer clear 5
timer on 5
forval i = 0/23 {
    mata: st_addvar("double", "v`i'")
}
timer off 5
timer list 5
drop v0-v23

* ============================================================================
* Test 6: Append empty dataset trick
* ============================================================================
di _n "Test 6: Append empty dataset with vars (0 obs)..."

* First create the empty template dataset
preserve
clear
gen double v0 = .
gen double v1 = .
gen double v2 = .
gen double v3 = .
gen double v4 = .
gen double v5 = .
gen double v6 = .
gen double v7 = .
gen double v8 = .
gen double v9 = .
gen double v10 = .
gen double v11 = .
gen double v12 = .
gen double v13 = .
gen double v14 = .
gen double v15 = .
gen double v16 = .
gen double v17 = .
gen double v18 = .
gen double v19 = .
gen double v20 = .
gen double v21 = .
gen double v22 = .
gen double v23 = .
drop if 1  // Make it 0 obs
tempfile empty_vars
save `empty_vars'
restore

timer clear 6
timer on 6
append using `empty_vars'
timer off 6
timer list 6
drop v0-v23

* ============================================================================
* Test 7: Append empty dataset - measure including creation
* ============================================================================
di _n "Test 7: Append empty dataset (including tempfile creation)..."
timer clear 7
timer on 7
preserve
clear
forval i = 0/23 {
    gen double v`i' = .
}
drop if 1
tempfile empty_vars2
save `empty_vars2'
restore
append using `empty_vars2'
timer off 7
timer list 7
drop v0-v23

* ============================================================================
* Summary
* ============================================================================
di _n "========================================"
di "SUMMARY"
di "========================================"
timer list 1
timer list 2
timer list 3
timer list 4
timer list 5
timer list 6
timer list 7
