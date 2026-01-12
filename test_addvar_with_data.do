* Test addvar performance with existing data
clear all
set more off

local classpath "/Users/Mike/Documents/GitHub/stata-ctools"

* Create dataset WITH data (like cmerge scenario) - 100M obs, 20 existing vars
set obs 100000000
forval i = 1/20 {
    gen x`i' = rnormal()
}
gen g = ceil(runiform()*100000)

di "Dataset has " c(k) " variables, " _N " observations"

* 24 new variables to create
local types "long+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+double+byte"
local names "v0+v1+v2+v3+v4+v5+v6+v7+v8+v9+v10+v11+v12+v13+v14+v15+v16+v17+v18+v19+v20+v21+v22+v23"

* Test 1: Mata _st_addvar
di _n "Test 1: Mata _st_addvar (24 vars)..."
timer clear 1
timer on 1
mata: _st_addvar(("long","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","double","byte"), ("v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23"))
timer off 1
timer list 1
drop v0-v23

* Test 2: Java batch (sequential)
di _n "Test 2: Java AddVarBatch sequential (24 vars)..."
timer clear 2
timer on 2
javacall AddVarBatch create, classpath(`classpath') args(`types' `names')
timer off 2
timer list 2
drop v0-v23

* Test 3: Java parallel
di _n "Test 3: Java AddVarParallel (24 vars, " c(processors_max) " threads)..."
timer clear 3
timer on 3
javacall AddVarParallel create, classpath(`classpath') args(`types' `names')
timer off 3
timer list 3
drop v0-v23

di _n "Done!"
