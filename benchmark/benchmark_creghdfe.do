********************************************************************************
* benchmark_creghdfe.do
* Speed benchmark comparing reghdfe and creghdfe performance
* 50M observations, 2 fixed effects
********************************************************************************

clear all
set more off

* Add ctools build directory to adopath
adopath ++ "./build"

********************************************************************************
* Configuration
********************************************************************************

local N = 50000000   // 50 million observations

di as txt _n _n
di as txt "================================================================================"
di as txt "                    REGHDFE vs CREGHDFE SPEED BENCHMARK"
di as txt "================================================================================"
di as txt "Date: `c(current_date)' `c(current_time)'"
di as txt "Stata: `c(stata_version)' `c(processors)'-core"
di as txt "Observations: " %12.0fc `N'
di as txt "================================================================================"

********************************************************************************
* Generate simulated data
********************************************************************************

di as txt _n "{hline 80}"
di as txt "Generating simulated data with N = " %12.0fc `N'
di as txt "{hline 80}"

set seed 12345
set obs `N'

* Create outcome and regressors
gen double y = rnormal()
gen double x1 = rnormal()
gen double x2 = rnormal()
gen double x3 = rnormal()
gen double x4 = rnormal()
gen double x5 = rnormal()

* Create fixed effect variables
* FE1: ~50,000 groups (like firm)
local n_fe1 = max(50000, floor(`N'/1000))
gen long fe1 = ceil(runiform() * `n_fe1')

* FE2: ~100 groups (like year)
gen int fe2 = ceil(runiform() * 100)

* Create cluster variable (same as fe1 for clustered tests)
gen long cluster_id = fe1

* Add FE effects to outcome
tempvar fe1_effect fe2_effect
bysort fe1: gen double `fe1_effect' = rnormal() if _n == 1
bysort fe1: replace `fe1_effect' = `fe1_effect'[1]
bysort fe2: gen double `fe2_effect' = rnormal() if _n == 1
bysort fe2: replace `fe2_effect' = `fe2_effect'[1]
replace y = y + 0.5*x1 + 0.3*x2 - 0.2*x3 + 0.1*x4 - 0.4*x5 + `fe1_effect' + `fe2_effect'

compress

di as txt "Data generation complete."
di as txt "  Observations: " %12.0fc _N
di as txt "  FE1 groups:   " %12.0fc `n_fe1'
di as txt "  FE2 groups:   100"

********************************************************************************
* Initialize results storage
********************************************************************************

* Store results: test_name, time_reghdfe, time_creghdfe, speedup
tempname results
matrix `results' = J(4, 4, .)
matrix colnames `results' = test time_reghdfe time_creghdfe speedup
local row = 0

********************************************************************************
* Test 1: 2 FEs, unadjusted VCE
********************************************************************************

local ++row
di as txt _n "{hline 80}"
di as txt "Test 1: 2 FEs (fe1 fe2), unadjusted VCE"
di as txt "{hline 80}"

timer clear 1
timer clear 2

* Run reghdfe
di as txt "Running reghdfe..."
timer on 1
qui reghdfe y x1 x2 x3 x4 x5, absorb(fe1 fe2)
timer off 1

* Store reghdfe results for comparison
tempname b_reghdfe V_reghdfe
matrix `b_reghdfe' = e(b)
matrix `V_reghdfe' = e(V)
local N_reghdfe = e(N)
local r2_reghdfe = e(r2)
local rmse_reghdfe = e(rmse)

* Run creghdfe
di as txt "Running creghdfe..."
timer on 2
qui creghdfe y x1 x2 x3 x4 x5, absorb(fe1 fe2)
timer off 2

* Store creghdfe results for comparison
tempname b_creghdfe V_creghdfe
matrix `b_creghdfe' = e(b)
matrix `V_creghdfe' = e(V)
local N_creghdfe = e(N)
local r2_creghdfe = e(r2)
local rmse_creghdfe = e(rmse)

qui timer list
local t1 = r(t1)
local t2 = r(t2)
local speedup = `t1' / `t2'

matrix `results'[`row', 1] = 1
matrix `results'[`row', 2] = `t1'
matrix `results'[`row', 3] = `t2'
matrix `results'[`row', 4] = `speedup'

di as txt _n "Timing:"
di as txt "  reghdfe:   " %10.2f `t1' " sec"
di as txt "  creghdfe:  " %10.2f `t2' " sec"
di as txt "  Speedup:   " %10.2f `speedup' "x"

di as txt _n "Output Comparison:"
di as txt "  N:         reghdfe=" %12.0fc `N_reghdfe' "  creghdfe=" %12.0fc `N_creghdfe' cond(`N_reghdfe'==`N_creghdfe', "  [MATCH]", "  [DIFFER]")
di as txt "  R-squared: reghdfe=" %12.6f `r2_reghdfe' "  creghdfe=" %12.6f `r2_creghdfe' cond(abs(`r2_reghdfe'-`r2_creghdfe')<1e-6, "  [MATCH]", "  [DIFFER]")
di as txt "  RMSE:      reghdfe=" %12.6f `rmse_reghdfe' "  creghdfe=" %12.6f `rmse_creghdfe' cond(abs(`rmse_reghdfe'-`rmse_creghdfe')<1e-6, "  [MATCH]", "  [DIFFER]")

di as txt _n "Coefficient Comparison:"
local coef_match = 1
forval i = 1/5 {
    local b_hdfe = `b_reghdfe'[1, `i']
    local b_chdfe = `b_creghdfe'[1, `i']
    local se_hdfe = sqrt(`V_reghdfe'[`i', `i'])
    local se_chdfe = sqrt(`V_creghdfe'[`i', `i'])
    local b_diff = abs(`b_hdfe' - `b_chdfe')
    local se_diff = abs(`se_hdfe' - `se_chdfe')
    if `b_diff' > 1e-6 | `se_diff' > 1e-6 {
        local coef_match = 0
    }
    di as txt "  x`i': b=" %10.6f `b_hdfe' " vs " %10.6f `b_chdfe' "  se=" %10.6f `se_hdfe' " vs " %10.6f `se_chdfe'
}
di as txt "  Coefficients/SEs: " cond(`coef_match'==1, "[ALL MATCH]", "[DIFFERENCES DETECTED]")

********************************************************************************
* Test 2: 2 FEs, clustered VCE
********************************************************************************

local ++row
di as txt _n "{hline 80}"
di as txt "Test 2: 2 FEs (fe1 fe2), clustered VCE (cluster_id)"
di as txt "{hline 80}"

timer clear 1
timer clear 2

* Run reghdfe
di as txt "Running reghdfe..."
timer on 1
qui reghdfe y x1 x2 x3 x4 x5, absorb(fe1 fe2) cluster(cluster_id)
timer off 1

* Store reghdfe results
matrix `b_reghdfe' = e(b)
matrix `V_reghdfe' = e(V)
local N_reghdfe = e(N)
local r2_reghdfe = e(r2)

* Run creghdfe
di as txt "Running creghdfe..."
timer on 2
qui creghdfe y x1 x2 x3 x4 x5, absorb(fe1 fe2) vce(cluster cluster_id)
timer off 2

* Store creghdfe results
matrix `b_creghdfe' = e(b)
matrix `V_creghdfe' = e(V)
local N_creghdfe = e(N)
local r2_creghdfe = e(r2)

qui timer list
local t1 = r(t1)
local t2 = r(t2)
local speedup = `t1' / `t2'

matrix `results'[`row', 1] = 2
matrix `results'[`row', 2] = `t1'
matrix `results'[`row', 3] = `t2'
matrix `results'[`row', 4] = `speedup'

di as txt _n "Timing:"
di as txt "  reghdfe:   " %10.2f `t1' " sec"
di as txt "  creghdfe:  " %10.2f `t2' " sec"
di as txt "  Speedup:   " %10.2f `speedup' "x"

di as txt _n "Output Comparison:"
di as txt "  N:         reghdfe=" %12.0fc `N_reghdfe' "  creghdfe=" %12.0fc `N_creghdfe' cond(`N_reghdfe'==`N_creghdfe', "  [MATCH]", "  [DIFFER]")
di as txt "  R-squared: reghdfe=" %12.6f `r2_reghdfe' "  creghdfe=" %12.6f `r2_creghdfe' cond(abs(`r2_reghdfe'-`r2_creghdfe')<1e-6, "  [MATCH]", "  [DIFFER]")

di as txt _n "Coefficient Comparison (with clustered SEs):"
local coef_match = 1
forval i = 1/5 {
    local b_hdfe = `b_reghdfe'[1, `i']
    local b_chdfe = `b_creghdfe'[1, `i']
    local se_hdfe = sqrt(`V_reghdfe'[`i', `i'])
    local se_chdfe = sqrt(`V_creghdfe'[`i', `i'])
    local b_diff = abs(`b_hdfe' - `b_chdfe')
    local se_diff = abs(`se_hdfe' - `se_chdfe')
    if `b_diff' > 1e-6 | `se_diff' > 1e-6 {
        local coef_match = 0
    }
    di as txt "  x`i': b=" %10.6f `b_hdfe' " vs " %10.6f `b_chdfe' "  se=" %10.6f `se_hdfe' " vs " %10.6f `se_chdfe'
}
di as txt "  Coefficients/SEs: " cond(`coef_match'==1, "[ALL MATCH]", "[DIFFERENCES DETECTED]")

********************************************************************************
* Test 3: 1 FE only, unadjusted VCE
********************************************************************************

local ++row
di as txt _n "{hline 80}"
di as txt "Test 3: 1 FE (fe1 only), unadjusted VCE"
di as txt "{hline 80}"

timer clear 1
timer clear 2

* Run reghdfe
di as txt "Running reghdfe..."
timer on 1
qui reghdfe y x1 x2 x3 x4 x5, absorb(fe1)
timer off 1

matrix `b_reghdfe' = e(b)
matrix `V_reghdfe' = e(V)
local N_reghdfe = e(N)
local r2_reghdfe = e(r2)

* Run creghdfe
di as txt "Running creghdfe..."
timer on 2
qui creghdfe y x1 x2 x3 x4 x5, absorb(fe1)
timer off 2

matrix `b_creghdfe' = e(b)
matrix `V_creghdfe' = e(V)
local N_creghdfe = e(N)
local r2_creghdfe = e(r2)

qui timer list
local t1 = r(t1)
local t2 = r(t2)
local speedup = `t1' / `t2'

matrix `results'[`row', 1] = 3
matrix `results'[`row', 2] = `t1'
matrix `results'[`row', 3] = `t2'
matrix `results'[`row', 4] = `speedup'

di as txt _n "Timing:"
di as txt "  reghdfe:   " %10.2f `t1' " sec"
di as txt "  creghdfe:  " %10.2f `t2' " sec"
di as txt "  Speedup:   " %10.2f `speedup' "x"

di as txt _n "Output Comparison:"
di as txt "  N:         reghdfe=" %12.0fc `N_reghdfe' "  creghdfe=" %12.0fc `N_creghdfe' cond(`N_reghdfe'==`N_creghdfe', "  [MATCH]", "  [DIFFER]")
di as txt "  R-squared: reghdfe=" %12.6f `r2_reghdfe' "  creghdfe=" %12.6f `r2_creghdfe' cond(abs(`r2_reghdfe'-`r2_creghdfe')<1e-6, "  [MATCH]", "  [DIFFER]")

di as txt _n "Coefficient Comparison:"
local coef_match = 1
forval i = 1/5 {
    local b_hdfe = `b_reghdfe'[1, `i']
    local b_chdfe = `b_creghdfe'[1, `i']
    local se_hdfe = sqrt(`V_reghdfe'[`i', `i'])
    local se_chdfe = sqrt(`V_creghdfe'[`i', `i'])
    local b_diff = abs(`b_hdfe' - `b_chdfe')
    local se_diff = abs(`se_hdfe' - `se_chdfe')
    if `b_diff' > 1e-6 | `se_diff' > 1e-6 {
        local coef_match = 0
    }
    di as txt "  x`i': b=" %10.6f `b_hdfe' " vs " %10.6f `b_chdfe' "  se=" %10.6f `se_hdfe' " vs " %10.6f `se_chdfe'
}
di as txt "  Coefficients/SEs: " cond(`coef_match'==1, "[ALL MATCH]", "[DIFFERENCES DETECTED]")

********************************************************************************
* Test 4: 1 FE, clustered VCE
********************************************************************************

local ++row
di as txt _n "{hline 80}"
di as txt "Test 4: 1 FE (fe1 only), clustered VCE (cluster_id)"
di as txt "{hline 80}"

timer clear 1
timer clear 2

* Run reghdfe
di as txt "Running reghdfe..."
timer on 1
qui reghdfe y x1 x2 x3 x4 x5, absorb(fe1) cluster(cluster_id)
timer off 1

matrix `b_reghdfe' = e(b)
matrix `V_reghdfe' = e(V)
local N_reghdfe = e(N)
local r2_reghdfe = e(r2)

* Run creghdfe
di as txt "Running creghdfe..."
timer on 2
qui creghdfe y x1 x2 x3 x4 x5, absorb(fe1) vce(cluster cluster_id)
timer off 2

matrix `b_creghdfe' = e(b)
matrix `V_creghdfe' = e(V)
local N_creghdfe = e(N)
local r2_creghdfe = e(r2)

qui timer list
local t1 = r(t1)
local t2 = r(t2)
local speedup = `t1' / `t2'

matrix `results'[`row', 1] = 4
matrix `results'[`row', 2] = `t1'
matrix `results'[`row', 3] = `t2'
matrix `results'[`row', 4] = `speedup'

di as txt _n "Timing:"
di as txt "  reghdfe:   " %10.2f `t1' " sec"
di as txt "  creghdfe:  " %10.2f `t2' " sec"
di as txt "  Speedup:   " %10.2f `speedup' "x"

di as txt _n "Output Comparison:"
di as txt "  N:         reghdfe=" %12.0fc `N_reghdfe' "  creghdfe=" %12.0fc `N_creghdfe' cond(`N_reghdfe'==`N_creghdfe', "  [MATCH]", "  [DIFFER]")
di as txt "  R-squared: reghdfe=" %12.6f `r2_reghdfe' "  creghdfe=" %12.6f `r2_creghdfe' cond(abs(`r2_reghdfe'-`r2_creghdfe')<1e-6, "  [MATCH]", "  [DIFFER]")

di as txt _n "Coefficient Comparison (with clustered SEs):"
local coef_match = 1
forval i = 1/5 {
    local b_hdfe = `b_reghdfe'[1, `i']
    local b_chdfe = `b_creghdfe'[1, `i']
    local se_hdfe = sqrt(`V_reghdfe'[`i', `i'])
    local se_chdfe = sqrt(`V_creghdfe'[`i', `i'])
    local b_diff = abs(`b_hdfe' - `b_chdfe')
    local se_diff = abs(`se_hdfe' - `se_chdfe')
    if `b_diff' > 1e-6 | `se_diff' > 1e-6 {
        local coef_match = 0
    }
    di as txt "  x`i': b=" %10.6f `b_hdfe' " vs " %10.6f `b_chdfe' "  se=" %10.6f `se_hdfe' " vs " %10.6f `se_chdfe'
}
di as txt "  Coefficients/SEs: " cond(`coef_match'==1, "[ALL MATCH]", "[DIFFERENCES DETECTED]")

********************************************************************************
* Summary
********************************************************************************

di as txt _n _n
di as txt "================================================================================"
di as txt "                           SPEED BENCHMARK SUMMARY"
di as txt "================================================================================"
di as txt _n

di as txt "{hline 76}"
di as txt %30s "Test" " | " %12s "reghdfe" " | " %12s "creghdfe" " | " %10s "Speedup"
di as txt "{hline 76}"

local test_names `""2 FEs, unadjusted" "2 FEs, clustered" "1 FE, unadjusted" "1 FE, clustered""'

forval i = 1/`row' {
    local t_hdfe = `results'[`i', 2]
    local t_chdfe = `results'[`i', 3]
    local speedup = `results'[`i', 4]
    local tname : word `i' of `test_names'

    di as txt %30s "`tname'" " | " as result %10.2f `t_hdfe' "s" as txt " | " as result %10.2f `t_chdfe' "s" as txt " | " as result %9.2f `speedup' "x"
}

di as txt "{hline 76}"

* Calculate average speedup
mata: st_local("avg_speedup", strofreal(mean(st_matrix("`results'")[,4])))
di as txt _n "Average speedup: " as result %6.2f `avg_speedup' "x"

di as txt _n "================================================================================"
di as txt "Note: Speedup = reghdfe time / creghdfe time (higher is better for creghdfe)"
di as txt "      All output comparisons should show [MATCH] for identical results"
di as txt "================================================================================"
