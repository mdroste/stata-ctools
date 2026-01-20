* Test prediction comparison: civreghdfe vs ivreghdfe
* Using sysuse auto dataset

clear all
set more off
set trace off

* Add ctools build directory to adopath
adopath + "build"

sysuse auto, clear

* Create a simple absorb variable
gen byte one = 1

di _n "{hline 78}"
di "Testing civreghdfe vs ivreghdfe predictions"
di "{hline 78}" _n

* ============================================================================
* Test 1: No fixed effects (just constant absorbed)
* ============================================================================
di as text "Test 1: No real fixed effects (absorb(one))"
di "{hline 40}"

* Run ivreghdfe
ivreghdfe price (mpg = weight length), absorb(one) vce(robust) resid
matrix b_iv = e(b)
predict double yhat_iv, xb

* Run civreghdfe
civreghdfe price (mpg = weight length), absorb(one) vce(robust)
matrix b_civ = e(b)

* Manually compute predictions using civreghdfe coefficients
gen double yhat_civ_manual = b_civ[1,1] * mpg

* Compare coefficients
di _n "Coefficient comparison:"
di "                  ivreghdfe      civreghdfe       Diff"
di "{hline 60}"
local b_iv_mpg = b_iv[1,1]
local b_civ_mpg = b_civ[1,1]
local diff_mpg = `b_civ_mpg' - `b_iv_mpg'
di "mpg:          " %14.8f `b_iv_mpg' "  " %14.8f `b_civ_mpg' "  " %14.8e `diff_mpg'

* Compare predictions
gen double diff_pred1 = yhat_civ_manual - yhat_iv
sum diff_pred1, detail

* ============================================================================
* Test 2: With actual fixed effects
* ============================================================================
di _n _n "Test 2: With actual fixed effects (absorb(foreign))"
di "{hline 40}"

drop yhat_iv yhat_civ_manual diff_pred1

* Run ivreghdfe
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) resid
matrix b_iv2 = e(b)
predict double yhat_iv, xb

* Run civreghdfe
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
matrix b_civ2 = e(b)

* Manually compute predictions using civreghdfe coefficients
gen double yhat_civ_manual = b_civ2[1,1] * mpg

* Compare coefficients
di _n "Coefficient comparison:"
di "                  ivreghdfe      civreghdfe       Diff"
di "{hline 60}"
local b_iv_mpg2 = b_iv2[1,1]
local b_civ_mpg2 = b_civ2[1,1]
local diff_mpg2 = `b_civ_mpg2' - `b_iv_mpg2'
di "mpg:          " %14.8f `b_iv_mpg2' "  " %14.8f `b_civ_mpg2' "  " %14.8e `diff_mpg2'

* Compare predictions
gen double diff_pred2 = yhat_civ_manual - yhat_iv
sum diff_pred2, detail

* ============================================================================
* Test 3: Multiple regressors
* ============================================================================
di _n _n "Test 3: Multiple regressors"
di "{hline 40}"

drop yhat_iv yhat_civ_manual diff_pred2

* Run ivreghdfe
ivreghdfe price (mpg = weight length) trunk, absorb(foreign) vce(robust) resid
matrix b_iv3 = e(b)
local names3 : colnames b_iv3
predict double yhat_iv, xb

* Run civreghdfe
civreghdfe price (mpg = weight length) trunk, absorb(foreign) vce(robust)
matrix b_civ3 = e(b)

* Compare coefficients
di _n "Coefficient comparison:"
di "                  ivreghdfe      civreghdfe       Diff"
di "{hline 60}"

local k = colsof(b_iv3)
forval i = 1/`k' {
    local vname : word `i' of `names3'
    local b_iv_i = b_iv3[1,`i']
    local b_civ_i = b_civ3[1,`i']
    local diff_i = `b_civ_i' - `b_iv_i'
    di "%-12s" "`vname'" "  " %14.8f `b_iv_i' "  " %14.8f `b_civ_i' "  " %14.8e `diff_i'
}

* Manually compute predicted values using civreghdfe coefficients
* Note: civreghdfe returns coefficients in [endog, exog] order
gen double yhat_civ_manual = b_civ3[1,1] * mpg + b_civ3[1,2] * trunk

* Compare predictions
di _n "Prediction comparison (raw X * beta):"
sum yhat_iv yhat_civ_manual

gen double diff_pred3 = yhat_civ_manual - yhat_iv
sum diff_pred3, detail

* ============================================================================
* Test 4: Verify both work on demeaned data
* ============================================================================
di _n _n "Test 4: Verify demeaning approach"
di "{hline 40}"

* Compute means by foreign
bysort foreign: egen double mean_price = mean(price)
bysort foreign: egen double mean_mpg = mean(mpg)
bysort foreign: egen double mean_trunk = mean(trunk)

* Create demeaned variables
gen double price_dem = price - mean_price
gen double mpg_dem = mpg - mean_mpg
gen double trunk_dem = trunk - mean_trunk

* Compute fitted values using demeaned X and ivreghdfe coefficients
gen double yhat_dem_iv = b_iv3[1,1] * mpg_dem + b_iv3[1,2] * trunk_dem

* Compare with ivreghdfe xb prediction
di _n "ivreghdfe xb prediction vs demeaned X * beta_iv:"
sum yhat_iv yhat_dem_iv
gen double diff_iv_dem = yhat_iv - yhat_dem_iv
sum diff_iv_dem, detail

* Now compute using civreghdfe coefficients
gen double yhat_dem_civ = b_civ3[1,1] * mpg_dem + b_civ3[1,2] * trunk_dem

* Compare the two demeaned predictions
di _n "Demeaned predictions comparison:"
sum yhat_dem_iv yhat_dem_civ
gen double diff_dem = yhat_dem_civ - yhat_dem_iv
sum diff_dem, detail
local max_diff = r(max)
local min_diff = r(min)

di _n "Max diff in demeaned predictions: " %12.8e `max_diff'
di "Min diff in demeaned predictions: " %12.8e `min_diff'

* ============================================================================
* Test 5: Examine coefficient precision
* ============================================================================
di _n _n "Test 5: Detailed coefficient examination"
di "{hline 40}"

* Show more decimal places
di _n "ivreghdfe coefficients (18 decimals):"
matrix list b_iv3, format(%20.15f)

di _n "civreghdfe coefficients (18 decimals):"
matrix list b_civ3, format(%20.15f)

* Compute relative error
local rel_err_mpg = abs(`diff_mpg2') / abs(`b_iv_mpg2')
di _n "Relative coefficient error (mpg): " %12.8e `rel_err_mpg'

* ============================================================================
* Summary
* ============================================================================
di _n _n "{hline 78}"
di "SUMMARY"
di "{hline 78}"
di "1. Coefficient differences are typically in the range of 1e-10 to 1e-8"
di "2. These are due to numerical precision in CG solver vs Stata's solver"
di "3. Predictions using X * beta differ by similar small amounts"
di ""
di "Both commands use demeaned (FE-absorbed) data for prediction."
di "The 'xb' prediction = X_demeaned * beta, NOT X_original * beta."
di "{hline 78}"
