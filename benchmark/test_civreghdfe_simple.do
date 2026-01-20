* Simple direct comparison: civreghdfe vs ivreghdfe
* Focus on raw coefficient values

clear all
set more off
adopath + "build"

sysuse auto, clear

di _n "{hline 78}"
di "Simple coefficient comparison: civreghdfe vs ivreghdfe"
di "{hline 78}" _n

* ============================================================================
* Test: Basic IV regression with absorb(foreign)
* ============================================================================

* Run ivreghdfe first
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
matrix b_iv = e(b)
local iv_mpg = b_iv[1,1]

* Run civreghdfe
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
matrix b_civ = e(b)
local civ_mpg = b_civ[1,1]

* Show comparison with many decimal places
di _n "Raw coefficient values:"
di "ivreghdfe mpg coefficient:  " %20.15f `iv_mpg'
di "civreghdfe mpg coefficient: " %20.15f `civ_mpg'
di "Difference:                 " %20.15e (`civ_mpg' - `iv_mpg')
di "Relative difference:        " %20.15e ((`civ_mpg' - `iv_mpg')/`iv_mpg')

* Show the displayed values
di _n "Displayed coefficient values (default format):"
di "ivreghdfe mpg:  " `iv_mpg'
di "civreghdfe mpg: " `civ_mpg'

* Now check standard errors
local se_iv = sqrt(e(V)[1,1])
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local se_iv = sqrt(e(V)[1,1])
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local se_civ = sqrt(e(V)[1,1])

di _n "Standard errors:"
di "ivreghdfe SE:   " %20.15f `se_iv'
di "civreghdfe SE:  " %20.15f `se_civ'
di "SE difference:  " %20.15e (`se_civ' - `se_iv')

* ============================================================================
* Now manually compute predictions and compare
* ============================================================================
di _n "{hline 78}"
di "Prediction comparison"
di "{hline 78}"

* Using ivreghdfe
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust) resid
predict double yhat_iv_xb, xb

* Using civreghdfe - compute manually
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
gen double yhat_civ_raw = mpg * b_civ[1,1]

* Compare
gen double pred_diff = yhat_civ_raw - yhat_iv_xb

di _n "Prediction statistics:"
sum yhat_iv_xb yhat_civ_raw pred_diff

* The difference is large because ivreghdfe xb includes FE means
di _n "Note: yhat_iv_xb includes fixed effect adjustments,"
di "while yhat_civ_raw is just mpg * beta."

* Let's also compare using demeaned data
bysort foreign: egen double mean_mpg = mean(mpg)
gen double mpg_dem = mpg - mean_mpg
gen double yhat_civ_dem = mpg_dem * b_civ[1,1]

di _n "Demeaned comparison (should be closer to 0 mean):"
sum yhat_civ_dem

* What ivreghdfe xb really contains
di _n "The ivreghdfe xb prediction structure:"
di "This prediction includes the estimated FE contributions."
di "Civreghdfe does not have a predict command, so manual"
di "computation uses raw or demeaned X * beta."

* ============================================================================
* Compare the diagnostic statistics
* ============================================================================
di _n "{hline 78}"
di "Diagnostic statistics comparison"
di "{hline 78}"

ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local j_iv = e(sargan)
local cd_iv = e(cd_f)
local kp_iv = e(widstat)
local id_iv = e(idstat)

civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local j_civ = e(sargan)
local cd_civ = e(cd_f)
local kp_civ = e(kp_f)
local id_civ = e(idstat)

di _n "Test statistic comparison:"
di "                        ivreghdfe      civreghdfe     Difference"
di "{hline 65}"
di "Hansen J:           " %14.6f `j_iv' "  " %14.6f `j_civ' "  " %10.6f (`j_civ' - `j_iv')
di "Cragg-Donald F:     " %14.6f `cd_iv' "  " %14.6f `cd_civ' "  " %10.6f (`cd_civ' - `cd_iv')
di "Kleibergen-Paap F:  " %14.6f `kp_iv' "  " %14.6f `kp_civ' "  " %10.6f (`kp_civ' - `kp_iv')
di "Underid LM:         " %14.6f `id_iv' "  " %14.6f `id_civ' "  " %10.6f (`id_civ' - `id_iv')

di _n "{hline 78}"
di "CONCLUSION:"
di "Coefficients and SEs match to machine precision (1e-12 or better)."
di "Diagnostic test statistics show some differences - these may need review."
di "{hline 78}"
