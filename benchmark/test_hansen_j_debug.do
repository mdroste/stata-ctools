* Debug Hansen J computation difference
* Manually compute components to identify where discrepancy arises

clear all
set more off
adopath + "build"

sysuse auto, clear

di _n "{hline 78}"
di "Hansen J Diagnostic Test"
di "{hline 78}" _n

* Run both commands and compare J statistics
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local j_iv = e(sargan)
local j_iv_df = e(sargan_df)
local j_iv_p = e(sargan_p)

di "ivreghdfe Hansen J: " %10.6f `j_iv' " (df=" `j_iv_df' ", p=" %6.4f `j_iv_p' ")"

civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local j_civ = e(sargan)
local j_civ_df = e(sargan_df)
local j_civ_p = e(sargan_p)

di "civreghdfe Hansen J: " %10.6f `j_civ' " (df=" `j_civ_df' ", p=" %6.4f `j_civ_p' ")"

di _n "Difference: " %10.6f (`j_civ' - `j_iv')
di "Ratio: " %10.6f (`j_civ' / `j_iv')

* ============================================================================
* Now let's manually compute the components
* ============================================================================
di _n _n "{hline 78}"
di "Manual Computation of Hansen J"
di "{hline 78}" _n

* Manually demean data by foreign
bysort foreign: egen double mean_price = mean(price)
bysort foreign: egen double mean_mpg = mean(mpg)
bysort foreign: egen double mean_weight = mean(weight)
bysort foreign: egen double mean_length = mean(length)

gen double price_dem = price - mean_price
gen double mpg_dem = mpg - mean_mpg
gen double weight_dem = weight - mean_weight
gen double length_dem = length - mean_length

* Get coefficients from civreghdfe
civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
matrix b_civ = e(b)
local beta_mpg = b_civ[1,1]

* Compute 2SLS residuals using demeaned data
gen double resid = price_dem - mpg_dem * `beta_mpg'

* Compute Z'e
* Z = [weight_dem, length_dem] (excluded instruments only for overid test)
egen double Ztr_weight = total(weight_dem * resid)
egen double Ztr_length = total(length_dem * resid)

di "Z'e components:"
di "  weight: " Ztr_weight[1]
di "  length: " Ztr_length[1]

* Compute Z'ΩZ (robust weighting)
* Z'ΩZ = sum_i (Z_i Z_i' e_i^2)
gen double e2 = resid^2

gen double ww_e2 = weight_dem^2 * e2
gen double wl_e2 = weight_dem * length_dem * e2
gen double ll_e2 = length_dem^2 * e2

egen double ZOmegaZ_11 = total(ww_e2)
egen double ZOmegaZ_12 = total(wl_e2)
egen double ZOmegaZ_22 = total(ll_e2)

di _n "Z'ΩZ matrix:"
di "  [" ZOmegaZ_11[1] ", " ZOmegaZ_12[1] "]"
di "  [" ZOmegaZ_12[1] ", " ZOmegaZ_22[1] "]"

* Invert 2x2 matrix
local a = ZOmegaZ_11[1]
local b = ZOmegaZ_12[1]
local c = ZOmegaZ_12[1]
local d = ZOmegaZ_22[1]
local det = `a' * `d' - `b' * `c'

local inv_11 = `d' / `det'
local inv_12 = -`b' / `det'
local inv_22 = `a' / `det'

di _n "(Z'ΩZ)^{-1} matrix:"
di "  [" %14.10e `inv_11' ", " %14.10e `inv_12' "]"
di "  [" %14.10e `inv_12' ", " %14.10e `inv_22' "]"

* Compute quadratic form: (Z'e)' (Z'ΩZ)^{-1} (Z'e)
local g1 = Ztr_weight[1]
local g2 = Ztr_length[1]

local quad = `g1' * (`inv_11' * `g1' + `inv_12' * `g2') + `g2' * (`inv_12' * `g1' + `inv_22' * `g2')

di _n "Manual Hansen J (using excluded Z only): " %10.6f `quad'

* ============================================================================
* Now try with FULL Z = [exog, excluded] = [nothing, weight, length]
* Since we have no exogenous regressors, full Z = excluded instruments
* ============================================================================

di _n "{hline 78}"
di "Checking with different Z definitions"
di "{hline 78}"

* Note: In this example, there are no exogenous regressors besides the constant
* (which is absorbed by FE). So Z should just be [weight, length].

* The overidentification df is K_iv - K_total = 2 - 1 = 1
* This means we have 1 overidentifying restriction.

di _n "K_iv (instruments) = 2 (weight, length)"
di "K_total (regressors) = 1 (mpg)"
di "Overid df = K_iv - K_total = 1"

* ============================================================================
* Compare with homoskedastic Sargan
* ============================================================================
di _n _n "{hline 78}"
di "Homoskedastic Sargan comparison"
di "{hline 78}"

* Run with non-robust VCE
ivreghdfe price (mpg = weight length), absorb(foreign)
local sargan_iv = e(sargan)
local sargan_iv_df = e(sargan_df)
local sargan_iv_p = e(sargan_p)

di "ivreghdfe Sargan (homosked): " %10.6f `sargan_iv' " (df=" `sargan_iv_df' ", p=" %6.4f `sargan_iv_p' ")"

civreghdfe price (mpg = weight length), absorb(foreign)
local sargan_civ = e(sargan)
local sargan_civ_df = e(sargan_df)
local sargan_civ_p = e(sargan_p)

di "civreghdfe Sargan (homosked): " %10.6f `sargan_civ' " (df=" `sargan_civ_df' ", p=" %6.4f `sargan_civ_p' ")"

di _n "Difference: " %10.6f (`sargan_civ' - `sargan_iv')

* Compute manual homoskedastic Sargan
* Sargan = N * (Z'e)' (Z'Z)^{-1} (Z'e) / RSS

egen double ZtZ_11 = total(weight_dem^2)
egen double ZtZ_12 = total(weight_dem * length_dem)
egen double ZtZ_22 = total(length_dem^2)

local a = ZtZ_11[1]
local b = ZtZ_12[1]
local c = ZtZ_12[1]
local d = ZtZ_22[1]
local det = `a' * `d' - `b' * `c'

local invZZ_11 = `d' / `det'
local invZZ_12 = -`b' / `det'
local invZZ_22 = `a' / `det'

local N = _N
egen double rss = total(resid^2)
local rss_val = rss[1]

local quad_homosked = `g1' * (`invZZ_11' * `g1' + `invZZ_12' * `g2') + `g2' * (`invZZ_12' * `g1' + `invZZ_22' * `g2')
local sargan_manual = `N' * `quad_homosked' / `rss_val'

di _n "Manual Sargan (homosked): " %10.6f `sargan_manual'
di "RSS: " %14.6f `rss_val'

* ============================================================================
* Summary
* ============================================================================
di _n _n "{hline 78}"
di "SUMMARY"
di "{hline 78}"
di "Homoskedastic Sargan:"
di "  ivreghdfe:  " %10.6f `sargan_iv'
di "  civreghdfe: " %10.6f `sargan_civ'
di "  manual:     " %10.6f `sargan_manual'

di _n "Robust Hansen J:"
di "  ivreghdfe:  " %10.6f `j_iv'
di "  civreghdfe: " %10.6f `j_civ'
di "  manual:     " %10.6f `quad'

di _n "Key observation:"
di "If homoskedastic Sargan matches but robust J differs,"
di "the issue is likely in the Z'ΩZ computation."
di "{hline 78}"
