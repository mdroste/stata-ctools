* Investigate the Hansen J ratio between ivreghdfe and civreghdfe

clear all
set more off
adopath + "build"

sysuse auto, clear

* Get values
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local j_iv = e(j)
local N = e(N)
local df_r = e(df_r)
local df_a = e(df_a)

civreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
local j_civ = e(sargan)

* Compute ratio
local ratio = `j_iv' / `j_civ'

di _n "{hline 78}"
di "Hansen J Analysis"
di "{hline 78}" _n

di "ivreghdfe J:   " %10.6f `j_iv'
di "civreghdfe J:  " %10.6f `j_civ'
di "Ratio (iv/civ):" %10.6f `ratio'

di _n "Sample info:"
di "  N:    " `N'
di "  df_r: " `df_r'
di "  df_a: " `df_a'

* Check if ratio matches any DOF adjustment
di _n "Possible DOF adjustments:"
di "  df_r / N = " %10.6f (`df_r' / `N')
di "  (N-1) / N = " %10.6f ((`N'-1) / `N')
di "  N / (N-1) = " %10.6f (`N' / (`N'-1))
di "  df_r / (df_r + 1) = " %10.6f (`df_r' / (`df_r' + 1))

* The ratio 3.969/5.204 = 0.763 ≈ something?
di _n "Ratio analysis:"
di "  Actual ratio:     " %10.6f `ratio'
di "  df_r / N (71/74): " %10.6f (71/74)
di "  df_r / (df_r+df_a) (71/73): " %10.6f (71/73)

* What if ivreghdfe divides by something extra?
* Check if civreghdfe's J * ratio = ivreghdfe's J
local implied_adj = `j_iv' / `j_civ'
di _n "To get ivreghdfe J from civreghdfe J, multiply by: " %10.6f `implied_adj'

* Check what scaling would give this
* If ivreghdfe uses J * (N-K-df_a)/N as DOF correction:
local K = 1
local dof_adj = (`N' - `K' - `df_a') / `N'
di "DOF adjustment (N-K-df_a)/N = " %10.6f `dof_adj'
di "J_civ * this adjustment = " %10.6f (`j_civ' * `dof_adj')

* Or if it's a different formula for robust VCE
* Hansen J for robust might use: (N-K)/(N-K-df_a) * adjusted_quadform
local adj2 = (`N' - `K') / (`N' - `K' - `df_a')
di _n "Alternative: (N-K)/(N-K-df_a) = " %10.6f `adj2'
di "J_civ / this = " %10.6f (`j_civ' / `adj2')

* Check e(W) and e(S) matrices from ivreghdfe
ivreghdfe price (mpg = weight length), absorb(foreign) vce(robust)
di _n "ivreghdfe e(W) matrix:"
matrix list e(W)
di _n "ivreghdfe e(S) matrix:"
matrix list e(S)

di _n "{hline 78}"
di "Conclusion: The ratio " %10.6f `ratio' " doesn't match simple DOF adjustments."
di "ivreghdfe may use a different robust Hansen J formula than the standard."
di "civreghdfe uses the textbook formula J = (Z'e)' (Z'ΩZ)^{-1} (Z'e)."
di "{hline 78}"
