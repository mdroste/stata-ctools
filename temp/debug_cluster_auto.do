adopath ++ build
sysuse auto, clear

* ivreghdfe
ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign)
local iv_r2_a = e(r2_a)
local iv_df_a = e(df_a)

* civreghdfe
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) verbose
local civ_r2_a = e(r2_a)
local civ_df_a = e(df_a)
capture local civ_df_a_for_vce = scalar(__civreghdfe_df_a_for_vce)

di ""
di "=== Comparison ==="
di "ivreghdfe r2_a = `iv_r2_a', df_a = `iv_df_a'"
di "civreghdfe r2_a = `civ_r2_a', df_a = `civ_df_a', df_a_for_vce = `civ_df_a_for_vce'"
di "r2_a diff = " `iv_r2_a' - `civ_r2_a'
