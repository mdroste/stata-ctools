adopath ++ build
webuse nlswork, clear
xtset idcode year

* ivreghdfe
di "=== ivreghdfe ==="
ivreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode)
local iv_r2_a = e(r2_a)
local iv_rmse = e(rmse)
local iv_df_a = e(df_a)

* civreghdfe
webuse nlswork, clear
xtset idcode year
di ""
di "=== civreghdfe ==="
civreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode) noid
local civ_r2_a = e(r2_a)
local civ_rmse = e(rmse)
local civ_df_a = e(df_a)

di ""
di "=== Comparison ==="
di "ivreghdfe r2_a = `iv_r2_a', rmse = `iv_rmse', df_a = `iv_df_a'"
di "civreghdfe r2_a = `civ_r2_a', rmse = `civ_rmse', df_a = `civ_df_a'"
di "r2_a diff = " `iv_r2_a' - `civ_r2_a'
di "rmse diff = " `iv_rmse' - `civ_rmse'
