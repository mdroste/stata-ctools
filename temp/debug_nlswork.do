adopath ++ build
webuse nlswork, clear
xtset idcode year

* ivreghdfe
di "=== ivreghdfe ==="
ivreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode)
local iv_r2 = e(r2)
local iv_r2_a = e(r2_a)
local iv_N = e(N)
local iv_df_r = e(df_r)
local iv_df_a = e(df_a)
local iv_df_m = e(df_m)
local iv_rmse = e(rmse)

* civreghdfe
webuse nlswork, clear
xtset idcode year
di ""
di "=== civreghdfe ==="
civreghdfe ln_wage tenure (ttl_exp = age), absorb(idcode) vce(cluster idcode)
local civ_r2 = e(r2)
local civ_r2_a = e(r2_a)
local civ_N = e(N)
local civ_df_r = e(df_r)
local civ_df_a = e(df_a)
local civ_df_m = e(df_m)
local civ_rmse = e(rmse)

* Compare
di ""
di "=== Comparison ==="
di "ivreghdfe: r2=`iv_r2' r2_a=`iv_r2_a' N=`iv_N' df_r=`iv_df_r' df_a=`iv_df_a' df_m=`iv_df_m' rmse=`iv_rmse'"
di "civreghdfe: r2=`civ_r2' r2_a=`civ_r2_a' N=`civ_N' df_r=`civ_df_r' df_a=`civ_df_a' df_m=`civ_df_m' rmse=`civ_rmse'"
di ""
di "r2 diff = " `iv_r2' - `civ_r2'
di "r2_a diff = " `iv_r2_a' - `civ_r2_a'
di "rmse diff = " `iv_rmse' - `civ_rmse'
