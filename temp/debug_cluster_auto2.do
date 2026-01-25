adopath ++ build
sysuse auto, clear
ivreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign)
local iv_r2 = e(r2)
local iv_r2_a = e(r2_a)
local iv_N = e(N)
local iv_df_m = e(df_m)
local iv_df_a = e(df_a)

sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) vce(cluster foreign) noid
local civ_r2 = e(r2)
local civ_r2_a = e(r2_a)
local civ_N = e(N)
local civ_df_m = e(df_m)
local civ_df_a = e(df_a)

di "ivreghdfe: r2=`iv_r2' r2_a=`iv_r2_a' N=`iv_N' df_m=`iv_df_m' df_a=`iv_df_a'"
di "civreghdfe: r2=`civ_r2' r2_a=`civ_r2_a' N=`civ_N' df_m=`civ_df_m' df_a=`civ_df_a'"

* Expected r2_a using ivreghdfe formula with df_a=0
local expected_r2_a = 1 - (1 - `iv_r2') * (`iv_N' - 1) / (`iv_N' - `iv_df_m' - 0 - 1)
di "Expected r2_a (df_a=0): `expected_r2_a'"

* civreghdfe r2_a using its own r2
local expected_civ_r2_a = 1 - (1 - `civ_r2') * (`civ_N' - 1) / (`civ_N' - `civ_df_m' - 0 - 1)
di "civreghdfe expected r2_a (df_a=0, its r2): `expected_civ_r2_a'"
