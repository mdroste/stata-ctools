adopath ++ build
sysuse auto, clear

* Run ivreghdfe
ivreghdfe price (mpg = weight), absorb(foreign)
local iv_r2 = e(r2)
local iv_r2_a = e(r2_a)
local iv_N = e(N)
local iv_df_r = e(df_r)
local iv_df_a = e(df_a)
local iv_df_m = e(df_m)

* Run civreghdfe
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign)
local civ_r2 = e(r2)
local civ_r2_a = e(r2_a)
local civ_N = e(N)
local civ_df_r = e(df_r)
local civ_df_a = e(df_a)
local civ_df_m = e(df_m)

* Compare
di ""
di "=== ivreghdfe ==="
di "r2 = `iv_r2'"
di "r2_a = `iv_r2_a'"
di "N = `iv_N'"
di "df_r = `iv_df_r'"
di "df_a = `iv_df_a'"
di "df_m = `iv_df_m'"

di ""
di "=== civreghdfe ==="
di "r2 = `civ_r2'"
di "r2_a = `civ_r2_a'"
di "N = `civ_N'"
di "df_r = `civ_df_r'"
di "df_a = `civ_df_a'"
di "df_m = `civ_df_m'"

* Check formula: r2_a = 1 - (1 - r2) * (N - 1) / df_r
local formula1 = 1 - (1 - `iv_r2') * (`iv_N' - 1) / `iv_df_r'
di ""
di "Formula 1 (N-1)/df_r: `formula1'"
local formula2 = 1 - (1 - `iv_r2') * (`iv_N' - 1) / (`iv_N' - `civ_df_m' - `civ_df_a')
di "Formula 2 (N-1)/(N-K-df_a): `formula2'"
local formula3 = 1 - (1 - `iv_r2') * (`iv_N' - 1) / (`iv_N' - `civ_df_m' - `civ_df_a' - 1)
di "Formula 3 (N-1)/(N-K-df_a-1): `formula3'"
