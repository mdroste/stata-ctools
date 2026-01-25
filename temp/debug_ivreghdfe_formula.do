adopath ++ build
sysuse auto, clear

* ivreghdfe - check internal values
ivreghdfe price (mpg = weight), absorb(foreign)
di "ivreghdfe: r2=" e(r2) ", r2_a=" e(r2_a) ", N=" e(N) ", df_m=" e(df_m)
di "ivreghdfe: dofminus=" e(dofminus) ", sdofminus=" e(sdofminus) ", df_a=" e(df_a)

* Calculate expected r2_a using ivreghdfe formula
local r2 = e(r2)
local N = e(N)
local rankxx = e(df_m)  // df_m = number of regressors
local dofminus = e(dofminus)
local sdofminus = e(sdofminus)
local expected_r2_a = 1 - (1 - `r2') * (`N' - 1) / (`N' - `rankxx' - `dofminus' - `sdofminus')
di "Expected r2_a = `expected_r2_a'"

* civreghdfe
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign) noid
di ""
di "civreghdfe: r2=" e(r2) ", r2_a=" e(r2_a) ", N=" e(N) ", df_m=" e(df_m) ", df_a=" e(df_a)

* What denom is civreghdfe using?
local civ_r2 = e(r2)
local civ_N = e(N)
local civ_K = e(df_m)
local civ_df_a = e(df_a)
local civ_denom = `civ_N' - `civ_K' - `civ_df_a' - 1
di "civreghdfe denom = N - K - df_a - 1 = `civ_N' - `civ_K' - `civ_df_a' - 1 = `civ_denom'"
local civ_expected = 1 - (1 - `civ_r2') * (`civ_N' - 1) / `civ_denom'
di "civreghdfe expected r2_a = `civ_expected'"
