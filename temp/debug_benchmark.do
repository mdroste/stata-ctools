adopath ++ build
sysuse auto, clear

* Run ivreghdfe
di "=== ivreghdfe ==="
ivreghdfe price (mpg = weight), absorb(foreign)
local ivreghdfe_N_hdfe = e(N_hdfe)
local ivreghdfe_r2_a = e(r2_a)

* Run civreghdfe
sysuse auto, clear
di ""
di "=== civreghdfe ==="
civreghdfe price (mpg = weight), absorb(foreign)
local civreghdfe_N_hdfe = e(N_hdfe)
local civreghdfe_r2_a = e(r2_a)

* Compare
di ""
di "=== Comparison ==="
di "ivreghdfe N_hdfe = `ivreghdfe_N_hdfe'"
di "civreghdfe N_hdfe = `civreghdfe_N_hdfe'"
if `ivreghdfe_N_hdfe' == `civreghdfe_N_hdfe' {
    di "N_hdfe: MATCH"
}
else {
    di "N_hdfe: MISMATCH"
}
di ""
di "ivreghdfe r2_a = `ivreghdfe_r2_a'"
di "civreghdfe r2_a = `civreghdfe_r2_a'"
local diff = `ivreghdfe_r2_a' - `civreghdfe_r2_a'
di "Difference = `diff'"
