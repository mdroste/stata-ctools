adopath ++ build
sysuse auto, clear
di "Running civreghdfe..."
civreghdfe price (mpg = weight), absorb(foreign)
di ""
di "=== After civreghdfe ==="
di "e(G) = " e(G)
di "e(N_hdfe) = " e(N_hdfe)
di "e(df_a) = " e(df_a)
local nhdfe_val = e(N_hdfe)
di "local nhdfe_val = `nhdfe_val'"
if `nhdfe_val' == 2 {
    di "SUCCESS: N_hdfe = 2"
}
else {
    di "FAILURE: N_hdfe = `nhdfe_val', expected 2"
}
