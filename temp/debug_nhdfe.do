adopath ++ build
sysuse auto, clear
civreghdfe price (mpg = weight), absorb(foreign)
di "G = " e(G)
di "N_hdfe = " e(N_hdfe)
capture confirm scalar __civreghdfe_num_levels_1
if _rc == 0 {
    di "Scalar exists, value = " scalar(__civreghdfe_num_levels_1)
}
else {
    di "Scalar does not exist"
}
scalar list
