sysuse auto, clear
di "=== ivreghdfe with absorb(foreign) ==="
ivreghdfe price (mpg = weight), absorb(foreign)
di "e(N_hdfe) = " e(N_hdfe)
di "e(G) = " e(G)
ereturn list
