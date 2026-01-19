* Test civreghdfe after optimization changes
* Tests: sortpreserve removed, coef reordering in C, nested FE detection in C

clear all

* Add build directory to ado path
adopath + build

set seed 12345
set obs 100000

* Generate data - standard IV setup (fewer FE levels to avoid near-singularity)
gen firm = ceil(_n/500)
gen year = mod(_n-1, 10) + 2000
gen industry = ceil(_n/1000)

* Exogenous variables
gen x1 = rnormal()
gen x2 = rnormal()

* Instruments
gen z1 = rnormal()
gen z2 = rnormal()

* Error term
gen u = rnormal()

* Endogenous variable (correlated with u, but z1/z2 are valid instruments)
gen endogx = 0.8*z1 + 0.4*z2 + 0.5*u + rnormal(0, 0.5)

* Outcome
gen y = 0.5*endogx + 0.3*x1 + 0.2*x2 + 0.8*u + rnormal()

* Load plugin
local plugin_name "ctools_mac_arm.plugin"
capture program ctools_plugin, plugin using("build/`plugin_name'")

di as text _n "============================================="
di as text "Test 1: Basic IV regression (no controls)"
di as text "=============================================" _n

timer clear 1
timer on 1
civreghdfe y (endogx = z1 z2), absorb(firm year) vce(cluster firm) verbose timeit
timer off 1
qui timer list 1
di as text _n "Total Stata-side elapsed time: " r(t1) " seconds" _n

* Check coefficients are in expected order
di as text "Coefficient vector:"
matrix list e(b)

di as text _n "============================================="
di as text "Test 2: Quick correctness check"
di as text "=============================================" _n

* Verify coefficient signs make sense
di as text "endogx coef (expected ~0.5): " _b[endogx]

* Check N and df
di as text _n "N = " e(N) " (expected ~100000)"
di as text "N_clust = " e(N_clust)

di as text _n "============================================="
di as text "Test 3: With exogenous controls"
di as text "=============================================" _n

civreghdfe y (endogx = z1 z2) x1 x2, absorb(firm year) vce(cluster firm) timeit verbose
di as text "Coefficients with controls:"
di as text "endogx coef (expected ~0.5): " _b[endogx]
di as text "x1 coef (expected ~0.3): " _b[x1]
di as text "x2 coef (expected ~0.2): " _b[x2]

di as text _n "============================================="
di as text "Test 4: Large dataset performance (1M obs)"
di as text "=============================================" _n

clear
set obs 1000000

* Generate large dataset (fewer FE levels)
gen firm = ceil(_n/500)
gen year = mod(_n-1, 10) + 2000
gen z1 = rnormal()
gen z2 = rnormal()
gen u = rnormal()
gen endogx = 0.8*z1 + 0.4*z2 + 0.5*u + rnormal(0, 0.5)
gen y = 0.5*endogx + 0.8*u + rnormal()

timer clear 2
timer on 2
civreghdfe y (endogx = z1 z2), absorb(firm year) vce(cluster firm) timeit
timer off 2
qui timer list 2
di as text _n "Large dataset elapsed time: " r(t2) " seconds"

di as text _n "============================================="
di as text "All tests completed successfully!"
di as text "============================================="
