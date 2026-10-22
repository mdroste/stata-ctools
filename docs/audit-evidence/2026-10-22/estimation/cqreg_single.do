clear all
set more off
capture log close _all
log using "/tmp/ctools-regression-audit/cqreg_single.log", text replace
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
set obs 300
set seed 845
gen double x=rnormal()
gen double y=.4*x+rnormal()
gen cl=1
capture noisily qreg y x, vce(excluster cl, fitted)
di "NATIVE SINGLE RC=" _rc
capture noisily cqreg y x, vce(cluster cl)
di "CTOOLS SINGLE RC=" _rc
matrix list e(V)
scalar v=e(V)[1,1]
di "CTOOLS V11=" v " MISSING=" missing(v)
capture noisily cqreg y x, vce(cluster cl) denmethod(kernel(gaussian))
di "CTOOLS KERNEL SINGLE RC=" _rc
log close
exit, clear
