clear all
set more off
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
log using "/tmp/ctools-data-audit/bins_metadata.log", text replace
set obs 20
gen x=_n
gen y=x^2+mod(_n,3)
gen g=mod(_n,2)
gen c=mod(_n,3)
cbinscatter y x, nq(2) by(g) controls(i.c) nograph
local recorded_by `"`e(by)'"'
local recorded_controls `"`e(controls)'"'
di "RECORDED_BY=`recorded_by'"
di "RECORDED_CONTROLS=`recorded_controls'"
capture confirm variable `recorded_by'
di "RECORDED_BY_EXISTS=" (_rc==0)
foreach v of local recorded_controls {
    capture confirm variable `v'
    di "CONTROL_`v'_EXISTS=" (_rc==0)
}
log close
