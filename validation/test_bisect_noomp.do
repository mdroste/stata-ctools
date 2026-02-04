* Test: civreghdfe + cdecode validation with OMP_NUM_THREADS=1
clear all
set more off
adopath ++ "build"
adopath ++ "../build"

di "Testing civreghdfe validation..."
capture noisily quietly do "validation/validate_civreghdfe.do"
di "civreghdfe complete"

di "Testing cdecode validation..."
do "validation/validate_cdecode.do"

di "Complete!"
