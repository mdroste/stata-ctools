* Bisect to find which validation causes heap corruption
clear all
set more off
adopath ++ "build"
adopath ++ "../build"

* Run ONLY civreghdfe
di "Testing civreghdfe..."
capture noisily quietly do "validation/validate_civreghdfe.do"
di "civreghdfe complete"

* Now run cdecode validation
di "Testing cdecode validation..."
do "validation/validate_cdecode.do"

di "civreghdfe + cdecode complete!"
