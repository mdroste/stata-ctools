* Run the actual civreghdfe validation, then cdecode to check for crash

clear all
set more off
adopath ++ "build"

di "=============================================="
di "Running full civreghdfe validation..."
di "=============================================="

quietly do "validation/validate_civreghdfe.do"

di ""
di "civreghdfe validation complete"
di ""

* Now run cdecode to test if heap is corrupted
di "Testing heap integrity with cdecode..."
clear
set obs 100
gen id = ceil(_n / 10)
gen id_label = "group_" + string(id)
encode id_label, gen(id_encoded)
cdecode id_encoded, gen(id_decoded)

di ""
di "If you see this, no crash occurred"
