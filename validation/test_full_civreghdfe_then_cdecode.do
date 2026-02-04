* Test: Full civreghdfe validation then cdecode
clear all
set more off

di "Step 1: Run full civreghdfe validation..."
quietly do "validation/validate_civreghdfe.do"
di "civreghdfe validation complete"

di "Step 2: Run cdecode..."
sysuse auto, clear
cdecode foreign, gen(foreign_str)

di "All tests completed successfully!"
