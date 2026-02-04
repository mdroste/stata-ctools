* Test: Find which command causes heap corruption
clear all
set more off

di "Testing csort..."
quietly do "validation/validate_csort.do"
di "csort complete"

di "Testing cmerge..."
quietly do "validation/validate_cmerge.do"
di "cmerge complete"

di "Testing cimport..."
quietly do "validation/validate_cimport.do"
di "cimport complete"

di "Testing cexport..."
quietly do "validation/validate_cexport.do"
di "cexport complete"

di "Testing creghdfe..."
quietly do "validation/validate_creghdfe.do"
di "creghdfe complete"

di "Testing cqreg..."
quietly do "validation/validate_cqreg.do"
di "cqreg complete"

di "Now testing cdecode (before civreghdfe)..."
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "cdecode passed!"

di "Now testing civreghdfe..."
quietly do "validation/validate_civreghdfe.do"
di "civreghdfe complete"

di "Now testing cdecode again (after civreghdfe)..."
sysuse auto, clear
cdecode foreign, gen(foreign_str2)
di "cdecode passed again!"

di "All tests completed successfully!"
