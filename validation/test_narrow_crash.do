* Test: Run validate_all steps one by one with cdecode check
clear all
set more off
adopath ++ "build"
adopath ++ "../build"

di "Test 1: cdecode before anything"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 2: After csort validation"
capture noisily quietly do "validation/validate_csort.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 3: After cmerge validation"
capture noisily quietly do "validation/validate_cmerge.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 4: After cimport validation"
capture noisily quietly do "validation/validate_cimport.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 5: After cexport validation"
capture noisily quietly do "validation/validate_cexport.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 6: After creghdfe validation"
capture noisily quietly do "validation/validate_creghdfe.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 7: After cqreg validation"
capture noisily quietly do "validation/validate_cqreg.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "Test 8: After civreghdfe validation"
capture noisily quietly do "validation/validate_civreghdfe.do"
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "  cdecode OK"

di "All tests completed successfully!"
