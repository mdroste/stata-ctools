* Test: Run validations then direct cdecode (not validate_cdecode.do)
clear all
set more off
adopath ++ "build"
adopath ++ "../build"

di "Testing csort..."
capture noisily quietly do "validation/validate_csort.do"
di "Testing cmerge..."
capture noisily quietly do "validation/validate_cmerge.do"
di "Testing cimport..."
capture noisily quietly do "validation/validate_cimport.do"
di "Testing cexport..."
capture noisily quietly do "validation/validate_cexport.do"
di "Testing creghdfe..."
capture noisily quietly do "validation/validate_creghdfe.do"
di "Testing cqreg..."
capture noisily quietly do "validation/validate_cqreg.do"
di "Testing civreghdfe..."
capture noisily quietly do "validation/validate_civreghdfe.do"

di "Now running direct cdecode tests (not validate_cdecode.do)..."

* Test 1: Basic cdecode
sysuse auto, clear
cdecode foreign, gen(foreign_str)
di "Test 1 PASSED"

* Test 2: Multiple cdecodes
sysuse auto, clear
cdecode foreign, gen(foreign_str)
cdecode rep78, gen(rep78_str)
di "Test 2 PASSED"

* Test 3: Large dataset
clear
set obs 100000
gen id = mod(_n, 100) + 1
label define idlab 1 "A" 2 "B" 3 "C" 4 "D" 5 "E" 10 "J" 50 "Middle" 100 "End"
label values id idlab
cdecode id, gen(id_str)
di "Test 3 PASSED"

di "All direct cdecode tests completed successfully!"
