capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

clear
set obs 3
gen byte first = 1
gen int source = _n
label define lbl 1 "A" 2 "B" 3 "C"
label values source lbl
gen byte third = 3

di "Dataset before cdecode:"
describe
list

cdecode source, generate(dest)
di ""
di "After cdecode:"
list
