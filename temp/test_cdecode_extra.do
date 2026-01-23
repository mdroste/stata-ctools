capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

clear
set obs 3

* Create an extra variable before our test vars
gen byte extra1 = 1

gen int source = _n
label define lbl 1 "A" 2 "B" 3 "C"
label values source lbl

di "Dataset:"
describe
list

cdecode source, generate(dest)

di ""
di "After cdecode:"
list
