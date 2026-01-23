capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

clear
set obs 3
gen str10 x = "A"
replace x = "B" in 2  
replace x = "C" in 3

di "Before cencode:"
describe
list

cencode x, generate(y) verbose

di ""
di "After cencode:"
list
summarize y
