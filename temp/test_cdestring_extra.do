capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

clear
set obs 3

* Create an extra variable before our test vars
gen byte extra1 = 1

gen str5 x = "1"
replace x = "2" in 2
replace x = "3" in 3

di "Dataset:"
describe
list

cdestring x, generate(y)

di ""
di "After cdestring:"
list
summarize y
