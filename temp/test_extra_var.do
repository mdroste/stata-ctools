capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

clear
set obs 3

* Create an extra variable before our test vars
gen byte extra1 = 1

gen str10 x = "A"
replace x = "B" in 2  
replace x = "C" in 3

di "Dataset:"
describe
list

* cencode should work - we're only passing x and y to the plugin
cencode x, generate(y)

di ""
di "After cencode:"
list
summarize y
