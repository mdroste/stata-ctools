capture do "validation/validate_setup.do"

clear
set obs 3
gen str10 x = "A"
replace x = "B" in 2
replace x = "C" in 3

di "Before cencode:"
describe

* Manually check what the plugin call would look like
local varlist = "x"
local generate = "y"
local ifcond = ""
local incond = ""

di "varlist = <`varlist'>"
di "generate = <`generate'>"
di "if = <`ifcond'>"
di "in = <`incond'>"

* Now run cencode
cencode x, generate(y)
list
