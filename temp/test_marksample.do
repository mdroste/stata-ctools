capture do "validation/validate_setup.do"

clear
set obs 3
gen str10 x = "A"

di "Before marksample:"
describe

* Simulate what cencode does
local varlist = "x"
marksample touse, strok

di ""
di "After marksample:"
describe

* Now try the plugin call
quietly gen long y = .
di ""
di "Plugin call vars: `varlist' y"
plugin call ctools_plugin `varlist' y, "cencode label=y"
list
