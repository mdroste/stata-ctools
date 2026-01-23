capture do "validation/validate_setup.do"
if _rc != 0 do "validate_setup.do"

webuse lifeexp, clear
di "Before cencode:"
describe
di ""

capture noisily cencode country, generate(test) verbose
di "rc = " _rc
list country test in 1/5
