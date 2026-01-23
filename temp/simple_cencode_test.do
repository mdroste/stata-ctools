* Load setup first
capture do "validation/validate_setup.do"

* Minimal test
clear
set obs 3
gen str10 x = ""
replace x = "A" in 1
replace x = "B" in 2
replace x = "C" in 3
list

* Try cencode
capture noisily cencode x, generate(y) verbose
di "rc = " _rc
list

* Check the return values
di "n_unique = " _cencode_n_unique
