* Debug long string key merge
clear
adopath + "build"

* Very long string keys (200+ characters)
set obs 5  // Smaller for debugging
gen str244 long_key = "prefix_" + substr("a" * 200, 1, 200) + "_" + string(_n)
gen value1 = _n * 10
tempfile str_long_master
save `str_long_master'

di "Master data:"
list

clear
set obs 5
gen str244 long_key = "prefix_" + substr("a" * 200, 1, 200) + "_" + string(_n)
gen value2 = _n * 100
tempfile str_long_using
save `str_long_using'

di "Using data:"
list

* Stata merge
use `str_long_master', clear
merge 1:1 long_key using `str_long_using', nogenerate
di "Stata merge result:"
list, abbreviate(20)
tempfile stata_result
save `stata_result'

* cmerge
use `str_long_master', clear
cmerge 1:1 long_key using `str_long_using', nogenerate noreport
di "cmerge result:"
list, abbreviate(20)

* Compare
cf _all using `stata_result'
