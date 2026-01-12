* Quick test of cmerge with Java variable creation
clear all
set more off
adopath + build

* Create master dataset
set obs 25000000
gen id = _n
gen x = runiform()
tempfile master
save `master'

* Create using dataset
clear
set obs 25000000
gen id = _n
gen y = runiform()
gen z = rnormal()
tempfile using
save `using'

* Run cmerge
use `master', clear
di "Running cmerge m:1..."
cmerge m:1 id using `using', verbose

* Check results
di _n "Checking results..."
assert _N == 25000000
assert _merge == 3
sum x y z
di "Success!"
