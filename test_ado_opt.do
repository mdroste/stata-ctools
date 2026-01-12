* Test ado-file optimizations for cmerge
clear all
set more off

* Create master dataset
clear
set obs 5000000
gen id = _n
gen x = runiform()
gen y = rnormal()
gen str20 name = "master" + string(_n)
tempfile master
save `master'

* Create using dataset
clear
set obs 5000000
gen id = _n
gen z = runiform()
gen w = rnormal()
gen str20 label = "using" + string(_n)
tempfile using
save `using'

* Run m:1 merge with verbose timing
use `master', clear
di as text _n "Running cmerge m:1 with 5M obs..."
cmerge m:1 id using `using', verbose

* Check results
assert _N == 5000000
assert _merge == 3 if !mi(_merge)
di as text _n "All tests passed!"
