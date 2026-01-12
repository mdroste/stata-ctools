* Test cmerge output table formatting
clear all
set more off
adopath + build

di _n "=============================================="
di "TEST 1: All matched - should NOT show detail rows"
di "=============================================="

sysuse auto, clear
gen id = _n
keep id make price
tempfile master
save `master'

sysuse auto, clear
gen id = _n
keep id weight length
tempfile using_data
save `using_data'

use `master', clear
cmerge 1:1 id using `using_data'

di _n "=============================================="
di "TEST 2: keep(using) - should show 0 for master and matched"
di "=============================================="

* Create partial overlap
clear
set obs 100
gen id = _n
gen x = runiform()
tempfile master2
save `master2'

clear
set obs 100
gen id = _n + 50  // IDs 51-150
gen y = runiform()
tempfile using2
save `using2'

use `master2', clear
di "Before merge: should show partial overlap"
cmerge 1:1 id using `using2'

use `master2', clear
di _n "With keep(using): should show 0 for master and matched"
cmerge 1:1 id using `using2', keep(using)

use `master2', clear
di _n "With keep(match): should show 0 for master and using (detail rows hidden)"
cmerge 1:1 id using `using2', keep(match)

use `master2', clear
di _n "With keep(master using): should show 0 for matched"
cmerge 1:1 id using `using2', keep(master using)

di _n "=============================================="
di "TEST 3: Compare with Stata merge output"
di "=============================================="

use `master2', clear
di _n "Stata merge with keep(using):"
merge 1:1 id using `using2', keep(using)
