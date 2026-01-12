* Debug 1:m merge failure
clear all
set more off
adopath + build

* Create master (4 regions)
clear
input byte region str20 division_name
1 "Northeast"
2 "Midwest"
3 "South"
4 "West"
end
tempfile master
save `master'
list

* Create using (50 states with regions)
sysuse census, clear
keep state region pop medage
tempfile using_data
save `using_data'
di "Using: " _N " obs"

* Run Stata merge
use `master', clear
merge 1:m region using `using_data'
tempfile stata_result
save `stata_result'
di ""
di "=== STATA MERGE RESULT ==="
di "N = " _N
tab _merge
list in 1/15

* Run cmerge
use `master', clear
cmerge 1:m region using `using_data', verbose
di ""
di "=== CMERGE RESULT ==="
di "N = " _N
tab _merge
list in 1/15

* Compare after sorting
di ""
di "=== COMPARISON (after sorting) ==="
sort region state
tempfile cmerge_sorted
save `cmerge_sorted'

use `stata_result', clear
sort region state
cf _all using `cmerge_sorted'
