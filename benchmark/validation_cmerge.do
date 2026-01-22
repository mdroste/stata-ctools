*! Validation tests for cmerge - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cmerge.log", replace text
validation_header "cmerge"

* === Setup test data ===
clear
set obs 100
gen id = _n
gen master_val = runiform()
gen shared_val = mod(_n, 10)
tempfile master
save `master'

clear
set obs 120
gen id = _n - 10
gen using_val = rnormal()
gen shared_val = mod(_n, 5)
tempfile using_data
save `using_data'

* Duplicates for m:1 tests
clear
set obs 200
gen id = mod(_n - 1, 100) + 1
gen val = runiform()
tempfile master_dup
save `master_dup'

* === Merge types ===
use `master', clear
capture drop _merge
run_test "1:1 merge" cmerge 1:1 id using `using_data'

use `master_dup', clear
capture drop _merge
run_test "m:1 merge" cmerge m:1 id using `master'

use `master', clear
capture drop _merge
run_test "1:m merge" cmerge 1:m id using `master_dup'

use `master_dup', clear
capture drop _merge
run_test "m:m merge" cmerge m:m id using `master_dup'

* === keep() option ===
foreach k in match master using {
    use `master', clear
    capture drop _merge
    run_test "keep(`k')" cmerge 1:1 id using `using_data', keep(`k')
}

use `master', clear
capture drop _merge
run_test "keep(1 3)" cmerge 1:1 id using `using_data', keep(1 3)

* === generate/nogenerate ===
use `master', clear
capture cmerge 1:1 id using `using_data', generate(my_merge)
capture confirm variable my_merge
local passed = (_rc == 0)
report_test "generate()" `passed'

use `master', clear
capture cmerge 1:1 id using `using_data', nogenerate
capture confirm variable _merge
local passed = (_rc != 0)
report_test "nogenerate" `passed'

* === Other options ===
use `master', clear
capture drop _merge
run_test "keepusing()" cmerge 1:1 id using `using_data', keepusing(using_val)

use `master', clear
capture drop _merge
run_test "verbose" cmerge 1:1 id using `using_data', verbose

use `master', clear
capture drop _merge
run_test "timeit" cmerge 1:1 id using `using_data', timeit

use `master', clear
capture drop _merge
run_test "noreport" cmerge 1:1 id using `using_data', noreport

use `master', clear
capture drop _merge
run_test "nolabel" cmerge 1:1 id using `using_data', nolabel

use `master', clear
capture drop _merge
run_test "nonotes" cmerge 1:1 id using `using_data', nonotes

use `master', clear
capture drop _merge
run_test "threads(2)" cmerge 1:1 id using `using_data', threads(2)

* === sorted option ===
use `master', clear
sort id
save `master', replace
use `using_data', clear
sort id
save `using_data', replace
use `master', clear
capture drop _merge
run_test "sorted" cmerge 1:1 id using `using_data', sorted

* === update/replace ===
clear
set obs 50
gen id = _n
gen value = runiform()
replace value = . in 1/10
tempfile master_miss
save `master_miss'

clear
set obs 50
gen id = _n
gen value = 999
tempfile using_nomiss
save `using_nomiss'

use `master_miss', clear
capture drop _merge
run_test "update" cmerge 1:1 id using `using_nomiss', update

use `master_miss', clear
capture drop _merge
run_test "replace" cmerge 1:1 id using `using_nomiss', replace

* === Multiple keys ===
clear
set obs 100
gen state = mod(_n - 1, 10) + 1
gen year = 2000 + mod(_n - 1, 10)
gen val1 = runiform()
tempfile compound
save `compound'

use `compound', clear
run_test "Multiple keys" cmerge 1:1 state year using `compound'

* === _n merge ===
use `master', clear
capture drop _merge
clear
set obs 100
gen using_nval = runiform()
tempfile using_n
save `using_n'
use `master', clear
capture drop _merge
run_test "1:1 _n merge" cmerge 1:1 _n using `using_n'

* === Return values ===
use `master', clear
capture drop _merge
capture cmerge 1:1 id using `using_data'
local passed = (r(N) > 0 & r(N_1) >= 0 & r(N_2) >= 0 & r(N_3) >= 0)
report_test "Return values" `passed'

* === Error cases ===
use `master', clear
capture drop _merge
run_test_error "Invalid merge type" 198 cmerge invalid id using `using_data'

use `master', clear
capture drop _merge
run_test_error "Missing file" 601 cmerge 1:1 id using "nonexistent_file.dta"

* === assert option ===
use `master', clear
capture drop _merge
capture cmerge 1:1 id using `using_data', assert(match)
local passed = (_rc == 9)
report_test "assert(match) fails" `passed'

* === Large merge ===
clear
set obs 50000
gen id = _n
gen val = runiform()
tempfile large
save `large'

clear
set obs 50000
gen id = _n
gen val2 = rnormal()
tempfile large_using
save `large_using'

use `large', clear
capture cmerge 1:1 id using `large_using'
local passed = (_rc == 0 & r(N_3) == 50000)
report_test "Large merge (50k)" `passed'

* === preserve_order ===
use `master', clear
capture drop _merge
run_test "preserve_order(1)" cmerge 1:1 id using `using_data', preserve_order(1)

validation_summary
log close
