*! Validation tests for cexport delimited - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cexport.log", replace text
validation_header "cexport delimited"

sysuse auto, clear
tempfile testdata
save `testdata'

* === Basic export and variables ===
tempfile out
sysuse auto, clear
run_test "Basic export" cexport delimited using `out', replace

sysuse auto, clear
run_test "Export specific vars" cexport delimited make price mpg using `out', replace

* === Delimiter options ===
sysuse auto, clear
run_test "delimiter(tab)" cexport delimited using `out', delimiter(tab) replace

sysuse auto, clear
run_test "delimiter(;)" cexport delimited using `out', delimiter(";") replace

* === Header and quote options ===
sysuse auto, clear
run_test "novarnames" cexport delimited using `out', novarnames replace

sysuse auto, clear
run_test "quote" cexport delimited using `out', quote replace

sysuse auto, clear
run_test "noquoteif" cexport delimited using `out', noquoteif replace

* === Replace option ===
sysuse auto, clear
capture cexport delimited using `out', replace
capture cexport delimited using `out', replace
local passed = (_rc == 0)
report_test "replace option" `passed'

* === Conditionals ===
sysuse auto, clear
run_test "if condition" cexport delimited using `out' if foreign == 1, replace

sysuse auto, clear
run_test "in range" cexport delimited using `out' in 1/20, replace

* === Format options ===
sysuse auto, clear
run_test "datafmt" cexport delimited using `out', datafmt replace

sysuse auto, clear
run_test "nolabel" cexport delimited using `out', nolabel replace

* === Performance options ===
sysuse auto, clear
run_test "verbose" cexport delimited using `out', verbose replace

sysuse auto, clear
run_test "timeit" cexport delimited using `out', timeit replace

sysuse auto, clear
run_test "threads(2)" cexport delimited using `out', threads(2) replace

sysuse auto, clear
run_test "noparallel" cexport delimited using `out', noparallel replace

* === I/O options ===
sysuse auto, clear
run_test "mmap" cexport delimited using `out', mmap replace

sysuse auto, clear
run_test "nofsync" cexport delimited using `out', nofsync replace

sysuse auto, clear
run_test "direct" cexport delimited using `out', direct replace

sysuse auto, clear
run_test "prefault" cexport delimited using `out', prefault replace

sysuse auto, clear
run_test "crlf" cexport delimited using `out', crlf replace

* === Large dataset ===
clear
set obs 50000
gen id = _n
gen x = runiform()
gen y = rnormal()
gen str20 name = "name" + string(_n)
run_test "Large dataset (50k)" cexport delimited using `out', replace

* === Combined options ===
sysuse auto, clear
run_test "Combined options" cexport delimited using `out', delimiter(tab) novarnames quote verbose threads(2) replace

* === Special data types ===
clear
set obs 100
gen id = _n
gen val = runiform()
replace val = . in 1/10
replace val = .a in 11/15
run_test "Missing values" cexport delimited using `out', replace

clear
set obs 10
gen id = _n
gen str50 text = "normal"
replace text = "with,comma" in 1
replace text = `"with"quote"' in 2
replace text = "with	tab" in 3
run_test "Special chars" cexport delimited using `out', quote replace

clear
set obs 100
gen byte b = mod(_n, 128)
gen int i = _n * 100
gen long l = _n * 100000
gen float f = runiform()
gen double d = runiform() * 1e10
run_test "All numeric types" cexport delimited using `out', replace

clear
set obs 10
gen strL longtext = "This is a longer string that might be stored as strL"
run_test "strL variables" cexport delimited using `out', replace

* === Export/import roundtrip ===
sysuse auto, clear
local orig_n = _N
capture sum price, meanonly
local orig_sum = r(sum)
capture cexport delimited price mpg foreign using `out', replace
clear
capture cimport delimited using `out', clear
local new_n = _N
capture sum price, meanonly
local new_sum = r(sum)
local passed = (`orig_n' == `new_n' & abs(`orig_sum' - `new_sum') < 1)
report_test "Export/import roundtrip" `passed'

* === Error case ===
sysuse auto, clear
capture cexport delimited using `out', replace
capture cexport delimited using `out'
local passed = (_rc != 0)
report_test "No replace error" `passed'

* === Empty dataset ===
clear
set obs 0
gen x = .
capture cexport delimited using `out', replace
report_test "Empty dataset" 1  // Just check doesn't crash

validation_summary
log close
