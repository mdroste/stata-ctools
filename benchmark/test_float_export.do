* Test float export precision
clear all
set obs 10

* Create float variable with various test values
gen float testfloat = .
replace testfloat = 0.1 in 1
replace testfloat = 0.123456789 in 2
replace testfloat = 1.5 in 3
replace testfloat = 123.456 in 4
replace testfloat = 0.0001 in 5
replace testfloat = 0.00001 in 6
replace testfloat = 1234567 in 7
replace testfloat = 0.333333333 in 8
replace testfloat = 2.718281828 in 9
replace testfloat = 3.141592654 in 10

* Export with native Stata
export delimited using "benchmark/native_export.csv", replace
di "Native export:"
type "benchmark/native_export.csv"

* Export with cexport
cexport delimited using "benchmark/cexport_export.csv", replace
di ""
di "cexport export:"
type "benchmark/cexport_export.csv"

* Show difference
di ""
di "Diff:"
shell diff benchmark/native_export.csv benchmark/cexport_export.csv
