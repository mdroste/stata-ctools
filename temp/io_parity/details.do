clear all
set more off
adopath ++ "build"
adopath ++ "temp/io_parity/build"
capture log close _all
log using "temp/io_parity/details.log", text replace
foreach fmt in %18.0g %17.0g %16.0g %12.0g %10.0g %9.0g {
 foreach x in 0.00001 0.000001 1e16 1e100 1.234567890123456 {
 di "`fmt' `x' = " string(`x',"`fmt'")
 }
}
foreach file in headers noheaders numbers {
 clear
 import delimited using "temp/io_parity/`file'.csv", clear encoding(utf-8)
 describe
 list
 foreach v of varlist _all {
 local lab : variable label `v'
 di "LABEL `v': `lab'"
 }
}
di "DETAILS_DONE"
log close
exit, clear
