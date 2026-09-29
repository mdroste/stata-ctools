clear all
set more off
capture clear
capture mata: st_addvar("byte", "_all")
di "ADDALL " _rc
describe
log close _all
log using "temp/io_parity/headerprobe.log", text replace
foreach file in stringdata numericheader mixedheader {
 import delimited using "temp/io_parity/`file'.csv", clear encoding(utf-8)
 describe
 list
}
foreach opts in asfloat asdouble {
 import delimited using "temp/io_parity/numbers.csv", clear encoding(utf-8) `opts'
 describe
}
clear
capture mata: st_addvar("byte", "_all")
di "ADDALL " _rc
describe
log close
exit, clear
