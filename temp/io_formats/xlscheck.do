clear all
set more off
adopath ++ "build"
log using "temp/io_formats/xlscheck.log", text replace
cimport excel using "temp/io_formats/native.xls", firstrow
describe
list, noobs
foreach v of varlist _all {
    local fm : format `v'
    di "XLS_META `v' `fm'"
}
use "temp/io_formats/source.dta", clear
cexport excel using "temp/io_formats/custom.xls", replace firstrow(variables)
clear
import excel using "temp/io_formats/custom.xls", firstrow
describe
list, noobs
log close
exit, clear
