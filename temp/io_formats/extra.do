clear all
set more off
adopath ++ "build"
log using "temp/io_formats/extra.log", text replace
di "FRED_KEY_CONFIGURED=" (strlen(c(fredkey))>0)
use "temp/io_formats/source.dta", clear
drop stamp
cexport dbase using "temp/io_formats/c.dbf", replace
import dbase using "temp/io_formats/c.dbf", clear
describe
list, noobs
cimport dbase using "temp/io_formats/native.dbf", clear
describe
list, noobs
use "temp/io_formats/source.dta", clear
export excel using "temp/io_formats/native.xls", replace firstrow(variables)
clear
import excel using "temp/io_formats/native.xls", firstrow
describe
list, noobs
char list
log close
exit, clear
