clear all
set more off
set linesize 255
adopath ++ "build"
log using "temp/io_formats/dbfprobe.log", text replace
use "temp/io_formats/source.dta", clear
drop stamp
export dbase using "temp/io_formats/native.dbf", replace
import dbase using "temp/io_formats/native.dbf", clear
describe
list, noobs
foreach v of varlist _all {
    local fm : format `v'
    local ty : type `v'
    di "DBFMETA `v' `ty' `fm'"
}
export dbase using "temp/io_formats/native4.dbf", version(4) datafmt replace
clear
import sasxport5 "temp/io_formats/formats.xpf", novallabels
describe
list, noobs
log close
exit, clear
