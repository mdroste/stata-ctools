clear all
set more off
set linesize 255
adopath ++ "build"
log using "temp/io_formats/optionprobe.log", text replace
foreach type in 28 31 {
 import shp "temp/io_formats/shape`type'.shp", clear
 capture noisily export shp "temp/io_formats/native_shape`type'.shp", replace shx
}
clear
set obs 3
gen byte mylongvariablename=_n
gen double mylongvariabletwo=_n+1.125
gen byte Abc=_n
gen byte ABC=_n
label define mylongvaluelabel 1 "Yes" 2 "No" 3 "Other"
label values mylongvariablename mylongvaluelabel
export sasxport5 * using "temp/io_formats/rename5.xpt", replace rename vallabfile(both)
import sasxport5 "temp/io_formats/rename5.xpt", clear
label dir
describe
label list
use "temp/io_formats/source.dta", clear
drop stamp
foreach fmt in %10.2f %10.2e %10.0g %12.2fc {
 format precise `fmt'
 capture noisily export dbase "temp/io_formats/native_`=substr("`fmt'",2,.)'.dbf", replace datafmt
 di "DBF `fmt' RC " _rc
}
use "temp/io_formats/source.dta", clear
drop stamp
export dbase "temp/io_formats/origdate.dbf", replace origdbfdate
cexport sasxport5 * using "temp/io_formats/custom5.xpt", replace
import sasxport5 "temp/io_formats/custom5.xpt", clear
describe
label list
import dbase "temp/io_formats/native.dbf", clear
char list
log close
exit, clear
