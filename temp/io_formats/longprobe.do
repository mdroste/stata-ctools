clear all
set more off
adopath ++ "build"
log using "temp/io_formats/longprobe.log", text replace
foreach file in fixture.zsav long.sas7bdat {
 local kind spss
 if "`file'"=="long.sas7bdat" local kind sas
 import `kind' "temp/io_allformats/`file'", clear
 describe
 capture noisily format text %2100s
 di "FORMAT 2100 " _rc
 capture noisily format text %2101s
 di "FORMAT 2101 " _rc
}
import sasxport8 "temp/io_allformats/native/output.v8xpt", clear
describe
capture noisily format text %2101s
di "XP FORMAT " _rc
log close
exit, clear
