clear all
set more off
log using "temp/io_formats/nativeparse.log", text replace
capture noisily import spss using temp/io_allformats/native/fixture.sav, clear
di "USING UNQUOTED " _rc
capture noisily import spss using "temp/io_allformats/native/fixture.sav", clear
di "USING QUOTED " _rc
capture noisily import spss using `"temp/io_allformats/native/fixture.sav"', clear
di "USING COMPOUND " _rc
capture noisily import spss "temp/io_allformats/native/fixture.sav", clear
di "POSITIONAL " _rc
import dbase "temp/io_formats/native.dbf", clear
char list
capture noisily export dbase "temp/io_formats/version4.dbf", replace version(iv)
di "VERSION IV " _rc
clear
set obs 2
gen strL text=2100*"a"
capture noisily export spss "temp/io_formats/long.sav", replace
di "LONG SPSS " _rc
capture noisily export sasxport8 "temp/io_formats/long.v8xpt", replace
di "LONG XP8 " _rc
log close
exit, clear
