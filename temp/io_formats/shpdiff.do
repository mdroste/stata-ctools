clear all
set more off
adopath ++ "build"
log using "temp/io_formats/shpdiff.log", text replace
foreach shape in shape28 shape31 {
 import shp "temp/io_formats/`shape'.shp", clear
 mata: a=st_sdata(.,"rec_header"); b=st_data(.,.)
 list, noobs
 mata: b
 cimport shp "temp/io_formats/`shape'.shp", clear
 list, noobs
 mata: st_data(.,.)
 mata: printf("NUM %g HDR %g LEN %g %g\n",all(b:==st_data(.,.)),all(a:==st_sdata(.,"rec_header")),strlen(a[1]),strlen(st_sdata(1,"rec_header")))
 mata: printf("NATIVE %s\nC %s\n",invtokens(strofreal(ascii(a[1]),"%03.0f")),invtokens(strofreal(ascii(st_sdata(1,"rec_header")),"%03.0f")))
}
use "temp/io_formats/source.dta", clear
export excel using "temp/io_formats/native.xls", firstrow(variables) replace
import excel using "temp/io_formats/native.xls", firstrow clear
foreach v of varlist _all {
 local fmt : format `v'
 di "FMT `v' `fmt'"
}
log close
exit, clear
