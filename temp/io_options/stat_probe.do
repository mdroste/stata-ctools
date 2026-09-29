capture log close _all
log using "temp/io_options/stat_probe.log", text replace
do "validation/io_test_helpers.do"
foreach pair in "spss sav" "sas sas7bdat" {
 gettoken kind ext : pair
 local ext=strtrim("`ext'")
 cio_import_test using "temp/io_allformats/formats.`ext'", kind(`kind') name("`kind' formats and user missing")
 import `kind' using "temp/io_allformats/formats.`ext'", clear
 describe
 list missing strmiss
}
clear
set obs 100
gen strL text=cond(_n<=80,200*"x","y")
compress text, nocoalesce
local type : type text
di "COMPRESS 200 80 TYPE `type'"
clear
set obs 100
gen strL text=cond(_n<=80,1000*"x","y")
compress text, nocoalesce
local type : type text
di "COMPRESS 1000 80 TYPE `type'"
log close
exit, clear
