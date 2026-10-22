clear all
set more off
set linesize 255
adopath ++ "/Users/mike/Documents/GitHub/stata-ctools/package/ado"
adopath ++ "/tmp/ctools-audit-oct22/build"
capture log close _all
log using "/tmp/ctools-data-audit/probes.log", text replace
set obs 4
gen x = _n
capture noisily cwinsor x, by(x) cuts(10 90) suffix(_w)
di "WINSOR_OVERLAP_RC=" _rc
list
clear
set obs 2
gen str10 s = "12©"
capture noisily destring s, generate(native) ignore("é¢")
di "NATIVE_UNICODE_RC=" _rc
capture confirm variable native
di "NATIVE_UNICODE_EXISTS=" (_rc == 0)
capture noisily cdestring s, generate(fast) ignore("é¢")
di "FAST_UNICODE_RC=" _rc
list
clear
set obs 2
gen str10 x = string(_n)
gen keep = 1
label variable x "Original variable label"
char x[custom] "custom metadata"
clonevar y = x
capture noisily destring x, replace
capture noisily cdestring y, replace
di "NATIVE_LABEL=|`: variable label x'| NATIVE_CHAR=|`: char x[custom]'|"
di "FAST_LABEL=|`: variable label y'| FAST_CHAR=|`: char y[custom]'|"
clear
set obs 5
gen double x = _n
gen double y = x
replace y=. in 2/4
capture noisily cipolate y x, generate(ya) by(x)
di "IPOLATE_BY_X_RC=" _rc
clear
set obs 2
gen byte wt = 0
gen s = "weight(wt)"
capture noisily cbsample if s == "weight(wt)"
di "BSAMPLE_STRING_IF_RC=" _rc
list
capture noisily cbsample, weight( wt )
di "BSAMPLE_WEIGHT_SPACES_RC=" _rc
log close
exit, clear
