adopath ++ build

sysuse bplong, clear
export delimited using temp/test_bplong.csv, replace

* First import with Stata
import delimited using temp/test_bplong.csv, clear
save temp/stata_bplong.dta, replace
describe

* Now import with cimport
cimport delimited using temp/test_bplong.csv, clear
describe

* Compare
cf _all using temp/stata_bplong.dta
