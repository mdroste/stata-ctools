* Debug webuse nlswork comparison

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* Load and export
di "=== Load webuse nlswork ==="
webuse nlswork, clear
desc, short
export delimited using "temp/nlswork.csv", replace

* Import with Stata
di _n "=== Import with Stata (varnames(1)) ==="
import delimited using "temp/nlswork.csv", varnames(1) clear
desc, short
list in 1/3

* Import with cimport
di _n "=== Import with cimport ==="
cimport delimited using "temp/nlswork.csv", clear
desc, short
list in 1/3

* Detailed comparison
di _n "=== Detailed comparison ==="
import delimited using "temp/nlswork.csv", varnames(1) clear
tempfile stata_data
save `stata_data', replace

cimport delimited using "temp/nlswork.csv", clear
tempfile cimport_data
save `cimport_data', replace

* Check types
di _n "Stata variable types:"
use `stata_data', clear
desc

di _n "cimport variable types:"
use `cimport_data', clear
desc

* Try cf with verbose
use `stata_data', clear
cf _all using `cimport_data', verbose
