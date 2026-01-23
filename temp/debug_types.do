* Debug type inference differences

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* Create test file
file open fh using "temp/pipe_delim.csv", write replace
file write fh "id|name|value" _n
file write fh "1|Alpha|100" _n
file write fh "2|Beta|200" _n
file write fh "3|Gamma|300" _n
file close fh

di "=== Stata import with varnames(1) ==="
import delimited using "temp/pipe_delim.csv", delimiters("|") varnames(1) clear
desc

di _n "=== cimport ==="
cimport delimited using "temp/pipe_delim.csv", delimiters("|") clear
desc

* Now try the comparison
di _n "=== Compare after saving ==="

import delimited using "temp/pipe_delim.csv", delimiters("|") varnames(1) clear
tempfile stata_data
save `stata_data', replace

cimport delimited using "temp/pipe_delim.csv", delimiters("|") clear
tempfile cimport_data
save `cimport_data', replace

di _n "Trying cf _all..."
use `stata_data', clear
capture cf _all using `cimport_data'
di "cf _rc = " _rc
