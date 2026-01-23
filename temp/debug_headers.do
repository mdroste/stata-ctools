* Debug header detection behavior

clear all
cd "/Users/Mike/Documents/GitHub/stata-ctools"
adopath ++ "build"

* ==== TEST A: Valid header names (timestamp,level,message) ====
di _n "=== TEST A: Valid string headers ==="
file open fh using "temp/test_valid_headers.csv", write replace
file write fh "timestamp,level,message" _n
file write fh "2023-01-15,INFO,Started" _n
file write fh "2023-01-16,WARN,Warning" _n
file close fh

di "File contents:"
type "temp/test_valid_headers.csv"

di _n "Stata import delimited:"
import delimited using "temp/test_valid_headers.csv", clear
desc, short
list

di _n "Stata import delimited with varnames(1):"
import delimited using "temp/test_valid_headers.csv", varnames(1) clear
desc, short
list

di _n "cimport (default):"
cimport delimited using "temp/test_valid_headers.csv", clear
desc, short
list

* ==== TEST B: All numeric first row ====
di _n "=== TEST B: All numeric first row ==="
file open fh using "temp/test_numeric_row1.csv", write replace
file write fh "1,2,3" _n
file write fh "10,20,30" _n
file write fh "11,21,31" _n
file close fh

di "File contents:"
type "temp/test_numeric_row1.csv"

di _n "Stata import delimited:"
import delimited using "temp/test_numeric_row1.csv", clear
desc, short
list

di _n "Stata import delimited with varnames(1):"
import delimited using "temp/test_numeric_row1.csv", varnames(1) clear
desc, short
list

di _n "cimport (default):"
cimport delimited using "temp/test_numeric_row1.csv", clear
desc, short
list

* ==== TEST C: Mixed first row - some valid varnames, some numeric ====
di _n "=== TEST C: Mixed first row ==="
file open fh using "temp/test_mixed_row1.csv", write replace
file write fh "id,2,value" _n
file write fh "1,A,100" _n
file write fh "2,B,200" _n
file close fh

di "File contents:"
type "temp/test_mixed_row1.csv"

di _n "Stata import delimited:"
import delimited using "temp/test_mixed_row1.csv", clear
desc, short
list

di _n "cimport (default):"
cimport delimited using "temp/test_mixed_row1.csv", clear
desc, short
list

* ==== TEST D: Check how benchmark_import calls both ====
* Check what the validation test actually does
di _n "=== TEST D: How validation script calls them ==="
type "temp/log_format.csv"
di _n "Stata:"
import delimited using "temp/log_format.csv", clear
desc, short
di _n "cimport:"
cimport delimited using "temp/log_format.csv", clear
desc, short
