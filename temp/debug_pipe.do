* Debug pipe delimiter type mismatch

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

di "File contents:"
type "temp/pipe_delim.csv"

di _n "=== Stata import ==="
import delimited using "temp/pipe_delim.csv", delimiters("|") varnames(1) clear
desc
list

di _n "=== cimport ==="
cimport delimited using "temp/pipe_delim.csv", delimiters("|") clear
desc
list
