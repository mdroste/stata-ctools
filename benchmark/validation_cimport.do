*! Validation tests for cimport - all options

clear all
set more off
capture log close _all
do "benchmark/validation_utils.do"
log using "benchmark/validation_cimport.log", replace text
validation_header "cimport delimited"

* === Create test CSV files ===
tempfile basic_csv noheader_csv tab_csv quoted_csv missing_csv large_csv decimal_csv

file open fh using `basic_csv', write text replace
file write fh "id,name,value,score" _n "1,Alice,100.5,85" _n "2,Bob,200.25,90" _n
file write fh "3,Charlie,150.75,78" _n "4,Diana,300.0,95" _n "5,Eve,175.5,88" _n
file close fh

file open fh using `noheader_csv', write text replace
file write fh "1,Alice,100.5" _n "2,Bob,200.25" _n "3,Charlie,150.75" _n
file close fh

file open fh using `tab_csv', write text replace
file write fh "id" _tab "name" _tab "value" _n "1" _tab "Alice" _tab "100.5" _n
file close fh

file open fh using `quoted_csv', write text replace
file write fh `"id,name,desc"' _n `"1,"John","A ""quoted"" string""' _n
file close fh

file open fh using `missing_csv', write text replace
file write fh "id,val1,val2" _n "1,100," _n "2,,200" _n "3,300,400" _n
file close fh

file open fh using `large_csv', write text replace
file write fh "id,x,y,z" _n
forvalues i = 1/1000 {
    file write fh "`i',`=runiform()',`=rnormal()',`=runiformint(1,100)'" _n
}
file close fh

file open fh using `decimal_csv', write text replace
file write fh "id;value" _n "1;100,50" _n "2;200,75" _n
file close fh

* === Basic import ===
clear
capture cimport delimited using `basic_csv', clear
local passed = (_rc == 0 & _N == 5)
report_test "Basic import" `passed'

* === varnames options ===
clear
capture cimport delimited using `basic_csv', clear varnames(1)
capture confirm variable id
local passed = (_rc == 0)
report_test "varnames(1)" `passed'

clear
capture cimport delimited using `noheader_csv', clear varnames(nonames)
capture confirm variable v1
local passed = (_rc == 0)
report_test "varnames(nonames)" `passed'

* === delimiters ===
clear
run_test "delimiters(tab)" cimport delimited using `tab_csv', clear delimiters(tab)

clear
run_test "delimiters(;)" cimport delimited using `decimal_csv', clear delimiters(";")

* === case options ===
foreach c in preserve lower upper {
    clear
    run_test "case(`c')" cimport delimited using `basic_csv', clear case(`c')
}

* === quote handling ===
clear
run_test "bindquotes(strict)" cimport delimited using `quoted_csv', clear bindquotes(strict)

clear
run_test "bindquotes(loose)" cimport delimited using `quoted_csv', clear bindquotes(loose)

clear
run_test "stripquotes" cimport delimited using `quoted_csv', clear stripquotes

* === rowrange ===
clear
capture cimport delimited using `basic_csv', clear rowrange(2:)
local passed = (_rc == 0 & _N == 4)
report_test "rowrange(2:)" `passed'

clear
run_test "rowrange(:3)" cimport delimited using `basic_csv', clear rowrange(:3)

clear
run_test "rowrange(2:4)" cimport delimited using `basic_csv', clear rowrange(2:4)

* === numeric options ===
clear
run_test "asfloat" cimport delimited using `basic_csv', clear asfloat

clear
run_test "asdouble" cimport delimited using `basic_csv', clear asdouble

clear
run_test "numericcols(1 3 4)" cimport delimited using `basic_csv', clear numericcols(1 3 4)

clear
run_test "stringcols(2)" cimport delimited using `basic_csv', clear stringcols(2)

* === other options ===
clear
run_test "verbose" cimport delimited using `basic_csv', clear verbose

clear
run_test "threads(2)" cimport delimited using `basic_csv', clear threads(2)

clear
run_test "decimalseparator(,)" cimport delimited using `decimal_csv', clear delimiters(";") decimalseparator(",")

clear
run_test "emptylines(skip)" cimport delimited using `basic_csv', clear emptylines(skip)

clear
run_test "emptylines(fill)" cimport delimited using `basic_csv', clear emptylines(fill)

clear
run_test "groupseparator(,)" cimport delimited using `basic_csv', clear groupseparator(",")

* === Large file ===
clear
capture cimport delimited using `large_csv', clear
local passed = (_rc == 0 & _N == 1000)
report_test "Large file (1000 rows)" `passed'

* === Missing values ===
clear
capture cimport delimited using `missing_csv', clear
capture count if missing(val1) | missing(val2)
local passed = (_rc == 0 & r(N) >= 2)
report_test "Missing values" `passed'

* === Combined options ===
clear
run_test "Combined options" cimport delimited using `basic_csv', clear varnames(1) case(lower) asdouble verbose threads(2)

* === Error case ===
clear
capture cimport delimited using "nonexistent_file.csv", clear
local passed = (_rc != 0)
report_test "File not found error" `passed'

* === Data integrity ===
clear
capture cimport delimited using `basic_csv', clear
local passed = (_rc == 0 & _N == 5 & id[1] == 1)
report_test "Data integrity" `passed'

validation_summary
log close
