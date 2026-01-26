adopath ++ "build"
sysuse auto, clear
keep make price mpg
keep in 1/3
replace price = . in 1
cexport excel using "test_verify.xlsx", cell(B5) missing("N/A") replace
