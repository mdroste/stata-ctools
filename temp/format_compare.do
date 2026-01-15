capture adopath ++ "build"
clear
set seed 12345
set obs 5
gen float x_float = runiform()
gen double x_double = runiform()

export delimited using "temp/stata_fmt.csv", replace
cexport delimited using "temp/cexport_fmt.csv", replace

di ""
di "=== Stata export delimited ==="
type "temp/stata_fmt.csv"
di ""
di "=== cexport delimited ==="
type "temp/cexport_fmt.csv"

erase "temp/stata_fmt.csv"
erase "temp/cexport_fmt.csv"
