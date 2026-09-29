clear all
set more off
adopath ++ "build"
adopath ++ "temp/io_parity/build"
capture log close _all
log using "temp/io_parity/probe.log", text replace
clear
set obs 18
gen double d = .
local i = 0
foreach n in 0.1 0.1234567890123456 1.234567890123456 123.4567890123456 0.0001 0.00001 1e16 1e20 1e-20 -0.1 1000000000 3.141592653589793 . .a .z 1e100 1e-100 9.999999999999999 {
    local ++i
    replace d = `n' in `i'
}
gen float f = d
gen str30 s = "hello world"
replace s = "" in 1
replace s = " abc " in 2
replace s = "a,b" in 3
replace s = char(34) + "hi" + char(34) in 4
format d %20.5f
export delimited using "temp/io_parity/native.csv", replace
cexport delimited using "temp/io_parity/custom.csv", replace
export delimited using "temp/io_parity/native_fmt.csv", replace datafmt
cexport delimited using "temp/io_parity/custom_fmt.csv", replace datafmt
foreach cmd in import cimport {
    `cmd' delimited using "temp/io_parity/native.csv", clear
    describe
    list in 1/4
    return list
}
di "PROBE_DONE"
log close
exit, clear
