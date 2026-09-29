clear all
set more off
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/build"
capture log close
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/csv_probe.log", text replace
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
foreach opt in `"delimiters("|;",collapse)"' `"delimiters("||",asstring)"' `"delimiters("|;")"' `"delimiters("|;") collapsedelimiters"' `"delimiters("|;",asstring collapse)"' {
    capture noisily import delimited using "`root'/multidelim.txt", varnames(nonames) `opt' clear
    local rc=_rc
    di `"OPT=`macval(opt)' RC=`rc'"'
    list
    describe
}
log close
exit, clear
