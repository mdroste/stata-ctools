clear all
set more off
set obs 3
gen x=_n
gen str8 s="old"
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
capture log close
log using "`root'/excel_probe.log", text replace
export excel using "`root'/native.xlsx", firstrow(variables) sheet(Original) replace
replace x=10+_n
replace s="new"
foreach opt in "sheet(New)" "sheet(Original)" "sheet(Original, modify)" "sheet(Original, replace)" "sheet(Other, modify)" "sheet(Other, replace)" "sheet(Original, modify) keepcellfmt" "sheet(Original) sheetmodify" "sheet(Original) sheetreplace" "sheet(Original, modify) replace" "keepcellfmt" "sheet(Original, replace) keepcellfmt" {
    copy "`root'/native.xlsx" "`root'/probe.xlsx", replace
    capture noisily export excel using "`root'/probe.xlsx", `opt' cell(B2)
    di "OPTIONS=`opt' RC=" _rc
}
foreach opt in "missing(99)" "missing(NA)" "locale(en_US)" "locale(xxx)" "firstrow(var)" "firstrow(varl)" {
    capture noisily export excel using "`root'/new", `opt' replace
    di "OPTIONS=`opt' RC=" _rc
}
log close
exit, clear
