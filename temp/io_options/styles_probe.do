clear all
set more off
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/build"
capture log close
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/styles_probe.log", text replace
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
set obs 4
gen double x=mdy(1,1,2020)+_n
gen str6 s="old"
format x %td
foreach ext in xlsx xls {
    export excel using "`root'/styled.`ext'", firstrow(variables) replace
}
replace x=_n
format x %9.0g
replace s="new"
replace x=. in 2
replace s="" in 2
foreach action in modify replace {
    foreach keep in "" keepcellfmt {
        copy "`root'/styled.xlsx" "`root'/ns-`action'-`keep'.xlsx", replace
        copy "`root'/styled.xlsx" "`root'/cs-`action'-`keep'.xlsx", replace
        export excel using "`root'/ns-`action'-`keep'.xlsx", sheet(,`action') `keep' cell(A2)
        cexport excel using "`root'/cs-`action'-`keep'.xlsx", sheet(,`action') `keep' cell(A2)
    }
}
foreach opt in "allstring" "allstring(%9.2f)" "allstring(%td)" "allstring(%21x)" "allstring(bad)" {
    import excel using "`root'/styled.xlsx", clear `opt'
    list
    describe
    di "OPT=`opt'"
}
log close
exit, clear
