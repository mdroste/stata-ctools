clear all
set more off
capture log close
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options/headerprobe.log", text replace
local root "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_options"
set obs 2
forvalues j=1/8 {
    gen v`j'=_n
}
label var v1 "A"
label var v2 "A"
label var v3 "123"
label var v4 "if"
label var v5 "A B"
label var v6 "C"
label var v7 "ΩABC"
label var v8 "ΩABC"
foreach ext in xls xlsx {
    export excel using "`root'/headers.`ext'", firstrow(varlabels) replace
    import excel using "`root'/headers.`ext'", firstrow clear
    describe
    clear
    set obs 2
    forvalues j=1/8 {
        gen v`j'=_n
    }
    label var v1 "A"
    label var v2 "A"
    label var v3 "123"
    label var v4 "if"
    label var v5 "A B"
    label var v6 "C"
    label var v7 "ΩABC"
    label var v8 "ΩABC"
}
foreach opt in 3 4 iii iv III IV "dBase III" "dBase IV" {
    capture noisily export dbase using "`root'/version.dbf", version(`opt') replace
    di "VERSION=`opt' RC=" _rc
}
foreach opt in "firstrow" "allstring allstring(%9.2f)" "allstring(%9s)" "allstring(%9.2f)" {
    capture noisily import excel using "`root'/numbers.xlsx", clear `opt'
    di "OPT=`opt' RC=" _rc
}
log close
exit, clear
