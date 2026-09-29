clear all
set more off
adopath ++ "build"
log using "temp/io_options/excel_selector_probe2.log", text replace
foreach ext in xls xlsx {
    foreach sel in "empty=C" "empty=D" "empty=E" "empty=F" "x=D y=D" "x=A empty=D" "x=A empty=E" "empty=E x=A" {
        capture noisily import excel `sel' using "temp/io_excel_options/numbers.`ext'", clear
        local nr=_rc
        if !`nr' {
            describe, fullnames
            list in 1/2, noobs
        }
        capture noisily cimport excel `sel' using "temp/io_excel_options/numbers.`ext'", clear
        local cr=_rc
        if !`cr' {
            describe, fullnames
            list in 1/2, noobs
        }
        di "SELECTOR `ext' `sel' NATIVE_RC=`nr' C_RC=`cr'"
    }
}
log close
exit, clear
