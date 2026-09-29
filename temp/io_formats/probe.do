clear all
set more off
set linesize 255
log using "temp/io_formats/probe.log", text replace
set obs 5
gen byte small = _n
gen double precise = 1.123456789012345 + _n
replace precise = .a in 4
replace precise = .z in 5
gen str12 text = "abc"
replace text = "" in 2
replace text = " abc " in 3
gen double day = td(01jan2020) + _n
format day %td
gen double stamp = clock("01jan2020 01:02:03", "DMYhms") + _n
format stamp %tc
label variable precise "Precision label"
label define yesno 1 "Yes" 2 "No" 3 "Other"
label values small yesno
label data "Test data label"
save "temp/io_formats/source.dta", replace
foreach kind in sas spss sasxport5 sasxport8 dbase {
    use "temp/io_formats/source.dta", clear
    di "EXPORT_KIND=`kind'"
    capture noisily export `kind' using "temp/io_formats/native_`kind'", replace
    local rc = _rc
    di "EXPORT_RC=`rc'"
    if !`rc' {
        clear
        capture noisily import `kind' using "temp/io_formats/native_`kind'"
        di "IMPORT_RC=" _rc
        describe
        list, noobs
        return list
        save "temp/io_formats/imported_`kind'.dta", replace
    }
}
log close
exit, clear
