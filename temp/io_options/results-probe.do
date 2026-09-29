clear all
set more off
set linesize 255
adopath ++ "build"
capture log close _all
log using "temp/io_options/results-probe.log", text replace
foreach pair in "sas temp/io_allformats/result/fixture.sas7bdat" "spss temp/io_allformats/native/fixture.sav" "sasxport5 temp/io_allformats/native/fixture.xpt" "sasxport8 temp/io_allformats/native/fixture.v8xpt" "dbase temp/io_allformats/native/fixture.dbf" "excel temp/io_excel_options/numbers.xls" "excel temp/io_excel_options/numbers.xlsx" {
    gettoken kind path : pair
    local path=strtrim("`path'")
    foreach command in import cimport {
        clear
        local using using
        if "`kind'"=="sasxport5" & "`command'"=="import" local using ""
        capture noisily `command' `kind' `using' "`path'", clear
        di "RESULTS: `command' `kind' `path' rc=" _rc
        return list
    }
}
log close
exit, clear
