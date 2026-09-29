clear all
set more off
set linesize 255
capture log close _all
log using "temp/io_options/state-probe.log", text replace
adopath ++ "build"
foreach pair in "delimited temp/io_delimited_options/range.csv" "excel temp/io_allformats/native/fixture.xlsx" "excel temp/io_allformats/native/fixture.xls" "spss temp/io_allformats/native/fixture.sav" "sas temp/io_allformats/result/fixture.sas7bdat" "sasxport5 temp/io_allformats/native/fixture.xpt" "sasxport8 temp/io_allformats/native/fixture.v8xpt" "dbase temp/io_allformats/native/fixture.dbf" "shp temp/io_allformats/shape1.shp" {
    gettoken kind path : pair
    local path=strtrim("`path'")
    local using using
    if "`kind'"=="sasxport5" local using ""
    foreach state in saved changed zeroobs zerovars empty {
        foreach prefix in "" "c" {
            clear
            if "`state'"!="empty" {
                quietly set obs 2
                quietly gen long sentinel=99
                save "temp/io_options/state-source.dta", replace
                if "`state'"=="changed" quietly replace sentinel=100 in 1
                if "`state'"=="zeroobs" quietly drop in 1/2
                if "`state'"=="zerovars" drop sentinel
            }
            local before=c(changed)
            capture noisily `prefix'import `kind' `using' "`path'"
            local rc=_rc
            di "STATE: `prefix'import `kind' `state' before=`before' rc=`rc' N=" _N " k=" c(k) " changed=" c(changed)
        }
    }
}
clear
set obs 2
gen byte sentinel=99
capture noisily import sasxport8 using "temp/io_allformats/native/fixture.v8xpt", clear vlabfile("temp/io_allformats/native/fixture.v8xpt")
di "COMPANION: rc=" _rc " changed=" c(changed)
describe
list
return list
log close
exit, clear
