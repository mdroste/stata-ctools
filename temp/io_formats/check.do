clear all
set more off
set linesize 255
adopath ++ "build"
log using "temp/io_formats/check.log", text replace
mata:
string matrix snapshot() {
    real scalar j
    string matrix v
    v=J(st_nobs(),st_nvar(),"")
    for(j=1;j<=st_nvar();j++) {
        if(st_isnumvar(j)) v[,j]=strofreal(st_data(.,j),"%21x")
        else v[,j]=st_sdata(.,j)
    }
    return(v)
}
end
program testimport
    args kind path
    clear
    import `kind' "`path'"
    mata: ref=snapshot()
    local n=_N
    local k=c(k)
    unab names : _all
    forvalues j=1/`k' {
        local v: word `j' of `names'
        local name`j' "`v'"
        local type`j': type `v'
        local format`j': format `v'
        local label`j': variable label `v'
        local vl`j': value label `v'
        di "NATIVE_META `j' `v' `type`j'' `format`j'' `vl`j''"
    }
    mata: ref
    clear
    capture noisily cimport `kind' "`path'"
    local rc=_rc
    di "CIMPORT_RC `kind' `rc'"
    if `rc' exit
    describe
    mata: got=snapshot(); got
    mata: st_local("same",strofreal(all(ref:==got)))
    di "VALUES_SAME `kind' `same'"
    unab names : _all
    forvalues j=1/`k' {
        local v: word `j' of `names'
        local ty: type `v'
        local fm: format `v'
        local vl: value label `v'
        di "C_META `j' `v' `ty' `fm' `vl'"
    }
end
use "temp/io_formats/source.dta", clear
replace text = strtrim(text)
capture noisily export sasxport5 * using "temp/io_formats/native5", replace
capture noisily testimport sasxport5 temp/io_formats/native5.xpt
testimport spss temp/io_formats/native_spss.sav
testimport sasxport8 temp/io_formats/native_sasxport8.v8xpt
testimport sas temp/io_formats/auto.sas7bdat
foreach kind in sas spss sasxport5 sasxport8 {
    use "temp/io_formats/source.dta", clear
    di "CEXPORT_KIND `kind'"
    capture noisily cexport `kind' * using "temp/io_formats/c`kind'", replace
    di "CEXPORT_RC `kind' " _rc
    capture noisily testimport `kind' temp/io_formats/c`kind'
}
log close
exit, clear
