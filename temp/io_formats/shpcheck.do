clear all
set more off
adopath ++ "build"
log using "temp/io_formats/shpcheck.log", text replace
mata:
string matrix shpsnapshot() {
    string matrix values
    real scalar j
    values=J(st_nobs(),st_nvar(),"")
    for(j=1;j<=st_nvar();j++) {
        if(st_isnumvar(j))values[,j]=strofreal(st_data(.,j),"%21x")
        else values[,j]=st_sdata(.,j)
    }
    return(values)
}
end
foreach shape in points lines polygon multi shape11 shape21 shape13 shape23 shape15 shape25 shape18 shape28 shape31 {
    import shp "temp/io_formats/`shape'.shp", clear
    mata: ref=shpsnapshot()
    capture noisily cimport shp "temp/io_formats/`shape'.shp", clear
    local rc=_rc
    di "SHPCIMPORT `shape' `rc'"
    if !`rc' {
        mata: st_local("same",strofreal(all(ref:==shpsnapshot())))
        di "SHPSAME `shape' `same'"
        describe
        char list
        capture noisily cexport shp "temp/io_formats/c_`shape'.shp", shx replace
        local rc=_rc
        di "SHPCEXPORT `shape' `rc'"
        if !`rc' {
            import shp "temp/io_formats/c_`shape'.shp", clear
            mata: st_local("same",strofreal(all(ref:==shpsnapshot())))
            di "SHPWRITESAME `shape' `same'"
        }
    }
}
log close
exit, clear
