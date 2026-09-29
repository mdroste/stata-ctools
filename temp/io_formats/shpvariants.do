clear all
set more off
set linesize 255
log using "temp/io_formats/shpvariants.log", text replace
foreach kind in 11 21 13 23 15 25 18 28 31 {
    di "SHAPE_KIND=`kind'"
    capture noisily import shp "temp/io_formats/shape`kind'.shp", clear
    di "SHAPE_RC=" _rc
    if !_rc {
        describe
        capture list _ID _X _Y shape_order, noobs
        foreach v of varlist _all {
            local ty: type `v'
            local fm: format `v'
            di "SHAPE_META `v' `ty' `fm'"
        }
        mata: st_varname(1..st_nvar())
        capture list _Z _M, noobs
        mata: strlen(st_sdata(.,"rec_header"))
    }
}
log close
exit, clear
