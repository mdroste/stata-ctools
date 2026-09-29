clear all
set more off
adopath ++ "build"
log using "temp/io_formats/shpprobe.log", text replace
foreach f in points lines polygon multi {
    import shp "temp/io_formats/`f'.shp", clear
    describe
    list, noobs
    char list
    return list
    export shp "temp/io_formats/native_`f'.shp", shx replace
}
log close
exit, clear
