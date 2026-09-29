clear all
set more off
adopath ++ "build"
log using "temp/io_formats/shptrace.log", text replace
set trace on
capture noisily cimport shp "temp/io_formats/points.shp", clear
set trace off
list
char list
mata: st_sdata(.,"rec_header")
log close
exit, clear
