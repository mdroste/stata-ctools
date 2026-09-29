clear all
adopath ++ "build"
log using "temp/io_formats/trace.log", text replace
set trace on
capture noisily cimport spss "temp/io_formats/native_spss.sav"
set trace off
log close
exit, clear
