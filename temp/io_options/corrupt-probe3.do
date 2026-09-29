clear all
capture log close _all
log using "temp/io_options/corrupt-probe3.log", text replace
do "validation/io_test_helpers.do"
cio_import_state_test using "temp/io_allformats/result/fixture.sas7bdat", kind(sas) state(changed) opts(clear encoding(bogus))
cio_import_state_test using "temp/io_allformats/native/fixture.sav", kind(spss) state(changed) opts(clear encoding(bogus))
cio_import_state_test using "temp/io_allformats/corrupt.sav", kind(spss) state(changed) opts(clear)
cio_import_state_test using "temp/io_allformats/corrupt.dbf", kind(dbase) state(changed) opts(clear)
cio_import_state_test using "temp/io_allformats/corrupt.shp", kind(shp) state(changed) opts(clear)
cio_import_state_test using "temp/io_allformats/corrupt.xls", kind(excel) state(changed) opts(clear)
log close
exit, clear
