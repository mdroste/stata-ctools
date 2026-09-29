clear all
capture log close _all
log using "temp/io_options/corrupt-probe.log", text replace
do "validation/io_test_helpers.do"
cio_import_state_test using "temp/io_options/corrupt/sas_invalid.sas7bdat", kind(sas) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sas_truncated.sas7bdat", kind(sas) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/spss_invalid.sav", kind(spss) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/spss_truncated.sav", kind(spss) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sasxport5_invalid.xpt", kind(sasxport5) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sasxport5_truncated.xpt", kind(sasxport5) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sasxport8_invalid.v8xpt", kind(sasxport8) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sasxport8_truncated.v8xpt", kind(sasxport8) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/dbase_invalid.dbf", kind(dbase) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/dbase_truncated.dbf", kind(dbase) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/xls_invalid.xls", kind(excel) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/xls_truncated.xls", kind(excel) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/xlsx_invalid.xlsx", kind(excel) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/xlsx_truncated.xlsx", kind(excel) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/shp_invalid.shp", kind(shp) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/shp_truncated.shp", kind(shp) state(changed) opts(clear)
log close
exit, clear
