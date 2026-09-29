clear all
capture log close _all
log using "temp/io_options/corrupt-probe2.log", text replace
do "validation/io_test_helpers.do"
cio_import_state_test using "temp/io_options/corrupt/sas_longinvalid.sas7bdat", kind(sas) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/spss_longinvalid.sav", kind(spss) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sasxport5_longinvalid.xpt", kind(sasxport5) state(changed) opts(clear)
cio_import_state_test using "temp/io_options/corrupt/sasxport8_longinvalid.v8xpt", kind(sasxport8) state(changed) opts(clear)
cio_import_state_test using "temp/io_allformats/result/fixture.sas7bdat", kind(sas) state(changed) opts("clear encoding(bogus)")
cio_import_state_test using "temp/io_allformats/native/fixture.sav", kind(spss) state(changed) opts("clear encoding(bogus)")
log close
exit, clear
