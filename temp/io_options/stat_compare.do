capture log close _all
log using "temp/io_options/stat_compare.log", text replace
do "validation/io_test_helpers.do"
foreach pair in "spss sav" "sas sas7bdat" {
 gettoken kind ext : pair
 local ext=strtrim("`ext'")
 cio_import_test using "temp/io_allformats/formats.`ext'", kind(`kind') name("`kind' formats and user missing")
}
log close
exit, clear
