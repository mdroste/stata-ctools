clear all
log using "temp/io_options/excel_formats_probe.log", text replace
capture noisily do "validation/io_test_helpers.do"
foreach opt in "firstrow" "firstrow allstring" "firstrow allstring(%9.2f)" {
    cio_import_test using "temp/io_excel_options/display_formats.xlsx", kind(excel) name("XLSX all date/display formats `opt'") opts(`opt')
}
di "FORMAT_PROBE_PASSED=$IO_FORMATS_PASSED FAILED=$IO_FORMATS_FAILED"
log close
exit, clear
