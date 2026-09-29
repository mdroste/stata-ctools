clear all
set more off
capture log close _all
log using "temp/io_parity/legacy.log", text replace
capture noisily do validation/validate_cimport.do
di "IMPORT_LEGACY_RC=" _rc " FAILED=$TESTS_FAILED"
capture noisily do validation/validate_cexport.do
di "EXPORT_LEGACY_RC=" _rc " FAILED=$TESTS_FAILED"
log close
exit, clear
