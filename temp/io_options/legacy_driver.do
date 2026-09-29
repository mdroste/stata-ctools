capture log close _all
log using "temp/io_options/legacy-components.log", text replace
capture noisily do "validation/validate_cimport.do"
local irc=_rc
di "LEGACY_IMPORT_RC=`irc' PASSED=$TESTS_PASSED FAILED=$TESTS_FAILED"
capture noisily do "validation/validate_cexport.do"
local erc=_rc
di "LEGACY_EXPORT_RC=`erc' PASSED=$TESTS_PASSED FAILED=$TESTS_FAILED"
log close
exit, clear
