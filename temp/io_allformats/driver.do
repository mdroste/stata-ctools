clear all
capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_allformats/regressions.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_io_formats.do" ""
local rc = _rc
di "IO_FORMATS_DRIVER_RC=`rc'"
log close
exit, clear
