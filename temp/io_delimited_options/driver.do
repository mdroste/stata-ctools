clear all
capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_delimited_options/regressions.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_io_delimited_options.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_delimited_options"
local rc = _rc
di "CSV_OPTIONS_DRIVER_RC=`rc'"
log close
exit, clear
