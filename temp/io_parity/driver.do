clear all
capture log close _all
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/io_parity/regressions.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/validation/validate_io_parity.do" ""
local rc = _rc
di "IO_PARITY_DRIVER_RC=`rc'"
log close
exit, clear
