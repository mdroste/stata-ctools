clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/production_candidate/acceptance.log", text replace
capture noisily do validation/validate_all.do
local rc = _rc
di "TRANSPORT_ACCEPTANCE_COMPLETE RC=`rc'"
log close
exit, clear
