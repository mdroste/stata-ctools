clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/profile_final_candidate/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/candidate.plugin")
program main
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/profile_final_candidate/case.do"
end
capture noisily main
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
