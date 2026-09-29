clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation/b0_baseline/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation/baseline.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation/b0_baseline/wide_million.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation/b0_baseline/maxstr_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation/b0_baseline/numeric_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_allocation/b0_baseline/strl_control.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
