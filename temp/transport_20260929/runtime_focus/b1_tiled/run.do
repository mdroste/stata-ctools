clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/tiled.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/n200k_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/n200k_k64_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/n200k_k64_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/n200k_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/n200k_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/n200k_k128_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/runtime_focus/b1_tiled/mixed20_str32_control.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
