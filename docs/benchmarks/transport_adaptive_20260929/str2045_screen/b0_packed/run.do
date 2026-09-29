clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/packed.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l8_k1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l8_k2.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l8_k4.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l64_k1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l64_k2.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l64_k4.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l512_k1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l512_k2.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l512_k4.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l2045_k1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l2045_k2.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_l2045_k4.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_empty_k4.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_full_t1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/s2045_full_filtered.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/str2045_screen/b0_packed/strl_nochange.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
