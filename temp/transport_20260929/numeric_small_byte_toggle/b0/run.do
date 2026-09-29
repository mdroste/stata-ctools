clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/toggle.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/n15625_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/n50000_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/n100000_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/n4096_k512_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/n50000_k256_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_small_byte_toggle/b0/n50000_k128_double.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
