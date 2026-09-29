clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/toggle.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k64_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k64_int.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k64_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k64_float.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k64_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k128_int.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k128_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k128_float.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k128_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k256_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k256_int.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k256_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k256_float.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n200000_k256_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n50000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_scheduler_toggle/b1/n10000_k2000_cycle.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
