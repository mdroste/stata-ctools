clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/toggle.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n1000000_k1_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n1000000_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n1000000_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n10000000_k1_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n10000000_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n10000000_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n1000000_k20_mixed8.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_order_toggle/b0/n1000000_k2_long8.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
