clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/dual.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/mixed_w8.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/mixed_w32.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/mixed_w244.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/numeric_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/large_sorted_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/dual_impl_control/b0/wide_byte.do"
end
capture noisily run_campaign
local rc = _rc
di "DUAL_TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
