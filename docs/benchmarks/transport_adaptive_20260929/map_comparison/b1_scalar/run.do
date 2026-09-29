clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison/b1_scalar/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison/scalar.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison/b1_scalar/dense_k1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison/b1_scalar/dense_k20.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison/b1_scalar/dense_k128.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/map_comparison/b1_scalar/sparse_k128.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
