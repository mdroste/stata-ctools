clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/tile64.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k256_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k512_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k64_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k512_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k128_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep/b0_tile64/k512_double.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
