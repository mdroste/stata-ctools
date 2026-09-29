clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/cells.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k256_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k512_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k64_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k512_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k128_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/k512_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/n10000_k200_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/n10000_k2000_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_tile_sweep_full/b1_cells/n50000_k128_cycle.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
