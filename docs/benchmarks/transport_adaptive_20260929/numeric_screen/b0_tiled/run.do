clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/tiled.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k32_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k32_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k32_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k64_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k64_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k128_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k256_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k256_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n200k_k256_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n50000_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n100000_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n1000000_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n50000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n100000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/n1000000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/sparse_wide.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/dense_wide.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/subset8_host256.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/wide_t1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b0_tiled/wide_t12.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
