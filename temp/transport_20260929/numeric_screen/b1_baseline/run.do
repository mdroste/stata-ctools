clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/baseline.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k32_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k32_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k32_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k64_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k64_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k128_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k128_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k256_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k256_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n200k_k256_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n50000_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n100000_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n1000000_k64_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n50000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n100000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/n1000000_k128_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/sparse_wide.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/dense_wide.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/subset8_host256.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/wide_t1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/numeric_screen/b1_baseline/wide_t12.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
