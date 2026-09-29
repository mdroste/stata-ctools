clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/interleave.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k8_numeric8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k8_numeric8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k8_mixed8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k8_mixed8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k8_mixed32_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k8_mixed32_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k13_numeric8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k13_numeric8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k13_mixed8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k13_mixed8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k13_mixed32_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k13_mixed32_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k64_numeric8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k64_numeric8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k64_mixed8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k64_mixed8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k64_mixed32_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n200000_k64_mixed32_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k20_numeric8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k20_numeric8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k20_mixed8_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k20_mixed8_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k20_mixed32_identity.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/store_screen/b0_interleave/n1000000_k20_mixed32_sorted.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
