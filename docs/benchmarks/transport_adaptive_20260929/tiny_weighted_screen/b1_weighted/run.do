clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/weighted.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k2_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k20_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k200_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k200_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1_k200_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k2_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k20_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k200_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k200_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n100_k200_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k2_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k20_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k200_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k200_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n1000_k200_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k2_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k20_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k200_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k200_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/n10000_k200_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n32767_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n32767_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n32769_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n32769_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n3275_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n3275_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n3277_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n3277_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n326_k200_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n326_k200_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n328_k200_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/tiny_weighted_screen/b1_weighted/boundary_n328_k200_mixed.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
