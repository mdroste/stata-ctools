clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/baseline.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k1_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k1_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k2_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k20_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k20_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n1000000_k20_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n10000000_k1_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n10000000_k1_long.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n10000000_k2_numeric.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n10000000_k2_mixed.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/no_sort_order_focused/b1_baseline/n10000000_k20_numeric.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
