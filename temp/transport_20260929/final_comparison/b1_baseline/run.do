clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/baseline.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/tiny_n1_k2.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/tiny_n100_k20.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/tiny_n1000_k20.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/tiny_mixed_n100_k200.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/tiny_strings_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/tiny_mixed_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_million.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_k2000.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/large_sorted_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/mixed_w8.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/mixed_w32.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/mixed_w244.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/maxstr_full.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/maxstr_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/maxstr_filtered.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/numeric_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/filtered_dense.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/filtered_sparse.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/unknown_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/strl_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/checked_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/empty.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_t1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_comparison/b1_baseline/wide_t12.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
