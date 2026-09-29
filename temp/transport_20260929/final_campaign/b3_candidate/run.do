clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/candidate.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/tiny_n1_k2.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/tiny_n100_k20.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/tiny_n1000_k20.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/tiny_mixed_n100_k200.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/tiny_strings_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/tiny_mixed_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_double.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_million.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_k2000.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_sorted.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/large_sorted_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/mixed_w8.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/mixed_w32.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/mixed_w244.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/maxstr_full.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/maxstr_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/maxstr_filtered.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/numeric_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/filtered_dense.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/filtered_sparse.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/unknown_strings.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/strl_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/checked_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/empty.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_t1.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_campaign/b3_candidate/wide_t12.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
