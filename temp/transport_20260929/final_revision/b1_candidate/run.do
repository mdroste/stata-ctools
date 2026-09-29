clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision"
program transport_io, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/candidate.plugin")
program run_campaign
version 16
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/tiny_n100_k20.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/wide_byte.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/wide_byte128.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/wide_cycle.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/wide_million.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/wide_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/large_sorted_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/mixed_w32.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/maxstr_full.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/maxstr_short.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/numeric_control.do"
do "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/final_revision/b1_candidate/filtered_dense.do"
end
capture noisily run_campaign
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
