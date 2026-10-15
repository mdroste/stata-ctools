clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_winsor/c1.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_winsor/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_candidate" c1 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_winsor" "w200 w1000 w_wide w_single rs1000" 4 "8" 1000000 "2"
local rc=_rc
di "GROUP_BENCHMARK_RC=`rc'"
log close
exit, clear
