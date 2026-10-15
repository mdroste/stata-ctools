clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/winsor_1m/c1.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/winsor_1m/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/combined" c1 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/winsor_1m" "w20 w200 w1000 w_single w_wide rs1000" 4 "1 8" 1000000 "2"
local rc=_rc
di "GROUP_BENCHMARK_RC=`rc'"
log close
exit, clear
