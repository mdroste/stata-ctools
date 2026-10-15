clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_3m/c2.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_3m/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/combined" c2 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_3m" "e1000 e10000 e_filtered e_existing" 4 "8" 3000000 "2"
local rc=_rc
di "ENCODE_BENCHMARK_RC=`rc'"
log close
exit, clear
