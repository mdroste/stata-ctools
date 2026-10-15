clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_1m/b1.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_1m/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/baseline" b1 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_1m" "" 4 "1 8" 1000000 "2"
local rc=_rc
di "ENCODE_BENCHMARK_RC=`rc'"
log close
exit, clear
