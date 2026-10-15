clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_encode/b1.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_encode/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_baseline" b1 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/integrated_encode" "e1000 e10000 e50000 e_long e_existing" 4 "8" 1000000 "2"
local rc=_rc
di "ENCODE_BENCHMARK_RC=`rc'"
log close
exit, clear
