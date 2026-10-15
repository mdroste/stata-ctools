clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_small/b1.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_small/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/baseline" b1 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/encode_small" "e10 e1000 e_existing" 8 "1 8" 1000 "4"
local rc=_rc
di "ENCODE_BENCHMARK_RC=`rc'"
log close
exit, clear
