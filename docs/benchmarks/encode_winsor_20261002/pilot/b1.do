clear all
set more off
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/pilot/b1.log", text replace
capture noisily do "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/pilot/benchmark.do" "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/baseline" b1 "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/pilot" "e10 e1000 e10000 e50000 e_unicode e_filtered e_existing" 2 "8" 200000 "2"
local rc=_rc
di "ENCODE_BENCHMARK_RC=`rc'"
log close
exit, clear
