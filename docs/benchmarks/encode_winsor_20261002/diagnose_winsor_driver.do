clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/diagnose_winsor_transport.log", text replace
cd "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/tests"
capture noisily do /Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/diagnose_winsor_transport.do
local rc=_rc
if $TESTS_FAILED > 0 local rc=9
di "TEST_RC=`rc'"
log close
exit, clear
