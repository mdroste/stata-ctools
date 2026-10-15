clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/test_winsor_transport.log", text replace
cd "/Users/Mike/Documents/GitHub/stata-ctools/temp/encode_performance_20261002/tests"
capture noisily do validation/validate_winsor_transport.do
local rc=_rc
if $TESTS_FAILED > 0 local rc=9
di "TEST_RC=`rc'"
log close
exit, clear
