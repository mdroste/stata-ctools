clear all
set more off
set linesize 220
capture log close _all
log using /tmp/ctools-audit-oct22/components.log, text replace
foreach component in csort cmerge csample cbsample cencode cdecode cdestring cipolate csplit crangejoin {
 capture noisily do validation/validate_`component'.do
 local rc = _rc
 noi di "AUDIT_COMPONENT=`component' RC=`rc' PASSED=$TESTS_PASSED FAILED=$TESTS_FAILED SKIPPED=$TESTS_SKIPPED COMPLETE=$CTOOLS_COMPONENT_COMPLETE"
}
log close
exit, clear
