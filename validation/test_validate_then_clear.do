/*******************************************************************************
 * test_validate_then_clear.do - Run civreghdfe validation then clear all
 ******************************************************************************/

* Run the full validation
do "validation/validate_civreghdfe.do"

* Now clear all - this should trigger the crash if the bug exists
di as text _n "Now running: clear all"
clear all

di as text "SUCCESS: No crash occurred!"
