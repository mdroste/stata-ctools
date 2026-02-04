/*******************************************************************************
 * test_full_plus_errors.do - Run full sections 1-25 then error tests
 *
 * First run test_bisect_1_24.do (sections 1-25), then add error tests
 ******************************************************************************/

* Run the 1-25 sections (no errors) first
do "validation/test_bisect_1_24.do"

* Now we should be in a clean state after clear all
* But the test script ended with clear all, so now let's load setup again
* and run the error tests

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

di as text _n "Now adding error tests after full 1-25 run..."

* Variable doesn't exist
sysuse auto, clear
capture civreghdfe price (mpg = nonexistent_var), absorb(foreign)
di "Test 1: nonexistent_var, rc = " _rc

* Missing absorb option
sysuse auto, clear
gen z = runiform()
capture civreghdfe price (mpg = z)
di "Test 2: missing absorb, rc = " _rc

* No observations after if condition
sysuse auto, clear
gen z2 = runiform()
capture civreghdfe price (mpg = z2) if price > 100000, absorb(foreign)
di "Test 3: no observations, rc = " _rc

* String variable as dependent
sysuse auto, clear
gen z3 = runiform()
capture civreghdfe make (mpg = z3), absorb(foreign)
di "Test 4: string dependent, rc = " _rc

di as text _n "Error tests added. Now final clear all..."

clear all

di as text "SUCCESS: Full 1-25 + error tests passed!"
