/*******************************************************************************
 * test_error_only.do - Test ONLY the error conditions section
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

di as text "Testing error conditions only then clear all..."

* Variable doesn't exist
sysuse auto, clear
capture civreghdfe price (mpg = nonexistent_var), absorb(foreign)
di "Test 1: nonexistent_var, rc = " _rc

* Missing absorb option
sysuse auto, clear
gen z = runiform()
capture civreghdfe price (mpg = z)
di "Test 2: missing absorb, rc = " _rc

* No instruments (underidentified)
sysuse auto, clear
capture civreghdfe price mpg, absorb(foreign)
di "Test 3: no instruments, rc = " _rc

* No observations after if condition
sysuse auto, clear
gen z2 = runiform()
capture civreghdfe price (mpg = z2) if price > 100000, absorb(foreign)
di "Test 4: no observations, rc = " _rc

* String variable as dependent
sysuse auto, clear
gen z3 = runiform()
capture civreghdfe make (mpg = z3), absorb(foreign)
di "Test 5: string dependent, rc = " _rc

di as text _n "Error tests complete. Now trying clear all..."

clear all

di as text "SUCCESS: No crash after error tests!"
