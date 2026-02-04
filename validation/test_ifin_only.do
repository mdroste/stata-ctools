/*******************************************************************************
 * test_ifin_only.do - Minimal test for if/in filtering bug
 ******************************************************************************/

capture do "validation/validate_setup.do"
if _rc != 0 {
    do "validate_setup.do"
}

di as text "Testing if/in conditions only..."

* Test 1: if condition
sysuse auto, clear
di "Test 1: if condition"
civreghdfe price (mpg = weight) if price > 5000, absorb(foreign)
di "  N = " e(N)

* Test 2: in condition
sysuse auto, clear
di "Test 2: in condition"
civreghdfe price (mpg = weight) in 1/50, absorb(foreign)
di "  N = " e(N)

* Test 3: Another if condition
sysuse auto, clear
di "Test 3: if foreign==1"
civreghdfe price (mpg = weight) if foreign==1, absorb(rep78)
di "  N = " e(N)

di as text "All if/in tests completed."
