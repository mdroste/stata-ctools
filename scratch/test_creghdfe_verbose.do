* Test script for creghdfe debugging - verbose mode
* Created: 2026-01-06

clear all
set more off

* Add build directory to adopath so Stata can find ctools commands
adopath + "../build"

* Load test data
sysuse auto

* ===========================================================================
* Test 1: Single FE (foreign)
* ===========================================================================
di as text _n "{hline 70}"
di as text "TEST 1: Single FE - absorb(foreign)"
di as text "{hline 70}"

di as text _n "=== reghdfe ==="
reghdfe price mpg, absorb(foreign)

di as text _n "=== creghdfe ==="
creghdfe price mpg, absorb(foreign)

* ===========================================================================
* Test 2: Two FEs (foreign trunk) - will have singletons
* ===========================================================================
di as text _n "{hline 70}"
di as text "TEST 2: Two FEs - absorb(foreign trunk)"
di as text "{hline 70}"

di as text _n "=== reghdfe ==="
reghdfe price mpg, absorb(foreign trunk)

di as text _n "=== creghdfe ==="
creghdfe price mpg, absorb(foreign trunk)

di as text _n "Test complete"
