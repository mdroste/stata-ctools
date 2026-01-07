* Test script for creghdfe debugging
* Created: 2026-01-06

clear all
set more off

* Add build directory to adopath so Stata can find ctools commands
adopath + "../build"

* Load test data
sysuse auto

* Display dataset info
describe
summarize price mpg foreign

* Try the failing command
di as text "Running: creghdfe price mpg, absorb(foreign)"
creghdfe price mpg, absorb(foreign)

di as text "Test complete"
