* Test IPM implementation with verbose output
clear all
set more off
adopath + "../build"

sysuse auto, clear

di as text "{hline 70}"
di as text "Testing IPM with verbose output"
di as text "{hline 70}"

* First show qreg result for reference
di _newline
di as result "=== Stata's qreg (target) ===" _newline
qreg price mpg weight, nolog
scalar qreg_obj = e(sum_adev)
di "Sum of abs deviations: " qreg_obj

* Run cqreg with verbose
di _newline
di as result "=== ctools cqreg (IPM with debug) ===" _newline
cqreg price mpg weight, verbose

di as text "{hline 70}"
