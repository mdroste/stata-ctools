* Debug IPM step by step
clear all
set more off
adopath + "../build"

sysuse auto, clear
keep price mpg

* Run both and show verbose output
di _newline "{hline 78}"
di "QREG"
di "{hline 78}"
qreg price mpg

di _newline "{hline 78}"
di "CQREG (verbose)"
di "{hline 78}"
cqreg price mpg, verbose

* Show the results
matrix list e(b)
