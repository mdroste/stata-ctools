* Debug residual computation
clear all
set more off

adopath + "../build"

sysuse auto, clear

di "Number of variables before: " c(k)
describe, simple

* Run creghdfe with resid and verbose
creghdfe price mpg weight, absorb(foreign trunk) resid verbose

di _n "Number of variables after: " c(k)
describe, simple

di _n "After creghdfe:"
capture describe _reghdfe_resid
summarize _reghdfe_resid
list _reghdfe_resid in 1/10

di _n "e(resid) = " "`e(resid)'"
