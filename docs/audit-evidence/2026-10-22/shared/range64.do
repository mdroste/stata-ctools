clear all
set more off
adopath ++ /Users/mike/Documents/GitHub/stata-ctools/package/ado
adopath ++ /tmp/ctools-audit-oct22/build
capture log close _all
log using /tmp/ctools-audit-oct22/range64.log, text replace
set obs 64
gen double x = _n-1
replace x=1e16 in 1
gen key = _n
crangestat (mean) cm=x (sum) cs=x, interval(key 0 0)
format x cm cs %22.15g
list in 1/10, noobs
count if cm != x
di "WRONG_MEANS=" r(N)
count if cs != x
di "WRONG_SUMS=" r(N)
log close
exit, clear
